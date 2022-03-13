function [events] = createStimEvents_v3a(sessLogFile, subject, sessionName, sessionNum,eegDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from stimMapping session. %%%%
%%%%
% update 1/10/16
% to run ensure stimMapping folder is filled and sessions are populated
% with DC10updown. If not, see below. Run this function to create events,
% out of DC10- there will be errors so after running- correct annotations
% so that they always begin with the correct format-- pretty easy, each
% annoation needs to start with a lowercase 'i ','b ','c ','o ' or 'ccep '.
% Run again until no errors exist.


%%%%    how to use this code:
%%%%     1) use eegStimPrep to create trigUpDown and annotation files...
%%%%     2) create behavioral/stimMapping/session_X folder, copy updown file there
%%%%     3) rename behavioral updown file with stim on/off pulses "session.log"
%%%%     4) hand-edit session.log file so annotations occur within pulse (or before?), change ANNOTATION to ELECTRODES where applicable
%%%%         (Cocjin's annotations (>=NIH024) easy to tweak... do find/replace  "ANNOTATE 	 i" --> "STIM_LEVEL     ",   "ANNOTATE    c" --> "ELECTRODES     "
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear
% subject     = 'NIH035';
%
% %rootEEGdir = '/Volumes/Macintosh HD 2/STIM_MAP_DATA'; %local
% rootEEGdir = '/Volumes/Shares/FRNU/dataWorking/eeg/';  %dataWorking
%
% sessionName = 'session_3'; sessionNum = 3;
% %
% sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/stimMapping',sessionName);
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fid = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    
    fprintf('\n\n **** Session.log not found: %s  ****', sessLogFile);
    
    %- not sure what this was intended to do...
    %chkdc11 = strrep(sessLogFile,'session.log','11.txt');
    %if exist(chkdc11,'file')
    
    tryFileList = { 'annotation.txt', 'trigDC11.updown.txt', 'trigDC10.updown.txt'};  dc10filename = '';
    for iF=1:length(tryFileList),
        if exist(fullfileEEG(fileparts(sessLogFile),tryFileList{iF}),'file'),
            dc10filename = fullfileEEG(fileparts(sessLogFile),tryFileList{iF});
        end
    end
    if isempty(dc10filename),
        fprintf('\n And no trigDC10.updown.txt file available to create one! Rerun eegStimPrep or create session.log from another trigDC file.\n Exiting\n\n');
        events = [];
        return ;
    end
    
    
    %- create session.log by copying the contents of trigDC10.  If John C was annotating, the file will be automatically modified to facilitiate extraction
    fin  = fopen(dc10filename,'r'); %if loops check to see that sessFile exists
    fout = fopen(sessLogFile,'w+');
    modFile=0;
    while ~feof(fin)
        s = fgetl(fin);
        sMod = strrep(s, 'ANNOTATE 	 c ',   'ELECTRODES     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 i',    'STIM_LEVEL     ');  if ~strcmpi(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 b',    'BEHAVIORAL     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 o',    'STIM_RESPONSE  ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 ccep', 'CCEP_HERE      ');  if ~strcmpi(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 CCEP', 'CCEP_HERE      ');  if ~strcmpi(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 s AD', 'AFTER_DISCHARGE ');  if ~strcmpi(s,sMod), s=sMod; modFile=modFile+1; end
        sMod = strrep(s, 'ANNOTATE 	 s sz', 'SEIZURE        ');  if ~strcmpi(s,sMod), s=sMod; modFile=modFile+1; end
        fprintf(fout,'%s\n',s);
    end
    fclose(fin);
    fclose(fout);
    
    
    fprintf('\n  session.log created from a copy of %s; \n %d modifications made based on John Cs annotation style \n ',dc10filename, modFile);
    if modFile==0,
        fprintf('\n  *** Zero automatic modifications were made: must do it by hand ***');
        fprintf('\n  *** Open new session.log and hand-edit so annotations occur within (or before?) pulses. ***');
        fprintf('\n  *** if subject>= NIH024, do find/replace:\n   "ANNOTATE 	 i" --> "STIM_LEVEL     "\n   "ANNOTATE 	 c" --> "ELECTRODES     "\n   "ANNOTATE     b" --> "BEHAVIORAL     "');
        fprintf('\n  *** paused in createStimEvents... edit file and continue to run for first pass extraction.\n\n');
        fclose('all');
        
    end
    
    
    fid = fopen(sessLogFile,'r'); %- open the newly minted session.log file and attempt to process it below
end


%- create and/or open the "fixLog"... a text file that saves information about potential errors in the events.mat
errorLog = fullfileEEG(fileparts(sessLogFile),'events.FOUND_ERRORS.txt');
fid_error = fopen(errorLog,'w+');
fprintf(fid_error,'\n ERROR LOG for stimMapping events.mat created from manually entered annotations');
fprintf(fid_error,'\n     events created from %s/session.log',sessionName);
fprintf(fid_error,'\n     using createStimEvents_v3a.m on %s\n\n',datestr(now));




%-initialize variables
serverDataPath = '/Volumes/Shares/FRNU/data/eeg/';      %- always point eegfile to this location
%serverDataPath = '/Volumes/Macintosh HD 2/STIM_MAP_DATA'; warning('Directing Data to local version!');

experiment    = 'stimMapping';
subject       = subject   ;
sessionName   = sessionName ;
sessionNum    = sessionNum ;
pulseFilename = '?';
fixCount      = 0;

%-params that vary
type          = '';
mstime        = -1;
eegfile       = '';

%-pulse params
isPulseEvent  = -1;
pulseStartTime = -1;
pulseStopTime  = -1;
pulseDuration  = -1;
electrodePair = '?';
stimulusLevel = NaN;
stimResponse  = '';
annotation    = '';

%-
inPulse       = 0;
isCCEP        = 0;
isICTAL = 0;



%- Read session.log line-by-line and convert to events structure
events      = struct([]);
index       = 0;
while true
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); break; end
    
    
    %- Generic text scan to get time, offset, and type
    [xTOT, pos]         = textscan(thisLine,'%d %s',1);
    mstime              = xTOT{1}(1);   %- must add (1) because numbers after the string in the above line cause overflow to first %f
    type                = xTOT{2}{1};
    info                = strtrim(thisLine(pos+1:end)); %trim leading/trailing spaces
    info                = info(find(info));             %trim empty (null) entries from end
    
    %-
    isPulseEvent        = 0; %-set to 0 every iteration... only PULSE_LO will set high (and increment index so event is created)
    isChannelChange     = 0; %-set to 0 every iteration... only PULSE_LO will set high (and increment index so event is created)
    
    %- default Parameters (details will be filled out/altered based on type)
    switch type
        case 'FILENAME'
            pulseFilename  = info;
            
        case 'EEGSTEM'
            eegfile        = fullfileEEG(serverDataPath,subject,'eeg.noreref',info);
            type  = 'SESS_START' ;
            index = index+1;        %-create session start event once the eegfile is logged
            
            
        case 'PULSE_HI'
            if inPulse==1, fprintf('\nERROR: pulse up, but was already up'); keyboard; end
            inPulse = 1;
            pulseStartTime = mstime;
            pulseStopTime  = -1;
            stimResponse   = '';    %-reset this param at begining of pulse
            annotation     = '';
            
        case 'PULSE_LO'
            if inPulse==0, fprintf('\nERROR: pulse down, but wasnt regsitered as being up'); keyboard; end
            inPulse = 0;
            pulseStopTime  = mstime;
            isPulseEvent   = 1;
            index = index+1;        %-pulse complete... create info that has all its info... after all events created will double events to create ups and down
            
        case 'ELECTRODES'
            electrodePair = info;
            annotation    = info;
            isChannelChange = 1;
            index = index+1;
            
            
        case 'STIM_LEVEL'
            stimulusLevel = sscanf(info,'%d');
            annotation    = info;
            
        case 'STIM_RESPONSE'
            stimResponse  = info;
            annotation    = info;
            
            
        case 'BEHAVIORAL'
            stimResponse  = info;
            annotation    = info;
            
        case 'ANNOTATE'
            %- user should replace ANNOTATION with ELECTRODES, STIM_LEVEL, RESPONSE... or make it automated here?
            annotation    = info;
            if inPulse==0, index=index+1; 
                fprintf(fid_error,'\nWarning: annotation at %d not within stim pulse: "%s" ', mstime, annotation);
                fprintf(          '\nWarning: annotation at %d not within stim pulse: "%s" ', mstime, annotation); 
                fixCount = fixCount +1; 
            end
            
        case 'SESS_END'
            %- user should replace ANNOTATION with ELECTRODES, STIM_LEVEL, RESPONSE... or make it automated here?
            index=index+1;
            
        case 'CCEP_HERE'
            annotation = type;
            stimulusLevel = sscanf(info,'%d');
            isCCEP = 1;
            
        case 'AFTER_DISCHARGE'
            annotation = 'AD';
            isICTAL = 1;
        case 'AFTER_DISCHARGE'
            annotation = 'SZ';
            isICTAL = 1;
        otherwise
            fprintf('\nWARNING: %s not expected type',type);
            
    end
    
    if isempty(stimulusLevel), stimulusLevel = NaN; end
    %- asign values to events array
    if index>length(events),
        
        clear thisEvent
        %-params that are fixed for the session
        thisEvent.experiment        = experiment  ;
        thisEvent.subject           = subject     ;
        thisEvent.sessionName       = sessionName ;
        thisEvent.sessionNum        = sessionNum  ;
        thisEvent.pulseFilename     = pulseFilename ;
        
        %-params that vary
        thisEvent.type              = type        ;
        thisEvent.mstime            = mstime      ;
        thisEvent.msduration        = nan         ;
        thisEvent.eegoffset         = mstime      ;
        thisEvent.eegfile           = eegfile     ;
        
        %-pulse params
        thisEvent.isPulseEvent      = isPulseEvent  ; %0-no, 1=stim start, 2=stim stop
        thisEvent.pulseStartTime    = pulseStartTime;
        thisEvent.pulseStopTime     = pulseStopTime ;
        thisEvent.pulseDuration     = pulseStopTime-pulseStartTime;
        thisEvent.isChannelChange   = isChannelChange;
        thisEvent.electrodePair     = electrodePair ;
        thisEvent.stimulusLevel     = stimulusLevel ;
        thisEvent.stimResponse      = stimResponse  ;
        thisEvent.annotation        = annotation    ;
        thisEvent.isCCEP            = isCCEP        ;
        thisEvent.isICTAL           = isICTAL       ;
        isCCEP = 0;
        isICTAL = 0;
        %-save to events structure array
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
        
        
        %-make stim start and stim stop events
        if isPulseEvent,
            %             thisEvent.type          = 'PULSE_LO';
            %             thisEvent.isPulseEvent  = 1;
            %             thisEvent.mstime        = pulseStartTime;
            %             thisEvent.msduration    = pulseStopTime-pulseStartTime;
            %             thisEvent.eegoffset     = pulseStartTime;
            %             events(index)           = thisEvent;
            
            thisEvent.type          = 'PULSE_HI';
            thisEvent.isPulseEvent  = 2;
            thisEvent.mstime        = pulseStopTime;
            thisEvent.msduration    = pulseStopTime-pulseStartTime;
            thisEvent.eegoffset     = pulseStopTime;
            index=index+1;
            events(index)           = thisEvent;
        end
    end
    
end
fclose(fid);  % close session.log


%%- now cleanup event timing... session
iChanChange  =  find( [events.isChannelChange]==1 );
for ii=1:length(iChanChange)-1,
    events(iChanChange(ii)).msduration     =  events(iChanChange(ii+1)).mstime - events(iChanChange(ii)).mstime;
    events(iChanChange(ii)).pulseStartTime =  nan;
    events(iChanChange(ii)).pulseStopTime  =  nan;
    events(iChanChange(ii)).pulseDuration  =  nan;
end
try
    events(iChanChange(end)).msduration     =  events(end).mstime - events(iChanChange(end)).mstime;
    events(iChanChange(end)).pulseStartTime =  nan;
    events(iChanChange(end)).pulseStopTime  =  nan;
    events(iChanChange(end)).pulseDuration  =  nan;
catch e
    % keyboard
end

%%- SESS_START duration should equal entire eeg timeseries
events(1).msduration = events(end).mstime;


%%- Create FAKE eeg.eeglog if one isn't there already
eegLogFile = fullfileEEG(fileparts(sessLogFile),'eeg.eeglog');
if ~exist(eegLogFile,'file'),
    %fprintf('\n\n **** eeg.eeglog not found: %s  ****', eegLogFile);
    %-
    fid = fopen(eegLogFile,'w+'); %
    [~,eegStem] = fileparts(eegfile);
    d = datetime(eegStem,'inputformat','yyMMdd_HHmm','TimeZone', 'local');
    oneMin = 0*1000;  
    fakeStartTime = posixtime(d)*1000 + oneMin;  %- tweak the start and end time to help prepAndAlign (looks like adding up to 2 min isnt helping..
    fakeEndTime   = fakeStartTime + double(events(1).msduration) - oneMin;
    if fakeEndTime<fakeStartTime, fprintf('\n uh oh... short session?'); keyboard; end
    fprintf(fid,'%d 0 fakeStartTime\n',fakeStartTime);
    fprintf(fid,'%d 1 CHANNEL_0_UP\n',fakeStartTime);
    fprintf(fid,'%d 0 fakeEndTime \n', fakeEndTime  );
    fclose(fid);
    
    fprintf('\n ****  eeg.eeglog created with fake start/stop times: %s  ****', eegLogFile);;
    
    %- sanity check
    %[mstimes]   = textread(eegLogFile,'%n%*[^\n]');   %read all ms time from the pulse file
    %dateMAT_act = datenum(epoch2date(mstimes(1)));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
    %dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    %dateMAT_actEnd = datenum(epoch2date(mstimes(end))); % *MST switch to epoch2date
    %dateStr_actEnd = datestr(dateMAT_actEnd,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
end


%- close the error file, rename if error free
fprintf(fid_error,'\n\nFix %d annotations in %s DC10 updown and rerun (this file will be automatically updated if error free)\n', fixCount, sessionName');
fclose(fid_error);
if fixCount,
    fprintf('Fix %d annotations in %s session.log (DC10 updown) and rerun\n', fixCount, sessionName'); 
else
    movefile(errorLog,fullfileEEG(fileparts(sessLogFile),'events.NO_ERRORS.txt'));
end



%%- only should be executed when running as script (not as function)
if exist('eventFile','var'),
    fprintf('running jwAttnTaskEvents directly: \n --> extracted %d events to %s\n', length(events), sessionDir);
    save(eventFile,'events');
end

