function [events] = createStimEvents_v4(sessLogFile, subject, sessionName, sessionNum,eegDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from stimMapping session. %%%%
%%%%
%% TS 12/2/16
% with new version, there is no reason include annoations with events file,
% instead this function will make events file out of DC10 file and thats it

%% JW 9/2017
%  I think tim used this for the "in-between" subjects before the GUI output got locked down, so NIH040 and possibly some others
%  renamed from "createStimEvents_v4" to "extractStimEventsFromOGstimLog_v10.m
%
%  REQUIREMENTS:   stimMapGUI/session_?/XXXX.stimLog   AND  trigDC10.upDown.txt
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%- Tims old method of telling the code below to look for DC11... might have to put this back
%chkdc11 = strrep(sessLogFile,'session.log','11.txt');
%if exist(chkdc11,'file')

tryFileList = { 'trigDC11.updown.txt', 'trigDC10.updown.txt'};  dc10filename = '';
for iF=1:length(tryFileList),
    if exist(fullfileEEG(fileparts(sessLogFile),tryFileList{iF}),'file'),
        dc10filename = fullfileEEG(fileparts(sessLogFile),tryFileList{iF});
    end
end
if isempty(dc10filename),
    fprintf('\n   And no trigDC10.updown.txt file available to create one! Rerun eegStimPrep or create session.log from another trigDC file.\n Exiting\n\n');
    events = [];
    return ;
end



%- create session.log by copying the contents of trigDC10.  If John C was annotating, the file will be automatically modified to facilitiate extraction
fin  = fopen(dc10filename,'r'); %if loops check to see that sessFile exists
fout = fopen(sessLogFile,'w+');
modFile=0;
while ~feof(fin)
    s = fgetl(fin);
    %         sMod = strrep(s, 'ANNOTATE 	 c ', 'ELECTRODES     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
    %         sMod = strrep(s, 'ANNOTATE 	 i', 'STIM_LEVEL     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
    %         sMod = strrep(s, 'ANNOTATE 	 b', 'BEHAVIORAL     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
    %         sMod = strrep(s, 'ANNOTATE 	 o', 'BEHAVIORAL     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
    %         sMod = strrep(s, 'ANNOTATE 	 ccep', 'CCEP_HERE     ');  if ~strcmp(s,sMod), s=sMod; modFile=modFile+1; end
    fprintf(fout,'%s\n',s);
    
end
fclose(fin);
fclose(fout);

fprintf('\n  session.log created from a copy of %s; \n %d modifications made based on John Cs annotation style \n ',dc10filename, modFile);
%     if modFile==0,
%         fprintf('\n  *** Zero automatic modifications were made: must do it by hand ***');
%         fprintf('\n  *** Open new session.log and hand-edit so annotations occur within (or before?) pulses. ***');
%         fprintf('\n  *** if subject>= NIH024, do find/replace:\n   "ANNOTATE 	 i" --> "STIM_LEVEL     "\n   "ANNOTATE 	 c" --> "ELECTRODES     "\n   "ANNOTATE     b" --> "BEHAVIORAL     "');
%         fprintf('\n  *** paused in createStimEvents... edit file and continue to run for first pass extraction.\n\n');
%         fclose('all');
%         %keyboard;
%     end

fid = fopen(sessLogFile,'r'); %- open the newly minted session.log file and attempt to process it below






%-initialize variables
serverDataPath = '/Volumes/Shares/FRNU/data/eeg/';      %- always point eegfile to this location
%serverDataPath = '/Volumes/Macintosh HD 2/STIM_MAP_DATA'; warning('Directing Data to local version!');

experiment    = 'stimMapping';
subject       = subject   ;
sessionName   = sessionName ;
sessionNum    = sessionNum ;
pulseFilename = '?';
fixCount = 0;

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
%- no longer necessary
% isCCEP        = 0;

warnCount=0;

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
            if inPulse==1, fprintf('ERROR: pulse up, but was already up'); keyboard; end
            inPulse = 1;
            pulseStartTime = mstime;
            pulseStopTime  = -1;
            stimResponse   = '';    %-reset this param at begining of pulse
            annotation     = '';
            
        case 'PULSE_LO'
            if inPulse==0, fprintf('ERROR: pulse down, but wasnt regsitered as being up'); keyboard; end
            inPulse = 0;
            pulseStopTime  = mstime;
            isPulseEvent   = 1;
            index = index+1;        %-pulse complete... create info that has all its info... after all events created will double events to create ups and down
            
            %         case 'ELECTRODES'
            %             error !
            %             electrodePair = info;
            %             annotation    = info;
            %             isChannelChange = 1;
            %             index = index+1;
            %
            %
            %         case 'STIM_LEVEL'
            %             error !
            %             stimulusLevel = sscanf(info,'%d');
            %             annotation    = info;
            %
            %         case 'STIM_RESPONSE'
            %             error !
            %             stimResponse  = info;
            %             annotation    = info;
            %             errror('Q');
            %
            %         case 'BEHAVIORAL'
            %             stimResponse  = info;
            %             annotation    = info;
            %
            %         case 'ANNOTATE'
            %             %- user should replace ANNOTATION with ELECTRODES, STIM_LEVEL, RESPONSE... or make it automated here?
            %             annotation    = info;
            %             if inPulse==0, index=index+1; fprintf('Warning: annotation at %d not within stim pulse: "%s" \n', mstime, annotation); fixCount = fixCount +1; end
            %
            
            %
            %         case 'CCEP_HERE'
            %             annotation = type;
            %             stimulusLevel = sscanf(info,'%d');
            %             isCCEP = 1;
        case 'SESS_END'
            %- user should replace ANNOTATION with ELECTRODES, STIM_LEVEL, RESPONSE... or make it automated here?
            index=index+1;
            
        otherwise
            if warnCount<5,
                fprintf('WARNING: %s not expected type\n',type);
            elseif warnCount==5,
                fprintf('WARNING: %s not expected type <supressing any further warnings about expected types>\n',type);
            end
            warnCount = warnCount+1;
            
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
        %thisEvent.stimResponse      = stimResponse  ;
        thisEvent.annotation        = annotation    ;
        %thisEvent.isCCEP            = isCCEP       ;
        %isCCEP = 0;
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
    
    fprintf('\n   eeg.eeglog created with fake start/stop times: %s ', eegLogFile);;
    
    %- sanity check
    %[mstimes]   = textread(eegLogFile,'%n%*[^\n]');   %read all ms time from the pulse file
    %dateMAT_act = datenum(epoch2date(mstimes(1)));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
    %dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    %dateMAT_actEnd = datenum(epoch2date(mstimes(end))); % *MST switch to epoch2date
    %dateStr_actEnd = datestr(dateMAT_actEnd,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
end

