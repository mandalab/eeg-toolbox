function [events] = extractStimEvents_v2(subject, sessName, behDir, eegDir)
%
%  eventsExtraction code for White Noise and Stim Mapping driven by Tim's GUI on the stim computer
%    stimlogs are generated on the stim comptuer, and they must be moved to subject/behavioral/stimMapping for this code to work
%
%    this code takes the stimLog and converts it into an events.mat?
%
%
%
% CS or TS created?
%
% JHW edited 7/12/2017 so it works with behavioral processing
%


%%- Thia uses the following to run locally
% rootDir = '/Volumes/Shares/FRNU/dataWorking/eeg/';
% funDir = '/Volumes/FRNU/people/steinhardt/Usefulmatfun';
% labfunDir = '/Users/steinhardtcr/Documents/eeg_toolbox_trunk';
% addpath(genpath(funDir));
% addpath(genpath(labfunDir));



dbl_DC = 1;
%typical annotations
annotation_cell = {'Counting','Reading','Naming','Pointing','Other','AD','Seizure'};


%%- Old way, use CD (yuk!)... new way, use paths
sessDir  =  [behDir '/' sessName];      %- JW tweak --- bad practice to use CD
fName    = dir([sessDir '/*stimlog']);
stimFile = [sessDir '/' fName.name];



%- open  the "stimlog" file and get the info
fid = fopen(stimFile,'r');
header = textscan(fid,'%s%s%s%s%s%s%.1f%s%s',1,'headerLines',1); %- get the the file header
fclose(fid);
GUI_ver = header{7};
[offset sites amp1 amp2 freq pulses type annotation]=textread(stimFile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t');  %- get the file contents


%- clean the output
sites=sites(4:end);
sitesbyName{size(sites,1)} = [];


%- replace function "getjckOrder" with its guts -- JW
jackMaster_file = fullfileEEG(behDir, '../../docs/jacksheetMaster.csv');
if ~exist(jackMaster_file,'file'), fprintf('\n ERROR: cant find jacksheet master at %s', jackMaster_file); events=[]; return; end
[jackChans, jackNames] = textread(jackMaster_file,'%s\t %s');
AllJackChanInfo=cellfun(@(x) strsplit(x,','),jackChans,'UniformOutput',0);

chNameId =find(strcmp(AllJackChanInfo{1},'chanName'));
cnt = 1;  AllTAGS={};
for aJs = 2:length( AllJackChanInfo)
    if strcmp(AllJackChanInfo{aJs}{chNameId}(1:2),'DC')
    else
        AllTAGS{cnt}=AllJackChanInfo{aJs}{chNameId};
        cnt =  cnt +1;
    end
end
%%- END GUTS OF "getjckOrder"


% get all the site pair names
sitepairs= unique(sites);
for n = 1:length(sitepairs)
    sitenums= strsplit(sitepairs{n},'-');
    stimNames=strcat([AllTAGS{str2num(sitenums{1})} '-' AllTAGS{str2num(sitenums{2})}]);
    [sitesbyName{find(strcmp(sites,sitepairs{n}))}] = deal(stimNames);
end
sitesbyName = sitesbyName(1:size(sites,1));
%% assume all NK uses channels 5 and 6
NKnames=strcat([AllTAGS{(5)} '-' AllTAGS{(6)}]);% pretty sure always use if it's not said



% make events structure
idx= 4;evcnt = 1; tot_pulses = 0;
events = struct([]); 
while idx < size(offset,1)
    events(evcnt).subject    = subject;
    events(evcnt).exptFolder = sessName;
    events(evcnt).amplitude  = str2num(amp1{idx});
    events(evcnt).frequency  = str2num(freq{idx});
    events(evcnt).stimType   = type{idx};
    events(evcnt).stimLocTag = sitesbyName{evcnt};
    events(evcnt).annotation = annotation{idx};
    events(evcnt).DCpulses   = 0;   % initialize as no pulses
    
    %- initialize these so all version of task (white noise, mapping, CCEP) have same event structure
    events(evcnt).tele        = [];
    events(evcnt).cumulative_pulses  = [];
    events(evcnt).NKreference = [];
    events(evcnt).pulses      = [];
    events(evcnt).offset      = [];
    events(evcnt).mstime      = [];
    events(evcnt).monopolar   = [];
    events(evcnt).isWN_REP    = [];
    events(evcnt).WN_REP_obj  = [];
    
    
    %% getting the dc pulse count
    if strncmp(type{idx},'SM',2)
        if GUI_ver > 5          % correction for extra DC pulses due to 1s max zap protection
            events(evcnt).DCpulses = floor(str2num(pulses{idx})/events(evcnt).frequency);
        else
            events(evcnt).DCpulses = 1* dbl_DC;
        end
        
    elseif strncmp(type{idx},'CC',2)
        events(evcnt).monopolar = str2double(type{idx}(3));
        events(evcnt).DCpulses = str2num(pulses{idx});
        
    elseif strncmp(type{idx},'WN_REP',6)
        %events(idx).stimType = 'WN';
        events(evcnt).DCpulses = 1;
        if length(type{idx})>6 %- annotation mode
            pieces = strsplit(type{idx},'_');
            isDone = strcmp(pieces{3},'DONE');
            WN_obj.start   = ~isDone;
            WN_obj.mod_num = str2double(pieces{3+isDone});
            WN_obj.tot_mod = str2double(pieces{4+isDone});
            WN_obj.type    = pieces{5+isDone};
            WN_obj.completed = isDone;
            events(evcnt).DCpulses = 0;
        end
        events(evcnt).isWN_REP = true;
        events(evcnt).WN_REP_obj = WN_obj;
        if str2num(amp1{idx}) < 0
            events(evcnt).monopolar = 0; % 0 is negative 1 is positie
        else
            events(evcnt).monopolar = 1;
        end
        
    elseif strncmp(type{idx},'WN',2)
        % events.type = 'WN';
        events(evcnt).DCpulses = 1;
        if str2num(amp1{idx}) < 0
            events(evcnt).monopolar = 0; % 0 is negative 1 is positie
        else
            events(evcnt).monopolar = 1;
        end
        
    elseif strncmp(type{idx},'Switching',9)
        events(evcnt).stimType = 'SWITCH';
        tmp = strsplit(type{idx},'_');
        elec_loc{1} = tmp{2};
        elec_loc{2} = tmp{3};
        events(evcnt).elec_loc = elec_loc;
        
    elseif any(strcmp(type{idx},annotation_cell))
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
        
    elseif strncmp(type{idx},'Telemark_',9)
        events(evcnt).stimType = 'INFO';
        telemark_mode = strncmp(fliplr(type{idx}),fliplr('_ON'),3);
        events(evcnt).tele = telemark_mode;
        
    else
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
    end
    
    tot_pulses = tot_pulses+events(evcnt).DCpulses;
    events(evcnt).cumulative_pulses = tot_pulses;
    %%
    
    
    curMatTime = str2num(offset{idx})/1e11;
    t          = addtodate(curMatTime, 4, 'hour');
    tmpdate    = datetime(t,'ConvertFrom','datenum');
    curTimePos = num2str(round(posixtime(tmpdate)*1000));
    
    events(evcnt).NKreference = NKnames;
    events(evcnt).pulses = str2num(pulses{idx}); %- JW converted this to num 8/2017
    events(evcnt).offset = str2num(offset{idx}); %- JW converted this to num 8/2017
    events(evcnt).mstime = str2num(curTimePos);  %- JW converted this to num 8/2017
    evcnt = evcnt+1;
    idx   = idx+1;
end



%%- EVENTS file should saved by calling function, behavioralProcessing() -- JW commented out the following two lines
%cd([rootDir subject '/behavioral/stimMapping/' sessName ]);
%save('eventstmp','events')


%%- check for 3 different files

if exist([sessDir '/eegDC09.eeglog'],'file') hasDC09=1; else hasDC09=0; end;
if exist([sessDir '/eegDC10.eeglog'],'file') hasDC10=1; else hasDC10=0; end;
if exist([sessDir '/eeg.eeglog'],'file')     hasELOG=1; else hasELOG=0; end;



%%- make the alignment file for later based on the dc and cumulative pulse count:
if hasDC10==0 & length(events)>0,
    pulseTimes = find(([events.DCpulses]>0)); % the events that aren't just information and accounting for multple frequencies of spiking
    eegLogFile = [sessDir '/eegDC10.eeglog'];
    if ~exist('eegLogFile','file')
        fileId =fopen(eegLogFile,'w');        % for now just replacing the eeg.eeglog.up file
        for n = 1:length(pulseTimes)
            fprintf(fileId,'%s\t%d\t%s\n',events(pulseTimes(n)).mstime,0,'CHANNEL_0_UP');
        end
        fclose(fileId);
    end
end

%%- while we are at it, look for pulse file in the same directory and convert that as well
if hasDC09==0,
    fNamePulse = dir([sessDir '/*pulselog']);
    if length(fNamePulse)>0,
        
        %- create an eeg.log with python/unix time for prepAndAlign
        eegLogFile = [sessDir '/eegDC09.eeglog'];
        fileId =fopen(eegLogFile,'w');      % for now just replacing the eeg.eeglog.up file
        
        %- open  the "pulse" file, convert to python/unit style times, and save out as eeg.eeglogDC09
        pulseFile  = [sessDir '/' fNamePulse.name];
        [timePCstr pulseStr]=textread(pulseFile,'%s\t %s', 'headerlines',1);  %- get the file contents
        
        %- now loop through all the "up" times and copy them to eeg.log with corrected mstime
        iUp = find(strcmp(pulseStr,'PULSE_HI'))';
        for ii=iUp,
            curMatTime = str2num(timePCstr{ii})/1e11;
            t          = addtodate(curMatTime, 4, 'hour');
            tmpdate    = datetime(t,'ConvertFrom','datenum');
            curTimePos = num2str(round(posixtime(tmpdate)*1000));
            
            fprintf(fileId,'%s\t%d\t%s\n',curTimePos,0,'CHANNEL_0_UP');
        end
        fclose(fileId);
        
        hasDC09 = 1;
    end
end %- hasDC09

%- now pick from the two logs to create a straight up "eeg.eeglog" for alignment... DC09 is the first choice
%-  if the proper log already exists dont mess with it
if hasELOG==0,
    if hasDC09,
        sourceStr = 'eegDC09.eeglog';
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile([sessDir '/' sourceStr],[sessDir '/eeg.eeglog'],'f');
        if SUCCESS, fprintf('\n  %s -->  eeg.eeglog is a copy of %s',sessDir,sourceStr);
        else        fprintf('\n   Uh oh'); keyboard; end
        
    elseif hasDC10,
        fprintf('\n WARNING: no DC09 pulse file found in %s, must use DC10 to align', sessDir);
        
        sourceStr = 'eegDC10.eeglog';
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile([sessDir '/' sourceStr],[sessDir '/eeg.eeglog'],'f');
        if SUCCESS, fprintf('\n  %s -->  eeg.eeglog is a copy of %s',sessDir,sourceStr);
        else        fprintf('\n   Uh oh'); keyboard; end
    
    else
        fprintf('\n Uh Oh: no log files... wont be able to align this!  Is it a valid session?');
        keyboard;
    end
end


%%- make a "fake" session log with entries matching stimlog so if there is an alignment issue we can figuer out which line is the culprit
%%- always remake the session.log... this can be used to debug alignment errors with events.mat
sessLog = [sessDir '/session.log'];
fileId = fopen(sessLog,'w'); %- always make a fresh copy of this
for iEv = 1:length(events),
    fprintf(fileId,'%d\t%d\t%d\t%s\n',events(iEv).mstime ,0, events(iEv).offset, events(iEv).stimLocTag);
end
fclose(fileId);
 


end  %- end main()
