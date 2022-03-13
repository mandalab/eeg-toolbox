function [events] = extractStimEventsFromStimLog(subject, sessName, behDir, eegDir)
%function [events] = extractStimEvents_v2(subject, sessName, behDir, eegDir)  %-- old name
%
%  8/2017 -- JW renamed to extractStimEventsFromStimLog_v10 to differentiate from old annotation-based stim mapping sessions
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
% JHW 9/2017 -- incorporate Thia's check for phys2chan and montage
%            -- ouptut helper files that inculde the date and unix time string for each line of the pulseLog and stimLog
% 12/2017 melkalliny 


%typical annotations
dbl_DC = 1;
annotation_cell = {'Counting','Reading','Naming','Pointing','Other','AD','Seizure'};



%%- Old way, use CD (yuk!)... new way, use paths
sessDir  =  [behDir '/' sessName];      %- JW tweak --- bad practice to use CD
fName    = dir([sessDir '/*stimlog']);
if isempty(fName)
    fName    = dir([sessDir '/*stimLog']);
end
fName.name = lower(fName.name);
stimFile = [sessDir '/' fName.name];
try % JW, take a look? trying to keep backwards compatibility here
    fDateStr = datetime(fName.name(1:strfind(fName.name,'.stimlog')-1),'InputFormat','yyyy_MM_dd_HH_mm');
catch
    try
        fDateStr = datetime(fName.name(1:strfind(fName.name,'.stimlog')-1),'InputFormat','yyyy_MM_dd_HHmm');   
    catch
        fDateStr = datetime(fName.name(1:strfind(fName.name,'.stimlog')-1),'InputFormat','dd-MMM-yyyy_HH_mm_SS');
    end
end
fDateNum = datenum(fDateStr);

%- open  the "stimlog" file and get the info
fid = fopen(stimFile,'r');
header = textscan(fid,'%s%s%s%s%s%s%.1f%s%s',1,'headerLines',1); %- get the the file header
fclose(fid);
GUI_ver = header{7};
[offset sites amp1 amp2 freq pulses type annotation]=textread(stimFile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t');  %- get the file contents
    
%%% Following crazy-ass hack is used to figure out if hour offset is required to match the machine times of the stim and EEG computer
%%%   use the date string from the stimLog file and compare that to the date that is compputed from eeg.eeglog by prepAndAlign.  
%%%   JW tried doing this more simply, but seems like the process of converting to posixtime and back again is important for figuring out the correct offset
if length(offset) > 4
    curMatTime = str2num(offset{5})/1e11; %- sometimes first time point (at offset{4}) is invalid... so look at the second
else
    curMatTime = str2num(offset{end})/1e11;
end
addHourOpt = [0 1];
t          = [addtodate(curMatTime, 4+addHourOpt(1), 'hour') addtodate(curMatTime, 4+addHourOpt(2), 'hour')]; %- try with and without additional hour
tmpdate    = datetime(t,'ConvertFrom','datenum');
curTimePos = num2str(round(posixtime(tmpdate)*1000));
dateMAT_act = datenum(epoch2date(str2num(curTimePos))); 
dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
[y i] = min(abs(dateMAT_act-fDateNum));
addHour = addHourOpt(i);
fprintf('\n note: stimLog filename = %s --- eeg.eeglog starttime = %s (addHour=%d)',datestr(fDateStr,'mm/dd/yy HH:MM PM'), dateStr_act(i,:),addHour);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- now make a "helper" stimLog and pulseLog file that is a copy of the originals, but with dateStrings and posix time strings for each row
fid_orig = fopen(stimFile,'r');
fid_help = fopen([stimFile '.helper'],'w');  %- create a version that can help with alignment issues
while 1,
    tline = fgetl(fid_orig);
    if ~ischar(tline), break, end
    timeStr    = sscanf(tline,'%s',1);
    curMatTime = str2num(timeStr)/1e11; %- sometimes first time point (at offset{4}) is invalid... so look at the second
    if ~isempty(curMatTime),
        tmpdate    = datetime(addtodate(curMatTime, 4+addHour, 'hour'),'ConvertFrom','datenum');
        curTimePos = num2str(round(posixtime(tmpdate)*1000));
        dateMAT_act = datenum(epoch2date(str2num(curTimePos)));
        dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    else
        dateStr_act = 'Date String:     no date';
        curTimePos  = 'UNIX/POSIX Time: no posix';
    end
    %- get the
    fprintf(fid_help,'<<< %s (addedHour=%d);  %s >>> \t %s\n',dateStr_act,addHour,curTimePos,tline);
end
fclose(fid_orig);
fclose(fid_help);


%- now make a "helper"  pulseLog file that is a copy of the originals, but with dateStrings and posix time strings for each row
fNamePulse = dir([sessDir '/*pulselog']);
if length(fNamePulse)>0
    pulseFile = [sessDir '/' fNamePulse.name];
    fid_orig = fopen(pulseFile,'r');
    fid_help = fopen([pulseFile '.helper'],'w');  %- create a version that can help with alignment issues
    while 1
        tline = fgetl(fid_orig);
        %disp(tline)
        if ~ischar(tline), break, end
        timeStr    = sscanf(tline,'%s',1);
        curMatTime = str2num(timeStr)/1e11; %- sometimes first time point (at offset{4}) is invalid... so look at the second
        if ~isempty(curMatTime)
            tmpdate    = datetime(addtodate(curMatTime, 4+addHour, 'hour'),'ConvertFrom','datenum');
            curTimePos = num2str(round(posixtime(tmpdate)*1000));
            dateMAT_act = datenum(epoch2date(str2num(curTimePos)));
            dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
        else
            dateStr_act = 'Date String:     no date';
            curTimePos  = 'UNIX/POSIX Time: no posix';
        end
        %- get the
        fprintf(fid_help,'<<< %s (addedHour=%d);  %s >>> \t %s\n',dateStr_act,addHour,curTimePos,tline);
    end
    fclose(fid_orig);
    fclose(fid_help);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- clean the output
sites=sites(4:end);
sitesbyName{size(sites,1)} = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rootEEGdir = behDir(1:strfind(behDir,subject)-2);
USE_MONTAGE = 0;
if USE_MONTAGE
    %%- new version from Thia... reads the montage AND uses phys2name to compair
    %   but this ended up having a lot of errors and wont work

    montFile = fullfileEEG(rootEEGdir, subject, '/docs/clinical_montage.xlsx');
    if exist(montFile,'file')==0,
        error('\nDidnt find docs/clinical_montage.xls.  Make sure file is there and named correctly.')
    end
    %- reading montage doesn't work for every subject... better to skip this requirement
    %----  old way ----%
    %[num,text] = xlsread(montFile,'Montage');
    [num,text] = xlsread(montFile,1); %- naming not consistent... best to say "first sheet"
    text = text(1:size(num,1), 3:4:end); % names
    num  = num(:,1:4:end);%numbers
    %---- New way  -----%
    [all_index, all_jacks, all_label] = readMontage(subject,montFile); %- this doesn't work for all subjects... e.g., gets stuck on NIH047
end

%- just grab any EEG file in the stim map folder for phys2name... not sure what the actual one is at this point
stim_map_folders = lsCell(fullfileEEG(rootEEGdir,subject,'raw/STIM_MAP'));
if isempty(stim_map_folders)
    fprintf('\n uh oh... nothing in STIM_MAP folder... looking for EEG file in raw to run phys2name');
    raw_folders = lsCell(fullfileEEG(rootEEGdir,subject,'raw'));     
    if isempty(raw_folders)
        
        fprintf('\n uh oh... no folders in raw folder.  Cant process stimMapGUI till there is at least one raw file grabbed');
        keyboard;
        return;
    end
end
eegFldr = stim_map_folders{1};
% I refactored the above stim folder code to be *more* correct, 
% but I left the logic that only grabs the first stim_map_folder
% because I'm not sure why it's done; it seems wrong to me. -Mike 03/18
    

% get all the site pair names
sitepairs= unique(sites);
for n = 1:length(sitepairs)
    sitenums= strsplit(sitepairs{n},'-');
    % phys2name(subj, rootEEGdir, timestamp, physNum, montageLabels, montageJack)
    physname{1} = phys2name(subject,rootEEGdir,eegFldr,str2num(sitenums{1}));
    physname{2} = phys2name(subject,rootEEGdir,eegFldr,str2num(sitenums{2}));
    
    stimNames = [physname{1} '-' physname{2}];
    if USE_MONTAGE
        [row1,col1]=ind2sub(size(num),str2num(sitenums{1}));% find specific electrode in montage
        [row2,col2]=ind2sub(size(num),str2num(sitenums{2}));
    
        if ~strcmp(text(row1,col1),physname{1}) | ~strcmp(text(row2,col2),physname{2})
        warning('\nPhys2name and montage do not match up! \nAlways trust the montage and look up the number on physical channel list and montage by eye to see if jacksheet, etc is funky or there is a typo! \nIf you can, fix it ! Call Thia if still confused!');
        stimNames = [text(row1,col1) '-' text(row2,col2)];
        end        
    end
    [sitesbyName{find(strcmp(sites,sitepairs{n}))}] = deal(stimNames);
end
sitesbyName = sitesbyName(1:size(sites,1));

%% assume all NK uses channels 5 and 6
if USE_MONTAGE
    %pretty sure always reference 5 and 6. Might need to find way to adjust for what happens where are switched...
    [rrow1,rcol1]=ind2sub(size(num),5);% find specific electrode in montage for references
    [rrow2,rcol2]=ind2sub(size(num),6);
    NKnames=strcat([text(rrow1,rcol1) '-' text(rrow2,rcol2)]);% pretty sure always use if it's not said
else
    NKnames=strcat([phys2name(subject,rootEEGdir,eegFldr,5) '-' phys2name(subject,rootEEGdir,eegFldr,6)]);% pretty sure always use if it's not said
end


% %- replace function "getjckOrder" with its guts -- JW
% jackMaster_file = fullfileEEG(behDir, '../../docs/jacksheetMaster.csv');
% if ~exist(jackMaster_file,'file'), fprintf('\n ERROR: cant find jacksheet master at %s', jackMaster_file); events=[]; return; end
% [jackChans, jackNames] = textread(jackMaster_file,'%s\t %s');
% AllJackChanInfo=cellfun(@(x) strsplit(x,','),jackChans,'UniformOutput',0);
% 
% chNameId =find(strcmp(AllJackChanInfo{1},'chanName'));
% cnt = 1;  AllTAGS={};
% for aJs = 2:length( AllJackChanInfo)
%     if strcmp(AllJackChanInfo{aJs}{chNameId}(1:2),'DC')
%     else
%         AllTAGS{cnt}=AllJackChanInfo{aJs}{chNameId};
%         cnt =  cnt +1;
%     end
% end
% %%- END GUTS OF "getjckOrder"
% 
% 
% % get all the site pair names
% sitepairs= unique(sites);
% for n = 1:length(sitepairs)
%     sitenums= strsplit(sitepairs{n},'-');
%     stimNames=strcat([AllTAGS{str2num(sitenums{1})} '-' AllTAGS{str2num(sitenums{2})}]);
%     [sitesbyName{find(strcmp(sites,sitepairs{n}))}] = deal(stimNames);
% end
% sitesbyName = sitesbyName(1:size(sites,1));
% %% assume all NK uses channels 5 and 6
% NKnames=strcat([AllTAGS{(5)} '-' AllTAGS{(6)}]);% pretty sure always use if it's not said



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
    events(evcnt).elec_loc    = {};
    
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
        WN_obj = struct;
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
    t          = addtodate(curMatTime, 4+addHour, 'hour');
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
%cd([rootEEGdir subject '/behavioral/stimMapping/' sessName ]);
%save('eventstmp','events')


%%- check for 3 different files

hasDC09=0; hasDC10=0; %- force an update of these each time, just like the session log
%if exist([sessDir '/eegDC09.eeglog'],'file') hasDC09=1; else hasDC09=0; end;
%if exist([sessDir '/eegDC10.eeglog'],'file') hasDC10=1; else hasDC10=0; end;
if exist([sessDir '/eeg.eeglog'],'file')     hasELOG=1; else hasELOG=0; end;



%%- make the alignment file for later based on the dc and cumulative pulse count:
if hasDC10==0 & length(events)>0,
    pulseTimes = find(([events.DCpulses]>0)); % the events that aren't just information and accounting for multple frequencies of spiking
   
    eegLogFile  = [sessDir '/eegDC10.eeglog'];
    fileId      = fopen(eegLogFile,'w');        % for now just replacing the eeg.eeglog.up file
    
    %eegLogFile2 = [sessDir '/eegDC10.eeglog.helper']; %- create a version that can help with alignment issues
    %fileId2     = fopen(eegLogFile2,'w');      % for now just replacing the eeg.eeglog.up file
       
    for n = 1:length(pulseTimes)
        fprintf(fileId,'%s\t%d\t%s\n',num2str(events(pulseTimes(n)).mstime),0,'CHANNEL_0_UP');
        %fprintf(fileId2,'%s\t%d\t%s\t<%s>\n',num2str(events(pulseTimes(n)).mstime),0,'CHANNEL_0_UP',num2str(events(pulseTimes(n)).offset));
    end
    fclose(fileId);
    %fclose(fileId2);
    hasDC10 = 1;
end

%%- while we are at it, look for pulse file in the same directory and convert that as well
if hasDC09==0,
    fNamePulse = dir([sessDir '/*pulselog']);
    if length(fNamePulse)>0,
        
        %- create an eeg.log with python/unix time for prepAndAlign
        eegLogFile = [sessDir '/eegDC09.eeglog'];
        fileId =fopen(eegLogFile,'w');      % for now just replacing the eeg.eeglog.up file
        
        %eegLogFile2 = [sessDir '/eegDC09.eeglog.helper']; %- create a version that can help with alignment issues
        %fileId2 =fopen(eegLogFile2,'w');      % for now just replacing the eeg.eeglog.up file
        
        %- open  the "pulse" file, convert to python/unit style times, and save out as eeg.eeglogDC09
        pulseFile  = [sessDir '/' fNamePulse.name];
        [timePCstr pulseStr]=textread(pulseFile,'%s\t %s', 'headerlines',1);  %- get the file contents
        
        %- now loop through all the "up" times and copy them to eeg.log with corrected mstime
        iUp = find(strcmp(pulseStr,'PULSE_HI'))';
        for ii=iUp,
            curMatTime = str2num(timePCstr{ii})/1e11;
            t          = addtodate(curMatTime, 4+addHour, 'hour');
            tmpdate    = datetime(t,'ConvertFrom','datenum');
            curTimePos = num2str(round(posixtime(tmpdate)*1000));
            
            fprintf(fileId,'%s\t%d\t%s\n',curTimePos,0,'CHANNEL_0_UP');
            %fprintf(fileId2,'%s\t%d\t%s\t<%s>\n',curTimePos,0,'CHANNEL_0_UP',timePCstr{ii});
        end
        fclose(fileId);
        %fclose(fileId2);
        
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
