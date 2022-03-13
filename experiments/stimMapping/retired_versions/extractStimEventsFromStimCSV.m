function [events] = extractStimEventsFromStimCSV(subject, sessName, behDir, eegDir)
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

% Number of Header Lines
numHeaderLines = 6;

%%- Old way, use CD (yuk!)... new way, use paths
sessDir  =  [behDir '/' sessName];      %- JW tweak --- bad practice to use CD
fName    = dir([sessDir '/*csv']);
fName.name = lower(fName.name);
stimFile = [sessDir '/' fName.name];
try % JW, take a look? trying to keep backwards compatibility here
    fDateStr = datetime(fName.name(1:strfind(fName.name,'.csv')-1),'InputFormat','yyyy_MM_dd_HH_mm');
catch
    try
        fDateStr = datetime(fName.name(1:strfind(fName.name,'.csv')-1),'InputFormat','yyyy_MM_dd_HHmm');   
    catch
        fDateStr = datetime(fName.name(1:strfind(fName.name,'.csv')-1),'InputFormat','dd-MMM-yyyy_HH_mm_SS');
    end
end
fDateNum = datenum(fDateStr);

% Load CSV
fileCSV = fullfile(sprintf('%s/%s',sessDir,fName.name));
stimCSV = readtable(fileCSV);

GUI_ver = 7; % Hard coding to 7 for conditional statement below

%- sanity check that posix time in the file matches the datestring on the session folder
curTimePos = stimCSV.offset(2); %- numHeaderLines above jumps over the header... so just grab 1st or 2nd row of real data for this step
dateMAT_act = datenum(epoch2date(curTimePos));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
delta1hour  = datenum('01/01/01 2:00')-datenum('01/01/01 1:00');

addHourOpt = [0 1];
t          = [dateMAT_act addtodate(dateMAT_act, addHourOpt(2), 'hour')]; %- try with and without additional hour
tmpdate    = datetime(t,'ConvertFrom','datenum');
dateMAT_act = datenum(tmpdate); 
dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
timeDiff = dateMAT_act-fDateNum;
y = min(timeDiff(timeDiff > 0));
i = find(timeDiff == y,1);
addHour = addHourOpt(i);
fprintf('\n note: stimLog filename = %s --- eeg.eeglog starttime = %s (addHour=%d)',datestr(fDateStr,'mm/dd/yy HH:MM PM'), dateStr_act(i,:),addHour);

if abs(dateMAT_act-fDateNum) > delta1hour/2 ,
    fprintf('\n uh oh... mismatch of more than 10 min in the date string in the session folder and the posix time in the stimLog.. check this');
    fprintf('\n session: %s    posix date = %s \n', sessName, dateStr_act);
     keyboard;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- clean the output
sitesbyName{size(stimCSV.site1,1)} = [];



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
sitepair1 = unique(stimCSV.site1);

% Get rid of NaN if it's there
sitepair1(isnan(sitepair1)) = '';

% micro id
microID = 16941;

% See if Micro, make all 16 if so
if stimCSV.id(2) == microID
    sitepair2 = zeros(length(sitepair1),1);
    sitepair2(1:length(sitepair1),1) = 16;
else
    sitepair2 = unique(stimCSV.site2);
    % Get rid of NaN if it's there
    sitepair2(isnan(sitepair2)) = '';
end


for n = 1:length(sitepair1)
    
    sitepairs{n} = cellstr(sprintf('%d-%d',sitepair1(n),sitepair2(n)));

    % Get name from montage if Macro
    if stimCSV.id(2) == microID
        physname{1} = cellstr(sprintf('utah%d',sitepair1(n)));
        physname{2} = cellstr(sprintf('utah%d',sitepair2(n)));
    else
        physname{1} = phys2name(subject,rootEEGdir,eegFldr,sitepair1(n));
        physname{2} = phys2name(subject,rootEEGdir,eegFldr,sitepair2(n));
    end
    
    stimNames = [physname{1} '-' physname{2}];
    
    if stimCSV.id(2) == microID
        sitepairs{n} = horzcat(stimNames{:});
    else
        sitepairs{n} = stimNames;
    end
    
    if USE_MONTAGE
        [row1,col1]=ind2sub(size(num),str2num(sitenums{1}));% find specific electrode in montage
        [row2,col2]=ind2sub(size(num),str2num(sitenums{2}));
    
        if ~strcmp(text(row1,col1),physname{1}) | ~strcmp(text(row2,col2),physname{2})
        warning('\nPhys2name and montage do not match up! \nAlways trust the montage and look up the number on physical channel list and montage by eye to see if jacksheet, etc is funky or there is a typo! \nIf you can, fix it ! Call Thia if still confused!');
        stimNames = [text(row1,col1) '-' text(row2,col2)];
        end        
    end
    [sitesbyName{find(stimCSV.site1==sitepair1(n))}] = deal(sitepairs{n});
end
sitesbyName = sitesbyName(1:size(stimCSV.site1,1));


%% assume all NK uses channels 5 and 6
if USE_MONTAGE
    %pretty sure always reference 5 and 6. Might need to find way to adjust for what happens where are switched...
    [rrow1,rcol1]=ind2sub(size(num),5);% find specific electrode in montage for references
    [rrow2,rcol2]=ind2sub(size(num),6);
    NKnames=strcat([text(rrow1,rcol1) '-' text(rrow2,rcol2)]);% pretty sure always use if it's not said
else
    NKnames=strcat([phys2name(subject,rootEEGdir,eegFldr,5) '-' phys2name(subject,rootEEGdir,eegFldr,6)]);% pretty sure always use if it's not said
end

% make events structure
idx= 1;evcnt = 1; tot_pulses = 0;
events = struct([]); 
while idx <= size(stimCSV.offset,1)
    
    if stimCSV.offset(idx) == 0
        idx = idx + 1;
        continue
    end
    
    events(evcnt).subject    = subject;
    events(evcnt).exptFolder = sessName;
    events(evcnt).amplitude  = stimCSV.amp1(idx);
    events(evcnt).frequency  = stimCSV.freq(idx);
    events(evcnt).stimType   = stimCSV.type(idx);
    events(evcnt).stimLocTag = sitesbyName{evcnt};
    events(evcnt).annotation = stimCSV.annotation(idx);
    events(evcnt).DCpulses   = 0;   % initialize as no pulses
    events(evcnt).cerestimID = stimCSV.id(idx);
    events(evcnt).NKreference = NKnames;
    events(evcnt).pulses = stimCSV.pulses(idx); %- JW converted this to num 8/2017
    events(evcnt).offset = stimCSV.offset(idx) + 3600000*addHour; %- JW converted this to num 8/2017
    events(evcnt).mstime = stimCSV.offset(idx) + 3600000*addHour;  %- JW converted this to num 8/2017 Was curTimePos{idx} previously - DY
    events(evcnt).mstime2 = stimCSV.offsetAfterPlay(idx) + 3600000*addHour;
    events(evcnt).mstime3 = stimCSV.offsetAfterDone(idx) + 3600000*addHour;
    
    %- initialize these so all version of task (white noise, mapping, CCEP) have same event structure
    events(evcnt).tele        = [];
    events(evcnt).cumulative_pulses  = [];
    events(evcnt).NKreference = [];
    events(evcnt).monopolar   = [];
    events(evcnt).isWN_REP    = [];
    events(evcnt).WN_REP_obj  = [];
    events(evcnt).elec_loc    = {};
    events(evcnt).MSeqElecList = {};
    events(evcnt).MSeqElecBin = [];
    
    %% getting the dc pulse count
    if strncmp(stimCSV.type(idx),'SM',2)
        if GUI_ver > 5          % correction for extra DC pulses due to 1s max zap protection
            events(evcnt).DCpulses = floor(stimCSV.pulses(idx)/events(evcnt).frequency);
        else
            events(evcnt).DCpulses = 1* dbl_DC;
        end
        
    elseif contains(stimCSV.type(idx),'MSeq') || contains(stimCSV.type(idx),'ClinicalGroupStim')
        % One event should coorrespond to a single time point of stim
        % List all elecs that were stimulated at once
        
        % Redefine it to MSeq since the list of elecs stimulated will be
        % defined elseware
        events(evcnt).stimLocTag = 'MSeqList';
        
        % Record all stim pairs
        events(evcnt).MSeqElecList = sitepairs;
        
        % Figure out how many electrodes were stimulated at once
        numAtThisTime = ismember(stimCSV.offset,stimCSV.offset(idx));
        
        % Number of electrodes stimulated at that offset
        numElecs = length(find(numAtThisTime == 1));
        
        % Map to sites by name to get the electrodes that were on
        elecsOn = sitesbyName(numAtThisTime);
        
        % Cycle through the electrode list to check binarize which were on
        % and off
        for iElec = 1:length(events(evcnt).MSeqElecList)
            events(evcnt).MSeqElecBin(iElec) = max(contains(elecsOn,events(evcnt).MSeqElecList(iElec)));
        end
        
        % Increment by the number of electrodes to proceed to the next time
        % point
        idx = idx + numElecs - 1;

    elseif strncmp(stimCSV.type(idx),'WN_REP',6)
        %events(idx).stimType = 'WN';
        events(evcnt).DCpulses = 1;
        WN_obj = struct;
        if length(stimCSV.type(idx))>6 %- annotation mode
            pieces = strsplit(stimCSV.type(idx),'_');
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
        
    elseif strncmp(stimCSV.type(idx),'WN',2)
        % events.type = 'WN';
        events(evcnt).DCpulses = 1;
        if stimCSV.amp1(idx) < 0
            events(evcnt).monopolar = 0; % 0 is negative 1 is positie
        else
            events(evcnt).monopolar = 1;
        end
        
    elseif any(strcmp(stimCSV.type(idx),annotation_cell))
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
        
    else
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
    end
    
    tot_pulses = tot_pulses+events(evcnt).DCpulses;
    events(evcnt).cumulative_pulses = tot_pulses;

    evcnt = evcnt+1;
    idx   = idx+1;
end

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
