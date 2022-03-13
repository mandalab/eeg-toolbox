function [events] = extractStimEventsFromAll(subject, sessName, behDir, eegDir)
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
% 9/2019 added mSeqGUI changes DY
% 2/12/2020 DY & SJ- changes to incorporate TRE and older subjects with different stimlog columns
%
%
% TO-DO: Need to figure out the stimlog.helper tool for older patients + TRE
%       *Currently, the stimlog has a non-posixtime 'mstime' ('offset' in this code), which will lead to
%       an error stating that the time exceeds the calendar. The .helper should convert these times into
%       usable times, but it doesn't currently run that way for older subjects - need to figure out a way
%       to identify which subjects need it?
%


%typical annotations
annotation_cell = {'Counting','Reading','Naming','Pointing','Other','AD','Seizure'};

% Number of Header Lines
numHeaderLines = 6;

%%- get the session directory
sessDir  =  [behDir '/' sessName];      %- JW tweak --- bad practice to use CD
subjNum = str2num(subject(4:6));
sessNum     = str2num(char(extractBetween(sessName,'session_','_')));


% Check if it has a CSV
if any(size(dir([sessDir '/*.csv' ]),1))
    EXTRACT_FROM_CSV = 1;
else
    EXTRACT_FROM_CSV = 0;
end

% Check if it's mSeqGUI
if contains(behDir,'mSeqGUI')
    EXTRACT_FROM_MSEQGUI_CSV = 1;
else
    EXTRACT_FROM_MSEQGUI_CSV = 0;
end

% Extract from CSV if it's there
if EXTRACT_FROM_CSV
    
    fNameExt = '.csv';
    fName    = dir([sessDir '/*csv']);
    if length(fName)>1,
        iKeep = find([fName.isdir]==0);
        if length(iKeep)~=1,
            fprintf('\n Should only be 1 CSV in stim behavioral folder.  Move original csvs to subfolder "original_CSV" and rerun extraction.');
            keyboard;
            error('\n fix it');
        end
        fName = fName(iKeep);
    end
    fName.name = lower(fName.name);
    stimFile = [sessDir '/' fName.name];
    
    % Load CSV
    fileCSV = fullfile(sprintf('%s/%s',sessDir,fName.name));
    stimCSV = readtable(fileCSV);
    
    %- early version of csv didn't include subsession, so add the column and set to 1's
    if ~strcmp('subSession',stimCSV.Properties.VariableNames)
        for iRow = 1:size(stimCSV,1)
            if stimCSV.offset(iRow) == 0
                stimCSV.subSession(iRow) = 0;
            else
                stimCSV.subSession(iRow) = 1;
            end
        end
    end
    
    %- convert empty annotation column gets read as NaNs... convert to empty strings
    if ~iscell(stimCSV.annotation)
        if sum(isnan(stimCSV.annotation))==size(stimCSV,1),
            annotation = {};
            for iRow = 1:size(stimCSV,1), annotation{iRow} = ''; end;
            stimCSV.annotation = annotation';
        end
    end
    
    %- sanity check, if there are only 2 "offset" values, then this file was probably opened and resaved in xls, which fucks everything
    if length(unique(stimCSV.offset))<=2,
        fprintf('\n SEVERE WARNING: looks like all the offset times are the same.  \n  Was this edited and resaved in xls?  \n  That fucks everything... go back to the original and edit with text.app \n');
        keyboard;
        error('Ask John if this doesnt make sense');
    end
    
    %- for NIH076 JW converted site1 and site2 to string... seems like maybe that was not a good call. correct here
    if iscell(stimCSV.site1), %- JW changed table output so its a string of numbers... seemed like that was the intent because 'N/A' was sometimes in the list, but causes problems beloew...
        siteList = [stimCSV.site1; stimCSV.site2]
        for iS=1:length(siteList),
            thisSite = siteList{iS};
            if isempty(thisSite),
                thisSiteNum = nan;
            elseif strcmp(thisSite,'N/A'),
                thisSiteNum = nan;
            else
                thisSiteNum = str2num(thisSite);
            end
            allsitesNum(iS) = thisSiteNum;
        end
        %- convert string to number
        stimCSV.site1 = [allsitesNum(1:length(stimCSV.site1))]';
        stimCSV.site2 = [allsitesNum(length(stimCSV.site1)+1:end)]';
    end
    
    
elseif ~EXTRACT_FROM_CSV
    % NO CSV, but this is subject 62 and above, so it's either one of the first few sessions when things were still rough (subj 62, sess<18)
    %   or a newer version that can be easily converted to a CSV on the fly
    
    fNameExt = '.stimlog';
    fName    = dir([sessDir '/*stimlog']);
    if isempty(fName)
        fName = dir([sessDir '/*stimLog']);
    end
    fName.name = lower(fName.name);
    stimFile = [sessDir '/' fName.name];
    
    
    %- old version of doing things before Dave made m-sequences and migrated all other tasks to the same format
    if subjNum == 62 && sessNum <= 18 && contains(subject,'NIH') % DY & SJ change
        
        [offset sites amp1 amp2 freq pulses id type annotation]=textread(stimFile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t','headerlines',numHeaderLines);  %- get the file contents
        
        for iSites = 1:length(sites)
            if strcmpi(sites(iSites),'-999-999')
                site1(iSites) = 0;
                site2(iSites) = 0;
                offsetVec(iSites) = 0;
                amp1Vec(iSites)  = 0;
                amp2Vec(iSites)  = 0;
                freqVec(iSites)  = 0;
                pulsesVec(iSites)  = 0;
                idVec(iSites)  = 0;
                subSession(iSites) = 0;
            else
                splitSites = strsplit(char(sites(iSites)),'-');
                site1(iSites) = str2num(splitSites{1});
                site2(iSites) = str2num(splitSites{2});
                offsetVec(iSites) = str2num(offset{iSites});
                amp1Vec(iSites)  = str2num(amp1{iSites});
                amp2Vec(iSites)  = str2num(amp2{iSites});
                freqVec(iSites)  = str2num(freq{iSites});
                pulsesVec(iSites)  = str2num(pulses{iSites});
                idVec(iSites)  = str2num(id{iSites});
                subSession(iSites) = 1;
            end
        end
        
        offsetAfterPlay = zeros(1,length(offset));
        offsetAfterDone = zeros(1,length(offset));
        if length(annotation) ~= length(offset)
            for iCell = 1:(length(offset) - length(annotation))
                annotation{end+1} = '';
            end
        end
        %- weird character getting in here
        for iCell = 1:length(annotation)
            if strcmp(annotation{iCell},'?') |  double(annotation{iCell})==13
                annotation{iCell} = '';
            end
        end
    
    % If old (before DY started in 7/2019)    
    elseif (subjNum < 62 && contains(subject,'NIH')) || (subjNum < 17 && contains(subject,'TRE')) % DY & SJ change (added to take out idVec and offset)
        
        [offset sites amp1 amp2 freq pulses type annotation]=textread(stimFile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t','headerlines',numHeaderLines);  %- get the file contents
        
        for iSites = 1:length(sites)
            if strcmpi(sites(iSites),'-999-999')
                site1(iSites) = 0;
                site2(iSites) = 0;
                offsetVec(iSites) = 0;
                amp1Vec(iSites)  = 0;
                amp2Vec(iSites)  = 0;
                freqVec(iSites)  = 0;
                pulsesVec(iSites)  = 0;
                subSession(iSites) = 0;
            else
                splitSites = strsplit(char(sites(iSites)),'-');
                site1(iSites) = str2num(splitSites{1});
                site2(iSites) = str2num(splitSites{2});
                offsetVec(iSites) = str2num(offset{iSites});
                amp1Vec(iSites)  = str2num(amp1{iSites});
                amp2Vec(iSites)  = str2num(amp2{iSites});
                freqVec(iSites)  = str2num(freq{iSites});
                pulsesVec(iSites)  = str2num(pulses{iSites});
                subSession(iSites) = 1;
            end
        end
        idVec = ones(1,length(offset)) * 65289; % Assume it was Macro 1 (DY & SJ change)
        offsetAfterPlay = zeros(1,length(offset));
        offsetAfterDone = zeros(1,length(offset));
        if length(annotation) ~= length(offset)
            for iCell = 1:(length(offset) - length(annotation))
                annotation{end+1} = '';
            end
        end
        %- weird character getting in here
        for iCell = 1:length(annotation)
            if strcmp(annotation{iCell},'?') |  double(annotation{iCell})==13
                annotation{iCell} = '';
            end
        end   
        
    else  % if subjNum == 62 & sessNum <= 18
        
        %%- Newest version of session log, but missing a CSV for some reason so make it on the fly
        [offset sites junk amp1 amp2 freq pulses id junk2 junk3 type annotation offset2 offset3 junk4]=textread(stimFile,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','headerlines',numHeaderLines,'delimiter','\t','whitespace','');  %- get the file contents
        
        %- could be
        if isempty(offset2{1}) && isempty(offset2{end}) && isempty(offset3{1}) && isempty(offset3{end})
            %fprintf('\n text read not matching file data... try again');
            [offset sites junk amp1 amp2 freq pulses id type junk2 offset2 offset3  ]=textread(stimFile,'%s %s %s %s %s %s %s %s %s %s %s %s','headerlines',numHeaderLines,'delimiter','\t','whitespace','');  %- get the file contents
            %[offset sites amp1 amp2 freq pulses id type offset2 offset3  ]=textread(stimFile,'%s%s%s%s%s%s%s%s%s%s','headerlines',numHeaderLines);  %- JW thinks automatic detection works better... until it doesnt. dont use this
            annotation = {};
            for iCell = 1:length(offset)
                annotation{end+1} = '';
            end
            annotation = annotation';
        end
        type = strtrim(type); %- get rid of the white spaces
        
        for iSites = 1:length(sites)
            if strcmpi(sites{iSites},'') || strcmpi(sites{iSites},'  ')
                site1(iSites)      = 0;
                site2(iSites)      = 0;
                offsetVec(iSites)  = 0;
                offsetAfterPlay(iSites) = 0;
                offsetAfterDone(iSites) = 0;
                amp1Vec(iSites)    = 0;
                amp2Vec(iSites)    = 0;
                freqVec(iSites)    = 0;
                pulsesVec(iSites)  = 0;
                idVec(iSites)      = 0;
                subSession(iSites) = 0;
                annotation{iSites} = junk2{iSites};
            else
                splitSites         = strsplit(char(sites(iSites)),'-');
                site1(iSites)      = str2num(splitSites{1});
                site2(iSites)      = str2num(splitSites{2});
                offsetVec(iSites)  = str2num(offset{iSites});
                offsetAfterPlay(iSites) = str2num(offset2{iSites});
                offsetAfterDone(iSites) = str2num(offset3{iSites});
                amp1Vec(iSites)    = str2num(amp1{iSites});
                amp2Vec(iSites)    = str2num(amp2{iSites});
                freqVec(iSites)    = str2num(freq{iSites});
                pulsesVec(iSites)  = str2num(pulses{iSites});
                idVec(iSites)      = str2num(id{iSites});
                subSession(iSites) = 1;
            end
        end
        
    end  %-  if subjNum == 62 & sessNum <= 18
    
    
    %- convert the extracted data from the session log into a table
    stimCSV = table(offsetVec',site1',site2',amp1Vec',amp2Vec',freqVec',pulsesVec',idVec',type,annotation,offsetAfterPlay',offsetAfterDone',subSession',...
        'VariableNames',{'offset' 'site1' 'site2' 'amp1' 'amp2' 'freq' 'pulses' 'id' 'type' 'annotation' 'offsetAfterPlay' 'offsetAfterDone' 'subSession'});
    
end


%- create a fake sub-session increment for stim mapping if/when the electrodes change (makes code below work better)
%   this will be done in stimMapGUI for patients after 68, but for now do it here
% Not necessary for mSeqGUI
if ~EXTRACT_FROM_MSEQGUI_CSV
    
    % DY 1/2020 Reincrement stimCSV to make new subSess number for each
    % unique run
    
    % Find subSessions that are zero
    zeroIdx = find(stimCSV.subSession == 0);
    
    % Loop through these to change each subSession number
    for iSubSess = 1:length(zeroIdx)-1
        
        % Change the events in between iSubSess and iSubSess + 1 indexes to
        % iSubSess
        stimCSV.subSession((zeroIdx(iSubSess)+1):(zeroIdx(iSubSess+1)-1)) = iSubSess;
        
        % Make sure the last one's are done too
        if iSubSess == length(zeroIdx)-1
            stimCSV.subSession((zeroIdx(iSubSess+1)+1):end) = iSubSess + 1;
        end
    end

    isStimMap  = any(contains(stimCSV.type,'SM'));
    multiChan  = length(unique([stimCSV.site1; stimCSV.site2]));
    numSubSess = length(unique([stimCSV.subSession]));
    if isStimMap & multiChan>2 & numSubSess==1,
        thisPair    = [stimCSV.site1(1) stimCSV.site2(1)];
        thisSubSess = stimCSV.subSession(1);
        for iRow = 2:size(stimCSV,1),
            if stimCSV.site1(iRow)~=thisPair(1) | stimCSV.site2(iRow)~=thisPair(2),
                thisPair = [stimCSV.site1(iRow) stimCSV.site2(iRow)];
                thisSubSess = thisSubSess+1;
            end
            stimCSV.subSession(iRow) = thisSubSess;
        end
    end
end

%%--------------------------------------------------------------------------------------------------------------------
%%   OLD VERSION OF STIM MAP SOMETIMES GOT DAYLIGHT SAVINGS WRONG
%%     --- following code uses THE DATESTRING FROM THE SESSION FOLDER TO CORRECT FOR DAYLIGHT SAVINGS ERRORS IN POSIX TIME
%%     --- JW thinks if no ".pulselog" then this correction shouldnt be applied because it wont be applied to the eeg.log
%%--------------------------------------------------------------------------------------------------------------------


%-
fNamePulse = dir([sessDir '/*pulselog']);
if ~isempty(fNamePulse), ADJUST_4_DST = 1;
else                     ADJUST_4_DST = 0; end %- code as written doesn't touch eeg.log if there is no pulse log, so fixing DST in session log will cause a mismatch with pulse log

    
    
%- grab the datestring from the session folder
try
    fDateStr = datetime(fName.name(1:strfind(fName.name,fNameExt)-1),'InputFormat','yyyy_MM_dd_HHmm');
catch
    fDateStr = datetime(fName.name(1:strfind(fName.name,fNameExt)-1),'InputFormat','dd-MMM-yyyy_HH_mm_SS');
end
fDateNum = datenum(fDateStr);


%- sanity check that posix time in the file matches the datestring on the session folder
curTimePos  = stimCSV.offset(2); %- numHeaderLines above jumps over the header... so just grab 1st or 2nd row of real data for this step
dateMAT_act = datenum(epoch2date(curTimePos));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
delta1hour  = datenum('01/01/01 2:00')-datenum('01/01/01 1:00');

% We shouldnt need this for mSeqGUI
if ~EXTRACT_FROM_MSEQGUI_CSV
    addHourOpt = [0 1 -1]; %- add and remove hour to correct for daylight savings (JW added -1 on 7/4/2019)
    t          = [dateMAT_act addtodate(dateMAT_act, addHourOpt(2), 'hour') addtodate(dateMAT_act, addHourOpt(3), 'hour')]; %- try with and without additional hour
    tmpdate    = datetime(t,'ConvertFrom','datenum');
    dateMAT_act = datenum(tmpdate);
    dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    timeDiff = dateMAT_act-fDateNum;
    [y    i]    = min(timeDiff(timeDiff > 0)); %- smallest positive difference... but could be VERY small negative difference
    [yAbs iAbs] = min(abs(timeDiff)); %- smallest positive difference... but could be VERY small negative difference
    
    addHour = addHourOpt(i);
    %fprintf('\n note: stimLog filename = %s --- eeg.eeglog starttime = %s (addHour=%d)',datestr(fDateStr,'mm/dd/yy HH:MM PM'), dateStr_act(i,:),addHour);
    
    
    if ADJUST_4_DST==0 & addHour~=0,
        fprintf('\n WARNING: CSV or session log timing mismatch from session folder date string, but no .pulseLog, so cant adjust for DST');
        knownOKsess = {'session_1b_MSeq_190503_1508'};
        if ~any(contains(knownOKsess, sessName)),
            keyboard
        end
        addHour = 0;
        i = 1;
    end

    if abs(dateMAT_act(i)-fDateNum) > delta1hour/2 ,
        timeDiffMin = 60 * abs(dateMAT_act(i)-fDateNum) / delta1hour;
        fprintf('\n uh oh... mismatch of %.1f min (usually <30 min) in the date string in the session folder and the posix time in the stimLog.. check this',timeDiffMin);
        fprintf('\n session: %s    posix date = %s \n', sessName, dateStr_act(i,:));
        %keyboard;
    end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- clean the output
sitesbyName{size(stimCSV.site1,1)} = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rootEEGdir = behDir(1:strfind(behDir,subject)-2);


%- phys2chan requires knowing a raw folder, but at this point task is unaligned, so we dont know that
%    if all raw folders had the same channel configuration we would be all good, but sometimes the system
%    switches from NK to blackrock... or electrodes could be moved around if a rearrangement.
%    so best to make sure the folder is at least close to correct (same day)
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
%- we have a list of stimMap folders, now find the closet match using datenum
stimDateList = datenum(datetime(stim_map_folders,'InputFormat','yyMMdd_HHmm'));
deltaTime    = fDateNum-stimDateList;
[minDelta, iClosest] = min(abs(deltaTime));
if abs(minDelta) > delta1hour*2 ,
    fprintf('\n uh oh... mismatch of more than 2 hour in the date string in the session folder and the closest STIM_MAP raw folder found');
    fprintf('\n this could mess up phys2chan if the recording channels changed (i.e., channel names changed or channels deactivated');
    fprintf('\n session: %s    closest stimMap = %s \n', sessName, stim_map_folders{iClosest});
    keyboard;
end
eegFldr = stim_map_folders{iClosest};  %- in theory this should be OK... even if not correct mapping recording device shouldn't change in short interval



% micro id
microID = 16941;

% Get number of subsessions
numSubSess = max(stimCSV.subSession);



% first grab all the phys2chan results so we are not looping over the same
% channels over and over again (phys2chan is a slow-ish step). do
% differently for mSeqGUI
if ~EXTRACT_FROM_MSEQGUI_CSV
    allsites = unique([stimCSV.site1 stimCSV.site2]);
    allsites = allsites(allsites>0 & allsites<=512 & ~isnan(allsites)); %- cut zero, >512, and nan
    phys2chanLookup = {};
    for iSites=1:length(allsites),
        phys2chanLookup{allsites(iSites)} = phys2name(subject,rootEEGdir,eegFldr,allsites(iSites));
    end
    
    % Loop through each subSession
    for iSub = 1:numSubSess
        
        % reset sitepair1 and sitepair2
        clear sitepair1 sitepair2
        
        % get all the site pair names
        sitepair1 = unique(stimCSV.site1(find(stimCSV.subSession == iSub)));
        
        % Get rid of NaN and zero if it's there
        sitepair1(isnan(sitepair1)) = '';
        sitepair1(~any(sitepair1,2)) = '';
        sitepair2 = []; %- initialize here in case first subSess has no pairs
        
        % Make a temporary sub table for next step
        tempTable = stimCSV(find(stimCSV.subSession == iSub),:);
        
        % See if Micro, make all cathodes=16 to represent the common ground
        if stimCSV.id(find(stimCSV.subSession == iSub)) == microID
            sitepair2 = zeros(length(sitepair1),1);
            sitepair2(1:length(sitepair1),1) = 16;
        else
            for iSite = 1:length(sitepair1)
                sitepair2(iSite) = unique(tempTable.site2(find(tempTable.site1 == sitepair1(iSite))));
            end
            
            % Get rid of NaN if it's there
            sitepair2(isnan(sitepair2)) = '';
            sitepair1(~any(sitepair2,2)) = '';
            
            if length(find(stimCSV.id==microID & stimCSV.subSession == iSub))>0,
                fprintf('\n probably should be counted as micro, but wasnt... look into it');
                keyboard;
            end
        end
        
        
        for n = 1:length(sitepair1)
            
            sitepairs{iSub}(n) = cellstr(sprintf('%d-%d',sitepair1(n),sitepair2(n)));
            
            %- make this into a 2xN numerical array... doesnt make sense to turn it into text
            physChans = [string(sitepair1(n)) '-' string(sitepair2(n))];
            %sitepairsPhys{iSub}(n) = cellstr(horzcat(physChans{:}));
            %sitepairsPhys{iSub}(1:2,n) = [sitepair1(n) sitepair2(n)];
            sitepairsPhys{iSub}{n} = [sitepair1(n) sitepair2(n)];
            
            % Get name of channels using phys2name if marco
            if stimCSV.id(find(stimCSV.subSession == iSub)) == microID
                physname{1} = cellstr(sprintf('utahS%d',sitepair1(n)));
                physname{2} = cellstr(sprintf('utahGND',sitepair2(n)));
            else
                %- old way... multiple calls to phys2name for same electrode (slow step)
                %physname{1} = phys2name(subject,rootEEGdir,eegFldr,sitepair1(n));
                %physname{2} = phys2name(subject,rootEEGdir,eegFldr,sitepair2(n));
                %- new way... make a look up table above
                physname{1} = phys2chanLookup{sitepair1(n)};
                physname{2} = phys2chanLookup{sitepair2(n)};
            end
            
            stimNames = [physname{1} '-' physname{2}];
            
            if stimCSV.id(find(stimCSV.subSession == iSub)) == microID
                sitepairs{iSub}(n) = cellstr(horzcat(stimNames{:}));
            else
                sitepairs{iSub}(n) = {stimNames};
            end
            
            [sitesbyName{find(stimCSV.site1==sitepair1(n) & stimCSV.subSession == iSub)}] = deal(sitepairs{iSub}(n));
        end
        
    end
    sitesbyName = sitesbyName(1:size(stimCSV.site1,1));
else
    allsites = unique([stimCSV.site1 stimCSV.site2]);
    allsites = allsites(allsites>0 & allsites<=1115 & ~isnan(allsites)); %- cut zero, >512, and nan
    phys2chanLookup = {};
    for iSites=1:length(allsites),
        if allsites(iSites) < 1000 % Implies macro
            phys2chanLookup{allsites(iSites)} = phys2name(subject,rootEEGdir,eegFldr,allsites(iSites));
        elseif allsites(iSites) < 1100 % Implies utah anode
            phys2chanLookup{allsites(iSites)} = sprintf('utahS%d',allsites(iSites)-1000);
        else
            phys2chanLookup{allsites(iSites)} = 'utahGND';
        end
    end
    
       % Loop through each subSession
    for iSub = 1:numSubSess
        
        % reset sitepair1 and sitepair2
        clear sitepair1 sitepair2
        
        % get all the site pair names
        sitepair1 = unique(stimCSV.site1(find(stimCSV.subSession == iSub)));
        
        % Get rid of NaN and zero if it's there
        sitepair1(isnan(sitepair1)) = '';
        sitepair1(~any(sitepair1,2)) = '';
        sitepair2 = []; %- initialize here in case first subSess has no pairs
        
        % Make a temporary sub table for next step
        tempTable = stimCSV(find(stimCSV.subSession == iSub),:);
        
        % Make sitepair2
        for iSite = 1:length(sitepair1)
            sitepair2(iSite) = unique(tempTable.site2(find(tempTable.site1 == sitepair1(iSite))));
        end
        
        % Get rid of NaN if it's there
        sitepair2(isnan(sitepair2)) = '';
        sitepair1(~any(sitepair2,2)) = '';

        for n = 1:length(sitepair1)
            
            sitepairs{iSub}(n) = cellstr(sprintf('%d-%d',sitepair1(n),sitepair2(n)));
            
            %- make this into a 2xN numerical array... doesnt make sense to turn it into text
            physChans = [string(sitepair1(n)) '-' string(sitepair2(n))];
            sitepairsPhys{iSub}{n} = [sitepair1(n) sitepair2(n)];
            
            % Get name of channels using phys2name if marco
            physname{1} = phys2chanLookup{sitepair1(n)};
            physname{2} = phys2chanLookup{sitepair2(n)};
            
            stimNames = [physname{1} '-' physname{2}];
            sitepairs{iSub}(n) = {stimNames};
            
            [sitesbyName{find(stimCSV.site1==sitepair1(n) & stimCSV.subSession == iSub)}] = deal(sitepairs{iSub}(n));
        end
        
    end
    sitesbyName = sitesbyName(1:size(stimCSV.site1,1));
    
end

%% assume NK always uses channels 5 and 6
%pretty sure always reference 5 and 6. Might need to find way to adjust for what happens where are switched...
NKnames=strcat([phys2name(subject,rootEEGdir,eegFldr,5) '-' phys2name(subject,rootEEGdir,eegFldr,6)]);% pretty sure always use if it's not said



% make events structure
idx= 1;evcnt = 1; tot_pulses = 0;
events = struct([]);

% prepare subSess tracker
subSess = 0;

while idx <= length(stimCSV.offset),
    
    
    if stimCSV.offset(idx) == 0 | stimCSV.offset(idx) == -999
        idx = idx + 1;
        continue
    end
    
    events(evcnt).subject      = subject;
    events(evcnt).exptFolder   = sessName;
    events(evcnt).subSession   = stimCSV.subSession(idx);
    events(evcnt).stimType     = stimCSV.type{idx};  %- want a character array
    %events(evcnt).stimType    = string(stimCSV.type(idx));
    events(evcnt).cerestimID   = stimCSV.id(idx);
    
    
    %- the the name of the stimulation electrodes.  If this is mSeq it will be overwritten further down
    stimLocTag = sitesbyName{idx};
    if iscell(stimLocTag),
        stimLocTag = stimLocTag{1};
    else
        stimLocTag = 'None';
    end
    
    events(evcnt).stimLocTag   = stimLocTag; %- cell array or string?
    events(evcnt).numStimLoc   = 1;                  %- assume one stimulation output (whether monopolar or biopolar)
    events(evcnt).amplitude    = stimCSV.amp1(idx);
    events(evcnt).frequency    = stimCSV.freq(idx);
    events(evcnt).pulses       = stimCSV.pulses(idx); %- JW converted this to num 8/2017
    if ~EXTRACT_FROM_MSEQGUI_CSV
        events(evcnt).mstime       = stimCSV.offset(idx) + 3600000*addHour; %- JW converted this to num 8/2017
        events(evcnt).mstimeBeforePlay = stimCSV.offset(idx) + 3600000*addHour;  %- JW converted this to num 8/2017 Was curTimePos{idx} previously - DY
        events(evcnt).mstimeAfterPlay  = stimCSV.offsetAfterPlay(idx) + 3600000*addHour;
        events(evcnt).mstimeAfterDone  = stimCSV.offsetAfterDone(idx) + 3600000*addHour;
    else
        events(evcnt).mstime       = stimCSV.offset(idx); %- JW converted this to num 8/2017
        events(evcnt).group        = stimCSV.group(idx);
    end
       
    
    %- initialize these so all version of task (white noise, mapping, CCEP) have same event structure
    %events(evcnt).DCpulses     = 0;   % initialize as no pulses... used to create eeg.eeglog10, but is that actually used for anything?
    %events(evcnt).cumulative_pulses = []; % counter of all DC10 pulses expected for this session
    %events(evcnt).tele        = [];
    %events(evcnt).cumulative_pulses  = [];
    events(evcnt).annotation   = stimCSV.annotation{idx};
    events(evcnt).NKreference  = NKnames;
    %events(evcnt).monopolar   = [];
    %events(evcnt).isWN_REP    = [];
    %events(evcnt).WN_REP_obj  = [];
    %events(evcnt).elec_loc    = {};
    events(evcnt).MSeqElecNumList = [];
    events(evcnt).MSeqElecList    = {};
    events(evcnt).MSeqElecBin     = [];
    events(evcnt).group           = NaN;

    DCpulses = 1; %- assume 1 pulse per event... will be adjusted below to zero if annotation, and more than one if prolonged stimMap
    
    
    %% getting the dc pulse count
    if strncmp(stimCSV.type(idx),'SM',2)
        % correction for extra DC pulses due to 1s max zap protection
        %events(evcnt).DCpulses = floor(stimCSV.pulses(idx)/events(evcnt).frequency);
        
        
    elseif contains(stimCSV.type(idx),'MSeq') || contains(stimCSV.type(idx),'ClinicalGroupStim') || contains(stimCSV.type(idx),'PulsesMaxOnly') || contains(stimCSV.type(idx),'Simple') || contains(stimCSV.type(idx),'GroupClear')
        % One event should coorrespond to a single time point of stim
        % List all elecs that were stimulated at once
        
        % Redefine it to MSeq since the list of elecs stimulated will be defined elseware
        events(evcnt).stimLocTag = 'MSeqList';
        
        % Record all stim pairs
        events(evcnt).MSeqElecNumList = sitepairsPhys{stimCSV.subSession(idx)};
        events(evcnt).MSeqElecList    = sitepairs{stimCSV.subSession(idx)};
        
        % Figure out how many electrodes were stimulated at once
        numAtThisTime = ismember(stimCSV.offset,stimCSV.offset(idx));
        
        % Number of electrodes stimulated at that offset
        numElecs = length(find(numAtThisTime == 1));
        events(evcnt).numStimLoc =numElecs;
        
        % Map to sites by name to get the electrodes that were on
        elecsOn = sitesbyName(numAtThisTime);
        
        % Cycle through the electrode list to check binarize which were on and off
        for iElec = 1:length(events(evcnt).MSeqElecList)
            events(evcnt).MSeqElecBin(iElec) = max(contains([elecsOn{:}],events(evcnt).MSeqElecList(iElec)));
        end
        
        % If subsess is same, keep the amplitude and cerestim vector the same
        if subSess ~= stimCSV.subSession(idx)
            
            % redefine subsess
            subSess = stimCSV.subSession(idx);
            
            % Make amplitude, group, and cerestimId a vector of the same size as MSeqElecBin
            events(evcnt).amplitude = zeros(size(events(evcnt).MSeqElecBin));
            events(evcnt).cerestimID = zeros(size(events(evcnt).MSeqElecBin));
            events(evcnt).group = zeros(size(events(evcnt).MSeqElecBin));
            
            % First grab rows of table we want
            tempTable = stimCSV(stimCSV.subSession == subSess,:);
            
            % Loop through to fill amplitude
            for iElec = 1:length(events(evcnt).MSeqElecList)
                numHolder = events(evcnt).MSeqElecNumList(iElec);
                elecIdx = find(tempTable.site1 == numHolder{1}(1));
                events(evcnt).amplitude(iElec) = tempTable.amp1(elecIdx(1));
                events(evcnt).cerestimID(iElec) = tempTable.id(elecIdx(1));
                if ismember('group',tempTable.Properties.VariableNames)
                    events(evcnt).group(iElec) = tempTable.group(elecIdx(1));
                else
                    events(evcnt).group(iElec) = NaN;
                end
            end
            
        else
            
            % Grab previous vectors
            events(evcnt).amplitude = events(evcnt-1).amplitude;
            events(evcnt).cerestimID = events(evcnt-1).cerestimID;
            events(evcnt).group = events(evcnt-1).group;
        end

        % Increment by the number of electrodes to proceed to the next time point
        idx = idx + numElecs - 1;
        
        
    elseif contains(stimCSV.type(idx),'Baseline')
        DCpulses = 0;
        events(evcnt).numStimLoc = 0;
        events(evcnt).annotation = stimCSV.annotation(idx);
        
        
    elseif any(strcmp(stimCSV.type(idx),annotation_cell))
        DCpulses = 0;
        events(evcnt).numStimLoc = 0;
        events(evcnt).stimType   = 'ANN';
        events(evcnt).annotation = stimCSV.annotation(idx);
        
        
    else
        DCpulses = 0;
        events(evcnt).numStimLoc = 0;
        events(evcnt).stimType   = 'ANN';
        events(evcnt).annotation = stimCSV.annotation(idx);
    end
    
    % events(evcnt).DCpulses = DCpulses;
    %tot_pulses = tot_pulses+DCpulses;
    %events(evcnt).cumulative_pulses = tot_pulses; %
    
    evcnt = evcnt+1;
    idx   = idx+1;
end


%- debugging step... show what the events are looking like
if 0,
    fprintf('\n');
    disp(events(end))
    %events(end).MSeqElecNumList
    disp(events(end).MSeqElecList)
end




%%% CREATE the EEG.EEGLOG if not already made

%%- 2 different versions of how an eeg.eeglog gets created
%     1) [NIH062-67] pulselog has DC09 pulses in matlab time, which needs to be converted to posix in a file called 'eeg.eeglog'
%     2) [NIH-68+]   eeg.eeglog is directly created... no pulseLog in the folder, no conversion to posix required

fNamePulse = dir([sessDir '/*pulselog']);
if length(fNamePulse)>0,
    
    %-----------------------------------------------------------------------------------------------------%
    %% First create a helper version of the pulselog to help with alignment issues
    %- now make a "helper"  pulseLog file that is a copy of the originals, but with dateStrings and posix time strings for each row
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
    
    
    
    %-----------------------------------------------------------------------------------------------------%
    %% Second convert from matlab to posix time as save as "eegDC09.eeglog"
    
    %- create an eeg.log with python/unix time for prepAndAlign
    eegLogFile = [sessDir '/eegDC09.eeglog'];
    fileId =fopen(eegLogFile,'w');      % for now just replacing the eeg.eeglog.up file
    
    %- open  the "pulse" file, convert to python/unit style times, and save out as eeg.eeglogDC09
    pulseFile  = [sessDir '/' fNamePulse.name];
    [timePCstr pulseStr]=textread(pulseFile,'%s\t %s', 'headerlines',1);  %- get the file contents
    
    %- check the first pulse... does it match up with the session folder date string?
    curMatTime = str2num(timePCstr{1})/1e11;
    t          = addtodate(curMatTime, 4+addHour, 'hour');
    tmpdate    = datetime(t,'ConvertFrom','datenum');
    curTimePos = round(posixtime(tmpdate)*1000);
    
    %- sanity check that posix time in the file matches the datestring on the session folder
    %curTimePos = str2num(offset{1}); %- numHeaderLines above jumps over the header... so just grab 1st or 2nd row of real data for this step
    dateMAT_act = datenum(epoch2date(curTimePos));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
    dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    if abs(dateMAT_act-fDateNum) > delta1hour,
        fprintf('\n uh oh... mismatch of more than 10 min in the date string in the session folder and the posix time in the stimLog.. check this');
        fprintf('\n session: %s    posix date = %s \n', sessName, dateStr_act);
        keyboard;
    end
    
    %- now loop through all the "up" times and copy them to eeg.log with corrected mstime
    iUp = find(strcmp(pulseStr,'PULSE_HI'))';
    for ii=iUp,
        curMatTime = str2num(timePCstr{ii})/1e11;
        t          = addtodate(curMatTime, 4+addHour, 'hour');
        tmpdate    = datetime(t,'ConvertFrom','datenum');
        curTimePos = num2str(round(posixtime(tmpdate)*1000));
        
        fprintf(fileId,'%s\t%d\t%s\n',curTimePos,0,'CHANNEL_0_UP');
    end
    fclose(fileId);
    
    
    
    %-----------------------------------------------------------------------------------------------------%
    %% 3rd step: delete any existing eeg.eeglog... this should be created fresh from the pulseLog -> eegDC10log that was just made
    fNameEEG = dir([sessDir '/eeg.eeglog*']);
    for iF=1:length(fNameEEG),
        thisFile = fullfileEEG(sessDir, fNameEEG(iF).name);
        delete(thisFile);
    end
    
    
    
    %-----------------------------------------------------------------------------------------------------%
    %% Last step: make eeg.eeglog a copy of eegDC09.eeglog
    sourceStr = 'eegDC09.eeglog';
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile([sessDir '/' sourceStr],[sessDir '/eeg.eeglog'],'f');
    if SUCCESS, %fprintf('\n  %s -->  eeg.eeglog is a copy of %s',sessDir,sourceStr);
    else        fprintf('\n   Uh oh'); keyboard; end
    
end




%%- make a "fake" session log with entries matching stimlog so if there is an alignment issue we can figure out which line is the culprit
%%- always remake the session.log... this can be used to debug alignment errors with events.mat
sessLog = [sessDir '/session.log'];
fileId = fopen(sessLog,'w'); %- always make a fresh copy of this
for iEv = 1:length(events),
    fprintf(fileId,'%d\t%d\t%s\n',events(iEv).mstime ,0, events(iEv).stimLocTag);
end
fclose(fileId);

%- little cleanup... if just making a new session.log then the align version is invalid anyway
sessLogAlign = [sessLog '.align'];
if exist(sessLogAlign,'file'),
    delete(sessLogAlign);
end



end  %- end main()
