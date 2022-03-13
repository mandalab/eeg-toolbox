function eegPrepAndAlign(subj, rootEEGdir, varargin)
%Function eegPrepAndAlign( subj, rootEEGdir)
%
%  Master function used to prepare a subject's data for alignment.
%  Rerun multiple times to guide alignment process.
%
%    --Preparation of EEG data--
%       (1) run to generate list of approximate eeg session times to guide Kareem's acquisition of RAW data
%       (2) run to convert single RAW eeg data to multiple channel files
%       (3) run to identify eeg sync channels and confirm existance of resultance eeg sync file
%       (4) preprocess signal: detrend and line noise removal --> eeg.processed folder
%       (5) find bad channels: high variance and persistant line noise --> bad_chans.csv
%       (6) rereference processed signal: all bipolar pairs, 2 average signals (all chans and all without bad_chans)
%
%    --Preparation of Behavioral data--
%       (1) run to extract behavioral events.mat for all tasks and sessions
%       (2) run to align behavioral sync file and eeg sync files
%
%  Outputs:
%    --Command line outputs describing status of RAW and Behavioral files
%    --text file saved to subject/behavior/alignmentSummary.txt
%    --if alignment possible, updated events.mat in each session and "master events.mat" found in task's root directory
%    --                       copies of session.log and eeg.eeglog with align info (session.log.align and eeg.eeglog.align)
%    --text file alignmentPairing_auto.txt that lists alignment pairs... this can be modified and saved as alignmentPairing_forced.txt to override auto pairing
%    --text file alignmentStats.txt that lists alignment fit info (sucess, R^2, max devitions, etc)
%    --csv file listing electrodes with high variance / line-noise
%
%  Inputs:
%    -- subj                        % subject string,        ex) 'NIH016'
%    -- rootEEGdir                  % data path up to "eeg", ex) '/Users/wittigj/DataJW/data/eeg'
%    -- FORCE_EVENT_REEXTRACTION    % 0 or 1: re-extract all behavioral events using "behavioralProcessing";
%                                             (events.mat always extracted if missing)
%    -- FORCE_MASTER_EVENTS_UPDATE  % 0 or 1: confirm all master events.mat up-to-date after alignment and update if need be
%                                             (master always updated when session aligned)
%  Name-Value Pair Inputs
%    -- freshJackandSplit=[true/false]   % False (default) to delete eeg.noreref, eeg.reref, eeg.processed, and jackSheetMaster.xls so they can all be freshly split
%    -- freshExtractAlign=[true/false]   % False (default) to re-extract all events.mats, even if they already exist, and delete docs/{aligmentStats,alignmentSummary,eventFilePaths}.
%    -- batch=[true/false]               % False (default) to prompt user with input. True to try skipping prompts (CAREFUL!)
%    -- disp=[true/false]                % False (default). Use true to have eeg_noise_metrics figures display when they are created. False to hide (they are always saved regardless)
%    -- skipProcessReref=[true/false]    % False (default). Use true to skip inquiry about processing and reref... useful when running BATCH on lots of subjects that will get processed (slow step) later
%    Note: "useStim" is no longer a valid input
%
%  Example call:  eegPrepAndAlign('NIH048', '/Volumes/JW24TB/data24TB/eeg', 'freshJackandSplit',1, 'freshExtractAlign',1, 'batch',0);%

%
%
% Revisions
%  9/10/2013 JHW created
%    07/2016 MST modified to use new file structure (element_info.csv/jacksheetMaster.csv)
%    05/2017 MST modified to create processed folder, rework reref folder, identify bad chans
%    06/2017 MST create processed folder with processAndReref (dungeon toolbox)
%    07/2020 SNJ fixed some end statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function Flow:
%       [1/5]  Searching RAW and/or EEG.NOREREF for Raw and/or Extracted Channel files
%               --> if raws missing, look for extracted channels and sync files and proceed (sometimes extraction present when raws aren't)
%               --> if raws present, confirm that each has been extracted and has a sync file;  extract if necessary;  rereference if necessary
%       [2/5] Searching behavioral directory for eeg.eeglog and events.mat files
%               --> extract events.mat and eeg.eeglogup if missing
%       [3/5] Preparing for alignment: identify pairings of eeg.eeglog <--> extracted sync pulses
%               --> try to identify pairs of behavioral and eeg files.  happens automatically, but can be overrode with text file
%       [4/5] Alignment: confirm all events.mat aligned; offer alignment if not:
%               --> steps through each alignment one-by-one;  update master events.mat.
%       [...] Alignment Summary output
%               --> graphs and text output indicating target dates of missing raw files, missing sync files, and unsuccessful alignments
%       [5/5] ProcessAndReref:
%               --> Use eeg.noreref files to create global/bipolar references in eeg.reref
%               --> Also create processed files and identify "good channels" in eeg.processed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- USER SETS THE FOLLOWING VARIABLES if running directly from M-file (instead of using as function... useful for debugging)

%clear all
%subj       = 'NIH019';
%rootEEGdir = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg';  %%- rootEEGdir is data path up to and including "EEG"
%rootEEGdir = '/Volumes/Kareem/data/eeg/';  %%- rootEEGdir is data path up to and including "EEG"
%rootEEGdir = 'C:/Users/jDub/DataJW/eeg';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- Global Variable Definition and Error check existance of rootEEGdir string
%%--
%- execution options
FORCE_RAW_REEXTRACTION      = 0 ;  % 0 or 1: re-extract all raw data to individual channels (Raw always extracted if matching channel data is missing)

FORCE_EVENT_REEXTRACTION    = 0 ;  % 0 or 1: re-extract all behavioral events using "behavioralProcessing"; if 0 only sessions missing events will be extracted
FORCE_MASTER_EVENTS_UPDATE  = 1 ;  % 0 or 1: confirm all master events.mat is up-to-date after alignment and update if need be (master always updated when session aligned)

FORCE_REALIGNMENT           = 0 ;  % 0 or 1: (re)run alignment on all event files.
BYPASS_ALIGNMENT_USER_QUERY = 0 ;  % 0 or 1: should ALWAYS be 0, unless you are sure that all alignments are OK and just need to realign a task.

CHECK_FILECLOSE = 1;   % a debugging step... trying to find some run away fopens

SKIP_EXTRACTION = 0;   % 0 or 1:  should ALWAYS be 0, unless you are doing a first pass at confirming behavior and physio line up... eeg.noreref and eeg.reref should be deleted after running this in setting 1
BATCH           = 0;   % 0 or 1: if 1, try not to prompt for user input

%- events eegfile entry all pointed to this root
serverDataPath = '/Volumes/Shares/FRNU/data/eeg/';
FORCE_ALIGNMENT_TO_SERVER_PATH = 1; % 0 or 1: if 1 (default), set aligned events eegfile to '/Volumes/Shares/FRNU/data/eeg/subj/eeg.reref/'; if 0, set aligned events to local data directory


% Input parsing
p = inputParser;
p.addParameter('freshJackandSplit', false);
p.addParameter('freshExtractAlign', false);
p.addParameter('disp', 0);
p.addParameter('batch', false);
p.addParameter('skipProcessReref', false);
p.KeepUnmatched = true;
parse(p, varargin{:});
if length(fieldnames(p.Unmatched))>0, fprintf('\n ERROR: your inputs to eegPrepAndAlign are not valid.  Please remove following argument and call again:'); disp(p.Unmatched); keyboard; return; end
freshJackandSplit  = p.Results.freshJackandSplit || FORCE_RAW_REEXTRACTION;
freshExtractAlign  = p.Results.freshExtractAlign || FORCE_EVENT_REEXTRACTION;
BATCH              = p.Results.batch || BATCH;
SKIP_PROCESS_REREF = p.Results.skipProcessReref;
disp_flag = p.Results.disp;



%- confirm subject directory exists
subjDir = fullfileEEG(rootEEGdir,subj);
if ~exist(subjDir,'dir'), error('Error: root eeg directory: %s \n       does not contain subject: %s',rootEEGdir, subj); end
fprintf('\n\n\n**************************************************************************************'); %%- Command line output
fprintf('\n*************************************** %s ***************************************', subj);     %%- Command line output
fprintf('\n**************************************************************************************');     %%- Command line output

if BYPASS_ALIGNMENT_USER_QUERY==1, fprintf('\n WARNING: BYPASS_ALIGNMENT_USER_QUERY set to 1... default state should be 0\n'); end
if SKIP_EXTRACTION==1,             fprintf('\n WARNING: SKIP_EXTRACTION set to 1... default state should be 0.  Raw Files will NOT be split.\n'); end



%%- delete all EEG split files + jacksheetMaster.csv for a fresh build
if freshJackandSplit,
    if ~BATCH,
        reply = input('\n\n  User selected "fresh jacksheet and RAW split". \n  Delete eeg.reref, eeg.noreref, eeg.proceeded, and jackSheetMaster.csv NOW?  Y/N [N]:','s');
        if isempty(reply)
            reply = 'N';
        end
    else
        reply = 'Y';
    end
    fprintf('\n');
    
    if ( reply(1)=='Y' || reply(1)=='y' )
        
        %- delete directories
        dir2delete = {'eeg.noreref' 'eeg.reref' 'eeg.processed' 'eeg.processedBP'};
        for ii = 1:length(dir2delete),
            thisDir = fullfileEEG(subjDir, dir2delete{ii});
            whtSpc  = ''; if length(thisDir)<80, whtSpc(80-length(thisDir))=' '; end %- just to make fprint output pretty
            if exist(thisDir,'dir'),
                [SUCCESS,MESSAGE,MESSAGEID] = rmdir(thisDir,'s');
                if SUCCESS, 
                    fprintf('\n %s %s --> deleted', thisDir, whtSpc);
                else
                    fprintf('\n %s %s --> error when deleting', thisDir, whtSpc); 
                    keyboard
                end
            else
                fprintf('\n %s %s --> dir not found, so not deleting it', thisDir, whtSpc);
            end
        end %- ii
        
        %- delete specific files
        file2delete = {'docs/jacksheetMaster.csv', 'docs/raw_info.csv', 'docs/raw_info.csv.temp.csv', 'docs/jacksheetMaster.csv.bak', 'docs/leads_bp.txt', './finalChecklist.csv', 'tal/leads_bp.txt'};
        for ii = 1:length(file2delete),
            thisFile = fullfileEEG(subjDir, file2delete{ii});
            whtSpc   = ''; if length(thisFile)<80, whtSpc(80-length(thisFile))=' '; end
            if exist(thisFile,'file')
                delete(thisFile);
                if ~exist(thisFile,'file')
                    fprintf('\n %s %s --> deleted', thisFile, whtSpc);
                else
                    fprintf('\n %s %s --> error when deleting', thisFile, whtSpc); 
                    keyboard
                end
            else
                fprintf('\n %s %s --> file not found, so not deleting it', thisFile, whtSpc); 
            end
        end %- ii
        
        %- search for and delete other stuff that shouldn't be in the subject folder
        fileRoots2delete = {'docs/electrodes.m', 'docs/jacksheet.txt', 'docs/jacksheetMaster.txt', 'docs/electrodes_*.png', 'docs/*.bak*', 'tal/*.bak*', 'behavioral/*.bak*'};
        for ii = 1:length(file2delete),
            thisFileList       =       dir(fullfileEEG(subjDir, fileRoots2delete{ii}));
            [thisFilePath,~,~] = fileparts(fullfileEEG(subjDir, fileRoots2delete{ii}));
            for iii = 1:length(thisFileList),
                thisFile = fullfileEEG(thisFilePath,thisFileList(iii).name);
                whtSpc   = ''; if length(thisFile)<80, whtSpc(80-length(thisFile))=' '; end
                if thisFileList(iii).isdir,
                    [SUCCESS,MESSAGE,MESSAGEID] = rmdir(thisFile,'s');
                else
                    delete(thisFile);
                    SUCCESS = ~exist(thisFile,'file');
                end
                if SUCCESS, 
                    fprintf('\n %s %s --> deleted', thisFile, whtSpc);
                else        fprintf('\n %s %s --> error when deleting', thisFile, whtSpc); 
                    keyboard
                end
            end
        end %- ii
        
        %- final checklist is a special case, because it gets notes in it. So rename finalChecklist.invalid.csv so it can be merged at the end
        thisFileSrc = fullfileEEG(subjDir,'docs/finalChecklist.csv');
        thisFileDst = fullfileEEG(subjDir,'docs/finalChecklist.invalid.csv');
        if exist(thisFileSrc,'file'),
            [SUCCESS,MESSAGE,MESSAGEID] = movefile(thisFileSrc,thisFileDst,'f');
            whtSpc   = ''; if length(thisFileSrc)<80, whtSpc(80-length(thisFileSrc))=' '; end
            if SUCCESS, 
                fprintf('\n %s %s --> renamed to docs/finalCHecklist.invalid.csv', thisFileSrc, whtSpc);
            else
                fprintf('\n %s %s --> error when moving', thisFileSrc, whtSpc); 
                keyboard
            end
        end
        
    end %- if reply(1)='Y'
end





%%- check for the jacksheetMaster.txt and create if not found.  This would normally happen during a call to nk_split, but add here for subjects that are already extracted
info     = getElementInfo(subj, rootEEGdir);
tagNames = info.tagName;
pulseTag = info{strcmpi(info.chanType,'SYNC'),'tagName'};

%fprintf('\n checking for jacksheetMaster.csv: ');
jackMaster_file = fullfileEEG(subjDir, 'docs', 'jacksheetMaster.csv');
leadsBP_file    = fullfileEEG(subjDir, 'docs', 'leads_bp.txt');


if ~exist(jackMaster_file, 'file')
    fprintf('\n\n No jacksheetMaster.csv found. Creating it now...');
    createMasterJack(subj, rootEEGdir);
    
elseif isOldJacksheet(subj, rootEEGdir)
    fprintf('\n\n jacksheetMaster.csv already exists but is old (element.info has been modified since it was created). \n Remaking jacksheet now...');
    createMasterJack(subj, rootEEGdir);
    
elseif ~exist(leadsBP_file, 'file'),
    fprintf('\n\n No leads_bp.txt found. Making fresh jackSheetMaster and leads_bp.txt now...');
    createMasterJack(subj, rootEEGdir);
end



%%%%%%- automatic selection of "useStim".
%  for now only split those RAWs if "stimMap" exists in behavioral, which started when Tim Sheehan create a stim GUI with the cerestim, circa NIH037 or so
%  subjects < NIH037 may have annotated stim mapping sessions, but not well controlled experimentally and not automatically extracted to subj/behavioral
taskListFound = dir(fullfileEEG(rootEEGdir,subj,'behavioral'));
subjNum       = str2num(subj(4:6)); %- NIHXXX or BEHXXX
if sum(strcmp({taskListFound.name},'stimMapAnn'))==0 & subjNum <= 45 & strcmp(subj(1:3),'NIH'),
    fprintf('\n  HEADS UP: behaviora/stimMapAnn NOT found, attempting to create from raw/STIM_MAP');
    fprintf('\n            attempting to create stimMapAnn directly from raw/STIM_MAP annotation files');
    createStimMapSessFromAnn_v01(rootEEGdir, subj);
end



if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>0,
        fprintf('\n uh oh: %s left %d files open:\n', subj, length(fIDs));
        for iF=1:length(fIDs), filename = fopen(fIDs(iF)); fprintf(' %d) %s\n',iF,filename); end
        %keyboard;
        fclose all;
    end
    numFIDS = length(fIDs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- search RAW directory
%%--        if raw data there extract the details (name/date)
%%--        confirm whether already extracted to noreref (extract if not)
%%--        look for associated sync files and pull TTL channel numbers if so
%%--
fprintf('\n\n[STEP 1/5] Searching RAW and/or EEG.NOREREF for Raw and/or Extracted Channel files:\n');



MISSING_RAW_AND_EXTRACTED = 0;
rawDir  = fullfileEEG(rootEEGdir,subj,'raw');
rawList = []; % struct list
extractedList = [];
rawTimes = NaN;

cRaw21s = getRaw21Es(rawDir);
stimMask = cellfun(@(x) ~isempty(strfind(x, 'STIM_MAP')), cRaw21s);
cRaw21s_stim = cRaw21s(stimMask);

%if ~useStim
%    % don't include in list to split to noreref
%    cRaw21s = cRaw21s(~stimMask);
%end

if (isempty(cRaw21s) || SKIP_EXTRACTION == 1)
    MISSING_ALL_RAW = 1;
    
    %%- Missing Raw, but possibly eeg.noreref (and sync files) already generated.  Check for that here:
    extractedList = [];
    for i = 1 : length(cRaw21s)
        raw21 = cRaw21s{i};
        [~,filename] = fileparts(raw21);
        peices = strsplit(filename, '_');
        extractFileroot = fullfileEEG(rootEEGdir, subj, strjoin(peices(1,3), '_')); % subj_date_time
        if findAnyExtractions(subj, rootEEGdir, extractFileroot)
            extractedList = [extractedList extractFileroot];
        end
    end
    
    if isempty(extractedList)
        MISSING_RAW_AND_EXTRACTED = 1;
        fprintf(' RAW EEG data not found... just output EEG times to grab based on session.log files:');
    else
        fprintf(' RAW EEG data not found, but %d EEG.NOREREF extractions found. Attempting to match to session.log files',length(extractedList));
    end
    
    MISSING_ALL_SYNC = 1; %assume no sync files at all, will be set to 0 if a single sync file found
    for iList=1:length(extractedList),
        
        %- read the tagName.txt list... last entry should specify the channel with pulses (either ECG, EKG, or DC)
        if (iList==1),
            if isempty(getTagNames(subj, rootEEGdir)) %- decommision "useOldSource"
                fprintf('\n WARNING: Could not get tag names. Guessing that EKG is the pulse channel.\n');
                pulseTag = {'EKG'}; %assume it is EKG...
            end
        end
        
        %-pull extracted root name and directory
        fName = extractedList(iList);
        extractRootPath = fullfileEEG(rootEEGdir, subj, 'eeg.noreref', fName);
        extractedList(iList).rootName = fName;
        extractedList(iList).extractRootPath = extractRootPath;
        
        %-pull extracted date/time for sorting relative to eeglog times later
        dateStr = fName(find(fName=='_', 1 )+1:end);
        %Cstr = textscan(fName,'%s%s%s','delimiter','_');
        extractedList(iList).dateN = datenum(dateStr,'yymmdd_HHMM',2001);  %--new version of date
        
        %-pull TTL channels from the jacksheet (used for creating sync file)
        jackFile        = fullfile(extractRootPath,'jacksheet.txt');
        if ~exist(jackFile,'file'), error(' ERROR: %s is missing\n', jackFile); end
        
        %-pull TTL channels from the jacksheet (used for creating sync file)
        [jackChans, jackNames] = textread2(jackFile,'%s %s');
        iTTLchan        = find(strncmp(jackNames,pulseTag,length(pulseTag)));  % last entry in tagNames.txt indicates sync channel
        TTLchan1        = str2double(jackChans{iTTLchan(1)});
        if length(iTTLchan)<2 || strcmp(pulseTag,'DC'),
            TTLchan2    = TTLchan1;
        else
            TTLchan2    = str2double(jackChans{iTTLchan(2)});
        end
        %fileTTLchan1    = sprintf('%03d', TTLchan1);
        %fileTTLchan2    = sprintf('%03d', TTLchan2);
        
        %-look for any sync file with the proper channel prefix... if not found try the predicted sync file names
        syncFileList = dir('*.trigDC09.sync.txt');
        if length(syncFileList)==1,
            syncFileFound = fullfileEEG(rootEEGdir,subj,'eeg.noreref',syncFileList(1).name);
        else
            if length(syncFileList)>1, fprintf('\n WARNING: more than 1 sync file found for %s. \n', extractRootPath); end
            syncFileFound = '';
        end
        
        %-look for expected specific sync file names
        syncFileExpect = fullfile(extractRootPath,'%03d.%03d.sync.txt', TTLchan1, TTLchan2);
        syncFile = syncFileExpect;
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'%03d.%03d.sync.txt', TTLchan2, TTLchan1); end
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'trigDC09.sync.txt'); end
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'trigDC10.syncStim.txt'); end
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'trigDC11.syncStim.txt'); end
        if (length(dir(syncFile))==0), syncFile = syncFileExpect; end %if don't find it try the two other methods, go back to original expectation for warning/error message
        
        
        %-does expected file match "found" file?
        if length(syncFileFound)>0 && ~strcmp(syncFile,syncFileFound),
            fprintf('\n WARNING: found [and will use] syncFile %s, but expected %s\n', syncFileList(1).name, syncFile);
            syncFile = syncFileFound;
        end
        
        %- confirm sync file exists and save the location
        if (length(dir(syncFile))>0)  hasSync = 1; syncStr = 'sync found';   MISSING_ALL_SYNC = 0; % allow for alignment even if just 1 sync present
        else                          hasSync = 0; syncStr = 'sync MISSING'; syncFile = '';  end
        
        %- save to the extractedList
        extractedList(iList).jackChanStr = sprintf('TTL chan: %d, %d', TTLchan1, TTLchan2);
        extractedList(iList).hasSync  = hasSync;
        extractedList(iList).syncStr  = syncStr;
        extractedList(iList).syncFile = syncFile;
        extractedList(iList).pairedToAlign = 0 ; % initialize to zero here... this will be a counter of the number of behavioral files paired with this raw
        
    end
    if (MISSING_ALL_SYNC==0)
        extrTimes = [extractedList.dateN];                                    %extracted times
        syncTimes = [extractedList([extractedList.hasSync]==1).dateN];  %sync times
    else
        fprintf('\n WARNING: at least 1 sync missing, but somewhat iffy code for estimating sync file name when RAW is missing... assumes pulses on "EKG".');
        extrTimes = [];
        syncTimes = [];
    end
    
else
    MISSING_ALL_RAW = 0;
    MISSING_EXTRACTION = 0;
    
    
    %%- create list of raw EEG times (to compare with eeglog list) (dir gets time).
    %%-         use .21E file creation date (.EEG creation date is not the same)
    for i = 1 : length(cRaw21s)
        rawFile = dir(cRaw21s{i});
        rawFile.rawPath = fullfileEEG(cRaw21s{i});
        parent = fileparts(cRaw21s{i});
        rawFile.rawDir  = fileparts(parent);
        rawFile.clnPath = rawFile.rawPath(length(rawFile.rawDir) + 2 : end); % rawPath = rawDir / clnPath
        if ~isempty(strfind(cRaw21s{i},'21E'));
            rawFile.eegPath = [rawFile.rawPath(1:end-3) 'EEG'];
        elseif ~isempty(strfind(cRaw21s{i},'ns3')) || ~isempty(strfind(cRaw21s{i},'nev'))
            rawFile.eegPath = [rawFile.rawPath(1:end-3) 'ns2'];
        elseif ~isempty(strfind(cRaw21s{i},'TRC'));
            rawFile.eegPath = [rawFile.rawPath(1:end-3) 'TRC'];
        end
        rawFile.isStim  = ismember(rawFile.rawPath, cRaw21s_stim); %- isStim specified by raw's coming from raw/STIM/XXXX_XXX folder
        rawList = [rawList rawFile];
    end
    
    %%- if any raw found, figure out the expected extraction file name
    for iList=1:length(rawList)
        
        EEG_file = rawList(iList).eegPath;  %switch the suffix from .21E to .EEG
        if ~exist(EEG_file,'file'), error('MISSING Raw .EEG file, but .21E file found: %s', cRaw21s(iList).rawPath); end; %% should never happen
        if ~isempty(strfind(EEG_file,'ns2'))
            % getting datetime manually
            NS2FileDir = EEG_file;
            temp = strfind(NS2FileDir,'/'); dateTime = NS2FileDir(temp(end-1)+1:temp(end)-1);
            extractRootName = dateTime;
            % getting datetime through header info
            %         NS2FileDir = EEG_file;
            %         NS2data = openNSx(NS2FileDir,'noread'); % load the header info for EEG
            %         T_year = num2str(NS2data.MetaTags.DateTimeRaw(1));
            %         T_year = T_year(3:4); T_year = str2double(T_year); % wtf
            %         T_month = NS2data.MetaTags.DateTimeRaw(2);
            %         T_day = NS2data.MetaTags.DateTimeRaw(4);
            %         T_hour = NS2data.MetaTags.DateTimeRaw(5);
            %         T_minute = NS2data.MetaTags.DateTimeRaw(6);
            %         extractRootName = sprintf('%02d%02d%02d_%02d%02d', T_year,T_month,T_day,T_hour,T_minute);
            %         clear NS2data
        
        elseif ~isempty(strfind(EEG_file,'TRC')),
            % JW code
            eeg_emptyset = struct('event',[]);
            % ---------------- Opening File------------------
            PARAM.filename=EEG_file;
            trcfile=PARAM.filename;
            fid=fopen(trcfile,'r');
            if fid==-1
                error('Can''t open *.trc file')
            end
            TRC=eeg_emptyset;
            %------------------reading patient & recording info----------
            fseek(fid,64,-1);
            surname = char(fread(fid,22,'char'))';
            name    = char(fread(fid,20,'char'))';
            if ~strcmp(surname(1:5),'*****'),
                fprintf('\n found patient info in file: %s,%s... wiping name now',surname,name);
                fclose(fid);
                fid=fopen(trcfile,'r+'); %- re-open for writing
                fseek(fid,64,-1);
                stars = '*************************';
                fwrite(fid,stars(1:22)); %- blank out surname
                fwrite(fid,stars(1:20)); %- blank out name
                fclose(fid);
                %- now confirm overwritten
                fid = fopen(trcfile,'r');
                fseek(fid,64,-1);
                surname = char(fread(fid,22,'char'))';
                name    = char(fread(fid,20,'char'))';
                fprintf(' \n new name/surname = %s,%s\n',surname,name);
            end
            %- recording date and time
            fseek(fid,128,-1);
            T_day   = fread(fid,1,'char');
            if length(num2str(T_day))<2; day = ['0' num2str(T_day)];
            else; day = num2str(T_day);; end
            T_month = fread(fid,1,'char');
            switch T_month
                case 1; month='JAN';
                case 2; month='FEB';
                case 3; month='MAR';
                case 4; month='APR';
                case 5; month='MAY';
                case 6; month='JUN';
                case 7; month='JUL';
                case 8; month='AUG';
                case 9; month='SEP';
                case 10; month='OCT';
                case 11; month='NOV';
                case 12; month='DEC';
            end
            T_year = fread(fid,1,'char')+1900-2000;
            year   = num2str(T_year);
            %- recording time
            T_hour   = fread(fid,1,'char');
            T_minute = fread(fid,1,'char');
            T_second = fread(fid,1,'char');
            fclose(fid);
            %fprintf(' Date of session: %d/%d/%d\n',T_month,T_day,T_year)
            %fprintf(' Time at start: %02d:%02d:%02d\n',T_hour,T_minute,T_second)
            strTime      = sprintf('%d/%d/%d %02d:%02d:%02d',T_month,T_day,T_year,T_hour,T_minute,T_second); %
            extractRootName = sprintf('%02d%02d%02d_%02d%02d',T_year,T_month,T_day,T_hour,T_minute);     % new version: file stem of extracted channels -- JHW 11/2013
        
        else %if nihon kohden
            %open EEG file, skip over initial info, then pull the date and time
            %   (see nk_split for original version of the following code)
            fid = fopen(EEG_file, 'r');
            
            %1) seek to EEG1 control block to get offset to waveform block
            offsetToEEG1 = 146 ;                    % skips device info (128 byte), skips block ID (1 byte), device type (16 byte), number of blocks (1 byte)
            fseek(fid,offsetToEEG1,'bof');          % fseek(fileID, offset, origin) moves to specified position in file. bof=beginning of file
            
            %2) seek to EEG2 waveform block to get offset of actual data
            offsetToEEG2 = fread(fid,1,'*int32');
            offsetToEEG2 = offsetToEEG2 + 18 ;      % skips block ID (1 byte), device type (16 byte), number of blocks (1 byte)
            fseek(fid,offsetToEEG2,'bof');
            
            %3) seek to actual data, skip over initial info then read date/time
            blockAddress = fread(fid,1,'*int32');
            blockAddress = blockAddress + 20 ;      % skips block ID (1 byte), device type (16 byte), number of blocks (1 byte), byte length of one data (1 byte), mark/event flag (1 byte)
            fseek(fid,blockAddress,'bof');          %
            
            %%- annonomous function to convert binary to decimal.  input is binary string created with dec2bin
            bcdConverter2 = @(strDec2bin)  10*bin2dec(strDec2bin(1:4)) + bin2dec(strDec2bin(5:8));
            
            % get the start time
            T_year   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
            T_month  = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
            T_day    = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
            T_hour   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
            T_minute = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
            fclose(fid);
            extractRootName = sprintf('%02d%02d%02d_%02d%02d', T_year,T_month,T_day,T_hour,T_minute);
        end
        
        % determine expected extraction root name and other extraction info (date, root, existance)
        extractRootPath = fullfileEEG(rootEEGdir,subj,'eeg.noreref',extractRootName);
        extractDateNum  = datenum(extractRootName,'yymmdd_HHMM',2001); % useful for sorting by absolute time difference
        
        rawList(iList).extractRootName = extractRootName ;  % ex) "70613_1057"
        rawList(iList).extractRootPath = extractRootPath ;  % ex) "/Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/170613_1057"
        rawList(iList).extractDateNum  = extractDateNum  ;  % ex) "7.3540e+05"
        
        % Assume one extraction means all channels extracted for a raw
        if ~FORCE_RAW_REEXTRACTION && findAnyExtractions(subj, rootEEGdir, extractRootName)
            rawList(iList).extracted = 1;
        else
            rawList(iList).extracted = 0;
            MISSING_EXTRACTION = 1;
        end
        
        rawList(iList).pairedToAlign   = 0 ; % initialize to zero here... this will be a counter of the number of behavioral files paired with this raw
    end
    
    
    %%- one or more raw files not extracted... extract now
    if (MISSING_EXTRACTION)
        needsExtr = find([rawList.extracted]==0);
        fprintf(' RAW EEG data found (%d files); but missing %d of %d EEG.NOREREF files: \n',  length(rawList), length(needsExtr), length(rawList) );
        
        %-selectively extract just the missing data
        MISSING_EXTRACTION = 0;
        for iList=needsExtr
            fprintf('\n>>>>>>>>>> extracting %s --> %s <<<<<<<<<<\n',rawList(iList).clnPath, rawList(iList).extractRootName);
            
            %- the workhorse... this will select the correct splitting code based on recording machine (NK, BK, CV)
            raw_split(subj, rootEEGdir, rawList(iList).eegPath, BATCH);
            
            %- Now confirm that data was extracted... drop error if not
            if findAnyExtractions(subj, rootEEGdir, rawList(iList).extractRootName)
                rawList(iList).extracted   = 1  ;
            else
                rawList(iList).extracted   = 0  ;
                MISSING_EXTRACTION = 1;
                fprintf('ERROR: %s not extracted to any channels\n',rawList(iList).clnPath);
                error('ERROR: just extracted but at least 1 RAW not extracted.'); % should never happen.. perhaps name mismatch
            end
        end
    else
        fprintf(' RAW EEG data found (%d files); all files already extrated to EEG.NOREREF', length(rawList))
    end
    
    %%- look for sync files (generated by using alignTool to convert TTL waveforms to pulse times)
    %%-    if missing, user needs to run alignTool to get pulse times from raw TTL waveforms
    MISSING_ALL_SYNC = 1;  % assume no sync files at all, will be set to 0 if a single sync file found
    for iList=1:length(rawList)
        
        
        %-all raw files are extracted at this point thanks to MISSING_EXTRACTION condition above... don't need to check whether extracted==1
        extractRootPath = rawList(iList).extractRootPath;
        jackFile        = fullfile(extractRootPath, 'jacksheet.txt');
        if ~exist(jackFile,'file'), error(' ERROR: %s is missing. You may need to delete eeg.noreref and rerun.\n', jackFile); end
        
        %-pull TTL channels from the jacksheet (used for creating sync file)
        fd = fopen(jackFile, 'r');
        jackCells = textscan(fd, '%s %s');
        fclose(fd);
        [jackChans, jackNames] = deal(jackCells{:});
        iTTLchan        = find(ismember(util_split_stringnum(jackNames), pulseTag));
        msg = ['\nPulse tag not found in raw jacksheet. You may need to delete jacksheetMaster.csv and eeg.noreref/ and rerun'...
            '.\nIf this fails, you may need to delete %s and associated files.'];
        assert(~isempty(iTTLchan), msg, jackFile);
        TTLchan1        = str2double(jackChans{iTTLchan(1)});
        if length(iTTLchan)<2 || strcmp(pulseTag,'DC')
            TTLchan2    = TTLchan1;
        else
            TTLchan2    = str2double(jackChans{iTTLchan(2)});
        end
        
        %-look for any sync file with the proper channel prefix... if not found try the predicted sync file names
        syncFileList = dir(fullfileEEG(extractRootPath,'*trigDC09.sync.txt'));
        if length(syncFileList)==1
            syncFileFound = fullfileEEG(extractRootPath,syncFileList(1).name);
        else
            if length(syncFileList)>1, fprintf('\n WARNING: more than 1 sync file found for %s. \n', extractRootPath); end
            syncFileFound = '';
        end
        
        %-look for expected specific sync file names
        %syncFileExpect = sprintf('%03d.%03d.sync.txt', TTLchan1, TTLchan2);
        %syncFileExpect = sprintf('%s.%s.sync.txt', jackNames{TTLchan1}, jackNames{TTLchan2});
        syncFileExpect = fullfile(extractRootPath,sprintf('%s.%s.sync.txt', jackNames{TTLchan1}, jackNames{TTLchan2}));
        syncFile = syncFileExpect;
        if (isempty(dir(syncFile))), syncFile = fullfile(extractRootPath,sprintf('%s.%s.sync.txt', jackNames{TTLchan2}, jackNames{TTLchan1})); end
        if (isempty(dir(syncFile))), syncFile = fullfile(extractRootPath, 'trigDC09.sync.txt'); end
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'trigDC10.syncStim.txt'); end
        if (length(dir(syncFile))==0), syncFile = fullfile(extractRootPath,'trigDC11.syncStim.txt'); end
        if (isempty(dir(syncFile))), syncFile = syncFileExpect; end %if don't find it try the two other methods, go back to original expectation for warning/error message
        
        %-does expected file match "found" file?
        if ~isempty(syncFileFound) && ~strcmp(syncFile,syncFileFound)
            fprintf('\n WARNING: found [and will use] syncFile %s, but expected %s\n', syncFileList(1).name, syncFile);
            syncFile = syncFileFound;
        end
        
        %- confirm sync file exists and save the location
        if (length(dir(syncFile))>0)  hasSync = 1; syncStr = 'sync found';   MISSING_ALL_SYNC = 0; % allow for alignment even if just 1 sync present
        else                          hasSync = 0; syncStr = 'sync MISSING'; syncFile = '';  end
        
        %- save to the rawList
        rawList(iList).jackChanStr = sprintf('TTL chan: %d, %d', TTLchan1, TTLchan2);
        rawList(iList).hasSync  = hasSync;      % ex) 1
        rawList(iList).syncStr  = syncStr;      % ex) sync found
        rawList(iList).syncFile = syncFile;     % ex) /Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/NIH016_280613_1411.083.084.sync.txt
    end
    rawTimes  = [rawList.extractDateNum];                              % used for plotting raw vs extracted vs behavior at end... define here to differentiate from case where raw is missing but extraction exists
    extrTimes = [rawList.extractDateNum];                              % extracted times
    syncTimes = [rawList([rawList.hasSync]==1).extractDateNum];  % sync times
    
    if (MISSING_ALL_SYNC || sum([rawList.hasSync])<length(rawList))
        %fprintf('\n\nAt least 1 sync file missing... USE alignTool to find TTL pulse times and create sync file\n\n');
        %- save pathdef.m for command-line instance of matlab (sans java) used to call alignTool for picking out peaks (not required... can call alignTool from graphical matlab command line)
        %savepath(fullfileEEG(rootEEGdir, subj, 'eeg.noreref/pathdef.m'))
        fprintf('; only %d sync files found in EEG.NOREREF', sum([rawList.hasSync]))
    else
        fprintf('; all sync files found in EEG.NOREREF', length(rawList))
    end
    
end % if (isempty(cRaw21s) || SKIP_EXTRACTION == 1)


if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open', length(fIDs)-numFIDS);
        keyboard;
        numFIDS = length(fIDs);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -- search "behavioral" directory
%%--        pull .eeglog times from all tasks, all sessions... make eeglog.up if not created already
%%--        confirm that behavioral events.mat already created... if not, create it
%%--
behDir   = fullfileEEG(rootEEGdir,subj,'behavioral');
fprintf('\n\n[STEP 2/5] Searching behavioral directory for eeg.eeglog and events.mat files:\n');


%%- generate list of tasks
dList = dir(behDir);
taskDir=[];
for i=1:length(dList),
    if dList(i).isdir == 1 && dList(i).name(1) ~= '.'
        taskDir = [taskDir dList(i)];
    end
end


%%- now exclude processsing of stimMapping folder for subj<40   %% everybody should be processed
%if useStim == 0,
%    taskDir = taskDir(~strcmp({taskDir.name},'stimMapping'));
%end


%- re-extract ALL events (unless exclued by "ONLY_THESE_TASKS".  This is a good step when doing a fresh rebuild of a subject's directory, as it will force aligment
if freshExtractAlign,
    
    if ~BATCH,
        reply = input('\n  User selected "fresh extract and align". Delete alignmentSummary, all events.mat, and re-extract all events now?  Y/N [Y]:','s');
        if isempty(reply)
            reply = 'Y';
        end
    else
        reply = 'Y';
    end
    fprintf('\n');
    
    
    if ( reply(1)=='Y' || reply(1)=='y' )
        %- remove a few key files so there is no confusion
        files2rm = {'alignmentPairing_auto.txt' 'alignmentStats.txt', 'alignmentSummary.txt', 'eventEegFilePaths.txt' 'readme_helper.txt'}; %- don't delete "alignment_forced.txt" if it exists
        for iF=1:length(files2rm),
            thisFile =  fullfileEEG(behDir,files2rm{iF});
            if exist(thisFile,'file'), delete(thisFile); end
        end
        
        %- delete all events.mat, session-level and aggregates
        delCnt=0;
        for iTask = [1:length(taskDir)],
            thisFile = fullfileEEG(behDir, taskDir(iTask).name, 'events.mat');
            if exist(thisFile,'file'), delete(thisFile); delCnt=delCnt+1; end
            
            sessions = dir(fullfileEEG(behDir, taskDir(iTask).name, 'session_*'));
            for iF=1:length(sessions),
                thisFile = fullfileEEG(behDir, taskDir(iTask).name, sessions(iF).name, 'events.mat');
                if exist(thisFile,'file'), delete(thisFile); delCnt=delCnt+1; end
            
                thisFile = fullfileEEG(behDir, taskDir(iTask).name, sessions(iF).name, 'eeg.eeglog.up');
                if exist(thisFile,'file'), delete(thisFile); delCnt=delCnt+1; end
            end
        end
        fprintf(' --> deleted %d events.mat and/or eeg.eeglog.up files\n', delCnt);
        
        %- now re-extract ALL events.mat
        for iTask = [1:length(taskDir)],
            behavioralProcessing(subj, rootEEGdir, taskDir(iTask).name, 'behavioral');         % extract iEEG tasks
        end
    end
end


%%- for each task generate list of session folders and extract the relevant info
allEEGlogs=[]; allSessionEvents=[]; taskStrAr=[]; numAlignedWithBlank=0;
for iTask = 1:length(taskDir)
    taskStr = taskDir(iTask).name;
    taskStrAr{iTask} = taskStr;
    
    
    sessDir = dir(fullfileEEG(behDir,taskStr,'session_*'));
    if length(sessDir)<1, fprintf(' WARNING -- task folder contains no session folders:  %s\n',taskStr); end
    
    
    %%- quick loop through session directories to correctly order double digit session numbers (else 10 comes after 1 and before 2)
    sessNumAr=[];
    for iSess = 1:length(sessDir)
        sessStr     = sessDir(iSess).name;
        strNumeric  = find( sessStr >= '0' & sessStr <= '9');
        sessNum     = str2double( sessStr(strNumeric) );  if isempty(sessNum), sessNum=iSess; end; %shouldn't need this catch...
        sessNumAr(iSess) = sessNum;
    end
    if length(sessDir)>0
        [~, sortInd] = sort(sessNumAr);
        sessDir = sessDir(sortInd);
    end
    
    
    for iSess = 1:length(sessDir)
        sessStr = sessDir(iSess).name ;
        
        eeglogfStr  = fullfileEEG(behDir,taskStr,sessStr,'eeg.eeglog');
               
        %- if stimMapping and no eeg.eeglog, try making it now
        if ~exist(eeglogfStr, 'file') & (strcmp(taskStr,'stimMapAnn') | strcmp(taskStr,'stimMapGUI')),
            %- for stimMapping, the eeg.eeglog is created during events extraction
            fprintf('\n HEADS UP: stimMap folder was missing eeg.eeglog... running eventExtraction code for stimMap now');
            behavioralProcessing(subj, rootEEGdir, taskStr, 'behavioral');
            fprintf('\n --------------------------------------------------------- \n\n');
            if ~exist(eeglogfStr, 'file')
                fprintf('\n EEG log not found. what up?');
                keyboard;
            end
        end
        
        %- error check... need eeg.eeglog (even just a fake one) to figure out when the session was run
        if ~exist(eeglogfStr, 'file')
            error('eeg.eeglog file not found: %s\n', eeglogfStr);
        end
        eeglogFile  = dir(eeglogfStr);
        dateStr     = eeglogFile.date;
        dateMAT_dir = eeglogFile.datenum;  % matlab datenum from directory listing (when was file created and/or modified
        dateStr_dir = datestr(dateMAT_dir,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
        
        %read mstime from eeg.eeglog to get pyepl start time (instead of just file creation/modification date)
        [mstimes]   = textread(eeglogfStr,'%n%*[^\n]');   %read all ms time from the pulse file
        dateMAT_act = datenum(epoch2date(mstimes(1)));     % *MST switch to epoch2date (JW used some java magic function... epoch2date is way better)
        dateStr_act = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
        dateMAT_actEnd = datenum(epoch2date(mstimes(end))); % *MST switch to epoch2date
        
        thisEEGlog.task     = taskStr;
        thisEEGlog.taskNum  = iTask;
        thisEEGlog.sess     = sessStr;
        thisEEGlog.taskDir  = fullfileEEG(behDir,taskStr);
        thisEEGlog.taskSess = fullfileEEG(taskStr,sessStr);
        thisEEGlog.dateMAT_dir = dateMAT_dir ;
        thisEEGlog.dateStr_dir = dateStr_dir ;
        thisEEGlog.dateMAT_act = dateMAT_act ;
        thisEEGlog.dateStr_act = dateStr_act ;
        thisEEGlog.dateMAT_actEnd = dateMAT_actEnd ;
        
        %- check for session lasting >90 min... probably means there is an error in the eeg.eeglog entries
        if strcmp(taskStr,'stimMapAnn') | strcmp(taskStr,'stimMapGUI'), addTimeDeltaSession = 25; else  addTimeDeltaSession=0; end;
        if abs(thisEEGlog.dateMAT_actEnd - thisEEGlog.dateMAT_act) > datenum('01/01/01 2:30')-datenum('01/01/01 1:00')+addTimeDeltaSession,
            fprintf(' SERIOUS WARNING: %s date discrepancy: end time >%d min different from start time %s vs %s\n',fullfileEEG(thisEEGlog.taskSess,'eeg.eeglog'),90+addTimeDeltaSession,datestr(dateMAT_act,'mm/dd/yy HH:MM PM'),datestr(dateMAT_actEnd,'mm/dd/yy HH:MM PM'));  %-- check conversion from dir to actual date
            if ~BATCH
                reply = input('               check session.log and eeg.eeglog file and remove spurious start times if appropriate.\n...  PRESS RETURN TO CONTINUE (break and rerun eegPrepAndAlign if eeg.eeglog has been modified)...\n');
            end
        end
        
        % check to see if eeglog.up file already created... if not, make it
        eeglogStr   = fullfileEEG(behDir,taskStr,sessStr,'eeg.eeglog');
        eeglogUpStr = fullfileEEG(behDir,taskStr,sessStr,'eeg.eeglog.up');
        if (isempty(dir(eeglogUpStr))),
            fixEEGLog(eeglogStr,eeglogUpStr);
        end;
        thisEEGlog.eeglogStr   = eeglogStr ;
        thisEEGlog.eeglogUpStr = eeglogUpStr ;
        
        % find the session log... alignment will created a modified version (session.log.align) that includes pointers to the eeg file
        sessionLogStr = fullfileEEG(behDir,taskStr,sessStr,'session.log');
        thisEEGlog.sessionLogStr = sessionLogStr;
        
        % load in the session log, and convert the first and last time entry into date numbers (instead of relying on directory date num) [possibly comment out following lines... no need to open/read file]
        extractSessionLogDates = 0;
        if extractSessionLogDates,
            [mstimes] = textread2(sessionLogStr,'%n%*[^\n]');
            sessStartDateMAT = datenum(char(cell(javaSDF.format(Date(mstimes(1))))));  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
            sessStartDateStr = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
            sessEndDateMAT   = datenum(char(cell(javaSDF.format(Date(mstimes(end))))));  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
            sessEndDateStr   = datestr(dateMAT_act,'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
            %fprintf('\n  session.log dates: %s (dir listing), %s (mstime conversion)',dateStr_dir,dateStr_act);
        end
        
        % look for events.mat... confirm whether alignment has already happened
        sessEventsStr = fullfileEEG(behDir,taskStr,sessStr,'events.mat');
        if (exist(sessEventsStr,'file')==0)
            eventsExists  = 0;
            eventsAligned = 0;
            eventsEEGfile = '';
            eventStr      = sprintf('events.mat MISSING');
        else
            eventsExists = 1;
            % events file exists... check to see whether already aligned
            events = [];
            load(sessEventsStr);
            if (isfield(events,'eegfile') && isfield(events,'eegoffset'))
                eventsAligned = 1;
                eventsEEGfile = events(1).eegfile;
                eventStr      = sprintf('events.mat aligned');
                numEmptyField = sum(strcmp({events.eegfile},''));
                if numEmptyField>0, eventStr = sprintf('%s [%d blank eegfile field]',eventStr,numEmptyField); numAlignedWithBlank = numAlignedWithBlank+1;end
            else
                eventsAligned = 0;
                eventsEEGfile = '';
                eventStr      = sprintf('events.mat NOT ALIGNED');
            end
        end
        thisEEGlog.eventsExists  = eventsExists;
        thisEEGlog.eventsAligned = eventsAligned;
        thisEEGlog.eventsEEGfile = eventsEEGfile;
        thisEEGlog.eventsFile    = sessEventsStr;
        thisEEGlog.eventsStr     = eventStr;
        
        % create list of all event logs
        allEEGlogs = [allEEGlogs thisEEGlog];
        
    end
end

%% - EVENTS MISSING: extract behavioral events now?
assert(~isempty(allEEGlogs), 'No eeg files found. Check that your behavioral and raw directories exist');
if ( sum([allEEGlogs.eventsExists])<length(allEEGlogs) )
    
    
    % create a list of potential tasks to extract
    iNeedEvents         = find( [allEEGlogs.eventsExists]==0 );
    iTaskList           = [allEEGlogs(iNeedEvents).taskNum];
    [uniqTask,iUniq,iC] = unique(iTaskList);  % iUniq is index into iNeedEvents that points to unique events (i think)
    
    fprintf(' EXTRACTION of behavioral events.mat. Missing %d total events.mat files from %d tasks: \n', length(iNeedEvents), length(uniqTask));
    for iList = 1:length(iUniq)
        fprintf('         %d events.mat missing from %s\n', length(find([allEEGlogs(iNeedEvents).taskNum]==allEEGlogs(iNeedEvents(iUniq(iList))).taskNum)), allEEGlogs(iNeedEvents(iUniq(iList))).task );
    end
    reply='Y';
    if (BYPASS_ALIGNMENT_USER_QUERY==0) && ~BATCH
        reply = input('Attempt to extract missing events.mat files from session.logs now? Y/N [Y]:','s');
        if isempty(reply)
            reply = 'Y';
        end
        fprintf('\n');
    end
    
    if ( reply(1)=='Y' || reply(1)=='y' )
        
        % behavioral processing does all sessions of task... only call once per task
        thisDir = pwd;
        iExtracted = [iNeedEvents];      % start with list of iNeedEvents
        for iList = 1:length(iUniq)
            thisEEGlog = allEEGlogs(iNeedEvents(iUniq(iList)));
            behavioralProcessing(subj,rootEEGdir,thisEEGlog.task);  % in eeg_toolbox/events: should create events for all sessions of the selected task...
            iExtracted = [iExtracted find(strcmp({allEEGlogs.task},thisEEGlog.task))];
        end
        chdir(thisDir);
        
        % now double check to see whether events were created
        iExtracted = unique(iExtracted);  % contains all iNeedEvents + any additional events that were extracted because they were in the same task as an iNeed
        eventExtStr = '\n';
        for iList = iExtracted
            thisEEGlog = allEEGlogs(iList);
            
            % look for events.mat... confirm whether alignment has already happened
            sessEventsStr = thisEEGlog.eventsFile;
            padTaskStr = sprintf('%s/%s',thisEEGlog.task, thisEEGlog.sess);
            padTaskStr(end+1:25) = ' ';
            if (exist(sessEventsStr,'file')==0)
                eventExtStr    = sprintf('%s WARNING: %s  NOT extracted\n', eventExtStr, padTaskStr);
            else
                %eventExtStr    = sprintf('%s%s  extracted\n',     eventExtStr, padTaskStr);
                eventStr       = sprintf('events.mat NOT ALIGNED');
                
                allEEGlogs(iList).eventsStr     = eventStr;
                allEEGlogs(iList).eventsExists  = 1;
                allEEGlogs(iList).eventsAligned = 0; % can't be aligned if just extracted
                
            end
        end
        fprintf(eventExtStr)
    end
end
fprintf(' BEHAVIORAL EEG.EEGLOG found (%d files); %d events.mat files found', length(allEEGlogs), sum([allEEGlogs.eventsExists]) );



if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open:\n ', length(fIDs)-numFIDS);
        for iF=1:length(fIDs), filename = fopen(fIDs(iF)); fprintf(' %d) %s\n',iF,filename); end
        keyboard;
        numFIDS = length(fIDs);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- PREP-ALIGNMENT: attempt to match eeglog times with raw times
%%--        - check to see whether "alignmentParing_force.txt" exists.. if so, use that to chose alignment pairs
%%--        - otherwise, automatically chose pairs by identifying gap > 1h between eeglog times and label those gaps as "isNewSession"
%%--         try to match raw file times with "NewSession" eeglogs
%%--
fprintf('\n\n[STEP 3/5] Preparing for alignment: identify pairings of eeg.eeglog <--> extracted sync pulses:\n');

%%- define delta time thresh for auto detection (will be ignored if alignment pairs are forced
deltaTimeThresh   = datenum('01/01/01 2:00')-datenum('01/01/01 1:00');  %any tasks>1h apart are likely from different sessions
[allDatesS,iLogs] = sort([allEEGlogs.dateMAT_act]);


%%- check for ALIGNMENT_PAIRS: read "alignmentPairing_force.txt" which lists attempted pairs. a precursor to this file is generated automatically if it's not found
fileForcedAlign = fullfileEEG(behDir,'alignmentPairing_forced.txt');
FORCE_PAIRS = 0;
if exist(fileForcedAlign,'file')
    fid = fopen(fileForcedAlign,'r');
    x = textscan(fid,'%s [ %[^]] ] <--> %s [ %[^]] ]');  %x is 4 cell arrays, one with beh_dir, beh_date, eeg_root, eeg_date
    fclose(fid);
    
    iListAr=1:find(strcmp(x{1},'***'))-1;  %look at all entries up to the '*****' row
    forcedPairs = [];
    for iList=iListAr
        this_behSession_dir  = x{1}{iList};
        this_eegExtract_root = x{3}{iList};
        
        this_behSession_dir(find(this_behSession_dir=='\')) = '/';  % shouldn't be necessary for any extraction using fullfileEEG, but older extractions from PC may have backslash
        
        iEegLog = find(strcmp({allEEGlogs.taskSess}, this_behSession_dir));
        
        if isempty(iEegLog)
            disp( {allEEGlogs.taskSess}' );
            disp( this_behSession_dir );
            error(' ERROR: alignmentPairing_forced.txt specifies raw file that is absent in raw directory\n');
        end
        
        if MISSING_ALL_RAW==0
            iRawOrExtrList = find(strcmp({rawList.extractRootName}, this_eegExtract_root));
        else
            iRawOrExtrList = find(strcmp({extractedList.rootName}, this_eegExtract_root));
        end
        % tricky, if there are multiple NSPs..
        forcedPairs(iList,1:2) = [iEegLog iRawOrExtrList(1)];
    end
    
    if size(forcedPairs,1)~=length(iLogs)
        error(' ERROR: alignmentPairing_forced.txt contains %d pairs, but %d eeg.eeglog files have been found! \nAll pairs must be specified\n', size(forcedPairs,1),length(iLogs));
    else
        fprintf(' Using FORCED Alignment Pairs (%d pairs specified in %s)\n', size(forcedPairs,1), fileForcedAlign);
    end
    FORCE_PAIRS = 1;
end


%%- loop through behavioral eeglogs (ordered by date, not task), and find the raw file with the closest date (unless alignmentPairing_force.txt overrides)
for iL=1:length(iLogs)
%     if (iL==1)
%         allEEGlogs(iLogs(iL)).deltaT = 999 ;
%     else
%         allEEGlogs(iLogs(iL)).deltaT = allEEGlogs(iLogs(iL)).dateMAT_act - allEEGlogs(iLogs(iL-1)).dateMAT_act;
%     end
%     if (allEEGlogs(iLogs(iL)).deltaT > deltaTimeThresh) isNewSess = 1;
%     else                                                isNewSess = 0; end;
%     allEEGlogs(iLogs(iL)).isNewSession = isNewSess;
    isNewSess = 1;
    allEEGlogs(iLogs(iL)).isNewSession = 1;
    
    %- if forced_pairs makes this raw different from the preceding raw, isNewSess should be triggered!
    if FORCE_PAIRS && iL>1
        thisForcedSess = forcedPairs(find(forcedPairs(:,1)==iLogs(iL)),2);
        lastForcedSess = forcedPairs(find(forcedPairs(:,1)==iLogs(iL-1)),2);
        if  thisForcedSess ~= lastForcedSess, isNewSess = 1;
        else                                  isNewSess = 0;  end;
    end
    
    rawStr = '';
    strAdd=sprintf('.................................................................................\n');
    if isNewSess && ~MISSING_ALL_RAW
        %-find the matching raw file
        diffList = allEEGlogs(iLogs(iL)).dateMAT_act - extrTimes;   %%
        
        %- normally grab the closest file in case laptop was started before
        %or after start/stop.  For stimMapGUI lots of short sessions
        %created so better to always make it positive definite 
        %11/2018 DY and JW
        if contains(allEEGlogs(iLogs(iL)).taskDir,'stimMapGUI'),
            delt= min(diffList(diffList > 0)); % finds the indices of the smallest POSITIVE value
            iList = find(diffList == delt,1);
        else
            iList = find(abs(diffList)<deltaTimeThresh);
        end
        
        if FORCE_PAIRS && ~isempty(find(forcedPairs(:,1)==iLogs(iL), 1))
            iList = forcedPairs(find(forcedPairs(:,1)==iLogs(iL)),2);
            assert(iList <= length(rawList), 'Something about forced alignment pairs is not quite right');
        end
        
        if (length(iList)==0)
            iList = find(abs(diffList)<deltaTimeThresh*2);
            fprintf(' WARNING: NO RAW FOUND WITHIN 1 HOUR OF EEGLOG %s/%s\n                  ... TRYING 2 HOURS --> %d found\n',allEEGlogs(iLogs(iL)).task,allEEGlogs(iLogs(iL)).sess, length(iList));
        end
        if (length(iList)==0)
            iList = find(abs(diffList)<deltaTimeThresh*3);
            if (length(iList)>0) strWarn = '[COULD BE WRONG RAW; CONFIRM HAS PULSES]';
            else                 strWarn = '[RAW APPEARS TO BE MISSING]'; end
            fprintf('                  ... TRYING 3 HOURS --> %d found %s\n',length(iList),strWarn);
        end
        if (length(iList)>1)
            if length(iList)==2 && strcmp(rawList(iList(1)).folder,rawList(iList(2)).folder),
                %- this is just from having 2 NSPs, no need to make a fuss
                iList = iList(1);
            else
                fprintf(' WARNING: MATCHED MULTIPLE RAW FILES TO EEGLOG %s/%s\n',allEEGlogs(iLogs(iL)).task,allEEGlogs(iLogs(iL)).sess);
                %[delt, iList] = min(abs(diffList));
                delt  = min(diffList(diffList > 0)); % Choose the smallest POSITIVE value to make sure the raw started first
                iList = find(diffList==delt,1,'first'); % JW found a bug 11/2018... iList was no longer indexing diffList... fixed with this line
            end
        end
        if (length(iList)>0)
            rawStr          = sprintf('%s>>%s  %s \t<<-- raw found, %s\n', strAdd, datestr(rawList(iList).extractDateNum,'mm/dd/yy HH:MM PM'), rawList(iList).clnPath, rawList(iList).syncStr);
            hasRaw          = 1 ;
            thisRawIndex    = iList;
            thisExtractRoot = rawList(iList).extractRootName;
            thisSyncFile    = rawList(iList).syncFile;  %can have a syncFile path even if it doesn't exist
            thisExtractDate = rawList(iList).extractDateNum;
            rawList(iList).pairedToAlign = rawList(iList).pairedToAlign+1;
            if (rawList(iList).hasSync==0)
                rawStr      = sprintf('%s--- need to run alignTool on %s, %s --- \n', rawStr, rawList(iList).extractRootName, rawList(iList).jackChanStr);
                hasSync     = 0 ;
            else
                hasSync     = 1 ;
            end
        else                        % matched raw not found
            rawStr = sprintf('%s>>****************************************** \t<<-- raw MISSING\n',strAdd) ;
            hasRaw          = 0 ;
            hasSync         = 0 ;
            thisRawIndex    = -1;
            thisExtractRoot = '';
            thisSyncFile    = '';
            thisExtractDate = nan;
        end
    elseif (isNewSess & MISSING_ALL_RAW==1)     % no raw's found
        
        rawStr = sprintf('%s>>****************************************** \t<<-- raw MISSING\n',strAdd) ;
        hasRaw          = 0 ;
        hasSync         = 0 ;
        thisRawIndex    = -1;
        thisExtractRoot = '';
        thisSyncFile    = '';
        thisExtractDate = nan;
        
        
        %-no raw file, but extracted files (and sync files) could exist
        if (MISSING_ALL_SYNC==0)
            diffList = allEEGlogs(iLogs(iL)).dateMAT_act - extrTimes;
            iList = find(abs(diffList)<deltaTimeThresh*2);  %give double the amount of time for extracted files
            
            if FORCE_PAIRS & ~isempty(find(forcedPairs(:,1)==iLogs(iL))), iList = forcedPairs(find(forcedPairs(:,1)==iLogs(iL)),2); end
            
            if (length(iList)>1)
                fprintf(' WARNING: MATCHED MULTIPLE EXTRACTED/SYNC FILES TO AN EEGLOG %d\n',iLogs(iL));
                [delt, iList] = min(diffList);
            end
            if (length(iList)>0)
                rawStr          = sprintf('%s>>%s  %s \t<<-- Extracted found; %s\n', rawStr, datestr(extrTimes(iList),'mm/dd/yy HH:MM PM'), extractedList(iList).rootName, extractedList(iList).syncStr);
                hasSync         = 1 ;
                thisExtractRoot = extractedList(iList).rootName;
                thisSyncFile    = extractedList(iList).syncFile;  %can have a syncFile path even if RAW doesn't exist
                thisExtractDate = extractedList(iList).dateN;
                extractedList(iList).pairedToAlign = extractedList(iList).pairedToAlign+1;
            end
        end
    end
    
    % associate raw data with task file: hasSync/SyncFile/ExtractFile only updated for new sessions
    allEEGlogs(iLogs(iL)).rawStr          = rawStr;
    allEEGlogs(iLogs(iL)).hasRaw          = hasRaw ;         % 1 or 0
    allEEGlogs(iLogs(iL)).hasSync         = hasSync ;        % 1 or 0
    allEEGlogs(iLogs(iL)).rawListIndex    = thisRawIndex ;   % -1 (if no raw), else index to rawList
    allEEGlogs(iLogs(iL)).extractRootName = thisExtractRoot; % ex) /Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/NIH016_130617_1057
    allEEGlogs(iLogs(iL)).syncFilePath    = thisSyncFile;    % ex) /Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/NIH016_130617_1057.trigDC10.sync.txt
    allEEGlogs(iLogs(iL)).extractDateN    = thisExtractDate; % matlab datenum, should match up with extracted channel file name root (e.g., NIH016_130617 --> 06/17/2013)
end



%%- ALIGNMENT_PAIRS: create "alignmentPairing_auto.txt" which lists attempted pairs.  this file can be altered and saved as "alignmentPairing_forced.txt" to override auto pairs
if FORCE_PAIRS==0, fileAlignmentPairs = sprintf('%s/alignmentPairing_auto.txt',behDir);
else               fileAlignmentPairs = sprintf('%s/alignmentPairing_forcedUsed.txt',behDir);
end
fid = fopen(fileAlignmentPairs,'w+');
for iList=1:length(allEEGlogs),
    thisEEGlog = allEEGlogs(iList);
    
    % output the behavioral directory and extracted file root... these can be rearranged to override auto pairs
    behSession_dir   = thisEEGlog.taskSess ;         % ex) playPass/session_7
    eegExtract_root  = thisEEGlog.extractRootName ;  % ex) NIH016_280613_1411
    
    
    behSession_dateStr = thisEEGlog.dateStr_act ;
    %should be able to get date number by searching rawList or extractedList... as long as raw pair was found!
    eegExtract_dateStr = '';
    if thisEEGlog.hasRaw, eegExtract_dateStr = datestr(thisEEGlog.extractDateN,'mm/dd/yy HH:MM PM'); end;
    
    
    % output to file and command line
    fprintf(fid, '%s [%s]\t <-->  %s [%s]\n', behSession_dir, behSession_dateStr,eegExtract_root,eegExtract_dateStr);
    %fprintf(     '%s [%s]\t <-->  %s [%s]\n', behSession_dir, behSession_dateStr,eegExtract_root,eegExtract_dateStr);
end
fprintf(fid, '\n*** completed list of extracted raw files [and dates] below... paste above then resave as "alignmentPairing_forced.txt" to force different pairs  ***\n');
%fprintf(     '\n*** completed list of extracted raw files (and dates) ***\n');
extractedString = {''};
for iList=1:length(rawList)+length(extractedList)
    if length(rawList)>0,
        thisExtractRoot = rawList(iList).extractRootName;
        thisExtractDate = datestr(rawList(iList).extractDateNum,'mm/dd/yy HH:MM PM');
    else
        thisExtractRoot = extractedList(iList).rootName;
        thisExtractDate = datestr(extractedList(iList).dateN,'mm/dd/yy HH:MM PM');
    end
    extractedString{iList} = sprintf('%s [%s]\n', thisExtractRoot,thisExtractDate);
    %fprintf(fid, '%s [%s]\n', thisExtractRoot,thisExtractDate);
    %fprintf(     '%s [%s]\n', thisExtractRoot,thisExtractDate);
end
orderedString = sort(extractedString);  %- this way they are listed in chronological order... easier to find targets
fprintf(fid,'%s',orderedString{:});
fclose(fid);
fprintf(' Alignment Pairing List saved to file %s\n   (modify and save as "alignmentPairing_forced.txt" to override auto-pairs)',fileAlignmentPairs);


if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open', length(fIDs)-numFIDS);
        keyboard;
        numFIDS = length(fIDs);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- ALIGNMENT and Updating Master EVENTS.mat
%%--         give user option to run alignment on any un-aligned pairs of behavioral & eeg files
%%--         propogate newly aligned events.mat files from each session to the aggregtaed events.mat for the task
%%--
fprintf('\n\n[STEP 4/5] Alignment: confirm all events.mat aligned; offer alignment if not:\n');

if FORCE_REALIGNMENT
    fprintf(' (forcing realignment of all alignable events)\n');
    iNeedAlign = find([allEEGlogs.eventsExists]==1 & [allEEGlogs.hasSync]==1 & [allEEGlogs.eventsAligned]==1);
    for iList = iNeedAlign, allEEGlogs(iList).eventsAligned = 0; end
end


%%- ALIGNMENT CHECK: are there any files that could/should be aligned?  If so, ask user if that is desirable
iNeedAlign   = find([allEEGlogs.eventsExists]==1 & [allEEGlogs.hasSync]==1 & [allEEGlogs.eventsAligned]==0);
iJustAligned = [];
if ( ~isempty(iNeedAlign) )
    
    fprintf('\nALIGNMENT of behavior to physio clock. %d events.mat can be aligned: \n', length(iNeedAlign));
    for iList = iNeedAlign,
        fprintf('         %s\n', fullfileEEG(allEEGlogs(iList).task, allEEGlogs(iList).sess, 'events.mat') );
    end
    reply='Y';
    if (BYPASS_ALIGNMENT_USER_QUERY==0) && ~BATCH
        reply = input('      Align now? Y/N [Y]:','s');
        if isempty(reply)
            reply = 'Y';
        end
        fprintf('\n');
    end
    
    if (reply(1)=='Y' || reply(1)=='y')
        
        %%- ALIGNMENT: Loop through Tasks/Sessions and align any events files that still need alignment
        fid = fopen(sprintf('%s/alignmentStats.txt',behDir),'a+');
        fprintf(fid,'\n\n=================================================================================================================');
        fprintf(fid,'\n=========================== Alignment stats from eegPrepAndAlign on %s ===========================', datestr(now,'mm/dd/yy HH:MM PM'));
        fprintf(fid,'\n=================================================================================================================');
        
        alignStr = '\nAlignment Summary:\n';
        for iList = iNeedAlign
            
            thisEEGlog = allEEGlogs(iList);
            
            %%- if (events not aligned & events exist (means eeglog exists) & sync exist (means chan extracted)
            if ( thisEEGlog.eventsAligned==0  &&  thisEEGlog.eventsExists==1  &&  thisEEGlog.hasSync==1 )
                
                behSync_file = thisEEGlog.eeglogUpStr ;     % ex) /Users/wittigj/DataJW/data/eeg/NIH016/behavioral/playPass/session_7/eeg.eeglog.up
                eegSync_file = thisEEGlog.syncFilePath ;    % ex) /Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/NIH016_280613_1411.083.084.sync.txt
                [~, eegChan_file] = findAnyExtractions(subj, rootEEGdir, thisEEGlog.extractRootName) ;  % ex) /Users/wittigj/DataJW/data/eeg/NIH016/eeg.noreref/NIH016_280613_1411.TG1
                events_file  = thisEEGlog.eventsFile ;      % ex) /Users/wittigj/DataJW/data/eeg/NIH016/behavioral/playPass/session_7/events.mat
                
                %keyboard
                %- pass the session log and pulse log so alignment adds info to those too (helpful if necessary to split session across two raw files
                singleLog_file = {events_file};
                multiLog_files = {thisEEGlog.eventsFile thisEEGlog.sessionLogStr thisEEGlog.eeglogStr};
                cellOfLogs     = multiLog_files ; % singleLog_file  or  multiLog_files
                
                % get the samplerate
                sampleRate = GetRateAndFormat(eegChan_file);
                if isempty(sampleRate), error('ERROR: GetRateAndFormat did not return a sampleRate.  Is params.txt present?'); end
                fprintf(     '\n>>>>>>>>>> alignment of [%s] and [%s] <<<<<<<<<<\n', fullfileEEG(thisEEGlog.task, thisEEGlog.sess), eegSync_file(strfind(eegSync_file,'eeg.noreref')+12:end));
                fprintf(fid, '\n-------------------------------------\n...alignment of [%s] and [%s]...\n', fullfileEEG(thisEEGlog.task, thisEEGlog.sess), eegSync_file(strfind(eegSync_file,'eeg.noreref')+12:end));
                
                
                
                % run the alignment
                alignmentWarning = 0;
                try
                    moreAccurateAlign = 0 ; %much slower if set to 1
                    alignInfo = runAlign(sampleRate, {behSync_file}, {eegSync_file}, {eegChan_file}, cellOfLogs, 'mstime', 0, 0, moreAccurateAlign);
                    fprintf(fid, '%s\n%s\n%s\n', alignInfo.strPulseAlign, alignInfo.strStats, alignInfo.strWarning);
                    %if ~isempty(strfind(alignInfo.strWarning,'WARNING')), alignmentWarning = 1; end  % used to pause for any warning, now just the "important" ones
                    iWarnDontStop = strfind(alignInfo.strWarning,'WARNING: pulse window reduced');
                    iWarnAll      = strfind(alignInfo.strWarning,'WARNING');
                    if length(iWarnAll)>length(iWarnDontStop), alignmentWarning = 1; end % only pause for user input if NOT the pulse window reduced warning
                catch err
                    report = getReport(err, 'basic','hyperlinks','on');
                    fprintf(     'runAlign threw an error: %s\n', report);
                    fprintf(fid, 'runAlign threw an error: %s\n', report);
                    alignmentWarning = 1;
                end
                %                 moreAccurateAlign = 0 ; %much slower if set to 1
                %                 alignInfo = runAlign(sampleRate, {behSync_file}, {eegSync_file}, {eegChan_file}, cellOfLogs, 'mstime', 0, 0, moreAccurateAlign);
                %                 fprintf(fid, '%s\n%s\n%s\n', alignInfo.strPulseAlign, alignInfo.strStats, alignInfo.strWarning);
                %                 if ~isempty(strfind(alignInfo.strWarning,'WARNING')), alignmentWarning = 1; end
                % query the user so each alignment is confirmed if a warning came up
                if alignmentWarning==1 && BYPASS_ALIGNMENT_USER_QUERY==0 && ~BATCH
                    reply = input(sprintf('\nConfirm that alignment results look OK: \n Is Max. Dev. < 5 ms?   Is R^2 > 0.98?   Are all events aligned?\n If not, type "N" to break so you can take a closer look at %s/%s: [Y]', thisEEGlog.task, thisEEGlog.sess) ,'s');
                    if (isempty(reply)) reply = 'Y'; end
                    if (reply(1)=='N' || reply(1)=='n')
                        fclose(fid);
                        error('--force break from eegPrepAndAlign.m so failed alignment can be examined--');
                    end
                end
                fprintf('\n');
                
                % confirm events file modified
                events = [];
                load(events_file);%
                if (isfield(events,'eegfile') && isfield(events,'eegoffset'))
                    alignStrMod = '';%
                    if FORCE_ALIGNMENT_TO_SERVER_PATH
                        for iEv=1:length(events)
                            events(iEv).eegfile = regexprep(events(iEv).eegfile, fullfileEEG(rootEEGdir,subj,''), fullfileEEG(serverDataPath,subj,'')); %
                        end
                        alignStrMod = ' [to server]';
                        saveEvents(events,events_file);
                    end
                    allEEGlogs(iList).eventsAligned = 1;
                    allEEGlogs(iList).eventsEEGfile = events(1).eegfile;
                    allEEGlogs(iList).eventsStr     = sprintf('events.mat aligned%s',alignStrMod);
                    alignStr                        = sprintf('%s  -%s/%s aligned%s\n', alignStr, thisEEGlog.task, thisEEGlog.sess,alignStrMod);
                    iJustAligned = [iJustAligned iList];
                else
                    allEEGlogs(iList).eventsStr     = sprintf('events.mat FAILED ALIGNMENT');
                    alignStr                        = sprintf('%s  -%s/%s alignment FAILED!!!\n', alignStr, thisEEGlog.task, thisEEGlog.sess);
                    %error('ERROR: Event not aligned... not sure why:\n\n%s',alignStr);  %This shouldn't happen, unless perhaps above Catch found an error too
                end
            end
        end
        fprintf(alignStr);
        fclose(fid);
    end
    fprintf('\n');
end


%%- CHECK MASTER EVENTS.MAT: if not aligned, aggregate aligned events structures together
if (FORCE_MASTER_EVENTS_UPDATE || ~isempty(iJustAligned))
    fprintf('\nCONCATENATED EVENTS: confirming master events.mat exists and is aligned\n');
    
    if (FORCE_MASTER_EVENTS_UPDATE) taskList = unique([allEEGlogs.taskNum]);
    else                            taskList = unique([allEEGlogs(iJustAligned).taskNum]); end
    for iTask = taskList
        
        % find all events from a single task
        numSess = length( find([allEEGlogs.taskNum] == iTask) ) ;
        iListAr = find( [allEEGlogs.taskNum] == iTask  &  [allEEGlogs.eventsAligned] ) ;
        if (length(iListAr)>0)
            % aggregate aligned events together
            allEvents    = [];
            for iList = iListAr,
                eventFile = allEEGlogs(iList).eventsFile;
                events = [];
                load(eventFile);
                allEvents=[allEvents, events];
            end
            
            % even if master events exists overwrite because this alignment could be an improvement
            eventMaster = fullfileEEG(rootEEGdir,subj,'behavioral',allEEGlogs(iListAr(1)).task,'events.mat');
            events = allEvents;
            fid = fopen(eventMaster,'w');
            if fid==-1, error('ERROR: cannot save master events file %s', eventMaster); else fclose(fid); end
            save(eventMaster, 'events', '-v7'); %- version 7 file format is compact... only need version 7.3 if >2GB
            
            % output result of check to terminal... everything OK?
            if (length(iListAr)<numSess) strMiss='<< MISSING SESSION'; else strMiss=''; end
            padStr = sprintf('%s/events.mat', allEEGlogs(iListAr(1)).task); padStr(end+1:25)=' ';
            fprintf('  aligned: %s (%d events; from %d of %d sessions) %s\n', padStr, length(events), length(iListAr), numSess, strMiss);
            
        end
    end
    
    %-- changes eegfile in events to point at the reref directory [this dones't really make sense here... should be done later so files point to server]
    fprintf('\n');
    
end


if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open', length(fIDs)-numFIDS);
        keyboard;
        numFIDS = length(fIDs);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- Text output pretty list for Kareem to use for grabbing raw files
%%--
%%--
fprintf('\n\n[SUMMARY] Prep and Align Results:\n');

fid = fopen(sprintf('%s/alignmentSummary.txt',behDir),'w+');
fprintf(fid,'List of EEG times (grabbed from all tasks, all times):\n');


%%- sort all the date/times for easier-to-read output
fprintf(fid,'\n\nSorted by task/session:\n');
%fprintf(    '\n\nSorted by task/session:\n');
for iList=1:length(allEEGlogs),
    thisEEGlog = allEEGlogs(iList);
    padTaskStr = sprintf('%s/%s',thisEEGlog.task, thisEEGlog.sess);
    padTaskStr(end+1:21) = ' ';
    
    % output to file and command line
    fprintf(fid, '  %s \t %s   -%s\n', padTaskStr, thisEEGlog.dateStr_act, thisEEGlog.eventsStr);
    %fprintf(     '  %s \t %s   -%s\n', padTaskStr, thisEEGlog.dateStr_act, thisEEGlog.eventsStr);
end


%%- same as above, but sorted by time
fprintf(fid,'\n\nSorted by time (and grouped by EEG session):\n');
fprintf(    '\nSORTED BEHAVIORAL DATA (sorted by date/time and grouped by EEG session):\n');
[~,iLogs] = sort([allEEGlogs.dateMAT_act]);
for iList=iLogs,
    thisEEGlog = allEEGlogs(iList);
    strAdd     = thisEEGlog.rawStr;
    padTaskStr = fullfileEEG(thisEEGlog.task, thisEEGlog.sess);   padTaskStr(end+1:25)=' ';
    
    % output to file and command line
    fprintf(fid, '%s  %s  %s \t<<-- %s\n', strAdd, thisEEGlog.dateStr_act, padTaskStr, thisEEGlog.eventsStr);
    fprintf(     '%s  %s  %s \t<<-- %s\n', strAdd, thisEEGlog.dateStr_act, padTaskStr, thisEEGlog.eventsStr);
end


%%- indicate any RAW that exists but ISN'T USED
fprintf(fid, '\n---------------------------------------------------------------------------------\n');
fprintf(     '\n---------------------------------------------------------------------------------\n');
iListExtAr = [];  iListRawAr = [];  numUnpairedRaw = 0;
if MISSING_RAW_AND_EXTRACTED==0,
    if MISSING_ALL_RAW==1, iListExtAr = find([extractedList.pairedToAlign]==0);
    else        [~,keepNSP1,~] = unique({rawList.extractRootName});
        rawList = rawList(keepNSP1);
        iListRawAr = find([rawList.pairedToAlign]==0);  end
    if length(iListExtAr)+length(iListRawAr)==0,
        fprintf(fid, '>>All %d raw and/or extracted channel files paired with at least 1 eeg.eeglog', length(rawList)+length(extractedList));
        fprintf(     '>>All %d raw and/or extracted channel files paired with at least 1 eeg.eeglog', length(rawList)+length(extractedList));
    else
        fprintf(fid, '>>WARNING: %d of %d raw and/or extracted channel files NOT paired with any eeg.eeglog:\n', length(iListExtAr)+length(iListRawAr), length(rawList)+length(extractedList));
        fprintf(     '>>WARNING: %d of %d raw and/or extracted channel files NOT paired with any eeg.eeglog:\n', length(iListExtAr)+length(iListRawAr), length(rawList)+length(extractedList));
        for iList=iListExtAr,
            fprintf(fid, '   %s  %s\n', datestr(extractedList(iList).dateN,'mm/dd/yy HH:MM PM'),    extractedList(iList).rootName);
            fprintf(     '   %s  %s\n', datestr(extractedList(iList).dateN,'mm/dd/yy HH:MM PM'),    extractedList(iList).rootName);
        end
        for iList=iListRawAr,
            numUnpairedRaw = numUnpairedRaw+1;
            fprintf(fid, '   %s  %s\n', datestr(rawList(iList).extractDateNum,'mm/dd/yy HH:MM PM'), rawList(iList).clnPath);
            fprintf(     '   %s  %s\n', datestr(rawList(iList).extractDateNum,'mm/dd/yy HH:MM PM'), rawList(iList).clnPath);
        end
    end
end


%%- a list of the raw dates (and associated extractions), sorted by time
fprintf(fid,'\n\nRAW data already saved to server:\n');
%fprintf(    '\n\nRAW data already saved to server:\n');
for iList=1:length(rawList),
    outStr = sprintf('  %s  %s', datestr(rawList(iList).extractDateNum,'mm/dd/yy HH:MM PM'), rawList(iList).clnPath);
    outStr(end+1:46) = ' '; % pad so text aligns to right
    
    if ispc, outStr(find(outStr=='\'))='/'; end % avoid Warning with PC version of matlab that "Escape sequence 'B' is not valid
    
    fprintf(fid, outStr);   % output to eegTimes textfile
    %fprintf(     outStr);   % also output to command line
    
    %also output associated sync file names (or instruct how to create)
    if (rawList(iList).extracted)
        if (rawList(iList).hasSync), fullSyncStr = sprintf(' [%s: %s]', rawList(iList).syncStr, rawList(iList).syncFile(max(find(rawList(iList).syncFile=='/' | rawList(iList).syncFile=='\'))+1:end)) ;
        else,                        fullSyncStr = sprintf(' [%s!] \n   -->> generate sync with alignTool on %s, %s <<--', rawList(iList).syncStr, rawList(iList).extractRootName, rawList(iList).jackChanStr);  end
        
        fprintf(fid, '     -%s \n', fullSyncStr);
        %fprintf(     '     -%s \n', fullSyncStr);
    else
        error('Shouldnt be possible to get here... raw always extracted above...')
    end
end


%%- Summary Counts
numTaskSess = length(       allEEGlogs                                                   );
numAligned  = length( find([allEEGlogs.eventsAligned]==1                               ) );
numNeedEvnt = length( find([allEEGlogs.eventsAligned]==0 & [allEEGlogs.eventsExists]==0) );
numNeedRaw  = length( find([allEEGlogs.eventsAligned]==0 & [allEEGlogs.hasRaw]==0      ) );
numNeedSync = length( find([allEEGlogs.eventsAligned]==0 & [allEEGlogs.hasSync]==0     ) );
numCldAlign = length (find([allEEGlogs.eventsAligned]==0 & [allEEGlogs.eventsExists]==1 & [allEEGlogs.hasSync]==1) );
if numAlignedWithBlank>0, strAlignedWithBlank = sprintf(' [%d events.mat with at least 1 blank eegfile field]', numAlignedWithBlank); else strAlignedWithBlank = ''; end
fprintf(fid, '\n\nSUMMARY: %02d of %02d task-sessions aligned %s \n         %02d missing events.mat \n         %02d missing eeg RAW data \n         %02d missing eeg sync\n         %02d alignment possible \n         %02d unpaired eeg RAW data', numAligned, numTaskSess, strAlignedWithBlank, numNeedEvnt, numNeedRaw, numNeedSync, numCldAlign, numUnpairedRaw);
fprintf(     '\n\nSUMMARY: %02d of %02d task-sessions aligned %s \n         %02d missing events.mat \n         %02d missing eeg RAW data \n         %02d missing eeg sync\n         %02d alignment possible \n         %02d unpaired eeg RAW data', numAligned, numTaskSess, strAlignedWithBlank, numNeedEvnt, numNeedRaw, numNeedSync, numCldAlign, numUnpairedRaw);


fclose(fid);
fprintf('\n\nTimes saved to %s \n', fullfileEEG(behDir,'alignmentSummary.txt'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%  Also create text output that matches format of README.txt
MAKE_README_HELPER = 0;  %- this was useful for NIH001-NIH021. The output has been edited as saved as "docs/testing_notes.txt" for those subjects.
%                                      Starting at NIH022 we used a lab clipboard for testing notes, which is scanned and saved in docs
if MAKE_README_HELPER,
    fileHelper = sprintf('%s/readme_helper.txt',behDir);
    fid = fopen(fileHelper,'w+');
    if fid==-1, fprintf('ERROR: could not create and/or open %s.  Confirm that you have write access and rerun eegPrepAndAlign', fileHelper); error('ERROR: cant create/write file'); end
    fprintf(fid,'\nPREP and ALIGN OUTPUT (sorted by date/time): to easily import to subjects README.txt\n');
    [~,iLogs] = sort([allEEGlogs.dateMAT_act]);
    for iList=iLogs,
        thisEEGlog = allEEGlogs(iList);
        
        syncChan = '??';  rawFilePath = '???';
        if thisEEGlog.hasRaw,
            thisRawFile = rawList(thisEEGlog.rawListIndex);
            syncChan    = thisRawFile.jackChanStr;
            rawFilePath = sprintf('%s --> %s', fullfileEEG(subj,'raw',thisRawFile.clnPath), thisRawFile.extractRootName);
        end
        syncFile = '??';
        if thisEEGlog.hasSync
            [~, syncFile, syncExt] = fileparts(thisEEGlog.syncFilePath);
            syncFile = sprintf(' [%s%s]', syncFile, syncExt);
        end
        
        behavDataDir = fullfileEEG('/eeg', subj, 'behavioral', thisEEGlog.taskSess);
        
        fprintf(fid, '\n-------------------------------------------------------------------\n');
        fprintf(fid, '%s\n',                    datestr(thisEEGlog.dateMAT_act, 'mm/dd/yyyy'));  %
        fprintf(fid, 'Task: %s\n',              thisEEGlog.task );
        fprintf(fid, 'Session: %s\n',           thisEEGlog.sess );
        fprintf(fid, 'Start: %s\n',             datestr(thisEEGlog.dateMAT_act, 'HH:MM PM') );
        fprintf(fid, 'End: %s\n',               datestr(thisEEGlog.dateMAT_actEnd, 'HH:MM PM') );
        fprintf(fid, 'Sync: %s %s\n',           syncChan, syncFile );
        fprintf(fid, 'Behavioral data: %s\n',   behavDataDir );
        fprintf(fid, 'Clinical EEG data: %s\n', rawFilePath );
        fprintf(fid, 'Testers: ??? \n');
        fprintf(fid, 'Notes: ??? \n');
        
    end
    fclose(fid);
    fprintf('Times also saved to %s \n', fullfileEEG(behDir,'README_helper.txt'));
    fprintf('******************************** %s *****************************\n\n', subj);
end


%%%%% Prompt user to create pulse align figure if all raws are accounted for but pulse figure is missing
if numNeedRaw==0 && ~exist(fullfileEEG(rawDir,'align_PlotPulseChannels.png'),'file') && BYPASS_ALIGNMENT_USER_QUERY==0
    fprintf(' NOTE: raw directory is missing pulse screenshot ''align_PlotPulseChannels.png''\n');
    
    if strcmp(pulseTag,'DC'),  pulseChan = {'DC09' 'DC10'};
    else                       pulseChan = {sprintf('%s2',pulseTag{1}),sprintf('%s1',pulseTag{1})}; end
    
    fprintf('       guessing (based on last entry in tagNames that pulse channel(s) are: ');
    for iPC=1:length(pulseChan), fprintf('%s ',pulseChan{iPC}); end
    reply = 'Y';
    if ~BATCH,
        reply = input('\n       Run pulseVisualize now to make the screenshot using these pulse channel(s)? Y/N [N]:','s');
        if isempty(reply)
            reply(1) = 'N';
        end
    end
    if upper(reply(1))=='Y'
        try
            eegPulseVisualize(rawDir, pulseChan);
        catch
            fprintf('\n counldnt compelte pulse visualize');
        end
    end
    fprintf('\n');
end


%%%%% Make sure eventEegFilePaths.txt is created once everything is aligned.  If any events just aligned recreate
%if numAligned==numTaskSess &  FORCE_ALIGNMENT_TO_SERVER_PATH==0  &  length(iJustAligned)>0 || ~exist(fullfileEEG(behDir,'eventEegFilePaths.txt'),'file'),
if numAligned==numTaskSess &&  ~isempty(iJustAligned) || ~exist(fullfileEEG(behDir,'eventEegFilePaths.txt'),'file')
    if  ~exist(fullfileEEG(behDir,'eventEegFilePaths.txt'),'file')
        fprintf(' NOTE: behavioral directory is missing ''behavioral/eventEegFilePaths.txt''... will generate using changeEventsEegFile now:\n');
    else
        fprintf(' NOTE: one or more experiments just aligned; updating ''behavioral/eventEegFilePaths.txt'' using changeEventsEegFile now:\n');
    end
    
    [uniqueEegFilePaths] = changeEventsEegFile(subj, rootEEGdir, '', '');   % if eventEegFilePaths.txt is not created make it now
    if length(uniqueEegFilePaths)>1,
        fprintf(' NOTE: multiple file paths found in events.mat:  \n       FIX any BLANKS [empty fields] by modifying non-aligned mstimes in a session.log.  \n       Re-extract, align, and/or re-map (using eegPrepAndAlign and/or changeEventsEegFile) any experiments with innaccurate filepaths\n ');
    elseif length(uniqueEegFilePaths)==1,
        fprintf(' LOOKING GOOD: single file path ''%s'' can be passed as "oldphrase" \n   to changeEventsEegFile(subj, rootEEGdir, oldphrase, newphrase) if/when events rereferencing is required\n', uniqueEegFilePaths{1});
    end
    fprintf('\n');
end



%%%% Check for ictal, interictal, and resected figures.  
% first confirm elementInfo filled out properly, then only make the figures if they are missing
jack  = getJackTable(subj, rootEEGdir); % cache for speed
ictal = getLeads(subj, rootEEGdir, 'markedAs', 'ictal',      'jackTable',jack);
inter = getLeads(subj, rootEEGdir, 'markedAs', 'interictal', 'jackTable',jack);
resec = getLeads(subj, rootEEGdir, 'markedAs', 'resected',   'jackTable',jack);
mask_sh = contains(jack.chanName, '_sh');
isShiftSubj = any(mask_sh);
        
electrodeFigs = {fullfile(rootEEGdir, subj, 'docs', 'electrodes_resected.png') fullfile(rootEEGdir, subj, 'docs', 'electrodes_ictal.png') fullfile(rootEEGdir, subj, 'docs', 'electrodes_interictal.png')};


%- only make the figs if at least one column is filled out... otherwise just hasn't been done yet
if isempty({ictal{:} inter{:} resec{:}})
    fprintf('\n\nWARNING: no channels marked as ictal, interictal, or resected in %s elementInfo. \n  Probably means Kareem hasnt entered them yet. \n  -Confirm that element_info is up to date\n  -Compare with patient_info and possibly ask kareem to update: ', subj);
    
    if ~BATCH,
        elFilename    = fullfile(rootEEGdir, subj, 'docs', 'element_info.csv');
        system(sprintf('open %s', elFilename));
        
        %- usually missing input from kareem, so probably cant update "now", but make it easy if so
        reply = input('\n\n >>>> Electrode Figures will NOT be generated until this is resolved. \n >>>> Manually UPDATE element_info.csv NOW, then break and rerun prepAndAlign to generate electrode figures \n Press [RETURN] to continue WITHOUT making electrode figures (or ENTER ''q'' to keyboard) \n','s');
        if ~isempty(reply) && (reply=='Q' || reply=='q')  fprintf('\n keyboard here so you can see what up... ');  keyboard;  end
    end  %-batch
else
    
 
	hFig = makeIctalResectedFig(subj, rootEEGdir); % makes a brain figure with ictal, interictal, and resected electrodes labeled   MT and JW added 4/2018
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- Graphical version of list: histogram of session times
%%--
if disp_flag,
    figure(1001); clf;
    
    %- sort the data
    allDates     = [allEEGlogs.dateMAT_act];
    deltaTimeHr  = datenum('01/01/01 2:00')-datenum('01/01/01 1:00');
    deltaTimeBin = deltaTimeHr*3;
    dBins  = [min(allDates)-deltaTimeBin/2:deltaTimeBin:max(allDates)+deltaTimeBin/2]; % create 2-hour bins
    taskHist=[];
    for iTask=1:length(taskDir)
        taskDates = [allEEGlogs(find([allEEGlogs.taskNum]==iTask)).dateMAT_act];
        [oc, dBins] = hist(taskDates,dBins);
        taskHist(iTask,:)=oc;
    end
    
    %-histogram of task sessions/day
    hB = barh(dBins,taskHist', 'stacked'); hold on
    datetick('y','mm/dd');
    for iTask=1:length(taskDir), set(hB(iTask),'BarWidth',1.0,'Edgecolor','k'); end
    xMax = max(get(gca,'xlim'));
    set(gca,'xlim',[0 xMax+2])
    title(sprintf('Testing Dates: .%s. ',subj),'fontsize',20)
    xlabel('Sessions Collected','fontsize',18)
    set(gca,'fontsize',14,'YDir','reverse')
    
    
    %-add markers of EEG recording breaks
    iCount = 1;
    for iList=iLogs
        thisEEGlog = allEEGlogs(iList);
        if (thisEEGlog.isNewSession)
            %hLine = plot([0 xMax+.2],thisEEGlog.dateMAT_act*[1 1],'r--');
            %hText = text(xMax+.45,thisEEGlog.dateMAT_act,sprintf('Test Session %d',iCount),'fontsize',15);
            
            hPt = plot(xMax+0.3,thisEEGlog.dateMAT_act,'ko'); set(hPt,'markersize',14,'MarkerFaceColor','r')
            iCount = iCount+1;
        end
    end
    
    
    %-add markers of raw start times
    if (MISSING_ALL_RAW==1 && MISSING_ALL_SYNC==1)
        lText = taskStrAr;
        lText{end+1} = 'Raw MISSING';
        legend([hB hPt],lText,'Location','Best')
    else
        hRaw = plot((xMax+0.3)*ones(size(rawTimes)),rawTimes,'ko');
        set(hRaw,'markersize',14,'MarkerFaceColor','k')
        
        if ~isempty(syncTimes)
            hSync = plot((xMax+0.3)*ones(size(syncTimes)),syncTimes,'ko');
            legSync = 'Raw Synced';
            set(hSync,'markersize',6,'LineWidth',0.5,'MarkerFaceColor','c','MarkerEdgeColor','k')
        else
            hSync = plot(max(get(gca,'xlim'))*2, 0,'kx');  %plot outside range of data
            legSync = 'No Sync Found';
            set(hSync,'markersize',12,'LineWidth',1,'MarkerFaceColor','c','MarkerEdgeColor','k')
        end
        
        lText = taskStrAr;
        lText{end+1} = 'Raw MISSING';
        lText{end+1} = 'Raw found';
        lText{end+1} = legSync ;
        legend([hB hPt hRaw hSync],lText,'Location','Best')
    end
end %- if disp_flag


if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open', length(fIDs)-numFIDS);
        keyboard;
        numFIDS = length(fIDs);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n[STEP 5 / 5] Rereferencing and looking for bad channels\n');
fprintf('\t- All channels will be run through a detrending, desaturation, and line-noise removal algorithm\n');
fprintf('\t- All Channels will be bipolar-rereference\n');
fprintf('\t- Subdural channels will assed for high variance and high amplitude, and outliers will be marked "bad"\n');
fprintf('\t- Subdural channels (both all and good) will be averaged into a global signal average for easy rereferencing\n\n');

leads = getLeads(subj, rootEEGdir,'jackTable',jack);
%%- ATTEMPT TO RE-REFERENCE IF leads exist

if isempty(leads)
    fprintf('  WARNING: rereferencing NOT possible because no leads found.\n');
    empty_leads = 1;
else
    empty_leads = 0;
    
    %- do a quick check to see that all eeg.noreref splits are in eeg.reref and eeg.processed
    noRerefDir = fullfileEEG(subjDir,'eeg.noreref'); % where no-reref data is located
    procDir    = fullfileEEG(subjDir,'eeg.processed');
    %rerefDir   = fullfileEEG(subjDir,'eeg.reref'); % where reref data will be saved to
    rerefDir   = fullfileEEG(subjDir,'eeg.processedBP'); % where reref data will be saved to
    
    needProcessAndReref = 0;  %- default -- assume we dont need to process or reref
    if ~exist(rerefDir,'dir') || ~exist(procDir,'dir')
        needProcessAndReref = 1;
        fprintf(' missing eeg.processed or eeg.processedBP folder... need to reref and clean!');
        
    else
        
        % --- offer user option to delete folders in eeg.* that are not in raw --- %
        
        % get names/paths
        rawDirPaths = findRaws(rawDir);
        rawDirPaths = cellfun(@fileparts,rawDirPaths,'uniformOutput',0);    % path ending with session folder
        [~,rawDirNames] = cellfun(@fileparts,rawDirPaths,'uniformOutput',0);% just session folder
        noreref_dates = lsCell(noRerefDir);
        noreref_dates = noreref_dates(cellfun(@isdir, fullfile(noRerefDir,noreref_dates)));
        reref_dates = lsCell(noRerefDir);
        reref_dates = reref_dates(cellfun(@isdir, fullfile(noRerefDir,reref_dates)));
        proc_dates = lsCell(noRerefDir);
        proc_dates = proc_dates(cellfun(@isdir, fullfile(noRerefDir,proc_dates)));
        
        % candidates are those in eeg.* that are not in raw/
        candidates = unique([fullfile(noRerefDir, setdiff(noreref_dates, rawDirNames)) , ...
            fullfile(rerefDir, setdiff(noreref_dates, reref_dates)) , ...
            fullfile(procDir, setdiff(proc_dates, rawDirNames))]);
        
        if ~isempty(candidates)
            fprintf('\nThe following session folders may be old and may be no longer needed since they were');
            fprintf('found in an eeg.* folder but not in the raw/ folder:\n');
            disp(char(candidates));
            yes = BATCH || inputYN('For your convenience, do you want to delete all these folders now?');
            if yes
                for i = 1:numel(candidates)
                    try
                        rmdir(candidates{i},'s');
                    catch e
                        disp(e.message);
                    end
                end
            end
        end
        
        %- root folders exists... do the split folders match (not checking inside here)
        noreref_dates = lsCell(noRerefDir);
        noreref_dates = noreref_dates(cellfun(@isdir, fullfile(noRerefDir, noreref_dates)));
        reref_dates   = lsCell(rerefDir);
        reref_dates   = reref_dates(cellfun(@isdir, fullfile(noRerefDir, reref_dates)));
        proc_dates    = lsCell(procDir);
        proc_dates    = proc_dates(cellfun(@isdir, fullfile(noRerefDir, proc_dates)));
        
        numMissingReref = length(setdiff(noreref_dates,reref_dates));
        numMissingProc  = length(setdiff(noreref_dates,proc_dates));
        
        if numMissingReref+numMissingProc>0,
            needProcessAndReref = 1;
            fprintf(' missing %d of %d file stems in eeg.reref or eeg.processed... need to reref and clean!', max([numMissingReref numMissingProc]), length(noreref_dates) );
        else
            fprintf(' all %d file stems found in eeg.reref and eeg.processed... already rerefed and processed!', length(noreref_dates) );
        end
    end
    
end

%- only offer reref if necessary
if ~empty_leads && needProcessAndReref & ~SKIP_PROCESS_REREF,
    
    if ~BATCH
        reply = input('  >>>>>>   Reref and clean now? Y/N [N]:','s');
        if isempty(reply)
            reply = 'N';
        end
        fprintf('\n');
    else
        reply = 'Y';
    end
    
    if reply(1)=='Y' | reply(1)=='y',
        %%- call reref wrapper to create all rereferenced channels too. only creates files that are missing
        %       (assumes filestem presence means all rerefs have been created)
        %fprintf(' re-reference channel files if not done already:\n');
        [stemsChecked, stemsRerefed] = processAndReref(subj, rootEEGdir, 'disp', disp_flag);
        fprintf('\n %d of %d file stems were rereferenced', stemsRerefed, stemsChecked );
        
        fprintf('\n running subject directory check now...');
        subjectDirCheck(subj,rootEEGdir);
    end
end


%- check to see if atlas simple was made, if not, make it now
atlasDir = fullfileEEG(rootEEGdir,subj,'tal/atlas'); %- if tal/atlas directory doesn't exist, don't bother checking for simples!
if exist(atlasDir,'dir'),
    atlasSimpleMP = fullfileEEG(atlasDir,'atlas_monopolar_simple.csv');
    if ~exist(atlasSimpleMP,'file'),
        fprintf('\n\n Missing atlas_monopolar_simple.csv ... making it now:\n');
        atlas_simplify(subj, rootEEGdir, 'allow_ties', 0);
    end
    atlasSimpleBP = fullfileEEG(atlasDir,'atlas_bipolar_simple.csv');
    if ~exist(atlasSimpleBP,'file'),
        fprintf('\n\n Missing atlas_bipolar_simple.csv ... making it now:\n');
        atlas_simplify(subj,rootEEGdir,'allow_ties', 0, 'fname_in', fullfileEEG(atlasDir,'atlas_bipolar.csv'));
    end
end


%- create final check file if not already created
finalCheckFile = fullfileEEG(subjDir, 'docs/finalCheck.csv');

% if ~empty_leads && ~needProcessAndReref && ~exist(finalCheckFile,'file'),
% Bypassing requirement for proccessed and reref to be created - DY, JW
% 4/2019
    
if ~empty_leads && ~exist(finalCheckFile,'file'),
    subjectDirCheck(subj,rootEEGdir);
end

fprintf('\n[STEP 5 / 5] COMPLETE!\n');




if CHECK_FILECLOSE,
    fIDs = fopen('all');
    if length(fIDs)>numFIDS,
        fprintf('\n uh oh: additional %d files open', length(fIDs)-numFIDS);
        keyboard;
        numFIDS = length(fIDs);
    end
end

end %function eegPrepAndAlign




%% Private functions
function d21Es = getRaw21Es(rawDir)
d21Es = {};
roots = findRaws(rawDir, {});
for i = 1 : length(roots)
    thisRootPath = roots{i};
    [~,thisRootFile,~]=fileparts(thisRootPath); 
    if     contains(thisRootFile,'ieeg')     | contains(thisRootFile,'INST'),  d21E = [thisRootPath '.ns3'];
    elseif contains(thisRootFile,'cervello') | contains(thisRootFile,'EEG'),   d21E = [thisRootPath '.TRC'];
    else                                                                       d21E = [thisRootPath '.21E']; 
    end
    if exist(d21E, 'file')
        d21Es = [d21Es d21E];
    else
        d21E = [thisRootPath '.nev'];
        if exist(d21E, 'file')
            d21Es = [d21Es d21E];
        end
    end
end
d21Es = unique(d21Es);
end


function [isExtracted, firstExtractedFilename] = findAnyExtractions(subj, rootEEGdir, extractRootName)
% look for a SINGLE expected extraction file. If we find ONE, assume that
% extraction is complete
chans = getLeads(subj, rootEEGdir);
for i = 1 : length(chans)
    expectedFilename = fullfileEEG(rootEEGdir,subj,'eeg.noreref',extractRootName, chans{i});
    if exist(expectedFilename, 'file')
        isExtracted = true;
        firstExtractedFilename = expectedFilename;
        return;
    end
end
isExtracted = false;
end


function varargout = textread2(filename, format_str)
% quick update of out-dated textread to textscan
fd = fopen(filename, 'r');
c = textscan(fd, format_str);
fclose(fd);
varargout = deal(c{:});
if ~iscell(varargout), varargout = {varargout}; end
end


function [doRecreate,jackDateStr,mostRecentModStr] = isOldJacksheet(subj, rootEEGdir)
% true if any raw timestamp folder or elementInfo modified after jacksheet
dbg = 0;

rawdir = fullfile(rootEEGdir, subj, 'raw');
mostRecentMod = 0;

fname_info = fullfile(rootEEGdir, subj, 'docs/element_info.csv');
if exist(fname_info, 'file')
    infoModDate   = getFileModDate(fname_info);
    mostRecentMod = max(infoModDate ,mostRecentMod);
else
    doRecreate = 1;
    return;
end

fname_21Es = strcat(findRaws(rawdir), '.21E');
timestamps = cellfun(@fileparts, fname_21Es, 'uniformoutput',0);
for i=1:length(fname_21Es)
    modeDate21E   = getFileModDate(fname_21Es{i});
    mostRecentMod = max(modeDate21E, mostRecentMod);
    if dbg, fprintf('%s\n',fname_21Es{i}); end
end
for i=1:length(timestamps)
    modDateSessDir = getFileModDate(timestamps{i});
    mostRecentMod  = max(modDateSessDir, mostRecentMod);
    if dbg, fprintf('%s\n',timestamps{i}); end
end

fname_jacksheet  = fullfile(rootEEGdir, subj, 'docs/jacksheetMaster.csv');
jackDate         = getFileModDate(fname_jacksheet);

% modified after jacksheet?
doRecreate       = mostRecentMod > jackDate;

jackDateStr      = datestr(jackDate);
mostRecentModStr = datestr(mostRecentMod);

if dbg
    fprintf('\nmost recent: %s\n', mostRecentModStr);
    fprintf('jacksheet  : %s\n', jackDateStr);
    fprintf('doRecreate : %d\n', doRecreate);
end

end %- isOldJacksheet


function date = getFileModDate(filename)
% creation time of file in matlab datenum
[~,result] = unix(sprintf('stat -f %%m "%s"', filename)); %maybe switch to %%c ? ; SJ 9/19/19 - change %s to "%s"
result = strsplit(result);
date = epoch2date(1000 * str2double(result{1}));
end


%#ok<*ISMT>
%#ok<*DTXTRD>
%#ok<*NOCOL>
%#ok<*NASGU> % unused var
%#ok<*SEPEX>
%#ok<*UNRCH>
%#ok<*SEPEX>
%#ok<*AGROW> % preallocation
%#ok<*SEPEX> % if comma
%#ok<*CTPCT> % unreachable code