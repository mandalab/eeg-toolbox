function extractMicrophysiology(subj,overrideOutputType,overrideDataPath56)
%subj         = 'NIH069';
%compMemoryGB = 64;
%
% optional: overrideOutputType = {'minimum' 'standard' 'all'};  %- defaults to 'standard' and 'stim' for "isStim" sessions if no input given
% optional: overrideDataPath56
% overrideDataPath56:  if data resides locally, instead of on FRNU56, put the path to the subjects here
%                     e.g. '/Volumes/JW24TB/data24TB/localFRNU56/';
%     alternatively, if you want to work on FRNU72, pass "72" for the override path
%

%% OLD NOTES
% adapted from extractUtahEphys
%
% Using information from the utahInfo CSV files, extract each utah
% channel's ephys time series for LFP(1-500Hz) and spikes(300-3000Hz)
%
% Inputs:
% subj: 'NIH037'
% outputTypes: structure with the following fields:
%   sort_unfilt, sort_reref, sort_noreref, lfp_reref, lfp_noreref
%   Set field value to 1 if you want to extract, 0 otherwise.
% compMemory: computer memory in gigabytes
%

%% NEW NOTES
%  JW's take on how this function should work in the plexon extraction pipeline (circa 2am on 2/4/2019)
%     1) run updateUtahInfo to generate list of potential sessions for plexon extraction
%     2)    --> manually select those sessions marking "y" in that file, then
%     3) re-run updateUtahInfo to generate list of channels to extract, matching jacksheets to getMicroSubjInfo list of output devices
%           --> this will also create subfolder notes/pruneSpikeChansWithNplay/ with a separate csv for each session that can be marked based on nPlay
%     4) extractMicrophysiology - split out ALL channels for LFP (to .mat) and SPIKES (to set of .bins), doing a global-reref of the spikes
%     3b or 4b) manually input "spikeChan2Prune" into a csv made for each session. (this can be done BEFORE or AFTER extractMicroPhys)
%     5) run updateUtahInfo to automatically move pruned files to data_processed/spikes_g/_PRUNEDnoSpikesOnNPlay/*.bin
%     6) sort the remaining .bin
%     7) run "extractSpikeInfo"
%
%  So what extractMicroPhys should do
%     1) load sess2extract and chan2extract csv files to get (1) sessions to extract and (2) input/output of old vs new channel names to look for
%     2) loop over each session
%              loop over each device in chan2extract
%                  load in all the channels from that device to compute global ave at 30kS, save that to disk
%              loop over each channel one-by-one
%                  write line to processed/session/jacksheetMicro.txt: for each channel   oldChanName -> newChanName  ; device num; which nsX file
%                  low-pass and downsample for LFP
%                  subtract global for that channels device, bandpass,
%                  save spike chan to disk as targetChanName.bin, placing in _PRUNED if pruned channels already specified
%              grab the pulse time series from any NSP that contributed to micro
%              save all the LFPs together as a single noreref.mat
%              save another mat with with chanNames + info from jacksheetMirco +  pulse time series
%
%
%
% Written by: Anthony I. Jang
% Edited by Mostafa El-Kalliny, November 2018
% reworked by JW 1/2019
% 11/19/2019 - SJ: Fixed if statement for subjects with empty utahMap in subjectInfo
% 11/20/2019 - SJ: changed chanNotStim to chanIsStim = []; to fix error from not initializing
% 2/24/2020  - SJ: fixed path 'spike' to 'sort'
% 2/2021     - SJ: changed dummySort to chanNoise & globalNoise


fprintf('\n\n ------------------ Starting extractMicrophys ::: %s ---------------------\n', subj);

if nargin<3
    overrideDataPath56 = '';
end
if nargin<2
    overrideOutputType = '';
end

PUBpath = '/Volumes/56PUB/readWrite/micro_forSorting';


%- this was a free parameter before... now we know you need a computer with >=64GB ram to run this.. no more free parameter
[memUsedGB, memFreeGB, memTotalGB] = getMemoryUsage();
compMemoryGB = memTotalGB;
if memTotalGB<64, fprintf('\n heads up, system RAM is probably too low to run this script. should have >=64GB RAM'); keyboard; end



%- grab the subj info here because the spike sort helper CSV will indicate which channels were microStim channels (they have weird noise so its helpful to know when sorting)
subjectInfo = getMicroSubjInfo_v11(subj,overrideDataPath56);
if ~isempty(subjectInfo.utahMap) && ~any(isnan(subjectInfo.utahMap.stimPigtail2Chan(:))) % SJ - changed first from isfield(subjectInfo,'utahMap')
    microStimChan = subjectInfo.utahMap.stimPigtail2Chan(:,2); %- electrode channels. location in this vector corresponds to pigtail number
else
    microStimChan = nan;
end

%- use this flag to force an update of the sort helper files EVEN if everything else is done
PUSH_NEW_SORT_HELPER = 1;  % default (0): set to 1 if new column, etc should be added to all spike folders
PUSH_SPK_EXTFIX = 0; %Both must be set to 1 for this to happen


% standard parameter used for stimMask (blanking out stimulation for stim sessions)
maskStimMsGrowDC10pulse = 5;  %- add 5ms of blank period before and after each DC10 or DC11 pulse


% Get Filter parameters
samprate_br  = 30000;
samprate_LFP = 1000;
samprate_nk  = 1000;


%- convert .bin and lfps to uV, or keep as 4x uV.   Better to keep as 4x so data isn't lost when converting back to int16 for disk
%%%%%%%% NOTE: IF YOU CHANGE CONVERT_2_UV, YOU NEED TO CHANGE IT FROM ZERO
%%%%%%%% IN fixSPKextractionInfo
CONVERT_2_uV = 0;    %- do NOT apply gain when calling concatOpenNSx, thus timeSeries is all 4x bigger than uV
gain_bin2uV = 0.25;  %- blackrock eCog and utah saved at 4x uV, so multiply by 0.25 to convert to uV.  GainBin2uB fixed because we ALWAYS write to disk as 4x uV

butterOrder         = 2;

% SJ Check
if CONVERT_2_uV ~= 0
    fprintf('%s\n', 'YOU CHANGED CONVERT_2_UV FROM ZERO. YOU ALSO NEED TO CHANGE IT FROM ZERO WHEN IT IS CALLED IN fixSPKextractionInfo: NSxData = concatOpenNSx(NSxFileDir,0,1,chanGroups{grp}); ((ZERO HERE))');
    keyboard;
end

% LFP
nyquist_freq_LFP    = samprate_br/2;
lfp_highBound       = 500;    % filter to 500Hz
wo_LFP              = lfp_highBound/nyquist_freq_LFP; % lowpass to Nyquist*(1-.00f1)
[b_LFP, a_LFP]      = butter(butterOrder,wo_LFP,'low');
str_LFP             = sprintf('LFP Band: butterworth low pass at %d Hz; filtfilt for zero phase shift (effective order %d)',lfp_highBound,butterOrder*2); %- filtfilt doubles order

% Spikes
samprate_spikes     = samprate_br;
nyquist_freq_spikes = (samprate_spikes/2);
spikeBand           = [300 3000];
wo_spike=spikeBand/nyquist_freq_spikes; % bandpass 300Hz to 3000HZ %mountainsort: 600 to 6000 Hz
[b_spike,a_spike]   = butter(butterOrder,wo_spike);
str_spike           = sprintf('Spike Band: butterworth bandpass %d-%d Hz; filtfilt for zero phase shift (effective order %d)',spikeBand,butterOrder*2); %- filtfilt doubles order

% Spikes - unfiltered
unfiltBand          = [1,5000];
wo_spike_unfilt     = unfiltBand/nyquist_freq_spikes; % "unfiltered": 1Hz - 5000Hz
[b_spike_unfilt,a_spike_unfilt] = butter(butterOrder,wo_spike_unfilt);
str_spike_unfilt    = sprintf('"Unfilt" Spike: butterworth bandpass %d-%d Hz; filtfilt for zero phase shift (effective order %d)',unfiltBand,butterOrder*2); %- filtfilt doubles order


% Get subject info
subjectInfo           = getMicroSubjInfo_v11(subj,overrideDataPath56);


% Directories for the critical input csv files
fileDir_bookkeeping   = fullfile(subjectInfo.dataPath56 ,subj,'_extraction_notes','micro_PickSess2Extract.xlsx');
fileDir_spikeChans    = fullfile(subjectInfo.dataPath56 ,subj,'_extraction_notes','micro_RenamedChanInfo.csv');

if ~exist( fullfile(subjectInfo.dataPath56 ,subj), 'dir'),
    fprintf('\n\n In extract microphys: Subject directory doesnt exist where specified so skipping subject\n%s', fullfile(subjectInfo.dataPath56 ,subj));
    return;
end
if ~exist(fileDir_bookkeeping,'file') | ~exist(fileDir_spikeChans,'file'),
    fprintf('\n\n In extract microphys: expected bookkeeping files are missing... need to run updateUtahInfoCSV on this subject first. Skipping subject');
    %keyboard;
    return;
end

%- output feedback to the screen or to a text file
fidOut = 1;  %- just output to the screen, good for debugging, but doesn't help much if matlab totally crashes due to insufficient memory
extractFile = fullfile(subjectInfo.dataPath56 ,subj,'_extraction_notes','_extractMicrophysLog.txt');
fidOut = fopen(extractFile ,'a+');
if fidOut==1, fprintf('\n WRITING LOG to SCREEN'); else fprintf('\n WRITING LOG TO DISK  %s',extractFile); end

%
[memUsedGB, memFreeGB, memTotalGB] = getMemoryUsage();
fprintf(fidOut,'\n\n\n\n ---- Processsing %s w/ %d GB comp memory (%.1f GB actual)::  %s ---- \n', subj, compMemoryGB, memTotalGB, datestr(now));


% Load utahInfo CSV files:  get the "SESSIONS" and the "CHAN STRINGS"
tSessList = readtable(fileDir_bookkeeping);
tSessList = tSessList(strcmpi(tSessList.toExtract,'y'),:); %- only look at the "to extract" sessions.  Dont need the ns5 file, just the session folder name

%-
if ~exist(fileDir_spikeChans,'file'),
    fprintf('\n ERROR: extractMicrophysiology cant find %s... skipping extraction\n',fileDir_spikeChans);
    return;
end
spikeChans = csv2cell(fileDir_spikeChans);
tChanList  = cell2table(spikeChans(4:end,1:5),'variablenames',spikeChans(3,1:5));
devNum = []; for iH=1:height(tChanList), devNum(iH) = str2num(tChanList.MicroDevNum{iH}); end; tChanList.MicroDevNum = devNum'; %- convert MicroDevNum to numeric... easier/faster to work with below


tStartAll = tic;


%- Loop over each "to extract" session
for iS=1:height(tSessList),
    
    tStartSess = tic;
    thisSess   = tSessList.folderName{iS};
    
    isStimSess = strcmpi(tSessList.isStim(iS),'y');
    if isStimSess, outputType = 'stimulation';
    else           outputType = 'standard'; end
    if ~isempty(overrideOutputType), outputType = overrideOutputType; end
    
    
    %- select which type of outputs are produced.  automatically changes for stim vs non-stim session
    switch outputType,
        case 'justLFP',
            extractOutputs = {                                                              'lfp_noreref'             };   %- just LFP
        case 'justRawSpike',
            extractOutputs = {                                                 'sort_raw'                           };   %- split the ns5/6 without reref or filter, only used for debugging
        case 'minimum',
            extractOutputs = {'sort_reref'                                                'lfp_noreref'             };   %- minimum outputs... wont be able to make a complete spikeInfo without sort_unfilt
        case 'standard'
            extractOutputs = {'sort_reref'                   'sort_unfilt'              'lfp_noreref'             };   %- typical outputs
        case 'stimulation',
            %extractOutputs= {'sort_reref'  'sort_noreref' 'sort_unfilt'              'lfp_noreref'             };   %- originally had additional outputs for stim sessions
            extractOutputs = {'sort_reref'                   'sort_unfilt'              'lfp_noreref'             };   %- now stim session split is initially the same, until masks are applied
        case 'all',
            extractOutputs = {'sort_reref'  'sort_noreref' 'sort_unfilt' 'sort_raw' 'lfp_noreref'             };   %- everything possible (lfp_reref not supported anymore)
        otherwise
            fprintf('\n error... output types not specified');
            keyboard;
    end
    %- sort_raw     = 30kS time series,                                                           .  Output 1 .bin per channel + global per device + channel name mapping
    %- sort_noreref = 30kS time series,                                    band pass 300 - 3000 Hz.  Output 1 .bin per channel + global per device + channel name mapping
    %- sort_reref   = 30kS time series, common average subtracted and then band pass 300 - 3000 Hz.  Output 1 .bin per channel + global per device + channel name mapping
    %- sort_unfilt  = 30kS time series, common average subtracted and then band pass   1 - 5000 Hz.  Output 1 .bin per channel + global per device + channel name mapping
    %- lfp_noreref    =  1kS time series,                          low pass 500Hz, then downsampled.   Output 1 mat file for ALL channels + channel name mapping
    
    %- LFP-REREF is no longer an option to be extracted here... reref and/or processed versions of LFP should be created after downsampling (no substantial difference in reref if using 30k, so dont bother)
    
    
    
    
    [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
    fprintf(1,      '\n --- %d of %d) %s Session %s, start %s --- [%.1f GB free, %.1f GB used]',iS,height(tSessList),subj,thisSess,datestr(now),memFreeGB,memUsedGB);
    fprintf(fidOut, '\n --- %d of %d) %s Session %s, start %s --- [%.1f GB free, %.1f GB used, %.1f GB compresssed]',iS,height(tSessList),subj,thisSess,datestr(now),memFreeGB,memUsedGB,memCompressGB);
    
    
    rawDir          = fullfile(subjectInfo.dataPath56,subj,'data_raw',thisSess);
    processedDir    = fullfile(subjectInfo.dataPath56,subj,'processed_SPK',sprintf('%s_ready',thisSess));
    if ~exist(processedDir,'dir')
        processedDir    = fullfile(subjectInfo.dataPath56,subj,'processed_SPK',thisSess);
    end
    processedDirLFP = fullfile(subjectInfo.dataPath56,subj,'processed_LFP',thisSess); %- LFPs are self contained with pulse info, they can get their own root
    
    extraFilesDir   = '_extractionInfo'; %- sub-directory in each of the sort_XXref folders that has all the extra extraction/alignment info. Underscore so at top of directory in finder
    
    
    %- make the required output directories if not made yet
    needsExtraction = zeros(length(extractOutputs),1);
    for iOut=1:length(extractOutputs)
        if contains(extractOutputs{iOut},'lfp')
            extractDir = fullfile(processedDirLFP,extractOutputs{iOut});
        else
            extractDir = fullfile(processedDir,extractOutputs{iOut});
        end
        
        if ~exist(extractDir,'dir')
            mkdir(extractDir); needsExtraction(iOut)=1;
            if contains(extractOutputs{iOut},'sort') %-SJ 2/24/20: changed 'spike' to 'sort'
                mkdir(fullfile(extractDir,extraFilesDir)); %- every sort_XX folder will have an extractionInfo subfolder
            end
        else
            fList = dir(extractDir);
            fList = fList([fList.isdir]==0); %- cut out directories, which could be '.', '..', and '_prunedNoSpikesNPlay'
            if length(fList)==0, needsExtraction(iOut)=1; end  %- are there ANY files in here?  later this could be made to match the number of expected files
        end
    end
    
    %- check for the stimMask folder... create if everything is there EXCEPT for that
    if sum(needsExtraction)==0 && isStimSess
        stimMaskDir = dir(fullfile(processedDir,'*stimMask*'));
        if length(stimMaskDir)==0
            fprintf('\n Session completely extracted EXCEPT for stimMask outputs... making those now');
            testOnSomeChans = 0;  pickSessDir = thisSess; %- JW found that 5ms margin worked well with NIH069 and 66...
            maskStimFromSpikes_v2(subj,overrideDataPath56, maskStimMsGrowDC10pulse, testOnSomeChans, pickSessDir);
        end                   
    end
    
    
    %-  Offer a cheap way of jumping past this session
    if sum(needsExtraction)==0 & PUSH_NEW_SORT_HELPER==0

        fprintf('\n Looks like this session was already processed. Skipping it\n');
        continue
    end
    
    
    %- do we need to compute the commonAve at 30kS for this extraction?
    
    if ~any(contains(extractOutputs,{'sort_reref'})) 
        GET_COMMON_AVE = 0; 
    else
        GET_COMMON_AVE = 1;
        commonAveDir = fullfile(processedDir, 'sort_reref',extraFilesDir,'rerefStuff');  
        if ~exist(commonAveDir,'dir'), mkdir(commonAveDir); end
    end
    
    
    %- load complete jackTable.  Will probably want to use t
    thisJack  = fullfile(rawDir,'jacksheetBR_complete.csv');
    tJackAll  = readtable(thisJack);
    
    
    %- merge tChanList and tJackAll so jacksheet has the new-string & device num info into tJackAll
    flag_chansNeedRename = 0; %- flag... set to 1 if jackSheetBR_complete does NOT have up to date micro names and device number
    for iH=1:height(tChanList)
        %- only check for name changes on the 30k sampleRate... those are the names that will be split out with this extraction (not the 1kS version of the same channel)
        iJack = find(strcmp(tJackAll.NSPsuffix,tChanList.NSPsuffix{iH}) & strcmp(tJackAll.ChanName,tChanList.ChanName{iH}) & tJackAll.SampFreq==30000); %- completely unambiguous: NSPsuffix + ChanName + sampleRate.
        if length(iJack)==1
            tJackAll.SortChanName{iJack} = tChanList.SortChanName{iH}; %- grab the SortChanName from micro_RenamedChanInfo so it is only defined in one place
            if ~strcmp(tChanList.ChanNameNew{iH},tJackAll.ChanNameNew{iJack}) | tChanList.MicroDevNum(iH)~=tJackAll.MicroDevNum(iJack)
                flag_chansNeedRename = 1;
                tJackAll.ChanNameNew{iJack} = tChanList.ChanNameNew{iH};
                tJackAll.MicroDevNum(iJack) = tChanList.MicroDevNum(iH);
            end
        elseif length(iJack)==0
            %- this could happen if a channel was on for some recordings and off for others.  Or if a channel switched to a different NSP.
            % tChanList has all possibilities, but jacksheet only has names from this particular session
        else
            fprintf('\n shouldnt happen');
            keyboard;
            error('\n figure out why this happened');
        end
    end
    
    if flag_chansNeedRename,
        fprintf('\n HEADS UP:  ChanNameNew and/or MicroDevNum not matching: \n         %s\n      vs %s.  \n      --> micro_Renamed values will be used in LFP.mat and split jacksheets\n',fileDir_spikeChans,thisJack);
        keyboard;
    end
    
    
    %- Note: Dont save a copy of the jacksheetBR_complete... makes things confusing that this top level copy *could* deviate from the individual sessions if done at a different time
    %- instead create a session and extraction specific "jackSheetBR_split", which conveys the naming scheme used in this particular split
    
    %- list of channels to extract as spikes and LFPs, this table also says what ns5 or ns6 file each channel belongs to
    tJackSplit = tJackAll(ismember(tJackAll.ChanName,tChanList.ChanName) & tJackAll.SampFreq==30000 & tJackAll.PhysicalChan<1000,:);
    tJackPulse = tJackAll(tJackAll.PhysicalChan>=1000,:);
    
    
    %- Confirm that "SortChanName" matches expectation.  Probably could remove this check...
    %   Name should be NSPoffset combined with the old channel name...  that is what the split channels are called. This column will be saved with all split versions of jacksheet
    for iChan=1:height(tJackSplit),
        if ~strcmp(tJackSplit.SortChanName{iChan}, sprintf('[%s]%s',tJackSplit.NSPsuffix{iChan},tJackSplit.ChanName{iChan}))
            fprintf('\n ERROR: sortChan name not matching expectation. Shouldnt ever happen. Investigate');
            keyboard
        end
    end
    
    
    %- loop over nsX files so we can be efficient with loading those data... but first confirm that a single device never spans multiple nsX files
    %- loop over the devices, devices should never span NSPs, or ns5 vs ns6 files... at least we dont need to code for that
    
    %- quick sanity check... device Nums should only source from one NSP (or one unique file)
    devNums = unique(tJackSplit.MicroDevNum);
    for iD = 1:length(devNums),
        fileNames = unique(tJackSplit.FileName(tJackSplit.MicroDevNum==devNums(iD)));
        if length(fileNames)>1, fprintf('\n Uh oh. Not expecting a device to span multiple ns5 or ns6 files. Not ready for this'); keyboard; error('\n crash and burn'); end
    end
    
    
    %- add the microstimChan to tJackSplit
    utahChan = mod(tJackSplit.PhysicalChan,128);
    utahChan(utahChan==0)=128;
    stimChan = zeros(size(utahChan));
    if ~isnan(microStimChan)
        stimChan(microStimChan) = [1:length(microStimChan)]; %- record the stim pigtail number
    end
    tJackSplit.microStimPigtail = stimChan;

    
    %- use tJackSplit to create a spike sorting helper function
    if any(contains(extractOutputs,'sort'))
        
        %- if processsing spikes create/update the spike helper here (so JW can push an update without processing spikes again)
        
        if PUSH_NEW_SORT_HELPER == 1
            writeSpikeHelperFiles(subj, tJackSplit, thisSess, isStimSess, processedDir, rawDir, subjectInfo, fidOut, compMemoryGB, PUSH_SPK_EXTFIX);
        
            if sum(needsExtraction)==0 %- and skip any further processing if just making an updated sort helper
                fprintf('\n Looks like this session was already processed, but PUSH_NEW_SPIKE_HELPER was flagged so did that then Skipping it\n');
                continue    
            end
        end
    end
    
    
    %- Identify which NSPs the sorted channels come from, and include a pulse structure (for subsequent alignemnt) for any NSP included
    fileNames = unique(tJackSplit.FileName,'stable'); %- stable to keep in order of jacksheetBR.  probably doesn't matter tho
    
    %- add a bit of info to the command line output:  #channels and recording duration
    fprintf(1,      '\n          -- Extracting %d micro chan; %d devNums; %d ns5/6 files; %.1fs min recording: ', height(tJackSplit), length(devNums), length(fileNames), tJackSplit.DurationMin(1));
    fprintf(fidOut, '\n          -- Extracting %d micro chan; %d devNums; %d ns5/6 files; %.1fs min recording: ', height(tJackSplit), length(devNums), length(fileNames), tJackSplit.DurationMin(1));
    
    
    %- create a single date/time string for this meta extraction... this can be used to link together all the extractions (LFP, noreref, reref)... 
    %    in the past I grabbed a fresh datestr when each extraction type was created, but this could lead to confusion about whether the same extraction code was used for all
    sessExtractDateStr = datestr(now,21); 
    
    
    %- Loop over ns5/6 files and extract
    pulsesStructCell   = {};
    postProcStructCell = {};
    fprintf(1, '(making pulse struct:'); ticPulse = tic;
    for iF = 1:length(fileNames)
        
        %- NSP suffix from jackTableBR... from a row with this fileName
        nspSuffix = tJackSplit.NSPsuffix{find(strcmp(tJackSplit.FileName,fileNames{iF}),1,'first')};
        
        if contains(fileNames{iF},'ns5')
            fileType = 'ns5';
        elseif contains(fileNames{iF},'ns6')
            fileType = 'ns6';
        end
        ns3File = strrep(fileNames{iF},fileType,'ns3');  ns3File = fullfile(rawDir,ns3File);
        nevFile = strrep(fileNames{iF},fileType,'nev');  nevFile = fullfile(rawDir,nevFile);
        if ~exist(ns3File); ns3File = strrep(ns3File,'ns3','ns4');  end
        
        
        %- extract the pulse up times and pulse time series, to append to LFP.mat for subsequent alignment
        [pulses] = makePulsesStruct('ns3_fpath',ns3File,'nev_fpath',nevFile);
        
        
        %- no longer saving the individual pulse structs... instead save the cell array of pulse structs in each extraction folder (e.g., noreref, reref, etc)

        %- this cell array will be attached to the LFP.mat below; and also save to processed_SPK/_extaction_notes/pulseStructCell.mat
        pulsesStructCell{iF,1} = nspSuffix;
        pulsesStructCell{iF,2} = fileNames{iF};
        pulsesStructCell{iF,3} = pulses;
    end
    fprintf(1,'%.1fs)',toc(ticPulse));
    
    
    %- now the big loop over nsX files.  will be >1 if using 2 NSPs, or for some reason mixed ns6 and ns5 data (did that with NIH050 microwires)
    for iF = 1:length(fileNames)
        
        %- these are the channels that we plan to split from this file
        tJackFileSplit = tJackSplit(strcmp(tJackSplit.FileName,fileNames{iF}),:);
        
        thisNSPsuffix  = tJackFileSplit.NSPsuffix{1};
        
        
        %- list of device number for the channels to split
        devNums  = unique(tJackFileSplit.MicroDevNum);
        devNames = {};
        for iiDev=1:length(devNums)
            iDev = find(subjectInfo.deviceNum==devNums(iiDev));
            devNames{end+1} = subjectInfo.newChanName{iDev,1}; %- create cell array of device names matching the device nums on this NSPs file
        end
        
        %- pull a "raw" jacksheet from the file of interest... best way to be 100% sure that channel names are matching channel numbers
        NSxFileDir     = fullfile(rawDir,fileNames{iF});
        dataFromFile   = openNSx(NSxFileDir,'noread');
        labelsFromFile = {dataFromFile.ElectrodesInfo.Label};
        labelsFromFile = cellfun(@deblank,labelsFromFile,'UniformOutput',false); % remove whitespace
        
        %- index of channels to pull from nsX
        iChanRead      = find(ismember(labelsFromFile,tJackFileSplit.ChanName));
        
        % determine how many channels to load at once based on computer memory
        [chanGroups, fileSizeGB] = getChanGroupsToLoad(NSxFileDir,iChanRead,compMemoryGB);
        if isempty(chanGroups), fprintf('\n uh oh... not enough memory to load a single channel'); keyboard; end
        
        %- extract group by group and compute the commonAve for each unique device in the file
        [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
        fprintf(fidOut,'\n LOADING ALL DATA TO COMPUTE COMMON AVERAGE BEFORE SPLITTING FILES  <current Memory: used %.1f GB, free %.1f GB, compressed %.1f GB>',memUsedGB,memFreeGB,memCompressGB);
        fprintf(fidOut,'\n [filesize %.0f GB; %d groups of <=%d chan to load]', fileSizeGB, length(chanGroups), length(chanGroups{1}));
        for grp = length(chanGroups):-1:1 %- start with last group because it can have fewer channels and we are gonna skip loading below (so want to end with most channels)
            
            %- load the data (slow step)
            fprintf(fidOut,'\n   loading group %d (%d chan)',grp,length(chanGroups{grp})); tLoad = tic;
            fprintf(1,'.');
            clear NSxData; %- make room because openNSx is going to use a lot of temp space?
            [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
            NSxData = concatOpenNSx(NSxFileDir,CONVERT_2_uV,1,chanGroups{grp});
            
            
            %- grab the concat NSx info... that way if things change in the future (i.e., gain; or whether we throw away "junk", etc... we'll know what was in this one)
            physio_nsx_postProc = NSxData.postProc;  %- get the "post process" string returned from concatOpenNSx and save it to LFP and spike structures
            postProcStructCell{iF,1} = thisNSPsuffix;
            postProcStructCell{iF,2} = fileNames{iF};
            postProcStructCell{iF,3} = NSxData.postProc; %- will be filled out below
            
            NSxData = NSxData.Data;
            [memUsedGB2, memFreeGB2, memTotalGB2, memCompressGB2] = getMemoryUsage();
            fprintf(fidOut,' [%.1fs; %.1fGB --> %.1f GB used;  %.1f -> %.1f GB free]',toc(tLoad),memUsedGB, memUsedGB2,memFreeGB,memFreeGB2);
            fprintf(1,'.');
            
            %- now that we know how long the time series is we can initailize the commonAve time series (1 per device)
            if grp==length(chanGroups)
                numSampBR = size(NSxData,2);
                commonAve = zeros(length(devNums),numSampBR);
                cAveCount = zeros(length(devNums),1);
                exampChan = commonAve; %- single trace from each device for the output figure
            end
            
            %- create a list of the channel names and device numbers that were just loaded
            chanStrLoad = labelsFromFile(chanGroups{grp});
            chanDevLoad = []; chanIsStim = []; %chanNotStim = []; - SJ changed to IsStim, chanNotStim not used anywhere else
            for iC=1:length(chanStrLoad)
                chanDevLoad(iC) = tJackFileSplit.MicroDevNum(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC})); 
                chanIsStim(iC)  = tJackFileSplit.microStimPigtail(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC}));  %- 0 for non-stim, >0 for stim
            end
            
            %- chanGroup might contain multiple devices, so have each channel contribute to the ave of its own device
            if GET_COMMON_AVE
                for iD=1:length(devNums)
                    iChanD = find(chanDevLoad==devNums(iD) & chanIsStim==0); %- all channels from this device EXCLUSING microStim channels, which have weird noise even if not stim session
                    if length(iChanD)>0
                        if cAveCount(iD)==0
                            exampChan(iD,:) = double(NSxData(iChanD(1),:)); %- save the first trace from each device to output in a figure below
                        end
                        commonAve(iD,:) = commonAve(iD,:)+sum(NSxData(iChanD,:),1); %- int16 can sum, dont need to convert to double here
                        cAveCount(iD)   = cAveCount(iD)+length(iChanD);
                    end
                end % for iD
            end
        end % for grp
        
        
        %- save out 30kS common average if outputing spikes or lfp_reref
        if GET_COMMON_AVE
            
            %- output the results
            fprintf(fidOut,'\n saving common average time series(s) '); tStart = tic;
            for iD=1:length(devNums)
                %- convert to a mean
                commonAve(iD,:) = commonAve(iD,:) / cAveCount(iD);
                
                % Save common average as .bin that could be loaded to plexon (only getting saved in processed, not processedLFP, because we only save noreref for LFP
                commonAvFileName_raw = sprintf('globalAvg_30kRaw_(%s)_[dev%d].bin', devNames{iD}, devNums(iD));
                fid = fopen(fullfile(commonAveDir,commonAvFileName_raw), 'w');
                fwrite(fid,commonAve(iD,:),'int16');  %Save for plexon... confirm this is saving in the correct format
                fclose(fid);
                
                % Save the reref figure for one channel... probably should do this separately for each device
                figSaveDir = fullfile(commonAveDir,sprintf('rerefPlot_(%s)_[dev%d]', devNames{iD}, devNums(iD)));
                getRerefFigs(subj,exampChan(iD,:),commonAve(iD,:),figSaveDir,'');
            end
            clear exampChan;
            fprintf(fidOut,' DONE [%.1fs] !\n', toc(tStart));
            fprintf(1,'.');
            
            %- sanity check figure...
            MAKE_GLOBAL_FIG = 0;
            if MAKE_GLOBAL_FIG
                figure(1);clf
                for i=1:size(commonAve,1)
                    hold on
                    subplot(1,size(commonAve,1),i)
                    plot(NSxData(:,1:1000:end)'); hold on;
                    plot(commonAve(1,1:1000:end),'k','linewidth',4); hold on;
                    tStr = sprintf('%s%s%s',processedDir,'-array ', subjectInfo.oldChanName{i,1});
                    tStr(tStr=='_')='-';
                    title(tStr);
                end
                keyboard
                % from the examplie I just looked it (NIH066), seems like it would make more sense to regress out the common signal rather than subtract from everything...
                %  most follow common but some do not.
            end
            
        end  %- if GET_COMMON_AVE
        
        
        
        %- Now loop over each channel, and output the requested data
        %
        
        %- create matricies to hold the LFP data... it will all be written at once after processing all channels
        if any(contains(extractOutputs,'lfp'))
            fprintf(fidOut,'\n Initalizing LFP output matricies.  [%.1f GB --> ', getMemoryUsage()); tInitLFP = tic;
            lfpDS = [1:(samprate_br/samprate_LFP):numSampBR];
            numSampLFP = length(lfpDS);
            if any(strcmp(extractOutputs,'lfp_noreref')), lfp_noreref = zeros(length(iChanRead),numSampLFP,'int16'); end
            clear tLFPdata; %- to become a table of lfp channels
            fprintf(fidOut,' %.1f GB;  %.1fs]', getMemoryUsage(), toc(tInitLFP));
        end
        
        %- time series of zeros used in place of common average for "noreref" outputs
        zeroCommon = zeros(1,numSampBR);
        
        [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
        fprintf(fidOut,'\n RE-LOADING DATA TO SPLIT FILE  <current Memory: used %.1f GB, free %.1f GB, compressed %.1f GB>',memUsedGB,memFreeGB,memCompressGB);
        fprintf(fidOut,'\n [filesize %.0f GB; %d groups of <=%d chan to load]', fileSizeGB, length(chanGroups), length(chanGroups{1}));
        for grp = 1:length(chanGroups)
            
            %- first group currently loaded and in memory, so start reload on group 2 (avoiding a reload of the same data)
            if grp > 1
                fprintf(fidOut,'\n   loading group %d (%d chan)',grp,length(chanGroups{grp})); tLoad = tic;
                clear NSxData; %- make room because openNSx is going to use a lot of temp space?
                [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
                NSxData = concatOpenNSx(NSxFileDir,CONVERT_2_uV,1,chanGroups{grp});
                NSxData = NSxData.Data;
                [memUsedGB2, memFreeGB2, memTotalGB2, memCompressGB2] = getMemoryUsage();
                fprintf(fidOut,' [%.1fs; %.1fGB --> %.1f GB used;  %.1f -> %.1f GB free]',toc(tLoad),memUsedGB, memUsedGB2,memFreeGB,memFreeGB2);
                fprintf(1,'.');
                
                %- create a list of the channel names and device numbers that were just loaded
                chanStrLoad = labelsFromFile(chanGroups{grp});
                chanDevLoad = []; for iC=1:length(chanStrLoad), chanDevLoad(iC) = tJackFileSplit.MicroDevNum(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC})); end;
            end
            
            
            %- now loop over each channel separately and create all the split offs
            fprintf(fidOut,'\n splitting channels:');
            for iC = 1:length(chanGroups{grp})
                
                tStartChan     = tic;
                chanDataRaw    = double(NSxData(iC,:));
                thisChanIndex  = chanGroups{grp}(iC);
                thisChanStrOld = chanStrLoad{iC};            %-channel name
                thisChanStrOut = tJackFileSplit.SortChanName{strcmp(tJackFileSplit.ChanName,thisChanStrOld)};  %- spike-sorting.bin filename [nspSuffix]ChanName.bin (use the new table column to make sure single source of this name)
                thisChanDevNum = chanDevLoad(iC);            %-actual device number, used for file name outputs
                thisChanDevInd = find(devNums==thisChanDevNum); %-inded into devNums, which is also index into commonAve
                
                fprintf(fidOut,'\n   %s (old string=%s; devNum=%d; indexThisNSP=%d)',thisChanStrOut,thisChanStrOld,thisChanDevNum,thisChanIndex);
                
                
                %- NSxData is 30kS raw data.
                for iOut=1:length(extractOutputs)
                    switch extractOutputs{iOut},
                        
                        case 'sort_unfilt',
                            %- sort_unfilt.  subtract the global mean, filter 1-5000, output the channel
                            numeratorCoeff   = b_spike_unfilt;
                            denominatorCoeff = a_spike_unfilt;
                            subCommon        = commonAve(thisChanDevInd,:);
                            strFilt          = str_spike_unfilt;
                            processedDirUse  = processedDir;
                            
                        case 'sort_reref',
                            %- sort_reref.  subtract the global mean, filter 300-3000, output the channel
                            numeratorCoeff   = b_spike;
                            denominatorCoeff = a_spike;
                            subCommon        = commonAve(thisChanDevInd,:);
                            strFilt          = str_spike;
                            processedDirUse  = processedDir;
                            
                        case 'sort_noreref',
                            %- sort_reref.  subtract the global mean, filter 300-3000, output the channel
                            numeratorCoeff   = b_spike;
                            denominatorCoeff = a_spike;
                            subCommon        = zeroCommon;
                            strFilt          = str_spike;
                            processedDirUse  = processedDir;
                            
                        case 'lfp_noreref',
                            %- lfp_noreref.   low pass filter and downsample
                            numeratorCoeff   = b_LFP;
                            denominatorCoeff = a_LFP;
                            subCommon        = zeroCommon;
                            strFilt          = str_LFP;
                            processedDirUse  = processedDirLFP;
                            
                        case 'sort_raw',
                            %- dont filter. just extract
                            numeratorCoeff   = nan;
                            denominatorCoeff = nan;
                            subCommon        = zeroCommon;
                            strFilt          = 'raw spike time series split from ns5 or ns6 without reref or filtering';
                            processedDirUse  = processedDir;
                    end
                    
                    fprintf(1,'.'); %- let the user know something is happening...
                    
                    %- possible to filter all channels at once?  one at a time takes 10x longer than just saving out the channel (30 vs 3sec per channel)
                    %keyboard
                    %fprintf('\n filering all channels at once (super memory intensive)'); tFiltStart = tic;
                    %filteredAll = [filtfilt(numeratorCoeff,denominatorCoeff,[double(NSxData)-(ones(length(chanGroups{grp}),1)*subCommon)]')]'; %- filtfilt operates along first non-zero dimension, so make it time x channels then convert back
                    %fprintf(' [%.1f]',toc(tFiltStart));
                    
                    %- apply zero-phase-shift filter...
                    if strcmp(extractOutputs{iOut}, 'sort_raw')
                        filtered = chanDataRaw; %- just split the channels without reref or filter
                    else
                        filtered = filtfilt(numeratorCoeff,denominatorCoeff,chanDataRaw-subCommon);
                    end
                    
                    %- different steps for LFP and SPIKES
                    if contains(extractOutputs{iOut},'lfp')
                        %- downsample the LFPs
                        tLFPdata(thisChanIndex,:) = tJackFileSplit(strcmp(tJackFileSplit.ChanName,thisChanStrOld),:); %- make a table matched to the LFP data
                        if     strcmp(extractOutputs{iOut},'lfp_noreref'), lfp_noreref(thisChanIndex,:) = int16(filtered(lfpDS));
                        elseif strcmp(extractOutputs{iOut},'lfp_reref'),   lfp_reref(thisChanIndex,:)   = int16(filtered(lfpDS)); end
                        
                    else
                        %- spike output .bin to disk
                        spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},sprintf('%s.bin',thisChanStrOut));
                        fchan = fopen(spkBinFile,'w','l');
                        fwrite(fchan,filtered,'int16'); %- convert to int16 before writing
                        fclose(fchan);
                    end
                    
                    %- save a copy of the jacksheet table in each output directory.  LFPs also get a copy embedded in .mat
                    if iF ==1 & grp ==1 & iC == 1    %- first channel dumped... only need to save the table once per outputDir
                        if contains(extractOutputs{iOut},'sort')
                            extraFilesDirUse = extraFilesDir; %- for spikes put all extraction/alignment stuff in a subfolder
                        else
                            extraFilesDirUse = ''; %- for LFPs there is not much in the folder so put them there
                        end
                        jackOut = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'jacksheetBR_justSplitChan.csv');
                        writetable(tJackSplit,jackOut); %- just the split channels contained within this folder..
                        
                        %- also make a little helper readme file with extraction date/time, computer name, filter settings
                        fidTXT = fopen( fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'_extractMicro_readme.txt'), 'a+');
                        fprintf(fidTXT, '\n\n\n**** extract MicroPhys on %s ***',  sessExtractDateStr);
                        fprintf(fidTXT, '\n  full directory path: %s ',       fullfile(processedDirUse,extractOutputs{iOut}) );
                        fprintf(fidTXT, '\n  run on computer %s, by user %s ',char(java.net.InetAddress.getLocalHost.getHostName),char(java.lang.System.getProperty('user.name')) );
                        fprintf(fidTXT, '\n  filtering:  %s ',          strFilt);
                        if all(subCommon==0)
                            fprintf(fidTXT, '\n  referencing: NO common average rerefrence subtracted during extraction');
                        else
                            fprintf(fidTXT, '\n  referencing: Common average rereference (by device) subtracted before filtering');
                            fprintf(fidTXT, '\n               %d stimulation micro channels removed from reference groups', sum(tJackFileSplit.microStimPigtail>0));
                            fprintf(fidTXT, '\n               total num channels per reference (device) group: ');
                            for iii=1:length(cAveCount), fprintf(fidTXT, '%d, ',cAveCount(iii)); end
                        end
                        fclose(fidTXT);
                        
                        %- if spike folder, save the pulseStruct and postProcStruct so it can be appended from this folder instead of requiring "other"
                        if contains(extractOutputs{iOut},'sort')
                            
                            %- output copy of filtered first channel
                            spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('%s.bin',thisChanStrOut));
                            fchan = fopen(spkBinFile,'w','l');
                            fwrite(fchan,filtered,'int16'); %- convert to int16 before writing
                            fclose(fchan);
                            
                            %- output copy of raw first channel
                            spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('%s_raw.bin',thisChanStrOut));
                            fchan = fopen(spkBinFile,'w','l');
                            fwrite(fchan,chanDataRaw,'int16'); %- convert to int16 before writing
                            fclose(fchan);
                            
                            %- and save a full copy of the pulseStructCell.mat, to be used with alignment
                            pulseCellOut  = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('pulseStructCell.mat',nspSuffix));
                            save( pulseCellOut, 'pulsesStructCell',  '-v7.3');
                       
                            %- and the postProc cell... THIS GIVES DEEPER INFO ABOUT THE NS5/6 SPLIT (openNSx)
                            postProcCellOut  = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('postProcStructCell.mat',nspSuffix));
                            save( postProcCellOut, 'postProcStructCell',  '-v7.3');
                        end

                    end %- copy files
                    
                end %- for extractionOut
                
                %- sanity check
                if 0 & iC==1
                    figure(1); clf;
                    subplot(211);  plot(double(lfp_noreref(2,1:1000)),'r'); title('noreref');
                    subplot(212);  plot(double(lfp_reref(2,1:1000)),'r');   title('reref');
                end
                
                fprintf(fidOut,' [%.1fs]', toc(tStartChan));
                
            end %- channels in this chanGroup
            
            clear NSxData;
            
        end %- for all chanGroups
        
        
        %- now save the LFPs to disk as a single mat file.  LFP matrix is already in int16 format above
        if any(contains(extractOutputs,'lfp'))
            
            %- first NSP
            if iF==1
                
                %- package the contents
                lfpStruct = struct;
                
                %- this is a copy of the mountainsort info string...
                infoStr = sprintf('this lfp.mat file, generated %s, contains the following fields:\n',datestr(now,21));
                infoStr = sprintf('%s\n createdDate         - a string indicating when this spikeInfo.mat was created',infoStr);
                infoStr = sprintf('%s\n chanIDperNSP        - a ( #NSPs used x 1 ) cell array. Each enty is a table corresponding to the channel data in lfp',infoStr);
                infoStr = sprintf('%s\n lfp                 - a ( #NSPs used x 1 ) cell array. Each cell entry is #channels x #time points matrix (int16). Time point dimension might be slightly different prior to NSP alignment',infoStr);
                infoStr = sprintf('%s\n gain_bin2uV         - multiplicative factor to convert lfp_noreref to uV',infoStr);
                infoStr = sprintf('%s\n samplingFreq        - samlple frequency of the data. typically 1000 Hz',infoStr);
                infoStr = sprintf('%s\n filterSettings      - a ( #NSPs used x 1 ) cell array.  Each cell entry is #channels x #time points matrix (int16). Time point dimension might be slightly different prior to NSP alignment',infoStr);
                infoStr = sprintf('%s\n sessStr             - the session name the spikes were extracted from (e.g., 190117_1336)',infoStr);
                infoStr = sprintf('%s\n sessDurMin          - session duration in minutes',infoStr);
                infoStr = sprintf('%s\n physio_nsx_postProc - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "nsx_postProc" struct containing information about how data was pulled from the NSx file (i.e., concatenate or junk additional segments)',infoStr);
                infoStr = sprintf('%s\n pulses              - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "pulse" struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries for that NSP',infoStr);
                infoStr = sprintf('%s\n alignedTo           - a string indicating which file the spikes have been aligned to',infoStr);
                infoStr = sprintf('%s\n alignmentChan       - a string indicating which channel in the pulses struct was used for alignment',infoStr);
                infoStr = sprintf('%s\n jackTableFull       - complete jacksheetBR table from this session (all channels) with device numbers and new channel names',infoStr);
                infoStr = sprintf('%s\n jackTableUsed       - just the jacksheet for the channels incorporated into this spikeInfo (combined across NSPs). this table can be used to go back and forth between original and new channel names',infoStr);
                
                %%-
                lfpStruct.readme              = infoStr;
                lfpStruct.createdDate         = datestr(now,21);
                lfpStruct.chanIDperNSP{1}     = tLFPdata; %- jackTable split by NSPsuffix, then ordered as the channels are placed into the lfp data matrix... so one-to-one correspondance of index with data
                lfpStruct.lfp{1}              = lfp_noreref;
                lfpStruct.rerefType           = extractOutputs{contains(extractOutputs,'lfp')}; %- no longer supporing reref, so only noreref should ever get here
                lfpStruct.gain_bin2uV         = gain_bin2uV;
                lfpStruct.samplingFreq        = samprate_LFP;
                lfpStruct.filterSettings      = str_LFP;
                lfpStruct.sessStr             = tJackSplit.RawDir{1};
                lfpStruct.sessDurMin          = tJackSplit.DurationMin(1);
                lfpStruct.physio_nsx_postProc = postProcStructCell; %- info about jumk segments
                lfpStruct.pulses              = pulsesStructCell;   %- will be used for alignment
                lfpStruct.alignedTo           = [];
                lfpStruct.alignmentChan       = [];
                lfpStruct.jackTableFull       = tJackAll;
                lfpStruct.jackTableUsed       = tJackSplit;
                
            elseif iF==2
                %- if two NSPs recording micros, then fill the cell array with the second jackFile and data time series matric
                lfpStruct.chanIDperNSP{2}   = tJackFileSplit;
                lfpStruct.lfp{2}            = lfp_noreref;
                
            end
            
            %- Done looping over fileNames, so dump the resultant matfile
            if iF==length(fileNames)
                
                fprintf(fidOut,'\n SAVING LFP.mat To DISK');
                fprintf(fidOut,'\n lfp_noreref) saving to disk'); tSave = tic;
                lfpMatFile  = fullfile(processedDirUse,'lfp_noreref','lfp.mat');
                save(lfpMatFile, 'lfpStruct', '-v7.3');
                fprintf(fidOut, ' [%.1f s]',toc(tSave));
                clear lfp_noreref;
            end
        end
        
    end %- for iF = 1:length(fileNames)   loop over file stems (i.e., INST0 and INST1)
    
    fprintf(1     ,'  Extraction finished for %s [%.1f min]\n',thisSess,toc(tStartSess)/60);
    fprintf(fidOut,'  Extraction finished for %s [%.1f min]\n',thisSess,toc(tStartSess)/60);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process the stim session masks using default parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- if a freshly extracted stim session, create a masked version of the spike time series for sorting.
    if isStimSess & any(contains(extractOutputs,'sort')) %- SJ 2/24/20: 'spike' changed to 'sort'
        %- now loop over extracted sessions and created masked version of splits for any stim session
        testOnSomeChans = 0;  pickSessDir = thisSess; %- JW found that 5ms margin worked well with NIH069 and 66...
        maskStimFromSpikes_v2(subj,overrideDataPath56, maskStimMsGrowDC10pulse, testOnSomeChans, pickSessDir);
    end
    
    
end %- loop over sessions to extract


fprintf(1     ,'\nExtraction finished for %d sessions [%.1f min]\n',height(tSessList),toc(tStartAll)/60);
fprintf(fidOut,'\nExtraction finished for %d sessions [%.1f min]\n',height(tSessList),toc(tStartAll)/60);



if fidOut~=1, fclose(fidOut); end

% SJ - Copy _extractMicrophysLog.txt to PUB 
PUBpathcomplete = [PUBpath filesep subj filesep '_extraction_notes'];
if exist(PUBpathcomplete,'dir')
    st1 = copyfile(extractFile,PUBpathcomplete);
    if st1 == 1
        fprintf('%s\n','_extractMicrophysLog.txt copied to PUB.')
    else
        fprintf('%s\n','ERROR!!! _extractMicrophysLog.txt was not copied to PUB!!!')
        keyboard
    end
else
    fprintf('\n%s\n','Subject not established in 56PUB yet. CSV will be copied in COPY_2_PUBforSorting step. Are you actually connected to 56PUB?');
    keyboard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the utahInfo csv file to show extracted sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
updateUtahInfoCSV_v10(subj,overrideDataPath56,true); % SJ - 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HELPER FUNCTION... make the sort helper files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeSpikeHelperFiles(subj, tJackSplit, thisSess, isStimSess, processedDirUse, rawDir, subjectInfo, fidOut, compMemoryGB, PUSH_SPK_EXTFIX)


%- make a sort helper if any spike channels are getting output
%- making a sort-helper csv



tJackSplit.stimChan = tJackSplit.microStimPigtail; %- for helper output change microStimPigtail name to "stimChan"
tJackSplit.chanNoise = cell(height(tJackSplit),1);  %- initialize here
tJackSplit.globalNoise = cell(height(tJackSplit),1);  %- initialize here

varKeep   = {'FileName'  'ChanName' 'ChanNameNew' 'stimChan' 'chanNoise' 'globalNoise' 'SortChanName'  };
tSortHelp = tJackSplit(:,varKeep);

% SJ changed dummySort to chanNoise and globalNoise
gradeStr  = '???';
%gradeStr  = 'need grade A-F';
%noteStr   = 'need comment on sort';
noteStr   = '-';
toSortStr = '-';
whichSrt  = 'reref'; %- always reref for non-stim session... exclude the column below
%if isStimSess, whichSrt = 'noreref(stimMask5)'; end
for iT=1:height(tSortHelp)
    tSortHelp.GradeNPlay1{iT} = gradeStr;
    tSortHelp.GradeNPlay2{iT} = gradeStr;
    tSortHelp.toSort{iT}      = toSortStr;
    %if isStimSess, tSortHelp.SortRefUsed{iT} = whichSrt; end %- only include this column for stim sesssions... non  stim always uses global
    tSortHelp.SortComment{iT} = noteStr;
    
    
    chanNoise = ' ';  if tSortHelp.stimChan(iT)>0, chanNoise = '???'; end
    globalNoise = ' ';
    tSortHelp.chanNoise{iT}   = chanNoise;
    tSortHelp.globalNoise{iT}   = globalNoise;
end

%varKeep2   = {'FileName'  'ChanName' 'ChanNameNew' 'stimChan' 'dummySort' 'SortChanName'  'GradeNPlay1' 'GradeNPlay2' 'toSort' 'SortRefUsed'};


%- see if a versio of this file exists, if so, dont make a new one
%sortRootName = fullfile(processedDirUse, sprintf('sortNotes(%s)_sortedBy', thisSess));
sortRootName = fullfile(processedDirUse, sprintf('sortNotes(%s)_sortedBy', thisSess));
dSortNotes = dir([sortRootName '*']);
if isempty(dSortNotes) || any(contains({dSortNotes.name},'??.xlsx'))  %- just the original "??" version... might as well overwrite
    notSortedYet = 1;
else
    notSortedYet = 0;
end


%- 
%OVERWRITE_NOTES = 1; %- overwrite all... if already a named version in there this will create a ?? version
%OVERWRITE_NOTES = isStimSess; %- pushing an update only to the stim sessions
%OVERWRITE_NOTES = isStimSess & notSortedYet;  %- this is part of the sortParty... just fix the ones that haven't been fixed
OVERWRITE_NOTES = notSortedYet;


if ~isempty(dSortNotes) & ~OVERWRITE_NOTES
    fprintf('\n found sort notes %s... not making a fresh one', dSortNotes(1).name);
else
    sortRootName = [sortRootName '??'];
    writetable(tSortHelp,[sortRootName '.xlsx']);
    fprintf('\n new sort notes file written to processed folder %s');
    
%     fid = fopen(fullfile(processedDirUse,'sortNotes_readme.txt'),'w+');
%     fprintf(fid,'\n sortNotes.xlsx readme');
%     fprintf(fid,'\n\n   1) rename the xlsx file so the initials of the sorter are included');
%     fprintf(fid,'\n\n   2) open the xlsx file and manually clear the "???" from the channels to be sorted (e.g., all utahs, but not micros).');
%     fprintf(fid,'\n\n   3) before manually sorting a session, use Windows/NPlay to view the ns5/ns6 and judge existence and quality of units');
%     fprintf(fid,  '\n     this judgement should be made twice per channel, once ~10% into the file (GradeNplay1), and once ~90% into the file (GradeNplay2)'); 
%     fprintf(fid,  '\n     for each channel, write a letter grade of that judgement:');
%     fprintf(fid,  '\n        "A" = definitely one or more units with great SNR that should be easy to sort/isolate');
%     fprintf(fid,  '\n        "B" = definitely one or more units with decent SNR. Definitely a unit, but might not be super easy to isolate');
%     fprintf(fid,  '\n        "C" = very likely one or more units. Best SNR is decent to poor.  Likely will get something in plexon, but poor isolation.');
%     fprintf(fid,  '\n        "D" = possibly one or more units. could be worth checking in plexon');
%     fprintf(fid,  '\n        "F" or Blank = definitely NO units, dont bother sorting. ');
%     fprintf(fid,'\n\n    4) when actually doing the sorts, make notes in the rightmost column, such as "3 clear units" "noise corrupted, cant sort", etc');
%     fprintf(fid,  '\n        for stimSession, update the "sortRefUsed" column if something other than the default (listed) reference folder was used for sorting');
%     fprintf(fid,  '\n        sort all the As, and then all the Bs, and then Cs.  If Cs were promising also sort Ds.  ');
%     fprintf(fid,  '\n        "toSort" column can be updated to help identify which ones to sort on each pass. Use the highst of the 2 grades for each channel.');
%     fprintf(fid,'\n\n        BE SURE TO NOTE in the comments section WHAT GRADE LETTER SORTING STOPPED ON. ');
%     fprintf(fid,  '\n           THIS IS HOW WE KEEP TRACK OF WHETHER MORE SORTS ARE POSSIBLE!!  ');
%     fprintf(fid,'\n\n    a copy of this table will be attached to the spikeInfo.mat when that is created');
%     fclose(fid);
    writeSortNotes_readme(fullfile(processedDirUse,'sortNotes_readme.txt')); %SJ added
    
end

if PUSH_SPK_EXTFIX == 1
    
    fixSPKextractionInfo(subj, tJackSplit, thisSess, isStimSess, processedDirUse, rawDir, subjectInfo, fidOut, compMemoryGB);
    fprintf('%s\n','Ran fixSPKextractionInfo');
    
end
end
end
