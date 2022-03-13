%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   JW's calls to eeg_toolbox utah processing and output functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

dbstop if error

% SJ CHANGE TEST

%- WHERE the micro data will be processed (i.e., split from ns5/6 to .bin; and where sorts are collected and converted to spikeInfo)
%      Intent is for this to be on 56PROC/micro_pristine, which is *not* shared to FRNU
%process_micro_dir = ''; %- blank means processing happens where the original data is stored, e.g., 56A/B/C/D.. we DONT want to do this anymore
process_micro_dir = '/Volumes/56PROC/micro_behavioral/micro_pristine'; %NORMAL
%process_micro_dir = '/Users/wittigj/Desktop/localData/micro_pristineTest';
%process_micro_dir = '/Volumes/Seagate Backup Plus Drive/local56/56PROC/micro_behavioral/micro_pristine';
%process_micro_dir = '/Volumes/56PROC/eeg_processing/.update/tempUpdate/sam_micro_align_update/micro_behavioral/micro_pristine';

rootMountainDir = '/Volumes/56PROC/micro_biowulfSorts/mountain_sorts';

%- WHERE sort.txt files are saved by plexon (i.e., intermingled with the .bin files).
%     Intent is for this to be on 56PUB/readWrite/micro_toSort, which is shared to FRNU
%sorted_micro_dir  = process_micro_dir;  %- use this line if sorts are in the same place as processed
sorted_micro_dir = '/Volumes/56PUB/readWrite/micro_forSorting'; %NORMAL
%sorted_micro_dir  = '/Users/wittigj/Desktop/localData/micro_forSorting/';  %- t
%sorted_micro_dir = '/Volumes/SeagateBackupPlusDrive/local56/56PUB/readWrite/micro_forSorting';


ecogDir = '/Volumes/56PROC/eeg_processing/.update/eeg'; %NORMAL %- potential location for offline processing... dont mess with this unless you know what you are doing %NORMAL
%ecogDir = '/Volumes/EEG_56STAGE-1/stage/working/';          %- use this for real-time processing
%ecogDir = '/Volumes/Seagate Backup Plus Drive/local56/56PROC/eeg_processing/.update/eeg';
%ecogDir = '/Volumes/56PROC/eeg_processing/.update/tempUpdate/sam_micro_align_update/eeg';

%- THINGS TO DO TO EACH SUBJECT: 
COPY_dataRaw_2_Micro = 0;  %- look at the rawDirList_NIHXXX.xlsx on /Volumes/56X/UTAH_X/subj/notes... copy any sessions with pulses to "process_micro_dir"
COPY_mnt_2_micro     = 0;
UPDATE_microCSV      = 1;  %- folders
EXTRACT_Microphys    = 0;  JUST_LFP = 0;
MAKE_LFP_PROCESSED   = 0;  %- SJ to do; 
COPY_2_PUBforSorting = 0;  %- ***for Sam***: copy any requested sessions from process_micro_dir to sorted_micro_dir (which includes some renaming/trimming of stuff)
COPY_sorts2proc      = 0;  %- just copy the sorts from "sorted_micro_dir" to "process_micro_dir"... dont necessarily make the spike info
EXTRACT_spikeInfo    = 0;  JUST_UPDATE_SUMMARY = 0;             %- copy any new sorts and make the spike infos, ***edit Sam***
COPY_lfpAndSorts2EEG = 0;  %- ***for Sam***: copy the LFPs and SpikeInfos to EEG/.update. Also copy sort figures back to sorted_micro_dir. Eventually this function will also copy from micro_biowulfSort folder. 
ALIGN_2_MACRO        = 0;  FORCE_ALIGN_ALL = true;  JUST_SPIKEINFO_SESS = 0; ONLY_TRC = true; %- run the macro-micro alignment script (w/in .update)

% No option to do JUST EXTRACT_spikeInfo -> always does COPY_sorts2proc
% first!
% noise redone: 37, 47, 57, 86, 88
%- WHICH SUBJECTS:
% NEED TO DO 51!!!!!!!!
% 37: 160129_1453_beh, 160209_1519_beh (let it try to align with 2556 drift), 160131_1232_beh, 160213_1006_beh
% 39: 160313_1250_beh, 160313_1417_beh 
% 42: 160627_1448_beh, 160613_1422, 160617_1526_beh, 160622_1329_beh, 160622_1529, 160623_1145_beh, 160626_1039_beh, 160627_1200_beh, 160618_1601_beh, 160618_1528
% 36: 151220_1754_beh, 151219_1509_beh, 151211_1524, 151219_1301_beh, 151221_1606_beh
% 34: 150925_1533_beh, 150926_1133_beh
% 60: 180628_1155: no alignment, 180628_1537: no alignment, 180630_1618: no alignment & max(sync2) = 0 for dc09 and dc12
% 62: stimMapGUI/session_24_CLIN_180924_1635: ERROR!!!! DC09 and DC12 are different lengths!!! Only calculating the duration for DC09! (no micro sessions in the end) (also sess 25, stimguess session 0_v5a)
% 64: 180930_1113_beh (very small additional zero needed maybe)
% 50: several sessions have din utah in addition to ainp micro & utah; 170625_1211_beh: alignment failed; 170712_1318_beh: drift is 104!! (usually like 30)
% Issue with zeroing for one in 52
%subjList = [50 54 57 59 62 64]; %RK process sorts **Wait for more sessions in 50 and 62
%subjList = [29:66 68:82];
%subjList = [64]; %54, 57, 59, 62 look messed up, 50 needs copy to eeg and aligned
%subjList = [50 54 57 59]; %RK batch #2 after corrupted files
%subjList = [62]; %RK Need to start from extract spike info after fixing 1 spike units
%subjList = [29 30 34 37 39 42 47]; %PAL sessions sorted by JZ, KS, NF
%subjList = [29 30 39 42 57 86 88]; %PAL sessions sorted by JZ, KS, NF, SJ (+88 for JZ); next: 37, 46, 50
%subjList = [37 46 50]; % next: 47, 86, 88 (JW's sorts)
%subjList = [47 86 88];
subjList = [57 60 62]; %UPDATE CSV!!!
%subjList = [60 62];

for iSubj = 1:length(subjList)
    subj = sprintf('NIH%03d',subjList(iSubj));
    fprintf('\n\n>>> CHECKING %s <<<<<',subj);
    
    subjectInfo = getMicroSubjInfo_v11(subj,process_micro_dir);
    if isempty(subjectInfo) || isempty(subjectInfo.dataPath56) % || ~exist(fullfile(overrideDataPath56,subj)),
        fprintf('  --> no subj info found, skipping');
        continue;
    else
        if isempty(process_micro_dir)
            process_micro_dir_use = subjectInfo.dataPath56; %- if no processing path was specified above, do the processing in the 56A/B/C/D folder
        else
            process_micro_dir_use = process_micro_dir; %- if a root path was specified above, use this.
        end
    end
    
    
    %- Copy files from 56A/B/C/E --> 56PUB for processing 
    % Do you mean --> 56PROC?? - Samantha
    if COPY_dataRaw_2_Micro
        rawPath56 = getMicroSubjInfo_v11(subj,'');
        rawPath56 = rawPath56.dataPath56;
        processEEGdir = fullfile(process_micro_dir,subj,'data_raw'); %- should be pointing to 56PUB/public/micro or something like that
        sessList = {};  JUST_CHECK = 0;
        getBR_grabFileList(subj,processEEGdir,rawPath56, sessList, JUST_CHECK);
    end
    
    
    %- last chance to jump out of this loop... does the subject folder exist?
    if ~exist(fullfile(process_micro_dir,subj),'dir')
        fprintf('\n Subject %s not found in %s; skipping to next subj', subj,process_micro_dir);
        continue;
    end
    
    if COPY_mnt_2_micro
        done_mnt = copyMountainSorts2EEG(subj,rootMountainDir,ecogDir);
    end
    
    
    %- Update/Create the UtahInfoCSV.  doesn't hurt to do this every time
    if UPDATE_microCSV
        SKIP_KEYBOARD    = 0; %- only set to 1 if you have already verified no errors
        Copy2UpdateandPUB = true;
        updateUtahInfoCSV_v10(subj, process_micro_dir_use, Copy2UpdateandPUB, SKIP_KEYBOARD); %- automatically marks sessions with pulses for extraction
    end
    
    
    %- the WORKHORSE... very slow step
    if EXTRACT_Microphys
        overrideOutputType   = ''; %- spike and LFP
        if JUST_LFP
            overrideOutputType = 'justLFP'; %- should be a lot faster if dont need to save all the spike files
        end
        extractMicrophysiology_v12b(subj,overrideOutputType,process_micro_dir_use);  %- automatically calls maskStimFromSpikes
    end
    
    
    %- MASK_SPIKES is automatically called by extractMicroPhys... use this separate call if you want to try a different msGrowDC10pulse or if a stimSession was not properly tagged
    %msGrowDC10pulse = 20;  testOnSomeChans = 10;  pickSessDir = '190712_0911';  pickSessDir = '';
    %maskStimFromSpikes_v2(subj,overrideDataPath56, msGrowDC10pulse, testOnSomeChans, pickSessDir);
    
    if COPY_2_PUBforSorting == 1
        done = copy2PUB_forSort(subj,process_micro_dir_use,sorted_micro_dir);
    end
    
    %- Create SpikeInfo... medium slow step (~10-20min)
    if EXTRACT_spikeInfo || COPY_sorts2proc
        %sessFolders  = {'190716_1618'}; %- specific session(s)
        sessFolders  = ''; %all sessions in data_processed
        grabSortedData_v2b(subj,process_micro_dir_use,sorted_micro_dir,sessFolders); %- probably should make this an automatic part of extractSpikeInfo...
        
        if EXTRACT_spikeInfo
            if JUST_UPDATE_SUMMARY == 1
                extractSpikeInfo_v3b(subj,process_micro_dir_use,sessFolders,1);
            else
                extractSpikeInfo_v3b(subj,process_micro_dir_use,sessFolders,0);
            end
        end
    end
    
    %- Copy lfp and sorts data too eeg update, then copy sort summary files
    %   to PUB
    if COPY_lfpAndSorts2EEG
        done = copy2EEG_LFP_Sorts(subj,process_micro_dir,ecogDir,sorted_micro_dir);
    end
    
    %- Incorporate LFP and Spikes with EEG database
    if ALIGN_2_MACRO
        
        %ecogDir = '/Volumes/EEG_56STAGE-1/stage/working/';          %- use this for real-time processing
        %ecogDir = '/Volumes/56PROC/eeg_processing/.update/eeg';  %- potential location for offline processing... dont mess with this unless you know what you are doing
        
        
        %- warning... only JW should do this step (until kinks are worked out or until this is pointed to a private resource, not the EEG backup
        if iSubj==1
            fprintf('\n\n WARNING: you are about to mess with the EEG directory: %s. \n For now only SJ/JW should do this.  Please break if not SJ or JW\n',ecogDir);
            keyboard;
        end
        
        
        %- check.. do we even see the EEG directory?
        if ~exist(fullfile(ecogDir,subj),'dir'), fprintf('\n cant find subject folder... check eCogDir path');
            keyboard;
            continue;
        end
        
        
        %- Now confirm the number of processed sesions we expect to find from 56E/micro directory...
        processEEGdir = fullfile(process_micro_dir_use,subj,'processed_SPK'); %- should be pointing to 56PUB/public/micro or something like that
        manSortProcessed = dir(processEEGdir);
        manSortProcessed = manSortProcessed(~strncmp({manSortProcessed.name},'.',1)); %- cut out '.' '..' '.DStore'
        if length(manSortProcessed)==0
            fprintf('\n WARNING: No processed folders in 56E/micro, probably shouldnt bother checking subject...');
        end
        
        
        %- Define the eCog path and confirm that micro/manualsorts folder is there
        manSortPath = fullfile(ecogDir,subj,'micro',filesep,'manualsort');
        manualSessFound = dir(manSortPath);
        manualSessFound = manualSessFound(~strncmp({manualSessFound.name},'.',1)); %- cut out '.' '..' '.DStore'
        if length(manualSessFound)==0
            fprintf('\n WARNING: no sessions found in %s \n skipping this subject');
            keyboard;
            continue;
        elseif length(manualSessFound) < length(manSortProcessed)
            fprintf('\n WARNING: missing some processed sessions in the EEG micro directory, should copy those over before aligning');
            keyboard;
            error('\n copy first!');  %- not really an error... in the future put the automatic copying code here
        end
        
        
        %- run alignment
        fprintf('\n FOUND %d sessions in %s. \n Attempting to align Micro and Macro now',length(manualSessFound),manSortPath);
        %keyboard %Switch the below statement to false after running test
        %[alignedFiles] = alignSpikesWithEcog_ManualORAuto(subj,ecogDir,true);
        [alignedFiles] = alignSpikesWithEcog_ManualORAuto(subj,ecogDir,FORCE_ALIGN_ALL,'onlyTRC',ONLY_TRC);
    end
    
end %- loop over subjects


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Alignment with eCog:
%     1)  if you dont even know which micro files to grab:  raw_fileList_findMatchingMicros.m.
%
%     2)  once the micros LFP/spikes are in a subject folder/micro, then run to align:
%              alignSpikesWithEcog_ManualORAuto
%              [alignedFiles] = alignSpikesWithEcog_ManualORAuto(subj,ecogDir)
%
%open alignSpikesWithEcog_ManualORAuto

if 0
    
    %- copy LFP and spikeInfo to eCog/micro directory
    
    subj = 'NIH076';
    subject = subj;
    ecogDir = '/Volumes/EEG_56STAGE-1/stage/working/'; %- root with subjects in it
    if ~exist(fullfile(ecogDir,subject),'dir'), fprintf('\n cant find subject folder... check eCogDir path'); keyboard; end
    [alignedFiles] = alignSpikesWithEcog_ManualORAuto(subj,ecogDir)
end

fprintf('\n\n%s\n','------------------------------------ ALL DONE ------------------------------------');

