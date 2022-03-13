function allTags = split_BR(subj, rootEEGdir, EEG_file, varargin)

% split_br - Splits a br .ns datafile into separate channels into
% the specified directory.
%
% FUNCTION:
%    split_br(subj,br_dir,EEG_file,...)
%
% INPUT ARGs:
% subj = 'TJ022'
% rootEEGdir = '/Users/damerasr/damerasr/data/eeg/'
% EEG_file = '/Users/damerasr/damerasr/data/eeg/160614_1008/x.EEG'
%
% Output
%	output is saved out to: output_dir = '/Users/damerasr/damerasr/data/eeg/[SUBJ]/eeg.noreref/' set in line 49
%
%
% Edited from previous version so that manual input of the ouput_dir would not be necessary
%
% 12/2017 Melkalliny - adapted the goals of nk_split this blackrock version
% 2/6/2020 SJ - fixed issue where align_nsps would use DC09 for NSP alignment but then pulsealign would
%               be given DC12 as the selected channel for alignment
%             - performance/readability improvement
%             - Added in 'end' for is useNEV == 0 where it was missing (was it running incorrectly
%               before...?)
% 3/2/2020 SJ - Implement check for jacksheetBR that ensures the channel names in ChanName and
%               ChanNameNew are correct, and gives user the option to automatically fix the naming issue.
% 3/9/2020 SJ - Get rid of call to transformSync- this is all already done in align_nsps anyway. Have
%               alignmentStats come from output of align_nsps instead of transformSync
% 5/11/2020 SJ- Finish fix from 3/2/20- make call to makejacksheetBR have overwrite = -1 (which is a fix
%               I just added to makejacksheetBR to prevent it from grabbing the local jacksheet already
%               in the folder)
% 7/2020   SJ - Changed readtable to readtableSafe
% 8/2020   SJ - removed signalFs from call to align_nsps (should be same as syncFs, also align_nsps only accepts 1 input now); Also check to make sure syncFs == signalFs
% 9/2020   SJ - Major update:
%                - Get rid of reliance on ChanNames.txt; if it still exists, checks against new process
%                (same with NSP_ChanCounts.txt)
% 9/14/2020 SJ- Put in fix for not grabbing ain17-20 from jacksheetBR
% 9/17/2020 SJ- Put in extra fix that solves jschans initialization issue
% 9/18/2020 SJ- (committed from Jiali's account)- Put in deblank in case jacksheetBR still has spaces
% 9/29/2020 SJ- Centralizing place where channel names are extracted from jacksheetBR/CV- now calling new function getChansFromJack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VERBOSE = 0 ; %[0,1] = output info about each channel's remapping

% load the jacksheetMaster, or create it, or give a warning that it can't be created...
subjDir             = fullfile(rootEEGdir, subj);
rawDir              = fullfile(subjDir, 'raw');
[rawSubDir,~,~]     = fileparts(EEG_file); %- path to the raw file being processed, so FRNU/data/eeg/NIHXXX/raw/121401_1415/
jackMaster_file_new = fullfile(subjDir, 'docs/jacksheetMaster.csv');

sessTitle           = EEG_file(strfind(EEG_file,subj):end);


% should we use analog or digital pulses?  Check jacksheetBR
jacksheetBR = dir(fullfile(rawSubDir,'jacksheetBR_local*.csv'));  %- why use the * in this... isn't name always the same?
if isempty(jacksheetBR)
    jacksheetBR = makeJacksheetBR(rawSubDir,'',0,1,0);
else %SJ - added the following to ensure that jacksheetBR has no changes, and if it does, the changes were made correctly.
    jacksheetBR_name = fullfile(rawSubDir,jacksheetBR(1).name);
    jacksheetBR = readtableSafe(fullfile(rawSubDir,jacksheetBR(1).name)); %SJ changed from readtable
    % Still need to make jacksheetBR to compare to old one
    jacksheetBR_new_comp = makeJacksheetBR(rawSubDir,'',-1,1,0); %Not writing out
    jacksheetBR_new_comp.ChanName = deblank(jacksheetBR_new_comp.ChanName); %SJ in case blanks still exist
    if ~isequaln(jacksheetBR,jacksheetBR_new_comp)
        % If they are not equal, this means someone changed something. Usually the 2 options are:
        % 1) Added channel to ChanNameNew (proper thing to do)
        % 2) Changed the name in ChanName instead of ChanNameNew - in this case, offer them the ability
        %    to automatically move the channel names to ChanNameNew and replace ChanName with the originals
        nonNameCols_old = jacksheetBR;
        nonNameCols_old.ChanName = []; nonNameCols_old.ChanNameNew = [];
        nonNameCols_new = jacksheetBR_new_comp;
        nonNameCols_new.ChanName = []; nonNameCols_new.ChanNameNew = [];
        
        if ~isequaln(nonNameCols_new,nonNameCols_old)
            fprintf('\n%s\n','ERROR!!! There have been changes to the non-channel name columns of jacksheetBR!!! WHY?? Ask SJ!');
            keyboard
        end
        
        %Ok, now just check the naming columns
        chanName_old = jacksheetBR.ChanName;
        chanNameNew_old = jacksheetBR.ChanNameNew;
        chanName_new = jacksheetBR_new_comp.ChanName;
        chanNameNew_new = jacksheetBR_new_comp.ChanNameNew;
        
        if isequaln(chanName_old,chanName_new) %These should always be the same!
            if isequaln(chanNameNew_old,chanNameNew_new)
                fprintf('\n%s\n','ERROR!!!! How do we get here?? Ask SJ!');
                keyboard
            else
                %This is okay to get to. This means the original chanName has stayed the same, but
                %chanNameNew has been updated (as it should)
                fprintf('\n%s\n','Detecting changes to the ChanNameNew column (which is how changes should be made). Keeping these changes.');
            end
        else %This means chanName was updated instead of just chanNameNew
            fprintf('\n%s\n','Uh oh. Looks like ChanName was edited, which is should not be!');
            diff_old = setdiff(chanName_old,chanName_new);
            diff_new = setdiff(chanName_new,chanName_old);
            fprintf('%s','Conflicting channels in current jacksheetBR.ChanName: ');
            fprintf('%s ',diff_old{:});
            fprintf('\n%s','What the channels should be in jacksheetBR.ChanName: ');
            fprintf('%s \n',diff_new{:});
            yn_in = input('Would you like to change jacksheetBR.ChanName back to the what they should be, and instead put the changes in the ChanNameNew column? (Y/N) ','s');
            if strcmpi(yn_in,'Y')
                idx_diff = find(ismember(jacksheetBR.ChanName,diff_old));
                if ~all(contains(jacksheetBR.ChanNameNew(idx_diff),'-'))
                    fprintf('%s\n','ERROR!!! ChanNameNew is not empty!!! What is going on? ChanNameNew contents: ');
                    fprintf('%s ',jacksheetBR.ChanNameNew{idx_diff});
                    fprintf('\n%s\n','Continue if you want to replace these anyway.');
                    keyboard
                end
                % All ChanNameNew spots are '-', so we can go ahead and replace!
                jacksheetBR.ChanName(idx_diff) = diff_new;
                jacksheetBR.ChanNameNew(idx_diff) = diff_old;
                
                if isequaln(jacksheetBR.ChanName,jacksheetBR_new_comp.ChanName)
                    fprintf('%s\n','Deleting old file and re-writing.')
                    delete(jacksheetBR_name)
                    writetable(jacksheetBR,jacksheetBR_name)
                else
                    fprintf('\n%s\n','ERROR!!! The ChanName Column still is not correct! Ask SJ!');
                    keyboard;
                end
            else
                fprintf('%s\n','Not correcting jacksheetBR, so still an issue!!!');
            end
        end
    end
end


%- ain1  = DC09 = ainp1 = PHYSICAL_CHAN 1001
%- ain1  = DC12 = ainp4 = PHYSICAL_CHAN 1004
%- ain17 = DC09 on nsp2 = ainp17 = PHYSCIAL_CHAN 1017
iDCchans = find(ismember(jacksheetBR.PhysicalChan,[1001 1004 1017 1020])); %- get row indicies for DC09 and DC12 (either NSP).
fileListDC = unique(jacksheetBR.FileName(iDCchans));
if sum(contains(fileListDC,'nev'))<length(fileListDC)
    thisFile    = fileListDC(~contains(fileListDC,'nev'));
    pulseString = thisFile{1}(end-2:end);
    useNEV = 0;
else
    thisFile    = fileListDC;
    pulseString = 'nev';
    useNEV      = 1;
    fprintf('\nusing digital (NEV) pulses - no nsX present!\n')
end
if ~exist(fullfile(rawSubDir,thisFile{1}),'file')
    fprintf('\n dont see the sync file we are expecting... perhaps because jacksheetBR is complete not local?');
    keyboard;
end


%- are there two NSPs?
doubleNSP = 0;
rawFileList = dir(fullfileEEG(rawSubDir,'*.ns2'));
if length(rawFileList)>1
    doubleNSP = 1;
else
    pulses_to_use = 1;
end


% create jacksheetMaster if it doesn't exist or hasn't looked at raws
if ~exist(jackMaster_file_new, 'file')
    createMasterJack(subj, rootEEGdir);
end

jacktable = getJackTable(subj, rootEEGdir);
jackMaster_names = jacktable.chanName;
jackMaster_chans = [1:length(jackMaster_names)];


% read the raw_info file to identify channels intentionally EXCLUDED from jacksheetmaster... no reason to give a warning if those are found below
raw_info_file = fullfile(subjDir, 'docs/raw_info.csv');
raw_info      = readtableSafe(raw_info_file);
chanDontSplit = raw_info.chanName(find(raw_info.in_jackSheetMaster==0));
chanDontSplit{end+1} = ''; %- dont report error or try to split channels with no name


% 1. define and create directories
outputDir = fullfile(rootEEGdir,'/',subj,'/eeg.noreref/'); %define noreref directory
if ~exist(outputDir)
    mkdir(outputDir); disp('eeg.noreref directory has been created');
end
folderOutDir = strrep(EEG_file,'raw','eeg.noreref'); % for split files
folderOutDir = strrep(folderOutDir,'STIM_MAP',''); % for raws from stim_map
temp = strfind(folderOutDir,'/'); folderOutDir(temp(end):end) = [];

if ~exist(folderOutDir); mkdir(folderOutDir);  end


% before extracting raw data, rearrange data as needed and run alignment
% between two files (IF there are two NSPs of iEEG data). if you're on
% the first NSP, use the output errorStr to summarize alignment in
% a .txt file in the raw folder.

if doubleNSP==1
    if contains(EEG_file,'ieeg1')
        signal1filepath = EEG_file;
        signal2filepath = strrep(EEG_file,'ieeg1','ieeg2');
        % ^ not valid anymore. could be ieeg2_utah
        sync2filepath = strrep(signal2filepath,'ns2',pulseString);
        sync1filepath = strrep(sync2filepath,'ieeg2','ieeg1');
        
        [sync, syncFs] = getBR_sync(sync1filepath,sync2filepath,signal1filepath);
        
        signal1 = concatOpenNSx(signal1filepath,0);  signalFs = signal1.MetaTags.SamplingFreq;
        signal1 = signal1.Data;
        signal2 = concatOpenNSx(signal2filepath,0);
        signal2 = signal2.Data;
        
        DC09or12 = 12; % always ask align_nsps to use DC12, when its blackrock
        
        if syncFs ~= signalFs %SJ
            fprintf('\n%s\n','STOP!! syncFs and signalFs are not equal, but they need to be because align_nsps only accepts one value!!');
            fprintf('\t%s\n',['syncFs = ' num2str(syncFs) ';  signalFs = ' num2str(signalFs)]);
            fprintf('%s\n','Go get SJ to come investigate!');
            keyboard
        end
            
        %[signal, pulses_to_use, syncStartIndex,errorStr, transforms] = align_nsps(sync, signal1, signal2, [syncFs signalFs], 'BR', DC09or12,rawSubDir,sessTitle,0);
        [signal, pulses_to_use, syncStartIndex,errorStr, transforms] = align_nsps(sync, signal1, signal2, syncFs, 'BR', DC09or12,rawSubDir,sessTitle,0); %SJ 8/28/2020: removed signalFs (should be same as syncFs, also align_nsps only accepts 1 input now)
        
        %- used to write alignment summary here...
        
        
    elseif contains(EEG_file,'ieeg2')
        signal2filepath = EEG_file;
        signal1filepath = strrep(EEG_file,'ieeg2','ieeg1');
        % ^ not valid anymore. could be ieeg2_utah
        sync1filepath = strrep(signal1filepath,'ns2',pulseString);
        sync2filepath = strrep(sync1filepath,'ieeg1','ieeg2');
        
        [sync, syncFs] = getBR_sync(sync1filepath,sync2filepath,signal1filepath);
        
        signal1 = concatOpenNSx(signal1filepath,0);  signalFs = signal1.MetaTags.SamplingFreq;
        signal1 = signal1.Data;
        signal2 = concatOpenNSx(signal2filepath,0);
        signal2 = signal2.Data;
        
        DC09or12 = 12; % always ask align_nsps to use DC12, when its blackrock
        
        if syncFs ~= signalFs %SJ
            fprintf('\n%s\n','STOP!! syncFs and signalFs are not equal, but they need to be because align_nsps only accepts one value!!');
            fprintf('\t%s\n',['syncFs = ' num2str(syncFs) ';  signalFs = ' num2str(signalFs)]);
            fprintf('%s\n','Go get SJ to come investigate!');
            keyboard
        end
        
        %[signal, pulses_to_use, syncStartIndex,~, ~] = align_nsps(sync, signal1, signal2, [syncFs signalFs], 'BR', DC09or12,rawSubDir,sessTitle,0);
        [signal, pulses_to_use, syncStartIndex,~, ~] = align_nsps(sync, signal1, signal2, syncFs, 'BR', DC09or12,rawSubDir,sessTitle,0); %SJ 8/28/2020: removed signalFs (should be same as syncFs, also align_nsps only accepts 1 input now)
        
    elseif contains(EEG_file,'INST1')
        signal2filepath = EEG_file;
        signal1filepath = strrep(EEG_file,'INST1','INST0');
        sync1filepath = strrep(signal1filepath,'ns2',pulseString);
        sync2filepath = strrep(sync1filepath,'INST0','INST1');
        
        [sync, syncFs] = getBR_sync(sync1filepath,sync2filepath,signal1filepath);
        
        signal1 = concatOpenNSx(signal1filepath,0);  signalFs = signal1.MetaTags.SamplingFreq;
        signal1 = signal1.Data;
        signal2 = concatOpenNSx(signal2filepath,0);
        signal2 = signal2.Data;
        
        DC09or12 = 12; % always ask align_nsps to use DC12, when its blackrock
        
        if syncFs ~= signalFs %SJ
            fprintf('\n%s\n','STOP!! syncFs and signalFs are not equal, but they need to be because align_nsps only accepts one value!!');
            fprintf('\t%s\n',['syncFs = ' num2str(syncFs) ';  signalFs = ' num2str(signalFs)]);
            fprintf('%s\n','Go get SJ to come investigate!');
            keyboard
        end
        
        %[signal, pulses_to_use, syncStartIndex, errorStr, transforms] = align_nsps(sync, signal1, signal2, [syncFs signalFs], 'BR', DC09or12,rawSubDir,sessTitle,0);
        [signal, pulses_to_use, syncStartIndex, errorStr, transforms] = align_nsps(sync, signal1, signal2, syncFs, 'BR', DC09or12,rawSubDir,sessTitle,0); %SJ 8/28/2020: removed signalFs (should be same as syncFs, also align_nsps only accepts 1 input now)
        
        
    elseif contains(EEG_file,'INST0')
        signal1filepath = EEG_file;
        signal2filepath = strrep(EEG_file,'INST0','INST1');
        sync2filepath = strrep(signal2filepath,'ns2',pulseString);
        sync1filepath = strrep(sync2filepath,'INST1','INST0');
        
        [sync, syncFs] = getBR_sync(sync1filepath,sync2filepath,signal1filepath);
        
        signal1 = concatOpenNSx(signal1filepath,0);  signalFs = signal1.MetaTags.SamplingFreq;
        signal1 = signal1.Data;
        signal2 = concatOpenNSx(signal2filepath,0);
        signal2 = signal2.Data;
        
        DC09or12 = 12; % always ask align_nsps to use DC12, when its blackrock
        
        if syncFs ~= signalFs %SJ
            fprintf('\n%s\n','STOP!! syncFs and signalFs are not equal, but they need to be because align_nsps only accepts one value!!');
            fprintf('\t%s\n',['syncFs = ' num2str(syncFs) ';  signalFs = ' num2str(signalFs)]);
            fprintf('%s\n','Go get SJ to come investigate!');
            keyboard
        end

        %[signal, pulses_to_use, syncStartIndex,errorStr, transforms] = align_nsps(sync, signal1, signal2, [syncFs signalFs], 'BR', DC09or12,rawSubDir,sessTitle,0);
        [signal, pulses_to_use, syncStartIndex,errorStr, transforms] = align_nsps(sync, signal1, signal2, syncFs, 'BR', DC09or12,rawSubDir,sessTitle,0); %SJ 8/28/2020: removed signalFs (should be same as syncFs, also align_nsps only accepts 1 input now)
        
        %- used to write alignment summary here...
        
    end
    
    
    % separate the output into the two NSPs of data again
    signal1 = signal{1,1}; signal2 = signal{1,2};
    
    checkPulses = 1;
    if checkPulses==1 && exist('errorStr') && exist('transforms') && errorStr~=1 && errorStr~=4%#ok<*EXIST>
        % transform original sync using transforms, then run get_triggers on it
        
        % SJ - reassign DC09or12 based on the errorStr of align_nsps (ie what was identified for
        % alignment, DC09 or DC12)
        oldDC09or12 = DC09or12;
        if errorStr == 2
            DC09or12 = 12;
        elseif errorStr == 3
            DC09or12 = 9;
        end
        if DC09or12 ~= oldDC09or12
            fprintf('%s\n','Stopped on a case where DC09or12 was assigned differently than the output of align_nsps errorStr. Check it out and then continue.')
            keyboard
        end
            
        %SJ: Comment out the following because this now happens in align_nsps
%         pulseType = DC09or12;
%         transformedSync = transformSync(sync,transforms);
%         
%         
%         split_filename = sprintf('%s/DC%02d_nsp1_transformed',folderOutDir,9); 
%         [fchan,msg] = fopen(split_filename,'w','l');
%         assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
%         fwrite(fchan,transformedSync{1,1},'int16');
%         fclose(fchan);
%          split_filename = sprintf('%s/DC%02d_nsp1_transformed',folderOutDir,12); 
%         [fchan,msg] = fopen(split_filename,'w','l');
%         assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
%         fwrite(fchan,transformedSync{1,2},'int16');
%         fclose(fchan);       
%         
%         % since transformSync only transforms the first input (micro data),
%         % lets flip things around so that it transforms the ecog, too
%         flipTransforms = transforms;
%         flipSync = sync; flipSync{1,4} = flipSync{1,3}; 
%         flipSync{1,2} = flipSync{1,1};
%         flipTransforms{1,1} = -flipTransforms{1,1};
%         flipTransforms{1,2} = -flipTransforms{1,2};
%         transformedSync_nsp2 = transformSync(flipSync,flipTransforms);
%         split_filename = sprintf('%s/DC%02d_nsp2_transformed',folderOutDir,9); 
%         [fchan,msg] = fopen(split_filename,'w','l');
%         assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
%         fwrite(fchan,transformedSync_nsp2{1,1},'int16');
%         fclose(fchan);
%          split_filename = sprintf('%s/DC%02d_nsp2_transformed',folderOutDir,12); 
%         [fchan,msg] = fopen(split_filename,'w','l');
%         assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
%         fwrite(fchan,transformedSync_nsp2{1,2},'int16');
%         fclose(fchan);    
% 
%         if pulseType==12
%             micro_ms = transformedSync{1,2};
%             ecog_ms = sync{1,1};
%         elseif pulseType==9
%             micro_ms = transformedSync{1,1};
%             ecog_ms = sync{1,3};
%         end
%         triggers_ecog = get_triggers(ecog_ms',1000);
%         triggers_micro = get_triggers(micro_ms',1000);
%         threshMS = 10; mywin = 200; samplerate = 1000;
%         % beh_ms and pulses should be vector of doubles
%         [ecog_pulses_use, micro_pulses_use] = pulsealign(triggers_ecog{1,1},triggers_micro{1,1},samplerate,threshMS,mywin,0,0,1);
%         ms_field = 'mstime';
%         [alignInfo] = logalign_microVsEcog({ecog_pulses_use},{micro_pulses_use},ms_field);
%         writeOutAlignInfo = sprintf('%s//alignmentStats.txt',rawSubDir);
%         fid  = fopen(writeOutAlignInfo,'w');
%         if fid==-1; fprintf('cant open file..\n'); keyboard; end
%         fprintf(fid,'NumPoints=%0.4f\n',alignInfo.reg_numPointFit);
%         fprintf(fid,'Intercept=%0.4f\n',alignInfo.reg_intercept);
%         fprintf(fid,'Slope=%0.4f\n',alignInfo.reg_slope);
%         fprintf(fid,'R^2=%0.4f\n',alignInfo.reg_Rsquare);
%         fprintf(fid,'MaxDev=%0.4f\n',alignInfo.reg_maxDev);
%         fprintf(fid,'MedianDev=%0.4f\n',alignInfo.reg_medianDev);
%         fprintf(fid,'PulsesUsed=%d\n',pulseType);
%         fclose(fid);
     end
                            
else
    syncStartIndex=1;
end



%% extract data from raw files & resample


% v1- pre-figuring out details of junks+splits
%NSxdata = openNSx(EEG_file); % load the data
%NSxdata.Data = getNSxData(NSxdata); % pass it through a check for a resync

%% v2 - working version - use signal processed through double-NSP alignment, or
% if single-NSP, use updated code

if doubleNSP==1
    if (contains(EEG_file,'ieeg2') || contains(EEG_file,'INST1'))
        NSxdata = signal2; clear signal1;
    elseif (contains(EEG_file,'ieeg1') || contains(EEG_file,'INST0'))
        NSxdata = signal1; clear signal2
    end
    
else
    NSxdata = concatOpenNSx(EEG_file,0); % pass it through a check for a resync
    NSxdata = NSxdata.Data;
end


%channellabels={NSxdata.ElectrodesInfo.Label}; % instead , use .txt later on
metaData = openNSx(EEG_file,'noread');
samprate_orig=metaData.MetaTags.SamplingFreq; % find sample rate


%% do not resample any more - 11/2018 - melkalliny
samprate_new=1000;
% [fsorig, fsres] = rat(samprate_orig/samprate_new); % resample
% chan=cell(length(NSxdata.ElectrodesInfo),1);
% for channel=1:length(NSxdata.ElectrodesInfo)
%     chan{channel}=resample(double(NSxdata.Data(channel,:)),fsres,fsorig);
% end

% instead, just write out at native sampling rate
chan=cell(size(NSxdata,1),1);
for channel=1:size(NSxdata,1)
    chan{channel}=double(NSxdata(channel,:));
end


%SJ replacing the following commented code with this (Centralized pulling of chanNames from
%jacksheet_local):
[jschans, ~, nums] = getChansFromJack(jacksheetBR);

% 
% % Get the chanNameNews, but also get rid of the double ains
% %jschans = cell(sum(~contains(jacksheetBR.FileName,'.nev')),1);
% jschans = cell(sum(~contains(jacksheetBR.FileName,'.nev') & cellfun(@isempty,regexp(jacksheetBR.ChanName,'ain.*(17|18|19|20)'))),1);
% %SJ: get nums now by searching for the 1st NSP but not the pulse channels,
% %which should hopefully be at the end...
% nums = sum(jacksheetBR.NSP == 1 & jacksheetBR.PhysicalChan < 1000);
% ii_row = 0;
% for ii = 1:numel(jacksheetBR.ChanName)
%     if ~contains(jacksheetBR.FileName{ii},'.nev') && isempty(regexp(jacksheetBR.ChanName{ii},'ain.*(17|18|19|20)')) %SJ added to take care of .ns2 ain17-20s
%         ii_row = ii_row + 1;
%         if strcmp(jacksheetBR.ChanNameNew{ii},'-')
%             jschans{ii_row} = jacksheetBR.ChanName{ii};
%         else
%             jschans{ii_row} = jacksheetBR.ChanNameNew{ii};
%         end
%     end
% end


%- also make a copy of jacksheetBR_local so eeg.noreref has low-level info about electrode numbers
jackCSV = fullfile(folderOutDir,'/','jacksheetBR_localUsed.csv'); %- call eeg.noreref copy "localUsed" to deliniate the source file in raw from the copy in eeg.noreref
writetable(jacksheetBR,jackCSV);


% cut raw channels that were not specified in master
% to do so, pull chanNames from .txt made while aligning w/ jacksheet (eegModifyAll21Es.m)
temp = strfind(EEG_file,'/');
channelTxtPath = strcat(EEG_file(1:temp(end)),'ChanNames.txt');
if exist(channelTxtPath,'file')
    fprintf('\n%s','ChanNames.txt still exists. Checking against jacksheetBR to make sure we get the same result ...');
    fid = fopen(channelTxtPath,'r');
    if fid < 0
        fprintf('Should not happen... how does it  exist then? \n');
        keyboard
    end
    names_txt = textscan(fid,'%s');
    %namesFromBRCell = {};
    %namesFromBRCell{1} = names_txt{1};
    fclose(fid);
    names_txt = names_txt{1};
    
    if ~isequal(names_txt,jschans)
        fprintf('\n%s\n','ERROR!! Not the same! Get SJ');
        keyboard
    else
        fprintf('%s\n',' we do! Deleting chanNames.txt...');
        delete(channelTxtPath)
    end
    
end    

%- open the "nums" file... channel counts for each NSP.
channelNumPath = strcat(EEG_file(1:temp(end)),'NSP_ChanCounts.txt');
if exist(channelNumPath,'file')
    fprintf('\n%s','NSP_ChanCounts.txt still exists. Checking against jacksheetBR to make sure we get the same result ...');
    fid = fopen(channelNumPath,'r');
    nums_txt = (textscan(fid,'%d')); 
    nums_txt = double(nums_txt{1});
    fclose(fid);
    if nums_txt ~= nums
        fprintf('\n%s\n','ERROR!! Not the same! Get SJ');
        keyboard
    else
        fprintf('%s\n',' we do! Deleting NSP_ChanCounts.txt...');
        delete(channelNumPath)
    end
end

% make jacksheet.txt using chans from both NSPs
%channellabels = namesFromBRCell{1}';
channellabels = jschans';

jackFile = fullfile(folderOutDir,'/','jacksheet.txt');
fileOut  = fopen(jackFile,'w','l');
if fileOut==-1; error('Jacksheet output directory is not found.'); end
for channel=1:length(channellabels)
        fprintf(fileOut,'%d %s\n',channel,channellabels{channel});
end
fclose(fileOut);

% only pull the chans recorded by the NSP from which raw was pulled
if doubleNSP==1
    if contains(EEG_file,'ieeg1') || contains(EEG_file,'INST0')
        channellabels = channellabels(1:nums);
        %namesFromBRCell = namesFromBRCell{1}(1:nums(1));
    elseif contains(EEG_file,'ieeg2') || contains(EEG_file,'INST1')
        channellabels = channellabels(nums+1:end);
        %namesFromBRCell = namesFromBRCell{1}(nums(1)+1:end);
    end
else
    %namesFromBRCell = namesFromBRCell{1};
end
%channellabels = namesFromBRCell';

% exclude raw chans that were recDontUse in jacksheetMaster
%[~,temp] = intersect(namesFromBRCell,chanDontSplit);
[~,temp] = intersect(channellabels',chanDontSplit);
chan(temp(:)) = []; channellabels(temp(:)) = [];
%namesFromBRCell(temp(:)) = [];
% exclude raw chans that were otherwise not present in jackMaster
%[mismatches, index] = setdiff(namesFromBRCell,jackMaster_names,'stable');
[mismatches, index] = setdiff(channellabels',jackMaster_names,'stable');
if ~isempty(mismatches)
    fprintf(['\n  HEADS UP: raw channel list contains chans that jacksheetMaster.csv does not, likely because of recNotUsed.'...
        '\n Channels in raws but not jacksheetMaster are:\n']) ;
    disp(mismatches);
    
    %- set the channel data to empty so its not split out
    chan(index(:)) = []; channellabels(index(:)) = [];
    %namesFromBRCell(index(:)) = [];
end

% for now, don't delete the data
% clear NSxdata

% what is the gain
ns2MaxAnalog = metaData.ElectrodesInfo(1).MaxAnalogValue;
ns2MaxDigi   = metaData.ElectrodesInfo(1).MaxDigiValue;
gainEEGchannels = 1/((double(ns2MaxDigi) / double(ns2MaxAnalog)));
if gainEEGchannels ~= 0.25
    fprintf('gain is a value never encountered before.  \n')
    keyboard
end


% make params.txt
paramsFile = fullfile(folderOutDir,'/','params.txt');
fileOut    = fopen(paramsFile,'w','l');
if fileOut==-1; error('params output directory is not found.'); end
fprintf(fileOut,'samplerate %0.11f\n',samprate_orig);
fprintf(fileOut,'dataformat ''int16''\n');
fprintf(fileOut,'gain %d\n',gainEEGchannels); % default
fclose(fileOut);


% make sourcetype.txt (which system recorded raw)
sourceFile=fullfile(folderOutDir,'/','sourcetype.txt');
fileOut = fopen(sourceFile,'w','l');
fprintf(fileOut,'blackrock raw\n');
fclose(fileOut);



% get analog DC, resample, and save analog files + trigDC09, trigDC12
pulseDir = EEG_file; pulseDir = strrep(pulseDir,'ns2',pulseString);
if useNEV==1
    pulsemetaData = openNEV(pulseDir,'nosave'); % For analog DC channel
    fprintf('\n getBlackrockPulsesDC has changed... will the line below even work?  Ask JW to take a look');
    pause(5);
    fprintf('\n moving forward... but best to look into this more later (see split_BR, useNEV==1)');
    pause(5);
    %keyboard
    %pulseData = getBlackRockPulsesDC(pulsemetaData,9);
    [pulseData, output_ts, IPIViolationFlag] = getBlackRockPulsesDC(pulsemetaData,9);
    pulseData = double(pulseData)./30;
else
    pulsemetaData = openNSx(pulseDir); % For analog DC channel
    pulseData = concatOpenNSx(pulseDir,0);
end


if useNEV==0
    %- do some checks on the ain channels here
    allChanNames = {pulsemetaData.ElectrodesInfo.Label}'; %- we expect "ain1, ain2, ain3, ain4... but if not under cervello controll could be ainp1 ainp2
    if length(allChanNames) < 4
        fprintf('\n uh oh... expected at least 4 ain channels, found %d.  look into this before moving forward (ask Mo or JW)', length(allChanNames));
        keyboard;
    end
    numDCchan = length(allChanNames);
else
    numDCchan = 4; %- assume 4 if using NEV
end


%- loop over all DC channels to output.  In the futuer this should be yolked to jacksheetBR and element_info... 
for DC=[1:numDCchan]
    
    if useNEV==1 && DC>1; continue;
    end
    
    if (contains(EEG_file,'ieeg2') || contains(EEG_file,'INST1'))
        DCtoPull = DC+16;
    else
        DCtoPull = DC;
    end
    
    if useNEV==0
        
        analogChanName = sprintf('ain%d',DCtoPull);
        if sum(contains(allChanNames,analogChanName))==0
            analogChanName = sprintf('ainp%d',DCtoPull);
            if sum(contains(allChanNames,analogChanName))==0
                fprintf('\n uh oh... expected ain%d or ainp%d but didnt see that in the recorded DC list. \nlook into this before moving forward (ask Mo or JW)', DC, DC);
                keyboard;
            end
        end
        
        analogDC     = pulseData.Data(find(not(cellfun('isempty', strfind(allChanNames, analogChanName)))),:);
        resamplerate = pulsemetaData.MetaTags.SamplingFreq/samprate_new;
        analogDC     = analogDC(floor(linspace(1,length(analogDC),length(analogDC)/resamplerate)));
        
        % multiply such that the gain to-be applied to DC is same as eeg channels (0.25)
        maxDigital   = pulsemetaData.ElectrodesInfo(1).MaxDigiValue;
        maxAnalog    = pulsemetaData.ElectrodesInfo(1).MaxAnalogValue;
        gainDCchannels = 1/(double(maxDigital)/double(maxAnalog));
        ns2MaxAnalog = metaData.ElectrodesInfo(1).MaxAnalogValue;
        ns2MaxDigi   = metaData.ElectrodesInfo(1).MaxDigiValue;
        gainEEGchannels = 1/((double(ns2MaxDigi) / double(ns2MaxAnalog)));
        analogDC     = analogDC .* (gainDCchannels/gainEEGchannels);
        % get timestamp of analog DC09 pulses
        [syncPulses] = get_triggers(double(analogDC'),samprate_new);
        
        
    elseif useNEV==1
        
        syncPulses = {pulseData};
        
    end
    
    
    % write out analog files
    
    if useNEV==0
        
        if (contains(EEG_file,'ieeg1') || contains(EEG_file,'INST0'))
            
            split_filename = sprintf('%s/DC%02d_nsp1',folderOutDir,DC+8);
            [fchan,msg] = fopen(split_filename,'w','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
            fwrite(fchan,analogDC,'int16');
            fclose(fchan);
            
            if pulses_to_use==1
                analogDC = analogDC(syncStartIndex:end);
                %analogDC = analogDC .* (gainDCchannels/gainEEGchannels); %- JW thinks already done above
                [syncPulses_aligned]=get_triggers(double(analogDC'),samprate_new);
                
                if DC==1; split_filename = sprintf('%s/DC0%s',folderOutDir,num2str(DC+8)); %easier syntax for this, I forgot
                else; split_filename = sprintf('%s/DC%s',folderOutDir,num2str(DC+8)); end
                [fchan,msg] = fopen(split_filename,'w','l');
                assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
                fwrite(fchan,analogDC,'int16');
                fclose(fchan);
            end
            
            
        elseif (contains(EEG_file,'ieeg2') || contains(EEG_file,'INST1'))
            
            split_filename = sprintf('%s/DC%02d_nsp2',folderOutDir,DC+8);
            [fchan,msg] = fopen(split_filename,'w','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
            fwrite(fchan,analogDC,'int16');
            fclose(fchan);
            
            if pulses_to_use==2
                analogDC = analogDC(syncStartIndex:end);
                %analogDC = analogDC .* (gainDCchannels/gainEEGchannels); %- JW thinks already done above
                [syncPulses_aligned]=get_triggers(double(analogDC'),samprate_new);
                
                if DC==1; split_filename = sprintf('%s/DC0%s',folderOutDir,num2str(DC+8)); %easier syntax for this, I forgot
                else; split_filename = sprintf('%s/DC%s',folderOutDir,num2str(DC+8)); end
                [fchan,msg] = fopen(split_filename,'w','l');
                assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
                fwrite(fchan,analogDC,'int16');
                fclose(fchan);
            end
            
            
        else
            %- single NSP, no alignment
            
            if DC==1; split_filename = sprintf('%s/DC0%s',folderOutDir,num2str(DC+8)); %easier syntax for this, I forgot
            else; split_filename = sprintf('%s/DC%s',folderOutDir,num2str(DC+8)); end
            [fchan,msg] = fopen(split_filename,'w','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
            fwrite(fchan,analogDC,'int16');
            fclose(fchan);
            
        end
        
    end%??????? end for if useNEV == 0??? - added by SJ
    
    % write out trig files
    % not writing out trigDC12.updown
    if DC==1 || DC==2 || DC==4

        if (contains(EEG_file,'ieeg1') || contains(EEG_file,'INST0'))
            if     DC==1
                DCsuffix='trigDC09_nsp1.sync.txt';
            elseif DC==2
                DCsuffix='trigDC10_nsp2.syncStim.txt';
            elseif DC==4
                DCsuffix='trigDC12_nsp2.syncNSP.txt'; 
            end
            fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
            fprintf(fchan,'%d\n',syncPulses{1}');
            fclose(fchan);

            if pulses_to_use==1
                if     DC==1
                    DCsuffix='trigDC09.sync.txt';
                elseif DC==2
                    DCsuffix='trigDC10.syncStim.txt';
                elseif DC==4
                    DCsuffix='trigDC12.syncNSP.txt'; 
                end
                fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
                fprintf(fchan,'%d\n',syncPulses_aligned{1}');
                fclose(fchan);
            end


        elseif (contains(EEG_file,'ieeg2') || contains(EEG_file,'INST1'))

            if     DC==1
                DCsuffix='trigDC09_nsp2.sync.txt';
            elseif DC==2
                DCsuffix='trigDC10_nsp2.syncStim.txt';
            elseif DC==4
                DCsuffix='trigDC12_nsp2.syncNSP.txt'; 
            end
            fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
            fprintf(fchan,'%d\n',syncPulses{1}');
            fclose(fchan);

            if pulses_to_use==2
                if     DC==1
                    DCsuffix='trigDC09.sync.txt';
                elseif DC==2
                    DCsuffix='trigDC10.syncStim.txt';
                elseif DC==4
                    DCsuffix='trigDC12.syncNSP.txt';
                end
                fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
                fprintf(fchan,'%d\n',syncPulses_aligned{1}');
                fclose(fchan);
            end

        else

            if     DC==1
                DCsuffix='trigDC09.sync.txt';
            elseif DC==2
                DCsuffix='trigDC10.syncStim.txt';
            elseif DC==4
                DCsuffix='trigDC12.syncNSP.txt'; 
            end

            fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
            if useNEV==1
                fprintf(fchan,'%0.0f\n',syncPulses{1}');
            else
                fprintf(fchan,'%d\n',syncPulses{1}');
            end
            fclose(fchan);

        end

    end %if DC=1 or 2
end  %for DC=[1:4]


% split and write all channels into noreref
fprintf('\nwriting files:')
ticker=0;
tick_inc=10;
nontrigchans=find(cellfun('isempty', regexp(channellabels, 'DC'))); % only for non-DC chans
for thisChan=nontrigchans
    if thisChan/length(nontrigchans)*100>=ticker
        fprintf(' %2.0f%%',ticker)
        ticker=ticker+tick_inc;
    end
    split_filename = sprintf('%s/%s',folderOutDir,channellabels{thisChan});
    [fchan,msg] = fopen(split_filename,'w','l');
    assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
    
    fwrite(fchan,NSxdata(thisChan,:),'int16');
    
    fclose(fchan);
    if ismac
        try
            % JHW - change files to "executable"... helps for sorting in mac finder
            % MST - avoid permission error (only owner can change permission)
            fileattrib(chanfile, '+x', 'a');
        catch
        end
    end
end


allTags = {};
end










