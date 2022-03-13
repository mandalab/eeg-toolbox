function updateUtahInfoCSV(subj, overrideDataPath56, Copy2UpdateandPUB, SKIP_KEYBOARD)
% getUtahInfoCSV(subj)
% Creates and updates the bookKeeping csv files with information about one subject's utah array data.
%
%  SUBJ, "NIH060"
%
%  OPTIONAL PARAMETERS:
%
% overrideDataPath56:  if data resides locally, instead of on FRNU56, put the path to the subjects here
%                     e.g. '/Volumes/JW24TB/data24TB/localFRNU56/';
%     alternatively, if you want to work on FRNU72, pass "72" for the override path
%
% Copy2UpdateandPUB:    (default true) If you want to copy the created CSV
% to PUB and eeg update. (You pretty much always want to do this)
%
%  SKIP_KEYBOARD:  (default 0)  should never skip the keyboards... unless you've already checked that they are OK
%          only skips warnings about channel names not matching... other keyboards/errors still in play
%
% Written by: Anthony I. Jang
% updated by JW 1/2019 -- removed "task" file
% updated by ME 1/2019 -- interacts with jacksheet/elec labels from raw data
% updated by SJ 10-11/2019 -- A few updates to this function:
%   Previous micro_PickSess2Extract.xlsx files were listed in a different
%   order and did not contain information on the task, notes, or sort
%   status. It is now updated to consider the old version and new version,
%   and update accordingly. Features added:
%       1. Recognize both old and new micro_PickSess2Extract.xlsx formats.
%          Update format to new if it is old.
%       2. Include information on whether sorts were attempted and
%          completed.
%       3. Include information from rawFileList on 'task' and 'notes'. This
%          uses a helper function called getRawFileData.m, which will iterate
%          through each rawFileList present in a subject's directory and
%          select the one with the most information.
%       4. Copy updated micro_PickSess2Extract.xlsx to 56PUB, as well as
%          56PROC eeg .update (Copy2UpdateandPUB flag)
% updated by SJ 12/2019 -- make sure output to micro_RenamedChanInfo.csv is
%           compatible with commas (change commas to semicolons)
% updated by SJ 1/9/2020 -- 
%       - take task from alignment
%       - if no alignment, write 'not_aligned'
%       - either case, put rawfilelist task & note in notes as: '(task)note'
% updated by SJ 2/13/2020 --
%       - check to make sure the paths to FRNU, eeg/.update, and FRNU are the typical paths
%
% 
% TO DO:
%   - Have microRenamedChanInfo use the current version of micro_PickSess2Extract instead of the one that
%     was previously written (for task)
%   - Include writeSortNotes_readme
% 


if nargin<4
    SKIP_KEYBOARD = 0;
end
if nargin<3
    Copy2UpdateandPUB = true;
elseif Copy2UpdateandPUB == 1
    Copy2UpdateandPUB = true;
end
if nargin<2
    overrideDataPath56 = '';
end

headerorder = {'folderName', 'toExtract', 'isStim', 'num_NS5or6', 'num_microChan', ...
    'pulses_beh', 'pulses_stim', 'br_duration', 'samplesAdded', 'maxMicroRng', 'extracted', ...
    'sort_attempted', 'sort_complete', 'task', 'notes'};


FRNUpath = '/Volumes/Shares/FRNU/data/eeg';
eegUpdatePath = '/Volumes/56PROC/eeg_processing/.update/eeg';
PUBpath = '/Volumes/56PUB/readWrite/micro_forSorting';

% FRNUpath = '/Volumes/Shares/FRNU/data/eeg';
% eegUpdatePath = '/Volumes/SeagateBackupPlusDrive/local56/56PROC/eeg_processing/.update/eeg';
% PUBpath = '/Volumes/SeagateBackupPlusDrive/local56/56PUB/readWrite/micro_forSorting';

if ~strcmp(eegUpdatePath,'/Volumes/56PROC/eeg_processing/.update/eeg') || ~strcmp(PUBpath,'/Volumes/56PUB/readWrite/micro_forSorting') || ~strcmp(FRNUpath,'/Volumes/Shares/FRNU/data/eeg')
    fprintf('%s\n','WAIT!!!! Your eegUpdate and PUBpath are not designated as the typical paths. Are you sure you want to continue??');
    keyboard
end


subjectInfo = getMicroSubjInfo_v11(subj,overrideDataPath56);
if ~isempty(overrideDataPath56)
    process_micro_dir = overrideDataPath56; 
else
    process_micro_dir = process_micro_dir; 
end


fprintf('\n\n ******** Updating %s Micro Session and Channel Selection CSVs in %s ******** \n', subj,process_micro_dir);

% Directories for the three csv files
fileDir_bookkeeping   = fullfile(process_micro_dir ,subj,'_extraction_notes','micro_PickSess2Extract.xlsx');
fileDir_spikeChans    = fullfile(process_micro_dir ,subj,'_extraction_notes','micro_RenamedChanInfo.csv');


%- define two fixed parameters for getBR_rawFileList
testRun = 0;  onlyReturnUtahInfo = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Create any utahInfo csv files that were not created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(fileDir_bookkeeping,'file')
    fprintf('No microPickSess2Extract.xlsx file found. Creating...\n');
    
    %- run through the raw_data folders, load jacksheetBR, and output those folders with ns5 or ns6 for potential extraction
    [sortFileStruct, sortFileCell] = getBR_rawFileList(subj,'',process_micro_dir ,testRun, onlyReturnUtahInfo);
    
    sortFileStruct = rmfield(sortFileStruct,'sorted');
    blanks_newstruct = repmat({'-'},numel(sortFileStruct),1);
    
    [sortFileStruct.sort_attempted] = blanks_newstruct{:};
    [sortFileStruct.sort_complete] = blanks_newstruct{:};
    
    sortFileStruct = orderfields(sortFileStruct,headerorder); %Reorder
    sortFileCell = cat(1,headerorder,squeeze(struct2cell(sortFileStruct))'); %Make sortFileCell match the struct 
    sortFileTable = struct2table(sortFileStruct);
    writetable(sortFileTable,fileDir_bookkeeping);
    
else 
    temptable = readtable(fileDir_bookkeeping);
    if ~isequal(temptable.Properties.VariableNames,headerorder) % If it already exists but doesn't match the current layout
        tempstruct = table2struct(temptable);
        tempstruct = rmfield(tempstruct,'sorted');
        blanks_newstruct = repmat({'-'},numel(tempstruct),1);
        [tempstruct.sort_attempted] = blanks_newstruct{:};
        [tempstruct.sort_complete] = blanks_newstruct{:};
        tempstruct = orderfields(tempstruct,headerorder); %Reorder
        temptable = struct2table(tempstruct);
    end
    sortFileTable = temptable;
    sortFileCell  = table2cell(sortFileTable); %- code below is looking for a cell matrix... 
    sortFileCell  = cat(1,sortFileTable.Properties.VariableNames,sortFileCell);
end



% is there a session marked for extraction? only if so, lets create csv for chan selection
createChanSelectionCSV = 0;
extractColIndex = find(strcmp(sortFileCell(1,:),'toExtract'));
sessToExtract   = find(strcmpi(sortFileCell(:,extractColIndex),'y'));
if size(sessToExtract,1) > 0
    createChanSelectionCSV = 1;
end


% create spike chans...  POINT of this file is 
%     (1) to centralize the definition newChanName to be used in LFP and spikeInfo creation. Names are grabbed from HERE
%     (2) So IF the newChanNames change in this file, that should trigger a rename of existing LFP and spikeInfo files [not implemented yet]
%     (3) Lastly, this file is meant to confirm that all sessions have the same set of NSP and oldChan name to be sure there is no ambiguity (i.e., change in channel names)  
if 1 || ~exist(fileDir_spikeChans,'file') && createChanSelectionCSV==1
    fprintf('\n Updating or Creating Channel Rename csv: %s',fileDir_spikeChans);
    
    if ~exist(fileDir_spikeChans,'file')
        fprintf('\n    CSV not found. Creating...',fileDir_spikeChans);
        microNameTableOld = table; %- initialize to empty 
    else
        fprintf('\n    CSV Found.  Confirming ChanName -> ChanNameNew is unchanged and sessions are up to date...',fileDir_spikeChans);
        spikeChansOld     = csv2cell(fileDir_spikeChans);
        microNameTableOld = cell2table(spikeChansOld(4:end,1:4),'VariableNames',spikeChansOld(3,1:4));
    end
    
    
    % load jacksheet and get names of all 30kHz channels present in all files marked for extraction
    tStartJackCheck = tic;
    microNameTableAll = table;
    allJackTables = {};
    for iSess=1:length(sessToExtract)
        jackSheetPath = fullfile(process_micro_dir ,subjectInfo.subj,'data_raw',sortFileCell{sessToExtract(iSess),1},'jacksheetBR_complete.csv');
        if ~exist(jackSheetPath,'file'),
            fprintf('\n jacksheet missing from %s... making it now',sortFileCell{sessToExtract(iSess),1});
            jackTable = makeJacksheetBR(fullfile(process_micro_dir ,subjectInfo.subj,'data_raw',sortFileCell{sessToExtract(iSess),1}),'',0,1);
        else
            jackTable = readtable(jackSheetPath);
        end
        
        allJackTables{iSess} = jackTable; %- will use this in next loop...
        
        %- JW thinks we can loop over microNameConvert for each session, concatinating outputs
        %   then confirm there are is only one new name per old+NSPsuffix
        %   then take the unique of the concatinated microNameTables
        [microNameTable, jackTable] = microNameConvert(subjectInfo,jackTable,SKIP_KEYBOARD);
        microNameTableAll = cat(1,microNameTableAll,microNameTable);
        
        microNameTableAll = unique(cat(1,microNameTableAll,microNameTable),'rows','stable');
        
    end %- loop over sess-to-extract to find complete set of 30 kS channel names
    
    
    %- SANITY CHECKS: now confirm that there are no more than one entry per channel (only possibility for this is if chan moved NSPs)
    testTable = microNameTableAll(:,{'NSPsuffix','ChanName'}); %- just nspSuffix and oldChanNAME
    if height(unique(microNameTableAll(:,{'NSPsuffix','ChanName'}),'rows'))<height(microNameTableAll)
        fprintf('\n WARNING: there is more than one copy of NSPsuffix and oldName for a particular channel.  Should never happen. Investigate if you get here');
        keybaord;
    end
    if height(unique(microNameTableAll(:,{'ChanNameNew'}),'rows'))<height(microNameTableAll) %- just newChanName
        if length(subjectInfo.nspSuffixAlt)==0% & length(subjectInfo.newChanName{:,1})==length(unique(subjectInfo.newChanName{:,1})),
            fprintf('\n HEADS UP: there is more than one newChanName for a particular channel, and no nspSuffixAlt specified. ');
            %- breaks here for NIH049 because JW used the same DevNum for two different chanNames (in begninning of recording called "chan", later LADM)
            %- real test would come from if in a single session the same name was used twice...
            %keyboard;
        end
    end
    
    
    %- Now see if an OLD table matches the new table.  This triggers a flag for updating spikeInfos and LFP.mat (not coded yet, so designed to throw an error below.
    if height(microNameTableOld)>0
        compTable = intersect(microNameTableOld(:,{'NSPsuffix','ChanName','ChanNameNew'}),microNameTableAll(:,{'NSPsuffix','ChanName','ChanNameNew'}));
        if height(compTable)~=height(microNameTableAll)
            fprintf('\n SERIOUS WARNING:  NEW CHANNEL NAMES DETECTED... SHOULD DELETE and RECREATE ALL existing SPIKE_INFO and LFP.mat NOW!\n Or write code that just updates names in existing SpikeInfo and LFP.mat');
            keyboard;
            %- first step should be to resave "old" csv so we dont loose the old transform
            oldFileName = sprintf('%s[newNames_%s].csv', fileDir_spikeChans(1:end-4),datestr(now,'YYmmDD_hhss'));
            cell2csv(spikeChansOld,oldFileName);
            
            error('\n Updated Channel Names: DELETE SPIKEINFO+LFP.mat+micro_RenamedChanInfo or Write code to fix');
        end
    end
    
    %- convert table contents to a cell array so we can have multiple headers:   filename/session folder/task
    spikeChanCell1 = cat(1,{'';'';'NSPsuffix'},         microNameTableAll.NSPsuffix);
    spikeChanCell2 = cat(1,{'';'';'ChanName'},          microNameTableAll.ChanName);
    spikeChanCell3 = cat(1,{'';'';'SortChanName'},      microNameTableAll.SortChanName);
    spikeChanCell4 = cat(1,{'';'';'ChanNameNew'},       microNameTableAll.ChanNameNew);
    spikeChanCell5 = {'';'';'MicroDevNum'};             for iR=1:height(microNameTableAll), spikeChanCell5{end+1} = num2str(microNameTableAll.MicroDevNum(iR)); end;
    spikeChanCell6 = {'fileName';'folderName';'task';}; for iR=1:height(microNameTableAll), spikeChanCell6{end+1} = '-'; end
    
    %- OUTPUT TABLE SHOULD CONTAIN:
    % nspSuff / oldChanName / newChanName / microDevNum...  then for each ns5 file: fileName, folderName, task
    spikeChans = {spikeChanCell1{:}; spikeChanCell2{:}; spikeChanCell3{:}; spikeChanCell4{:}; spikeChanCell5{:}; spikeChanCell6{:}}';
    
    
    %- Now push the table info to a cell array?  or just extend the table with a new column for each session?
    for iSess=1:length(sessToExtract)
        
        jackTable     = allJackTables{iSess};
        %row30kChans   = find(jackTable.SampFreq==30000 & ~contains(jackTable.FileName,'.nev'));
        row30kChans   = find(jackTable.MicroDevNum>0 & jackTable.SampFreq==30000);
        files30kChan  = unique(jackTable.FileName(row30kChans),'stable');  %- use 'stable' so returned in the order of NSPs
        
        task          = sortFileCell{strcmp(sortFileCell(:,1),jackTable.RawDir{1}),find(strcmp(sortFileCell(1,:),'task'))}; %- SJ - changed from #2
        if contains(task,',')
            task = strrep(task,',',';');
        end
        if sum(strcmp(sortFileCell(:,1),jackTable.RawDir{1}))==0
            fprintf('\n WARNING: flaw in detection of task name.  Figure it out and fix the error');
            keyboard;
        end
        
        
        %- loop over 30k files for this session
        for iF=1:length(files30kChan)
            
            thisFile30k = files30kChan{iF};
            
            %- add a column to the CSV and mark which channels were found from this file.  For dual NSP recordings (2 utah array, utah+micro), different files may contains different channels
            colToAdd  = size(spikeChans,2) + 1;
            
            spikeChans{1,colToAdd} = thisFile30k; %- ns5 or ns6 file
            spikeChans{2,colToAdd} = jackTable.RawDir{1}; %- session directory or foldername
            spikeChans{3,colToAdd} = task; %- task string
            
            for i=4:size(spikeChans,1)  %- row 4 is first electrode row
                matchWithRef = find(strcmp(jackTable.NSPsuffix,spikeChans{i,1}) & strcmp(jackTable.ChanName,spikeChans{i,2}) & strcmp(jackTable.FileName,thisFile30k)); 
                if length(matchWithRef)>1
                    fprintf('\n ERROR: shouldnt be possible.  how did we get here');
                    keyboard;
                elseif length(matchWithRef) == 0
                    % nope, elec not found
                    spikeChans{i,colToAdd} = 'N/A';
                else
                    %spikeChans{i,colToAdd} = '-';
                    spikeChans{i,colToAdd} = 'Found';
                end
            end
        end %- length(files30kChan)
    end
    
    %- save the result
    cell2csv(spikeChans,fileDir_spikeChans);
    fprintf(' [%.1fs]',toc(tStartJackCheck));
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Once all utahInfo csv files are created, update with latest info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(fileDir_bookkeeping,'file')
    fprintf('\n Updating bookkeeping CSV: %s',fileDir_bookkeeping);  tUpdateBook = tic;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% utahInfo_bookkeeping %%%%%%%%%%
    %
    %  1) scan data_raw to see if bookkeeping is missing any new files.  Offer to update (wihtout loosing mark up) if so.
    %  2) sort the bookkeeping rows with rows to extract separated at the top
    %  3) look in data_processed to see if each session has been extracted, and if so, sorted
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bookkeeping = sortFileCell; %- already loaded above and not manipulated.. so no need to reload here
    clear bkStruct;
    % Change back to structure so it's easier to work with
    for ff = 1:size(bookkeeping,2) %SJLOOK- will need to add cols
        for row = 2:size(bookkeeping,1)
            bookkeepValue = bookkeeping{row,ff};
            if sum(strcmp(bookkeeping{1,ff},{'toExtract' 'isStim' 'extracted' 'sort_attempted' 'sort_complete' 'task' 'notes'}))>0 && isempty(bookkeepValue)
                bookkeepValue = '-'; %- enforce that "blanks" are indicated by '-'... otherwise a potential crash below
            end
            sortFileStruct(row-1).(bookkeeping{1,ff}) = bookkeepValue;
        end
    end
    
    % First, see if new data files have been added since the bookkeeping csv was last updated.
    allRelevantDataFiles = {};
    allRelevantDataDir   = {};
    for rawDirFolder = getDirFileNames(fullfile(process_micro_dir ,subj,'data_raw'))
        sessFileNames = getDirFileNames(fullfile(process_micro_dir ,subj,'data_raw',char(rawDirFolder)));
        % Only look at utah array files
        toUse_i = false(length(sessFileNames),length(subjectInfo.utahChanCnt)); %- this was length(seubjectInfo.utahSuffix, but that field isn't used anymore
        for suf = 1:length(subjectInfo.utahChanCnt)  %- this was length(seubjectInfo.utahSuffix, but that field isn't used anymore
            toUse_i(:,suf) = (~cellfun('isempty',strfind(sessFileNames,'.ns5')) | ~cellfun('isempty',strfind(sessFileNames,'.ns6')))';
        end
        if isempty(sessFileNames(logical(sum(toUse_i,2)))'); continue; end
        allRelevantDataFiles = cat(1,allRelevantDataFiles,sessFileNames(logical(sum(toUse_i,2)))');
        allRelevantDataDir{end+1} = char(rawDirFolder);
    end
    
    % If any new data has been added, prompt the user whether they want to update the bookkeeping csv.
    if any(~ismember(allRelevantDataDir,{sortFileStruct.folderName}))
        fprintf('\nThe following sessions are not in %s:\n',fileDir_bookkeeping);
        for k = find(~ismember(allRelevantDataDir,{sortFileStruct.folderName})), fprintf('  %s\n',allRelevantDataDir{k}); end
        reRun_bookkeepingUpdate = input('Update bookkeeping csv with new sessions (y,n)? ','s');
        if strcmpi(reRun_bookkeepingUpdate,'y')
            
            [rawFile_struct, sortFileCell] = getBR_rawFileList(subj,'',process_micro_dir ,testRun,onlyReturnUtahInfo);
            
            i_append = length(sortFileStruct)+1;
            for ff = 1:length(rawFile_struct)
                if ~ismember(rawFile_struct(ff).folderName,{sortFileStruct.folderName})
                    % First, fill with '-'
                    for fn = fieldnames(sortFileStruct)', sortFileStruct(i_append).(char(fn)) = '-'; end
                    sortFileStruct(i_append).folderName    = rawFile_struct(ff).folderName;
                   
                    sortFileStruct(i_append).toExtract     = rawFile_struct(ff).toExtract; %- getBR_rawFileList take a guess as to whether a session should be extracted or isStim
                    sortFileStruct(i_append).isStim        = rawFile_struct(ff).isStim;
        
                    sortFileStruct(i_append).num_NS5or6    = rawFile_struct(ff).num_NS5or6;
                    sortFileStruct(i_append).num_microChan = rawFile_struct(ff).num_microChan;
   
                    sortFileStruct(i_append).pulses_beh    = rawFile_struct(ff).pulses_beh;
                    sortFileStruct(i_append).pulses_stim   = rawFile_struct(ff).pulses_stim;
                    sortFileStruct(i_append).br_duration   = rawFile_struct(ff).br_duration;
                    sortFileStruct(i_append).samplesAdded  = rawFile_struct(ff).samplesAdded;
                    sortFileStruct(i_append).maxMicroRng   = rawFile_struct(ff).maxMicroRng;
                    
                    i_append = i_append+1;
                end
            end
        end
    end
    
    % Sort bookkeeping rows
    LIST_EXTRACTED_FIRST = 0 ; %- this nice option by Anthony made sense when the list had EVERY session recorded... now its a pretty short list
    if LIST_EXTRACTED_FIRST
        % Show files to be extracted on top, sorted
        toExtract_i    = find(ismember({sortFileStruct.toExtract},'y'));
        [~,toExtract_sort_i] = sort({sortFileStruct(toExtract_i).folderName});
        toExtract_i    = toExtract_i(toExtract_sort_i);
        % Show all other files, sorted
        toNotExtract_i = find(~ismember({sortFileStruct.toExtract},'y'));
        [~,toNotExtract_sort_i] = sort({sortFileStruct(toNotExtract_i).folderName});
        toNotExtract_i = toNotExtract_i(toNotExtract_sort_i);
        sortFileStruct = sortFileStruct([toExtract_i,toNotExtract_i]);
    else
        % Sort everything by folername (date string)
        [~,sort_i] = sort({sortFileStruct(:).folderName});
        sortFileStruct = sortFileStruct([sort_i]);
    end
    
    % Go through each file and check whether it's been extracted and sorted
    % (SJ)
    processedDir = fullfile(process_micro_dir ,subj,'processed_SPK');
    sortedDir = fullfile(process_micro_dir ,subj,'sorts_manual');
    processedDirFiles = getDirNamesRegexp(processedDir,'^\d.*');
    sortedDirFiles = getDirNamesRegexp(sortedDir,'^\d.*');
    % Call getRawFileData to get tasks and notes
    table2use_raw = getRawFileData(subj, FRNUpath);
    table2use_alignment = getTaskFromAlignment(subj, eegUpdatePath);
    for ff = 1:length(sortFileStruct)
        % Is it aligned to DC09?
        fold1 = sortFileStruct(ff).folderName;
        if any(~cellfun('isempty',strfind(processedDirFiles,sprintf('%s',sortFileStruct(ff).folderName))))
            sortFileStruct(ff).extracted = 'y';
            % Are the spikes sorted and extracted (look for spikeInfo.mat)?
        end
        if any(contains(sortedDirFiles,fold1)) % Fill in sort_attempted
            sortFileStruct(ff).sort_attempted = 'y';
            sort_summ_file = char(getDirNamesRegexp([sortedDir filesep fold1 filesep 'reref*sortedBy*'],'^sorts_.*\.txt'));
            if contains(sort_summ_file,'(complete)') % Fill in sort_complete
                sortFileStruct(ff).sort_complete = 'y';
            end
        end
        % Check raw file list (SJ)
        %rawfiles = getDirNamesRegexp([FRNUpath filesep subj filesep 'raw'],'.*rawFile.*');
        %rawFiles = 1;
        if ~isempty(table2use_alignment)
            align_task = char(table2use_alignment.task(strcmp(table2use_alignment.folderName,fold1)));
            if isempty(align_task)
                align_task = 'not_aligned';
            end
        else
            align_task = 'not_aligned';
        end
        sortFileStruct(ff).task = align_task;
        
        if ~isempty(table2use_raw)
            raw_task = char(table2use_raw.task(contains(table2use_raw.folderName,fold1(1:11))));
            raw_notes = char(table2use_raw.notes(contains(table2use_raw.folderName,fold1(1:11))));
            % Now make sure you get rid of any cases where there are hidden
            % new lines
            if size(raw_task,1) > 1
                raw_task = raw_task(1,:);
            end            
            if size(raw_notes,1) > 1
                raw_notes = raw_notes(1,:);
            end
            if ~isempty(raw_task) && ~strcmp(raw_task,'-')
                %sortFileStruct(ff).task = raw_task;
                raw_task_str = ['(' raw_task ')'];
            else
                raw_task_str = '(-)';
            end
            if ~isempty(raw_notes) && ~strcmp(raw_notes,'-')
                %sortFileStruct(ff).notes = raw_notes;
                raw_notes_str = raw_notes;
            else
                raw_notes_str = '-';
            end
            notes_final_str = [raw_task_str raw_notes_str];
        else
            notes_final_str = '--';
        end
        sortFileStruct(ff).notes = notes_final_str;
        
    end
    
    %- convert to a table for output of xlsx using "writetable", and so following code can quickly compare notes and task fields
    sortFileStruct(end+1) = sortFileStruct(end); %- little trick to make table use cell arrays even if only one element in the structure... only really matters if length(sortFileStruct)==1
    sortFileTable = struct2table(sortFileStruct);
    sortFileTable = sortFileTable(1:end-1,:);
    
    % Save bookkeeping file
    writetable(sortFileTable,fileDir_bookkeeping); %- xlsx format
    
    
    % Copy to eeg .update and PUB after creation
    if Copy2UpdateandPUB
        % Copy bookkeeping files to eeg update to be pushed to FRNU
        st1 = copyfile(fileDir_bookkeeping,[eegUpdatePath filesep subj filesep 'micro']); 
        if st1 == 1 
            fprintf('\n%s\n','xlsx copied to EEG .update');
        else
            fprintf('\n%s\n','ERROR!!! xlsx not able to be copied to EEG .update!!');
            keyboard
        end

        % Copy bookkeeping files to 56PUB (replace previous)
        PUBpathcomplete = [PUBpath filesep subj filesep '_extraction_notes'];
        if exist(PUBpathcomplete,'dir')
            st2 = copyfile(fileDir_bookkeeping,PUBpathcomplete); 
            st22 = copyfile(fileDir_spikeChans,PUBpathcomplete); 
            if st2 == 1 && st22 == 1
                fprintf('%s\n','CSV and xlsx copied to 56PUB');
            else
                fprintf('\n%s\n','ERROR!!! CSV and/or xlsx not able to be copied to 56PUB!!');
                keyboard
            end
        else
            fprintf('\n%s\n','Subject not established in 56PUB yet. CSV will be copied in COPY_2_PUBforSorting step. Is that right? Are you connected to 56PUB?');
            keyboard
        end
    end
    
    fprintf(' [%.1fs]', toc(tUpdateBook)); 

    fprintf('\n all done \n');

end
end