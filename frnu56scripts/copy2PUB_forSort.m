function done = copy2PUB_forSort(subj, rootEEGdir, rootPUBdir)

% copy any new subjects/sessions from PROC to PUB (which includes some renaming/trimming of stuff)
%
% from 56PROC/micro_behavioral/micro_pristine/NIH0##/processed_SPK/sess ->
% 56PUB/readWrite/micro_forSorting/NIH0##/sess
% Inputs:
%   subj - string specifying the subject, ex) 'NIH060'
%   rootEEGdir - string specifying the root eeg directory that has your
%                   subject of interest in it 
%                   ex) '/Volumes/56PROC/micro_behavioral/micro_pristine'
%   rootPUBdir - string specifying the directory where you want to copy the
%                   sort information (sort figs, spikeinfo, sorts, sort log)
%                   ex) '/Volumes/56PUB/readWrite/micro_forSorting'
%
%
% 10/1/2019 Created by Samantha Jackson 
% 1/6/2020 SJ: Added infunctionality for if things weren't copied corectly the first time (also to fix issues from extraction info and stim folders not being created)
% 3/10/2020 SJ: Quick fix for stim sessions

% subj = 'NIH066';
% rootEEGdir = '/Volumes/56PROC/micro_behavioral/micro_pristine';
% rootPUBdir = '/Volumes/56PUB/readWrite/micro_forSorting';

%%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% %%%%%%%%%% 
stuff2copyfromRAW = {'jacksheetBR_complete.csv', ...
    '*.ns5', ...
    '*.ns6'};
stuff2copyfromSPK = {'sort_reref*', ...
    'sortNotes_readme.txt', ...
    'sortNotes(*)_sortedBy??.xlsx'};

stimMaskstr = 'stimMask5'; % In case this ever changes!!

FIX_STIM_AND_EXTRACTION_ISSUE = false;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s\n',['Preparing to copy from ' rootEEGdir ' to ' rootPUBdir]);
fprintf('\n%s\n',['Subject ' subj ': ']);

% if FIX_STIM_AND_EXTRACTION_ISSUE
%     fprintf('%s\n','WARNING! You have FIX_STIM_AND_EXTRACTION_ISSUE set to TRUE. Please stop and change this if you are just trying to run under default conditions.');
% end


SPKsessdir = dir([rootEEGdir filesep subj filesep 'processed_SPK' filesep '*_*']);
SPKsesslist = {SPKsessdir.name};
SPKsesslist = SPKsesslist(~ismember(SPKsesslist,'.DS_Store'));

% What is already in PUB? We only want to copy sessions that aren't already there!!
PUBsessdir = dir([rootPUBdir filesep subj filesep '*_*']);
PUBsesslist = {PUBsessdir.name};
PUBsesslist = PUBsesslist(~ismember(PUBsesslist,'.DS_Store'));
PUBsesslist = PUBsesslist(find(~cellfun(@isempty, regexp(PUBsesslist,'\d{6}_\d{4}'))));
PUBsesslist_clean = regexp(PUBsesslist,'\d{6}_\d{4}(_beh)?(_premie)?(?=.*)','match'); % remove 'grabbed' or 'ready'
PUBsesslist_clean = horzcat(PUBsesslist_clean{:});

newsess_list = setdiff(SPKsesslist,PUBsesslist_clean);

if isempty(newsess_list)
    if FIX_STIM_AND_EXTRACTION_ISSUE
        nsans = input('No new sessions for this subject, but FIXER activated! Would you like to continue anyway (C, default), or skip this subject (S)?','s');
        if strcmp(nsans,'S') || strcmp(nsans,'s')
            done = 1;
            fprintf('%s\n','Skipping ... ');
            return;
        else
            fprintf('%s\n','Continuing to copying...');
            readygrabbed_list = PUBsesslist(contains(PUBsesslist,'grabbed')|contains(PUBsesslist,'ready'));
            % If the session is already sorted, it doesn't need to be fixed
            if ~isempty(readygrabbed_list)
                readygrabbed_list_clean = regexp(readygrabbed_list,'\d{6}_\d{4}(_beh)?(_premie)?(?=.*)','match');% remove 'grabbed' or 'ready'
                readygrabbed_list_clean = horzcat(readygrabbed_list_clean{:});
                newsess_list = setdiff(SPKsesslist,readygrabbed_list_clean);
            else
                newsess_list = SPKsesslist; 
            end
        end
    else
        done = 1;
        fprintf('%s\n','No new sessions for this subject! Skipping ... ');
        return;
    end
elseif ~isempty(newsess_list) && FIX_STIM_AND_EXTRACTION_ISSUE
    fprintf('%s\n','ERROR!!!! newsess_list is not empty but FIX_STIM_AND_EXTRACTION_ISSUE is triggered!!!');
    keyboard
end
    
PROCsubpath = [rootEEGdir filesep subj];
PUBsubpath = [rootPUBdir filesep subj];
%pubt2 = split(pubt,'_',1)

copy_list = cell(1,2);

% Make Subject Directory
if exist(PROCsubpath,'dir') && ~exist(PUBsubpath,'dir')
    mdstatus = mkdir(PUBsubpath);
    if mdstatus == 1
        fprintf('\t%s\n',['CREATED DIRECTORY: ' PUBsubpath]);
    else
        fprintf('\n\n%s\n',['ERROR!!! could not make directory: ' PUBsubpath]);
        keyboard;
    end
elseif exist(PROCsubpath,'dir') && exist(PUBsubpath,'dir')
    fprintf('\t%s',['Subject folder already exists at ' PUBsubpath]);
else
    fprintf('\t%s',['Subject folder does not exist in PROC!! HUH?? (' PROCsubpath ')']);
end

% Copy _extraction_notes if new subject or not any on PUB yet
extraction_path = [PROCsubpath filesep '_extraction_notes'];
if exist(extraction_path,'dir')
    extraction_path_PUB = [PUBsubpath filesep '_extraction_notes'];
    if ~exist(extraction_path_PUB,'dir')
        fprintf('\t%s\n',[extraction_path_PUB ' does not exist ... will create it now.']);
        fprintf('\t%s',[extraction_path ' -> ' extraction_path_PUB]);
        status1 = copyfile(extraction_path, extraction_path_PUB);
        if status1 == 1
            fprintf('%s\n',' ... done.');
        else
            fprintf('\n\n%s\n',['ERROR!!! could not copy ' extraction_path ' -> ' extraction_path_PUB]);
            keyboard;
        end
    end
else
    fprintf('\t%s\n',[extraction_path ' does not exist! Why not?']);
    keyboard;
end
        

for ii = 1:numel(newsess_list)
    newsess = newsess_list{ii};
    
    % Check to see if session exists in PUB
    
    %nPlay_visualize folder with jacksheetcomplete and INST0.ns5/6 from
    %data_raw
    PROCrawpath = [PROCsubpath filesep 'data_raw' filesep newsess];
    PUBnPlaypath = [PUBsubpath filesep newsess filesep 'nPlay_visualize'];
    if ~exist(PUBnPlaypath,'dir')
        mdstatus = mkdir(PUBnPlaypath);
        if mdstatus == 1
            fprintf('\t%s\n',['CREATED DIRECTORY: ' PUBnPlaypath]);
        else
            fprintf('\n\n%s\n',['ERROR!!! could not make directory: ' PUBnPlaypath]);
            keyboard;
        end
    end
    for rr = 1:numel(stuff2copyfromRAW)
        rawstr = stuff2copyfromRAW{rr};
        rawstrdir = dir([PROCrawpath filesep rawstr]);

        if ~isempty(rawstrdir)
            if numel(rawstrdir) > 1 && ~contains(rawstr,'.ns')
                fprintf('\n\n%s\n',['ERROR!! multiple occurrences of ' PROCrawpath filesep rawstr '. Please check!']);
                keyboard
            end
            
            thisname_list = {rawstrdir.name};
            for mm = 1:numel(thisname_list)
                thisname = thisname_list{mm};
                raw_from = [PROCrawpath filesep thisname];
                raw_to = [PUBnPlaypath filesep thisname];
                if isempty(copy_list{1})
                    copy_list(1,:) = {raw_from, raw_to};
                else
                    copy_list(end+1,:) = {raw_from, raw_to};
                end
            end
        else
            fprintf('\t%s\n',['No occurrences of ' rawstr ' in ' newsess ' to copy. Should be okay if .ns6 or .ns5']);
        end
    end
    
    % sort_reref, sortNotes txt, sortNotes xlsx from processed_SPK
    % 
    PROCSPKpath = [PROCsubpath filesep 'processed_SPK' filesep newsess];
    PUBsesspath = [PUBsubpath filesep newsess];
    nplaypath = [PUBsesspath filesep 'nPlay_visualize'];
    for ss = 1:numel(stuff2copyfromSPK)
        SPKstr = stuff2copyfromSPK{ss};
        SPKstrdir = dir([PROCSPKpath filesep SPKstr]);
        if numel(SPKstrdir) > 1 % Only copy stimMask session if it is there
            SPKstrlist = {SPKstrdir.name};
            thisname = SPKstrlist{find(contains(SPKstrlist,stimMaskstr))};
            % Now, delete the non-stim one?
            stim2del = SPKstrlist{find(~contains(SPKstrlist,stimMaskstr))};
            if ~contains(stim2del,'sort_reref')
                fprintf('%s\n',['ERROR! Multiple cases of this, but it is not a stim sort reref. If you continue, the following file will be deleted: ' stim2del]);
                keyboard
            end
            if exist([PUBsesspath filesep stim2del],'dir')
                fprintf('\n%s\n','Non-stim mask folder found in PUB when there should be stim mask. Deleting and replacing now.');
                rmdir([PUBsesspath filesep stim2del],'s');
            end
        else
            thisname = SPKstrdir.name;
        end
        

        SPK_from = [PROCSPKpath filesep thisname];
        SPK_to = [PUBsesspath filesep thisname];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%% DELETE THIS AFTER FIXING EXTRACTION ISSUE %%%%%%%%%%%%%%%%%%%%
        if FIX_STIM_AND_EXTRACTION_ISSUE
            if contains(thisname,'sort_reref') && ~contains(thisname,stimMaskstr) %Just regular 

                del_files = strcat(SPK_to,filesep,getDirNamesRegexp(SPK_to,'(_extractMicro_readme|jacksheetBR_justSplitChan)'));
                if numel(del_files) ~= 2
                    fprintf('%s\n','ERROR!!!!! Should be deleting exactly 2 files for _extractionInfo!!!');
                    keyboard
                end
                for dd = 1:numel(del_files)
                    delete(del_files{dd});
                end

                % Now delete the _extractionInfo directory
                extractionInfo_path = [PUBsesspath filesep thisname filesep '_extractionInfo'];
                status_rmdir = rmdir(extractionInfo_path,'s');
                if status_rmdir ~= 1
                    fprintf('%s\n',['ERROR!!! Could not remove this directory: ' extractionInfo_path]);
                    keyboard
                end
                SPK_from = [SPK_from filesep '_extractionInfo'];
                SPK_to = extractionInfo_path;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(copy_list{1})
            copy_list(1,:) = {SPK_from, SPK_to};
        else
            copy_list(end+1,:) = {SPK_from, SPK_to};
        end
    end    

end

% Copy time
% There should be an error if the copy_to already exists!!! do not let
% overwrite

for jj = 1:size(copy_list,1)
    copy_from = copy_list{jj,1};
    copy_to = copy_list{jj,2};
    
    if exist(copy_from)
        if ~exist(copy_to)
            fprintf('\t%s',['COPY: ' copy_from ' -> ' copy_to]);
            status2 = copyfile(copy_from, copy_to);
            if status2 == 1
                fprintf('%s\n',' ... done.');
            else
                fprintf('\n\n%s\n',['ERROR!!! could not copy ' copy_from ' -> ' copy_to]);
                keyboard;
            end
        else
            fprintf('\n\n%s\n',['ERROR!!!! ' copy_to ' already exists!! Not copying.']);
        end
    else
        fprintf('\n\n%s\n',['ERROR!!! ' copy_from ' does not exist! You cannot copy from something that does not exist!']);
        keyboard
    end
end

fprintf('\n\n%s\n',['SUBJECT ' subj ' COMPLETE!']);
done = 1;

end

