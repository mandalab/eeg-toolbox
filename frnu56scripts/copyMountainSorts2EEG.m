function done = copyMountainSorts2EEG(subj,rootMountainDir,rootupdatedir)

% rootMountainDir = '/Volumes/56PROC/micro_biowulfSorts/mountain_sorts';
% rootupdatedir = '/Volumes/56PROC/eeg_processing/.update/eeg';
% subj = 'NIH029';


% NOTE! This function is not finished. It will only copy mountain sorts if
% there is no mountain sort folder for that subject in eeg/.update!!!

%sublist = getDirNamesRegexp(rootSortsDir,'NIH\d{3}');

% Created by Samantha N. Jackson //2019
% 
% 2/10/2020 - SJ: Added in conditions for NIH070, 71, 72, 74, and 76, which were processed by SB and have
%              a different file structure than the previous subjects
% 
% 
% 
% 
% For NIH070:76, this is the file structure:
%   ORIG:                         ACTUAL:
%   /raw/spikeInfo.mat            /spike/outputs/spikeInfo.mat          
%   /sorting/sortSummary.csv      /spike/outputs/sortSummary.csv
%   /sorting/sortFigs             /spike/outputs/sortFigs

up_micro = 'micro';
up_mountain_sorts = 'mountain_sorts';
micro_ext = 'micro_pickSess2Extract.xlsx';

dirs2copyfrom2update = {'raw'         , 'sorting'};
% Use regular expressions below
stuff2copy2update =    {{'.*_spikeInfo.*\.mat'},{'sortFigs', '.*sortSummary.*\.csv'}};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Subject ' subj ' Mountain Sorts preparing to be copied ... ']);

mt_subpath = [rootMountainDir filesep subj];

if ~exist(mt_subpath,'dir')
    fprintf('%s\n',['Directory ' mt_subpath ' does not exist! Continue to return and move to the next subject.']);
    keyboard
    done = 0;
    return
end

up_subpath = [rootupdatedir filesep subj];
if ~exist(up_subpath,'dir')
    fprintf('%s\n',['Directory ' up_subpath ' does not exist!']);
    keyboard
end

up_mountain_sorts_path = [up_subpath filesep up_micro filesep up_mountain_sorts];
up_mountain_sorts_micro_ext_path = [up_subpath filesep up_micro filesep micro_ext];
extTab = readtable(up_mountain_sorts_micro_ext_path);
%Find which sessions are to extract
ext_sess_list = extTab.folderName(strcmp(extTab.toExtract,'y'));
ext_sess_list_onlynums = regexp([ext_sess_list{:}],'\d{6}_\d{4}','match')';

if ~exist(up_mountain_sorts_path,'dir') %Make the 'mountain_sorts' folder in 'micro' folder in eeg/.update
    
    fprintf('%s\n',[up_mountain_sorts_path ' does not exist. Making it now and preparing to copy.']);
    mkdir(up_mountain_sorts_path);
    
    
else
    fprintf('%s\n',[up_mountain_sorts_path ' already exists!!! WHY? Function not able to copy properly yet if you continue.']);
    keyboard
end

mt_sesslist_orig = getDirNamesRegexp(mt_subpath,'\d{6}_.*')'; %ignore the '_swarms' folder in NIH070:76
mt_sesslist_use = mt_sesslist_orig(contains(mt_sesslist_orig,ext_sess_list_onlynums));
copy_list = cell(1,2);

for sess = 1:numel(mt_sesslist_use)
    mt_session_path = [mt_subpath filesep mt_sesslist_use{sess}];
    up_mountain_sorts_sess_path = [up_mountain_sorts_path filesep mt_sesslist_use{sess}];
    if ~exist(up_mountain_sorts_sess_path,'dir') %Make the session folder in eeg/.update
        fprintf('%s\n',[up_mountain_sorts_sess_path ' does not exist. Making it now and preparing to copy.']);
        mkdir(up_mountain_sorts_sess_path);
    else
        fprintf('%s\n',[up_mountain_sorts_sess_path ' already exists!!! WHY? Function not able to copy properly yet if you continue.']);
        %keyboard
    end    
    
    for dd = 1:numel(dirs2copyfrom2update)
        mt_dir_path = [mt_session_path filesep dirs2copyfrom2update{dd}];
        up_dir_path = [up_mountain_sorts_sess_path filesep dirs2copyfrom2update{dd}];
        if ~exist(mt_dir_path,'dir')
            fprintf('%s\n',[mt_dir_path ' does not exist. Continue to skip?']);
            keyboard     
            continue
        end
        if ~exist(up_dir_path,'dir')
            fprintf('%s\n',[up_dir_path ' does not exist. Making it now and preparing to copy.']);
            mkdir(up_dir_path);
        end   
        for ss = 1:numel(stuff2copy2update{dd})
            stuff_name = char(getDirNamesRegexp(mt_dir_path,stuff2copy2update{dd}{ss}));
            if isempty(stuff_name)
                fprintf('%s\n',[stuff_name ' does not exist. Continue to skip?']);
                keyboard     
                continue
            end
            mt_stuff_path = [mt_dir_path filesep stuff_name];
            up_stuff_path = [up_dir_path filesep stuff_name];
            if isempty(copy_list{1})
                copy_list(1,:) = {mt_stuff_path, up_stuff_path}; % Don't append if copy_list is empty
            else
                copy_list(end+1,:) = {mt_stuff_path, up_stuff_path}; % Append to list; can initialize better in future
            end
        end
    end
end

%Now copy

for cc = 1:size(copy_list,1)
    copy_from = copy_list{cc,1};
    copy_to = copy_list{cc,2};
    
    if exist(copy_from,'file')
        fprintf('%s\n',['Copying: ' copy_from ' -> ' copy_to]);
        if exist(copy_to)
            fprintf('%s\n','ERROR!!!! copy_to file/dir already exists!!!!');
            keyboard
        else 
            st1 = copyfile(copy_from,copy_to);

            if st1 == 1
                fprintf('\t%s\n','Success');
            else
                fprintf('%s\n','ERROR!!!! Copying was not successful!');
                keyboard                
            end
        end
    else
        fprintf('%s\n',['ERROR!!!! ' copy_from 'does not exist! How can it copy?']);
        keyboard
    end
    
end
done = 1;
end




