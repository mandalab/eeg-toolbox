function [done] = copy2EEG_LFP_Sorts(subj, rootEEGdir, rootupdatedir, rootPUBdir)

% Copy sort files of interest from PROC behavioral to PROC update eeg;
% also create sorts_log.log
%
% Inputs:
%   subj - string specifying the subject, ex) 'NIH060'
%   rootEEGdir - string specifying the root eeg directory that has your
%                   subject of interest in it 
%                   ex) '/Volumes/56PROC/micro_behavioral/micro_pristine'
%   rootupdatedir - string specifying the directory containing subjects where you will copy to,
%                   ex) '/Volumes/56PROC/eeg_processing/.update/eeg'
%   rootPUBdir - string specifying the directory where you want to copy the
%                   sort information (sort figs, spikeinfo, sorts, sort log)
%                   ex) '/Volumes/56PUB/readWrite/micro_forSorting'
% 
%
% Files that will be copied: (Otherwise update below)
%   To PROC (rootupdatedir):
%       lfp_noreref (whole directory) 
%       sortFigs
%       spikeInfo*.mat
%       sortNotes
%       sorts_*.txt
%       sort_log.log
%   To PUB (rootPUBdir):
%       sortFigs
%       spikeInfo*.mat       
%       sorts_*.txt
%       sort_log.log
% 
% Below are the current names of the folders this function will reference.
% They can be changed if these names change in the future. There is also a
% place to specify which spike sorting files/directories will get copied to PROC
% (rootupdatedir) vs PUB (rootPUBdir) which can also be changed if needed
% in the future. lfp_noreref will always get copied.
% 
% 2/13/2020 SJ: Updated to look at mod date on top of whether the file/folder already exists.
% 2/19/2020 SJ: Add option to overwrite newer files
%
%
%%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
up_micro = 'micro';
up_manualsort = 'manualsort';
up_lfp = 'lfp_noreref';
up_spk = 'spikes';

mp_sortsmanual = 'sorts_manual';

% Fill in a text string that will only be found in the file/directory you
% want copied
sortstuff2copy2PROC = {'sortFigs', ... % sortFigs
    '_spikeInfo*.mat', ... % ######_##_spikeInfo.mat
    'sortNotes', ... % sortNotes(######_####)_sortedByXX.xlsx
    'sorts_*.txt', ... % sorts_(complete).txt
    'sort_log.log'};  %sort_log.log

sortstuff2copy2PUB = {'sortFigs', ... % sortFigs
    '_spikeInfo*.mat', ... % ######_##_spikeInfo.mat
    'sorts_*.txt', ... % sorts_(complete).txt
    'sort_log.log'};  %sort_log.log

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Subject ' subj ' preparing to be copied ... ']);

if ~exist([rootupdatedir filesep subj],'dir')
    fprintf('%s\n',['Directory ' rootupdatedir filesep subj ' does not exist!']);
    keyboard
end

%cd([rootEEGdir filesep subj filesep mp_sortsmanual]) %Dont think we need this because all paths are from root

% Get lists of LFP sessions and sorts_manual sessions
LFPsessdir = dir([rootEEGdir filesep subj filesep 'processed_LFP' filesep '*_*']);
LFPsesslist = {LFPsessdir.name};
LFPsesslist = LFPsesslist(~ismember(LFPsesslist,'.DS_Store'));

SPKsessdir = dir([rootEEGdir filesep subj filesep 'processed_SPK' filesep '*_*']);
SPKsesslist = {SPKsessdir.name};
SPKsesslist = SPKsesslist(~ismember(SPKsesslist,'.DS_Store'));
mansortdir = dir([rootEEGdir filesep subj filesep mp_sortsmanual filesep '*_*']);
mansortlist = {mansortdir.name};
mansortlist =  mansortlist(~ismember(mansortlist,'.DS_Store'));

rerefdir = dir([rootEEGdir filesep subj filesep mp_sortsmanual filesep '*' filesep 'reref*_sortedBy*']);
rereflist = {rerefdir.name};
rereflist = rereflist(~ismember(rereflist,'.DS_Store'));
rereflist = strcat(mansortlist,filesep,rereflist);

% Return the union of LFP and sorts_manual sessions
ALLsesslist = union(LFPsesslist,mansortlist);

if numel(rereflist) ~= numel(mansortlist)
    fprintf('%s\n',[rootupdatedir filesep subj ' has an unequal number of sessions vs reref folders!']);
    keyboard
end

% Check update eeg for existing folders (spikes) and create them if not
spikes_path = strcat(rootupdatedir,filesep,subj,filesep,up_micro,filesep,up_manualsort,filesep, ...
    ALLsesslist,filesep,up_spk); % Should create all parent directories too

for ff = 1:numel(spikes_path)
    fold = spikes_path{ff};
    fold_sess = fold(regexp(fold,'(?<=\/)\d{6}_\d{4}.*(?=\/)'):regexp(fold,'(?<=\/)\d{6}_\d{4}.*(?=\/)','end'));
    if isempty(find(contains(SPKsesslist,fold_sess))) && isempty(find(contains(mansortlist,fold_sess))) % Session is in LFP but not SPK or mansort
        fold = fold(1:regexp(fold,'.*(?=\/spikes)','end'));
    end
    if ~exist(fold,'dir')
        fprintf('\t%s',[fold ' does not exist. Making it now ...']);
        status1 = mkdir(fold);
        if status1 == 1
            if (find(contains(SPKsesslist,fold_sess)) > 0) && isempty(find(contains(mansortlist,fold_sess))) % Session is in SPK list but hasn't been sorted yet (not in mansort)
                txtfile = fopen([fold filesep 'SPK_available_but_unsorted.txt'],'w');
                fclose(txtfile);
            end
            fprintf('%s\n',' done.');
        else
            fprintf('\n\t%s\n',' ERROR!!! Folder could not be made!!');
            keyboard
        end
    else
        if ~isempty(find(contains(mansortlist,fold_sess))) && exist([fold filesep 'SPK_available_but_unsorted.txt'],'file') % If there is stuff to copy in spikes but the txt file is already there
            delete([fold filesep 'SPK_available_but_unsorted.txt']);
        end
    end
end

% Make list of files/directories to copy from (col 1) and to (col 2)
copy_list = cell(1,2);
for ss = 1:numel(ALLsesslist)
    sess = ALLsesslist{ss};
    rootCopyTo = [rootupdatedir filesep subj filesep up_micro filesep up_manualsort];
    
    % LFP
    % Copy whole directory
    LFPpath = [rootEEGdir filesep subj filesep 'processed_LFP' filesep sess filesep up_lfp];
    if exist(LFPpath,'dir')
        LFPpath_to = [rootCopyTo filesep sess filesep up_lfp];
        if isempty(copy_list{1})
            copy_list(1,:) = {LFPpath, LFPpath_to}; % Don't append if copy_list is empty
        else
            copy_list(end+1,:) = {LFPpath, LFPpath_to}; % Append to list; can initialize better in future
        end
    else
        if sum(strcmp(ALLsesslist{ss},LFPsesslist))==1
            fprintf('\n\t%s\n',['ERROR!! ' LFPpath ' session is included but directory full path does not exist.']); 
            keyboard
        else
            fprintf('\n\t%s\n',[sess ' session exists for SPK but not for LFP.']); 
        end
    end
    
    % SPK
    SPKpath_proc = [rootEEGdir filesep subj filesep mp_sortsmanual filesep sess];
    if ~isempty(find(contains(rereflist,sess)))
        SPKpath_fold = [rootEEGdir filesep subj filesep mp_sortsmanual filesep rereflist{find(contains(rereflist,sess))}];
        [~,sortersess] = fileparts(SPKpath_fold);
        %sorter = char(regexp(rereflist{ss},'(?<=_sortedBy)\w+','match'));
        %cd(SPKpath)
        slog = fopen([SPKpath_fold filesep 'sort_log.log'],'a');
        num_sorts_dir = dir([SPKpath_fold filesep 'sort_txt' filesep '*.txt']);
        num_sorts = size({num_sorts_dir.name},2);
        compstatus_dir = dir([SPKpath_fold filesep 'sorts_*.txt']);
        compstatus = {compstatus_dir.name};
        if contains(compstatus,'INCOMPLETE')
            compstatus1 = 'INCOMPLETE';
        else
            compstatus1 = 'Complete';
        end
        fprintf(slog,'%s\t%s\n',datetime,['Processed with ' num2str(num_sorts) ' sorts. (' compstatus1 ')']);
        fclose(slog);
        
        % Fill in copy list with paths to PROC update
        for cc = 1:numel(sortstuff2copy2PROC)
            SPKtxt = sortstuff2copy2PROC{cc};
            SPKtxtdir = dir([SPKpath_fold filesep '*' SPKtxt '*']);
            if isempty(SPKtxtdir)
                SPKtxtfullname = 'dne';
            else
                SPKtxtfullname = SPKtxtdir.name;
            end
            SPKpath = [SPKpath_fold filesep SPKtxtfullname];
            if exist(SPKpath)
                SPKpath_to = [rootCopyTo filesep sess filesep up_spk filesep SPKtxtfullname];
                if isempty(copy_list{1})
                    copy_list(1,:) = {SPKpath, SPKpath_to};
                else
                    copy_list(end+1,:) = {SPKpath, SPKpath_to};
                end
            else
                if sum(strcmp(ALLsesslist{ss},SPKsesslist))==1
                    fprintf('\n\t%s\n',['ERROR!! ' SPKpath ' session is included but directory dne, SPKtxt2PROC = ' SPKtxt]); 
                    fprintf('\n\t%s\n','Continue if SPKtxt2PROC is sortFigs and spike info has [NoUnits] in its name.');
                    keyboard
                else
                    fprintf('\n\t%s\n',[sess ' session exists for LFP but not for SPK.']); 
                end
            end    
        end
        
        % Do the same for PUB
        if ~exist(rootPUBdir,'dir')
            fprintf('\n%s\n','LOOKS LIKE 56PUB IS NOT MOUNTED. PLEASE CONNECT TO 56PUB AND CONTINUE.');
            keyboard;
        end
        if ~exist([rootPUBdir filesep subj filesep sess '_(grabbed)' filesep sortersess],'dir')
            mkdir([rootPUBdir filesep subj filesep sess '_(grabbed)' filesep sortersess]);
        end
        for ccpub = 1:numel(sortstuff2copy2PUB)
            SPKtxt = sortstuff2copy2PUB{ccpub};
            SPKtxtdir = dir([SPKpath_fold filesep '*' SPKtxt '*']);
            if isempty(SPKtxtdir)
                SPKtxtfullname = 'dne';
            else
                SPKtxtfullname = SPKtxtdir.name;
            end
            SPKpath = [SPKpath_fold filesep SPKtxtfullname];
            if exist(SPKpath)
                SPKpath_to = [rootPUBdir filesep subj filesep sess '_(grabbed)' filesep sortersess filesep SPKtxtfullname];
                if isempty(copy_list{1})
                    copy_list(1,:) = {SPKpath, SPKpath_to};
                else
                    copy_list(end+1,:) = {SPKpath, SPKpath_to};
                end
            else
                if sum(strcmp(ALLsesslist{ss},SPKsesslist))==1
                    fprintf('\n\t%s\n',['ERROR!! ' SPKpath ' session is included but directory dne, SPKtxt2PUB = ' SPKtxt]); 
                    fprintf('\n\t%s\n','Continue if SPKtxt2PUB is sortFigs and spike info has [NoUnits] in its name.'); 
                    keyboard
                else
                    fprintf('\n\t%s\n',[sess ' session exists for LFP but not for SPK.']); 
                end
            end
        end
    else
        
    end
end

% Copy
owka = false; % Overwrite for all
ow = true; % Overwrite just once
if ~isempty(copy_list{1}) % If get BR file list grabbed non-micro data by accident
    for cp = 1:size(copy_list,1)
        copy_from = copy_list{cp,1};
        copy_to = copy_list{cp,2};
        % Check to see if summart txt file is incomplete vs complete
        if contains(copy_from,'sorts_') && contains(copy_from, '.txt')
            [~, summary_from_file, summary_from_ext] = fileparts(copy_from);
            summary_from = [summary_from_file summary_from_ext];
            [summary_to_path, ~] = fileparts(copy_to);
            summary_to_dir = dir([summary_to_path filesep 'sorts_*.txt']);
            if isempty(summary_to_dir) %If hasn't been run yet
                summary_to = summary_from;
            else
                summary_to = summary_to_dir.name; %If has already been run and you need to compare them
            end
            if ~strcmp(summary_from,summary_to) % If they don't have the same name - ie would write another instead of overwriting
                fprintf('\n\t%s\n',['Summary text has changed from ' summary_to ' to ' summary_from '.']);
                fprintf('\t%s\n',['Check ' summary_to_path filesep summary_to ' against ' copy_from ' and then continue if this is correct.']);
                keyboard
                fprintf('\t%s\n',['Deleting ' summary_to ' now. Will copy new summary text.']);
                to_delete = [summary_to_path filesep summary_to];
                delete(to_delete)
            end
        end

        %ow = true;
        if exist(copy_from)
            if exist(copy_to)
                if ~owka
                    if strcmp(newerFile(copy_to,copy_from),copy_from) %If your source file is newer than your destination (need to still copy)
                        fprintf('\t%s\n',[copy_to ' already exists but is older, overwriting.']);
                        ow = true;
                        owka = false;
                    elseif isempty(newerFile(copy_to,copy_from)) %Same mod date, do not overwrite
                        fprintf('\t%s\n',[copy_to ' has same mod date, not overwriting.']);
                        ow = false;
                        owka = false;
                    elseif contains(upper(copy_to),'LFP') %Not re-copying LFPs if they are already there % destination file is newer that the source?????
                        fprintf('\n\t%s\n',['LFP folder already exists. Not re-copying.']);
                        ow = false;
                        ow = false;
                    else
                        fprintf('\n\t%s\n\t',['ERROR!!!' copy_to ' already exists and is newer than the source file!!!. WHY????']);
                        oworno = input('Should you overwrite (o) or skip (s)? ','s');
                        if strcmp(upper(oworno),'O')
                            ow = true;
                        else
                            ow = false;
                        end
                        keyboard
%                         ow_ans = input('\nWould you like to overwrite this one time (O), overwrite all existing directories/files (A), keep all as is (K), or stop (S) [default]? \n\n','s');
%                         fprintf('\n\n');
%                         if strcmp(ow_ans,'O') || strcmp(ow_ans,'o')
%                             ow = true;
%                             owka = false;
%                             fprintf('\t%s\n',' ... will overwrite.');
%                         elseif strcmp(ow_ans,'A') || strcmp(ow_ans,'a')
%                             ow = true;
%                             owka = true;
%                             fprintf('\t%s\n',' ... will overwrite.');
%                         elseif strcmp(ow_ans,'K') || strcmp(ow_ans,'k')
%                             ow = false;
%                             owka = true;
%                         else
%                             return
%                         end
                    end
                else
                    %fprintf('\t%s\n',' ... will overwrite.');
                end
                %fprintf('\n\n');
            end

            if ~exist(copy_to) || ow
                fprintf('\t\t%s',[copy_from ' -> ' copy_to]);
                status2 = copyfile(copy_from, copy_to);
                if status2 == 1
                    fprintf('%s\n',' ... done.');
                else
                    fprintf('%s\n','ERROR!!!! Not copied!!!!');
                end
            end
        else
            fprintf('%s\n',['ERROR!!! ' copy_from ' does not exist!!!!!!']);
        end
    end
else
    fprintf('%s\n','No micro data! copy_list is empty. Verify that this is correct and that getBRrawfilelist grabbed only non-micro data.');
    keyboard;
end

fprintf('\n\n%s\n',['SUBJECT ' subj ' COMPLETE!']);
done = 1;

end