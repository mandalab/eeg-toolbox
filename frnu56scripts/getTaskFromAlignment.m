function [table2use] = getTaskFromAlignment(subj, rootupdatedir)
% 
% Used in updateUtahInfoCSV_v10.m to get the task associated with each micro session from the alignment
% data 
% 
% 
% 

% 1/10/2020 - Created by Samantha Jackson
% 2/18/2020 - SJ: added functionality for task session names that contain non-alpha-numeric symbols
% 3/9/2020 - SJ: changed to work with new alignmentToSpikes_manual.csv format
%
%
% COMMENT THIS OUT:
% rootupdatedir = '/Volumes/56PROC/eeg_processing/.update/eeg';
% subj = 'NIH050';

%
up_micro = 'micro';
up_manualsort = 'manualsort';



if nargin < 2
    rootupdatedir = '/Volumes/56PROC/eeg_processing/.update/eeg';
end


manualsort_path = fullfile(rootupdatedir,subj,up_micro,up_manualsort);
alignmentfile = getDirNamesRegexp(manualsort_path,'alignment.*Spikes.*');
if numel(alignmentfile) == 1
    alignmentfile_path = fullfile(manualsort_path,char(alignmentfile));
elseif numel(alignmentfile) > 1
    fprintf('\n%s',['ERROR!!! ' num2str(numel(alignmentfile)) ' alignment files. Should be exactly 1 in ' manualsort_path]);
    keyboard    
else
    table2use = [];
    return
end

alignment_table = readtable(alignmentfile_path);

micro_task_list = cell(1,2);
% Loop through each micro session that exists in the alignment table and return task name
% Therefore, only returns the micro sessions that were successfully aligned (keep that in mind)
ct = 0;

% SJ- for older forms of the alignment table, just leave it be and don't fill in alignment stuff
if ~any(strcmp('micro_sess',alignment_table.Properties.VariableNames))
    keboard
end
    
for ii = 1:numel(alignment_table.micro_sess)
    %micro_full_path = alignment_table.micro_sess{ii};
    %Note: in regexp, need greedy expression so it captures 'beh' and 'premie' without the rest of the
    %expression
    micro_sess = alignment_table.micro_sess{ii};
    %micro_sess = char(regexp(micro_full_path,['(?<=' up_manualsort '\' filesep ')\d{6}_\d{4}.*?(?=\' filesep ')'],'match'));
    % Now get task associated with this micro session
    task_full_path = alignment_table.beh_sess{ii};
    if strcmp(task_full_path,'-') 
        continue      
    elseif strcmp(micro_sess,'no match') || strcmp(micro_sess,'-')
        %fprintf('%s\n','Hey, wait a minute, the micro path no match or -, but task is not - ???');
        %keyboard
    elseif ~isempty(regexp(task_full_path,['[^a-zA-Z_0-9' filesep ']'],'once'))
        fprintf('%s\n%s\n','ERROR!!! The behavioral session path in the alignment table contains a non-alpanumeric + underscore character other than filesep!', ...
            'Please investigate this. Regexp will not grab the correct task name because of this otherwise.');
        task = char(regexp(task_full_path,['(?<=behavioral\' filesep ')\w*?(?=\' filesep ')'],'match'));
        fprintf('%s\n',['Parsed task is: ' task '. Continue if this is correct. Fix it somehow if not.']);
        keyboard
        ct = ct+1;
        micro_task_list{ct,1} = micro_sess;
        micro_task_list{ct,2} = task;    
    else
        ct = ct+1;
        task = char(regexp(task_full_path,['(?<=behavioral\' filesep ')\w*?(?=\' filesep ')'],'match'));
        micro_task_list{ct,1} = micro_sess;
        micro_task_list{ct,2} = task;        
    end
end

% micro_task_list contains repeats since it occurs for each session. Now, only take unique, but check if
% there are differing task names for different sessions (there can be). In this case, make the task be
% called "task1&task2"

if isempty(micro_task_list{1,1}) %all alignment failed
    table2use = [];
    return
end

[micro_list_unique,ia,ic] = unique(micro_task_list(:,1));
micro_task_list_unique = cell(1,2);

for jj = 1:numel(micro_list_unique)
    thissess = micro_list_unique{jj};
    thissess_tasks = strjoin(unique(micro_task_list(ic==jj,2)),'&'); %Get unique list of tasks associated with this session and join with &
    micro_task_list_unique{jj,1} = thissess;
    micro_task_list_unique{jj,2} = thissess_tasks;
end

table2use = cell2table(micro_task_list_unique,'VariableNames',{'folderName','task'});


end