function [table2use] = getRawFileData(subj, FRNUpath)
% 
% Get task and notes from rawFileList to be used in updateUtahInfoCSV_v10.m
% 
% INPUTS:
%       subj:       NIH000
%       FRNUpath:   Path to FRNU (where rawFileList will be read from)
%                   (default: '/Volumes/Shares/FRNU/data/eeg')
% OUTPUT:
%       table2use:  A #sessions-by-3 table with a subject's micro sessions,
%                   tasks from rawFileList, and notes from rawFileList
%
% This function will work with the older versions of rawFileList (.txt
% files), which do not have notes associated with them. getRawFileList will
% search through each rawFileList present in the subject's directory and 
% use the one with the most tasks, followed by notes.
% 
% 
% Created by Samantha Jackson
%
% 11/2019 - Created by SJ
% 

if nargin < 2
    FRNUpath = '/Volumes/Shares/FRNU/data/eeg';
end


rawFiles = getDirNamesRegexp([FRNUpath filesep subj filesep 'raw'],'^rawFile.*');

% filepath, tweak?, csv txt or xlsx, #task, #notes
rawCell = cell(numel(rawFiles),5);

for ii = 1:numel(rawFiles)
    rawfilepath = [FRNUpath filesep subj filesep 'raw' filesep rawFiles{ii}];
    % ---- #1: Filepath ----
    rawCell{ii,1} = rawfilepath;
    if contains(rawfilepath,'tweak')
        tw = 1;
    else
        tw = 0;
    end
    % ---- #2: Tweak? ----
    rawCell{ii,2} = tw;
    fileender = char(regexp(rawfilepath,'\..*','match')); %file ender
    % ---- #3: File ender ----
    rawCell{ii,3} = fileender;
    
    if tw == 1 || strcmp(fileender,'.txt')
        % If old version (text file), need to loop through each line, convert date
        % to folder name format, and find the task listed next to date
        if strcmp(fileender,'.txt')
            fid = fopen(rawfilepath);
            tline = '';
            T_cell = cell(1,3); % 1) session folder name, 2) task, 3) notes
            while ~contains(tline,'filename') && ~feof(fid) % Get past header of text file 
                tline = fgetl(fid);
            end
            ct = 0;
            while ~isempty(tline) && ~feof(fid) % Stop when you get to an empty line or end of the file
                tline = fgetl(fid);
                ct = ct + 1;
                %first get file name
                full_date = char(regexp(tline,'(?<=\[)\d+\/.*(?=\])','match'));
                day_str = char(regexp(full_date,'\d+\/\d+\/\d+','match'));
                day_cell = strsplit(day_str,'/');
                if isempty(tline) % Sometimes there are a few empty lines at the end of the file, just make sure they aren't elsewhere
                    fprintf('\n%s\n','Empty line found in rawFileList.txt (probably at end of file). Please investigate and continue if okay.');
                    keyboard
                    continue
                elseif strcmp(day_cell,'') % At least one case (NIH029) where there is one line at the end with date information but not in the correct format
                    fprintf('\n%s\n','Issue with finding date information for a line in rawFileList.txt. Please investigate and continue if okay.');
                    fprintf('%s\n',['Offending line: ' tline]);
                    keyboard;
                    T_cell{ct,1} = 'XXXXXX_XXXX';
                    T_cell{ct,3} = '-'; %for notes, never in .txt file
                    T_cell{ct,2} = '-';
                    continue
                end
                % Rearrange the day string to fit the session folder name
                % format
                day_cell{3} = day_cell{3}(end-1:end);
                day_cell_new = pad(day_cell,'left','0');
                day_str_new = [day_cell_new{3} day_cell_new{1} day_cell_new{2}];
                
                % Rearrange the time string to fit the session folder name
                % format
                time_str = char(regexp(full_date,'\d{2}\:\d{2}\:\d{2}','match'));
                time_cell = strsplit(time_str,':'); %no need to pad because already correct length
                time_str_new = [time_cell{1} time_cell{2}];
                
                % Put day and time together...
                sess_name = [day_str_new '_' time_str_new];
                T_cell{ct,1} = sess_name;
                
                % Now get task
                task_name = regexp(tline,'(?<=<< )\w*.*$','match');
                if ~isempty(task_name)
                    T_cell{ct,2} = task_name;
                else
                    T_cell{ct,2} = '-';
                end
                
                T_cell{ct,3} = '-'; %for notes, never in .txt file
                
            end
            fclose(fid);
            
            T_short = cell2table(T_cell,'VariableNames',{'folderName','task','notes'});
            num_notes = 0;
            
        else % Tweak file for .csv
            T = readtable(rawfilepath);
            T_short = [T(:,'folderName') T(:,'task') T(:,'notes')];
            num_notes = sum(~strcmp(T_short.notes,'-'));
        end
        num_tasks = sum(~strcmp(T_short.task,'-'));
        
        % ---- #4: Tasks ----
        rawCell{ii,4} = num_tasks;
        % ---- #5: Notes ----
        rawCell{ii,5} = num_notes;
        % ---- #6: Entire Table ----
        rawCell{ii,6} = T_short; % Save whole table in cell
    else
        rawCell{ii,4} = 0;
        rawCell{ii,5} = 0;
    end
end

[maxTask,imaxTask] = max([rawCell{:,4}]); % Which rawFileList has the most tasks?
if maxTask > 0
    table2use = rawCell{imaxTask,6};
else % Now take the one with the most notes
    [maxTask2,imaxTask2] = max([rawCell{:,5}]); % Which rawFileList has the most notes?
    if maxTask2 > 0
        table2use = rawCell{imaxTask2,6};
    else
        table2use = [];
    end  
end

end