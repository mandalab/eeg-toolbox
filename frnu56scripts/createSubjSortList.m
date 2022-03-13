function status = createSubjSortList(rootPUBdir,rootUPDATEdir)
% Create summary file for sorted data located in /Volumes/56PUB/readWrite/micro_forSorting
% 
% Creates the file sortBuild_######.xlsx, which provides a summary of
% complete vs incomplete sorts per subject and session
% Inputs:
%       rootPUBdir: PUB directory where the sorts are created and grabbed
%                   ex) rootPUBdir = '/Volumes/56PUB/readWrite/micro_forSorting';
% 
%       rootUPDATEdir: .update directory where the file should also be present
%                   ex) rootUPDATEdir = '/Volumes/56PROC/eeg_processing/.update/eeg';
% 
% Output:
%       status: status = 1 when everything is successful, 0 when not
% 
% Summary of what fields are present in sortBuild_######.xlsx:
%   Sheet ALL_SESS:
%       subject: NIHXXX or TREXXX
%       session: subject's session (either grabbed or not)
%       sort_attempted: if a sort has been attempted and grabbed
%       sort_complete: if the sort has been deemed complete by the sort summary text file
%       num_sorts_missing: number of individual sorts deemed incomplete by the sort summary text file
%       sorter: initials of the sorter
%       notes: notes, such as if there are no units
%   Sheet PER_SUBJ_SUMMARY:
%       subject: NIHXXX or TREXXX
%       num_sess_total: total number of sessions for that subject in PUB
%       num_sess_unsorted: number of sessions that have not been grabbed
%       num_sess_sorted: number of sessions that have been grabbed
%       num_sess_sorted_complete: number of sessions that have been grabbed and deemed complete by the sort summary text file
%       num_sess_sorted_incomplete: number of sessions that have been grabbed but deemed incomplete by the sort summary text file
%       task_aligned:
%       rawFile_task_note:
% 
% Created by Samantha Jackson
%   
% 10/2019 - Created by SJ
% 10/21/2019 - Added in columns: num_sorted_channels, num_units_total, and
%               median_units_per_sorted_chan
% 1/10/2019  - Added in columns: task_aligned & rawFile_task&note
% 2-3/2021   - Removed rawFile_task&note, added in alignment success
%
% Add in case where a session folder isn't either nothing or grabbed? (ready)
% 
% rootPUBdir = '/Volumes/Seagate Backup Plus Drive/local56/56PUB/readWrite/micro_forSorting';
% rootUPDATEdir = '/Volumes/Seagate Backup Plus Drive/local56/56PROC/eeg_processing/.update/eeg';
%rootPROCdir = '/Volumes/Seagate Backup Plus Drive/local56/56PROC/micro_behavioral/micro_pristine';

rootPROCdir = '/Volumes/56PROC/micro_behavioral/micro_pristine';

% Check first to make sure you have the right directory... 
if ~strcmp(rootPUBdir,'/Volumes/56PUB/readWrite/micro_forSorting')
    fprintf('%s\n',['CAUTION!! Your rootPUBdir is set to: ' rootPUBdir])
    fprintf('%s\n','In most cases it should be set to: /Volumes/56PUB/readWrite/micro_forSorting');
    fprintf('%s\n','Continue only if you are sure what you are doing.');
    keyboard;
end
if ~strcmp(rootUPDATEdir,'/Volumes/56PROC/eeg_processing/.update/eeg')
    fprintf('%s\n',['CAUTION!! Your rootUPDATEdir is set to: ' rootUPDATEdir])
    fprintf('%s\n','In most cases it should be set to: /Volumes/56PROC/eeg_processing/.update/eeg');
    fprintf('%s\n','Continue only if you are sure what you are doing.');
    keyboard;
end

outfile = ['sortBuild_' datestr(now,'yymmdd') '.xlsx'];
outpathPUB = [rootPUBdir filesep outfile];
outpathUPDATE = [rootUPDATEdir filesep outfile];

remem = false;
if exist(outpathPUB,'file')
    fprintf('%s\n', ['File already exists: ' outpathPUB '. Deleting it now.']);
    delete(outpathPUB);
elseif ~isempty(getDirNamesRegexp(rootPUBdir,'sortBuild_\d{6}.*'))
    fprintf('%s\n', 'File in PUB exists with an older date. Remember to delete the older version after the new one is created!');
    remem = true;
end
if exist(outpathUPDATE,'file')
    fprintf('%s\n', ['File already exists: ' outpathUPDATE '. Deleting it now.']);
    delete(outpathUPDATE);
elseif ~isempty(getDirNamesRegexp(rootUPDATEdir,'sortBuild_\d{6}.*'))
    fprintf('%s\n', 'File in UPDATE exists with an older date. Remember to delete the older version after the new one is created!');
    remem = true;
end
    
cd(rootPUBdir)

subjList = getDirNamesRegexp(rootPUBdir,'((NIH)|(TRE))\d{3}');
sessCT = 0;
for ss = 1:numel(subjList) % Iterate over subjects
    subj = subjList{ss};
    subpath = [rootPUBdir filesep subj];
    sesslist = getDirNamesRegexp(subpath,'\d{6}_\d{4}');
    fprintf('%s\n',['Subject ' subj ': (' num2str(ss) '/' num2str(numel(subjList)) ')']);
    fprintf('\t%s','Sessions processed: ');
    for ss2 = 1:numel(sesslist) % Iterate over session folders
        sessCT = sessCT + 1; % Increase count per session for new lines in excel file
        sess = sesslist{ss2};
        sesspath = [subpath filesep sess];
        % Initialize data values
        sort_attempted = 0;
        sorter = '';
        sort_complete = NaN;
        num_sorted_channels = NaN;
        num_sorts_missing = NaN;
        num_units_total = NaN;
        median_units_per_sorted_chan = NaN;
        notes = '';
        task_aligned = '';
        alignment_success = '';
        %rawFile_task_note = '';
        
        fprintf('%s',' ... ');
        
        % COLUMN 3: sort_attempt
        if contains(sess,'grabbed')
            sort_attempted = 1;
            rerefdir = dir([sesspath filesep 'reref*_sortedBy*']);
            if numel(rerefdir) == 1
                rerefname = rerefdir.name;
                rerefpath = [sesspath filesep rerefname];
                
                % COLUMN 4: sort_complete
                    % Sort Summary
                sortsumdir = dir([rerefpath filesep 'sorts_*.txt']);
                if numel(sortsumdir) ~= 1
                    fprintf('\n%s\n',['ERROR!!!! ' numel(sortsumdir) ' sort summary text files. Supposed to have 1.']);
                    keyboard;
                end
                sortsumname = sortsumdir.name;
                
                
                if contains(sortsumname,'INCOMPLETE')
                    sort_complete = 0;
                    % COLUMN 6: num_sorts_missing
                    % Read the sort summary text file and count the number
                    % of micro arrays listed under each grade
                    sortsumID = fopen(fullfile(sortsumdir.folder,sortsumname),'r');
                    sortsumtxt = textread(fullfile(sortsumdir.folder,sortsumname),'%s','delimiter','\n');
                    fclose(sortsumID);
                    inc_grade = find(~cellfun(@isempty, regexp(sortsumtxt,'Grade: \w')),1); % Find first instance of "Grade: " and look below
                    sortsumtxt2 = sortsumtxt(inc_grade:end);
                    num_sorts_missing = sum(~cellfun(@isempty, regexp(sortsumtxt2,'[^(Grade: \w)]'))); % Disclude any other occurrences of "Grade: "

                elseif contains(sortsumname,'(complete)')
                    num_sorts_missing = 0;
                    sort_complete = 1;
                else
                    fprintf('\n%s\n','ERROR!!!! Sort summary name does not contain INCOMPLETE or complete)!');
                end
                
                % Nagivate to 56PROC/micro_behavioral/micro_pristine to get
                % sort data
                % COLUMN 5: num_sorted_channels
                sess_PROC = char(regexp(sess,'.*(?=_\(grabbed.*)','match'));
                PROCtxtpath = [rootPROCdir filesep subj filesep 'sorts_manual' filesep sess_PROC filesep rerefname filesep 'sort_txt'];
                txtdir = dir([PROCtxtpath filesep '*.txt']);
                txtnames = {txtdir.name};
                num_sorted_channels = numel(txtnames);
                
                % COLUMN 7: num_units_total
                num_units_list = zeros(1,num_sorted_channels);
                for sc = 1:num_sorted_channels % Iterate through sort .txt files
                    thischan = txtnames{sc};
                    txttable = readtable([PROCtxtpath filesep thischan]);
                    if isempty(txttable) % for [NoUnit] cases
                        num_units = 0;
                    else
                        num_units = numel(unique(txttable.Var1))-1; % subtract noise (unit 0)
                    end
                    num_units_list(sc) = num_units;
                end
                num_units_total = sum(num_units_list);
                
                % COLUMN 8: median_units_per_sorted_chan
                ave_units_per_sorted_chan = mean(num_units_list); % Not reported now, can just do num_units_total/num_sorted_channels in excel file
                median_units_per_sorted_chan = median(num_units_list);
                
                % COLUMN 9: sorter
                sorter = char(regexp(rerefname,'(?<=.*_sortedBy)\w+','match'));
                
                % COLUMN 10: notes
                % So far only 1 note - for if no units
                    % Spike Info
                spkinfodir = dir([rerefpath filesep '*spikeInfo*.mat']);
                if numel(spkinfodir) ~= 1
                    fprintf('\n%s\n',['ERROR!!!! ' numel(spkinfodir) ' spike infos. Supposed to have 1.']);
                    keyboard;
                end                

                    % Notes column
                if contains(spkinfodir.name,'NoUnits')
                    notes = 'No Units sorted. There could be an issue with the array. Check sortNotes to verify.';
                end
                

                
                
            else
                fprintf('\n%s\n',['ERROR!!!! Reref dir has ' numel(rerefdir) ' elements. Supposed to have 1.']);
                keyboard;
            end
        else % (doesn't have _grabbed)
            sort_attempted = 0;
        end
        
        
        
%         path_to_micro_Pick = fullfile(rootUPDATEdir,subj,'micro','micro_PickSess2Extract.xlsx');
%         if ~exist(path_to_micro_Pick,'file')
%             fprintf('%s\n',['ERROR!!! Cannot locate: ' path_to_micro_Pick]);
%             keyboard
%             continue
%         end
%         ext_table = readtable(path_to_micro_Pick);
        
        path_to_alignment_file = fullfile(rootUPDATEdir,subj,'micro','manualsort','alignmentToSpikes_manual.xlsx');
        if ~exist(path_to_alignment_file,'file')
            fprintf('%s\n',['ERROR!! Cannot locate: ' path_to_alignment_file]);
            keyboard
            continue
        end
        align_table = readtableSafe(path_to_alignment_file);

        % Removed grabbed or ready from name to search in micropick2extract file
        if contains(sess,'grabbed') || contains(sess,'ready')
            sess_ext = char(regexp(sess,'\d{6}_\d{4}.*?(?=_(\(grabbed|ready))','match'));
        else
            sess_ext = sess;
        end
        % COLUMN 11: task_aligned
        %task_aligned = ext_table.task{strcmp(ext_table.folderName,sess_ext)};
        % SJ: edited to reflect task in alignment table instead
        matching_beh_sess = align_table.beh_sess(strcmp(align_table.micro_sess,sess_ext));
        task_aligned_list = regexp(matching_beh_sess(:),'(?<=behavioral\/).*(?=\/)','match');
        if numel(task_aligned_list) > 1
            task_aligned_list_unique = unique([task_aligned_list{:}]);
            if numel(task_aligned_list_unique) > 1
                task_aligned = strjoin(task_aligned_list_unique,'&');
            else
                task_aligned = char(task_aligned_list_unique);
            end
        elseif numel(task_aligned_list) < 1
            %keyboard % micro sess is not in alignment table, so no match, make sure this works
            task_aligned = 'no match';
        else
            task_aligned = char(task_aligned_list{:});
        end
        
        % COLUMN 12: rawFile_task_note
        %rawFile_task_note = ext_table.notes{strcmp(ext_table.folderName,sess_ext)};
        % SJ: edited to reflect alignment success in alignment table instead
        alignment_success_list = align_table.alignmentSuccess(strcmp(align_table.micro_sess,sess_ext));
        if numel(alignment_success_list) > 1
            alignment_success_list_unique = unique(alignment_success_list);
            if numel(alignment_success_list_unique) > 1
                keyboard %This can't happen...right??
                alignment_success = strjoin(alignment_success_list_unique,'&');
            else
                alignment_success = char(alignment_success_list_unique);
            end
        elseif numel(alignment_success_list) < 1
            %keyboard % micro sess is not in alignment table, so no match, make sure this works
            alignment_success = 'alignment not attempted';
        else
            alignment_success = char(alignment_success_list);
        end
        % Just get the whole row.......?
        
        % Combine data values for all sessions from all subjects
        subj_sessLIST{sessCT,1} = subj; % 1
        subj_sessLIST{sessCT,2} = sess; % 2
        
        sort_attempted_complete_missingLIST(sessCT,1) = sort_attempted; % 3
        sort_attempted_complete_missingLIST(sessCT,2) = sort_complete; % 4
        sort_attempted_complete_missingLIST(sessCT,3) = num_sorted_channels; % 5
        sort_attempted_complete_missingLIST(sessCT,4) = num_sorts_missing; % 6
        sort_attempted_complete_missingLIST(sessCT,5) = num_units_total; % 7
        sort_attempted_complete_missingLIST(sessCT,6) = median_units_per_sorted_chan; % 8
        
        sorter_notesLIST{sessCT,1} = sorter; % 9
        sorter_notesLIST{sessCT,2} = notes; % 10
        sorter_notesLIST{sessCT,3} = task_aligned; % 11
        sorter_notesLIST{sessCT,4} = alignment_success; % 12
        %sorter_notesLIST{sessCT,4} = rawFile_task_note; % 12
        
        fprintf('%s',[' [' num2str(ss2) '/' num2str(numel(sesslist)) '] ']);
    end
    fprintf('\n');
end

% Convert cell arrays to table (T) and then back into one (Tcell) 
% ** Go back and change this
T = [cell2table(subj_sessLIST,'VariableNames',{'subject', 'session'}), ...
    array2table(sort_attempted_complete_missingLIST,'VariableNames',{'sort_attempted', 'sort_complete', 'num_sorted_channels', 'num_sorts_missing', 'num_units_total', 'median_units_per_sorted_chan'}), ...
    cell2table(sorter_notesLIST,'VariableNames',{'sorter', 'notes', 'task_aligned', 'alignment_status'})];

% Initialize variables for sheet 2
num_sess_total = NaN(numel(subjList),1);
num_sess_unsorted = NaN(numel(subjList),1);
num_sess_sorted = NaN(numel(subjList),1);
num_sess_sorted_complete = NaN(numel(subjList),1);
num_sess_sorted_incomplete = NaN(numel(subjList),1);

% Use table from sheet 1 to compute sheet 2 values:
for ii = 1:numel(subjList)
    subj = subjList{ii};
    num_sess_total(ii) = nansum(strcmp(T.subject,subj));
    num_sess_unsorted(ii) = nansum(isnan(T.sort_complete(strcmp(T.subject,subj))));
    num_sess_sorted(ii) = nansum(T.sort_attempted(strcmp(T.subject,subj)));
    num_sess_sorted_complete(ii) = nansum(T.sort_complete(strcmp(T.subject,subj)));
    num_sess_sorted_incomplete(ii) = nansum(T.sort_complete(strcmp(T.subject,subj))==0);
    if (num_sess_sorted_incomplete(ii) + num_sess_sorted_complete(ii) ~= num_sess_sorted(ii)) || (num_sess_unsorted(ii) + num_sess_sorted(ii) ~= num_sess_total(ii))
        fprintf('\n%s\n',['ERROR!!!! Numbers do not add up for ' subj '!!! Check table: ']);
        disp(T(strcmp(T.subject,subj),:))
        keyboard;
    end
end
subject = subjList';
T2 = table(subject, num_sess_total, num_sess_unsorted, num_sess_sorted, num_sess_sorted_complete, num_sess_sorted_incomplete);

% Convert both sheet 1 and sheet 2 data back to cell array to write to excel file
Tcell = [T.Properties.VariableNames;table2cell(T)];
T2cell = [T2.Properties.VariableNames;table2cell(T2)];

% Name the sheets
sh1 = 'ALL_SESS';
sh2 = 'PER_SUBJ_SUMMARY';

% Write to excel file
% Note: don't use writetable if you have specific sheet names b/c it will
% append your sheets after blank Sheets 1-3
fprintf('%s\n', ['Now writing excel file to ' outpathPUB ' ... ']);
status1 = xlwrite(outpathPUB,Tcell,sh1,'A1');
status2 = xlwrite(outpathPUB,T2cell,sh2,'A1');
fprintf('%s\n','DONE.');
%writetable(T,outfile,'Sheet',sh1)
%first_line_date = {sprintf('Created %s', datestr(now))};
%header = {"subject", "session", "sort_attempted", "sort_complete", "num_sorts_missing", "sorter", "notes"};
% strt2 = [char('A'+size(subj_sessLIST,2)) '2'];
% strt3 = [char(strt2(1) + size(sort_attempted_complete_missingLIST,2)) '2'];
% status1 = xlwrite(outfile,header,sh1,'A1'); % Write header
% status2 = xlwrite(outfile,subj_sessLIST,sh1,'A2'); % Write subject and session
% status3 = xlwrite(outfile,sort_attempted_complete_missingLIST,sh1,strt2);
% status4 = xlwrite(outfile,sorter_notesLIST,sh1,strt3);

% Check to make sure sheets 1 and 2 wrote correctly and copied to UPDATE
if exist(outpathPUB,'file') && (status1+status2 == 2) % PUB file exists, both sheets written
    fprintf('\n%s\n', ['Success. ' outpathPUB ' exists and was written properly.']);
    fprintf('%s', ['Now copying to ' outpathUPDATE ' ... ']);
    copyfile(outpathPUB,outpathUPDATE); % Copy PUB file to UPDATE
    fprintf('%s\n','DONE.');
    
    if exist(outpathUPDATE,'file') % UPDATE file exists (copied from PUB)
        status = 1;
        fprintf('%s\n', ['Success. ' outpathUPDATE ' was copied correctly.']);
    else                    % UPDATE file does not exist (not copied from PUB correctly)
        status = 0;
        fprintf('\n%s\n', ['ERROR!!! ' outpathUPDATE ' does not exist. It was not copied over correctly.']);
        keyboard;
    end
    
elseif exist(outpathPUB,'file') && (status1+status2 ~= 2) % PUB file exists but either or both sheets were not written properly
    status = 0;
    fprintf('\n%s\n', ['ERROR!! ' num2str(status1+status2) ' write(s) to ' outpathPUB '. Not performed properly.']);
    keyboard;
else                        % PUB file does not exist! Never wrote anything
    status = 0;
    fprintf('\n%s\n', ['ERROR!! ' outpathPUB ' does not exist. It was not created properly.']);
    keyboard;    
end

if remem % Don't forget to delete older sortBuilds!
    fprintf('\n\n%s\n', 'Remember to delete the older versions of the file in:');
    fprintf('\t%s\n',[rootPUBdir ' , AND ']);
    fprintf('\t%s\n',rootUPDATEdir);
end
end

