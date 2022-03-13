function returnLogTable=subjectDirCheck(subj, rootEEGdir, varargin)

%%
%   looks for regex directory and file patterns in the subject directory
%   and gives warnings if things are missing or should be deleted.
%
%   steps:
%   1. recursively searches the subject directory and makes list of all files
%   2. deletes all hidden files starting with . or ~
%   3. deletes old events.mat files
%   4. looks for matches of patterns in directory_patterns.csv
%   5. prints warnings and info based on expected vs actual pattern matches
%   6. looks through alignmentStats.txt for warnings to be fixed
%   7. looks through alignmentSummary.txt for warnings to be fixed
%   8. looks through eventEegFilePaths.txt for correct last 2 lines
%
% -output
%   printing to STDOUT and writing to finalChecklist.txt in subj docs
%   directory
%
% -optional inputs
%    writeOut (1)   %- save the checklist to subject/docs/finalChecklist.csv
%    printOut (1)   %- present the checklist on the command line
%  returnTable (0)  %- return a table structure with the output
% 'copyInvalidNotes',(1)  %- copy any notes found in "invalid" to the new list
% 'deleteInvalid' (1)     %- delete the invalid after copying 
%
% created: CZ Aug/24/2017
% edited: CZ Dec/2017
%    4/4/18 -- JW did minor tweaks and added some comments
%   4/12/18 -- JW removes requirement for finalChecklistFields.txt... field list is now dynamically creatd from patterns.xls
%    5/2/18 -- JW imports notes from invalid to the new checklist, then deletes invalid (so a clean rebuild of a subj doesn't loose the notes or require manual copy)
%   6/19/20 -- SJ replace readtable with readtableSafe

p = inputParser;
p.addParameter('writeOut', 1, @isint);      %- save the checklist to subject/docs/finalChecklist.csv
p.addParameter('printOut', 1, @isint);      %- present the checklist on the command line
p.addParameter('returnTable', 0, @isint);
p.addParameter('copyInvalidNotes',1, @isint);  %- if invalid is present, take its notes (if no note for that entry)
p.addParameter('deleteInvalid',1, @isint);     %- if invalid is present, delete it
parse(p, varargin{:});

global printOut;
printOut = p.Results.printOut;
global writeOut;
writeOut = p.Results.writeOut;


checkListFile    = fullfileEEG(rootEEGdir,subj,'docs/finalChecklist.csv');           %- JW update, move to docs
checkListInvalid = fullfileEEG(rootEEGdir,subj,'docs/finalChecklist.invalid.csv');   %- JW update, move to docs


if ~exist(strcat(rootEEGdir,'/',subj))
    error('%s does not exist or cannot be accessed.\n',strcat(rootEEGdir,'/',subj));
end

patterns_filename = which('directory_patterns.csv');
if isempty(patterns_filename)
    error('%s directory_patterns.csv file does not exist or cannot be accessed.\n');
end
patterns = readtableSafe(patterns_filename,'ReadVariableNames',1, 'Delimiter', ','); %- read table now to genearte field list for open function JW 4/12/18; SJ: replace with readtableSafe 6/19/2020
fieldsInPattern = patterns.name(strcmp(patterns.expectedMatches,'-1')==0);


%- these are additional fields tested within this function. Define here so they can be chcked in "openLogFile"  -JW addition 4/12/18
fieldsNotInPattern = {'alignmentStats_warnings',
    'alignmentSummary_task_sessions'
    'alignmentSummary_missing_events.mat'
    'alignmentSummary_missing_eeg_raw'
    'alignmentSummary_missing_eeg_sync'
    'alignmentSummary_alignment_possible'
    'alignmentSummary_unpaired_eeg'
    'eventEEGFilePaths_2ndToLastLine'
    'eventEEGFilePaths_lastLine'
    'patientInfo_containsResection'
    'elementInfo_emptyField_resected'
    'elementInfo_emptyField_ictal'
    'elementInfo_emptyField_interictal'
    };


log_table = openLogFile(checkListFile, fieldsInPattern, fieldsNotInPattern);

if exist(checkListInvalid,'file'),
    if p.Results.copyInvalidNotes,
        log_table_invalid = openLogFile(checkListInvalid, fieldsInPattern, fieldsNotInPattern);
        for iRow=1:size(log_table,1),
            rowInvalid = find(strcmp(log_table.name{iRow},log_table_invalid.name));
            %- if invalid has a matching row, and that row has a note, and the valid one DOESNT have a note... copy the note from invalid to valid
            newlogEmpty=0; invalidhasNote=0;  %- initalize
            if isempty(log_table.notes{iRow}) || strcmp(log_table.notes{iRow},'--')  newlogEmpty=1; end
            if ~isempty(rowInvalid) &&  ~isempty(log_table_invalid.notes{rowInvalid}) && ~strcmp(log_table_invalid.notes{rowInvalid},'--') invalidhasNote=1; end
            if newlogEmpty & invalidhasNote,
                log_table.notes{iRow} = log_table_invalid.notes{rowInvalid};
                if printOut, fprintf('\n copying note from invalid to new checklist: %s --> ''%s''',log_table.name{iRow}, log_table.notes{iRow}); end
            end
        end
    end
    if p.Results.deleteInvalid,
        delete(checkListInvalid);
        if printOut, fprintf('\n invalid checklist found and deleted'); end
    end
end

if p.Results.writeOut
    fprintfWrapper(sprintf('Log is being written to %s\n\n', checkListFile));
end



fprintfWrapper(sprintf('gathering files in %s/\n', subj));

fprintfWrapper(sprintf('-----------------------------------------------------------\n'));

%if p.Results.writeOut -- dlete in all cases, not just "writeOut"
fprintfWrapper(sprintf('deleting hidden files beginning with . or ~\n'));
%end

%looping BFS to get all files in patient directory
search_paths = {strcat(rootEEGdir,'/',subj)};
all_files = '';
hidden_file_count = 0;

while ~isempty(search_paths)
    
    current_search_path = search_paths{1};
    
    %trim current path to prepend to filenames
    split_path = strsplit(current_search_path, '/');
    subj_match = find(strcmpi(subj, split_path),1);
    
    trimmed_search_path = '';
    %concat all parts of file path starting after subject ID
    if ~isempty(subj_match) && ~isequal(subj_match,length(split_path))
        for k = subj_match+1:length(split_path)
            trimmed_search_path = strcat(trimmed_search_path,'/',split_path{k});
        end
    end
        
    %get search path contents
    current_files = dir(current_search_path);
    current_files = struct2table(current_files);
    
    %find hidden files starting with . or ~
    current_hidden_files = current_files( find( cellfun(@(c) isequal(c(1),'.') || isequal(c(1),'~') ,current_files.name) ) , current_files.Properties.VariableNames);
    current_hidden_files = current_hidden_files( find( ~current_hidden_files.isdir ) ,'name');
    
    hidden_files = strcat(current_search_path,'/',current_hidden_files.name);
    hidden_file_count = hidden_file_count + length(hidden_files);
    
    %if p.Results.writeOut %- delete in all cases, not just "writeOut"  -JW
    cellfun(@delete, hidden_files);
    %end
    
    %add found files to all_files
    file_paths = strcat(trimmed_search_path,'/',current_files.name);
    
    if isequal(all_files,'')
        all_files = file_paths;
    else
        all_files = [all_files;file_paths];
    end
    
    %add found directories to search_paths
    
    %non hidden files
    non_hidden_files = current_files( find( cellfun(@(c) ~isequal(c(1),'.') && ~isequal(c(1),'~') ,current_files.name) ) , current_files.Properties.VariableNames);
    current_folders = non_hidden_files( find( non_hidden_files.isdir ) ,'name');
    
    %add folders to search path array, if any exist
    if size(current_folders,1) > 0
        full_paths_folders = strcat(current_search_path, '/', current_folders.name);
        search_paths = [search_paths;full_paths_folders];
    end
    
    %remove first element from search_paths
    search_paths = search_paths(2:end);
    
end



if p.Results.writeOut
    
    fprintfWrapper(sprintf('deleted %d hidden files.\n', hidden_file_count));
    
    %delete old events.mat files
    fprintfWrapper(sprintf('deleting old events.mat files (events.mat.old OR oldEvents.mat) in task and session directories\n'));
    delete_patterns = {'^\/behavioral\/.+\/events\.mat\.old$' '^\/behavioral\/.+\/oldEvents\.mat$' '^\/behavioral\/.+\/session_\d+\/events\.mat\.old$' '^\/behavioral\/.+\/session_\d+\/oldEvents\.mat$'};
    
    old_event_count = 0;
    for i = 1:length(delete_patterns)
        
        current_pattern = delete_patterns{i};
        pattern_present = ~cellfun(@isempty, regexp(all_files, current_pattern, 'match'));
        match = find(pattern_present);
        
        old_event_count = old_event_count + length(match);
        
        if ~isempty(match)
            %fprintfWrapper(sprintf('%s\n', full_names{match});
            for j = 1:length(match)
                file_idx = match(j);
                delete(strcat(rootEEGdir,'/',subj,'/',all_files{file_idx}));
                all_files{file_idx} = '';
            end
        end
    end
    
    fprintfWrapper(sprintf('deleted %d old events.mat files.\n', old_event_count));
end

fprintfWrapper(sprintf('-----------------------------------------------------------\n'));


% open patterns that are known to be in a subjects folder
% patterns = readtable(patterns_filename,'ReadVariableNames',1, 'Delimiter', ','); %- JW does this above now
patterns.target = cellfun(@(c) regexprep(c,'\[subject\]', subj), patterns.target, 'UniformOutput', 0);
patterns.regex  = cellfun(@(c) regexprep(c,'\[subject\]', subj), patterns.regex, 'UniformOutput', 0);
patterns.messageIfPresent = cellfun(@(c) regexprep(c,'\[subject\]', subj), patterns.messageIfPresent, 'UniformOutput', 0);
patterns.present = zeros(size(patterns,1), 1);
patterns.matches = cell(size(patterns,1), 1);

for pattern_idx = 1:size(patterns,1)
    current_pattern = patterns{pattern_idx, 'regex'};
    
    if isequal(0, patterns{pattern_idx, 'caseSensitiveMatching'})
        pattern_present = ~cellfun(@isempty, regexp(all_files, current_pattern, 'match', 'ignorecase'));
    else
        pattern_present = ~cellfun(@isempty, regexp(all_files, current_pattern, 'match'));
    end
    
    match = find(pattern_present);
    patterns.present(pattern_idx) = length(match);
    patterns.matches{pattern_idx} = match;
end


%if any values in expectedMatches are not numbers, look to replace the
%value with the other row it references in directory_patterns.csv
non_digit_match = cellfun(@isempty, regexp(patterns.expectedMatches, '^-?\d+$', 'match'));
non_digit_rows = find(non_digit_match);
for j = 1:length(non_digit_rows)
    current_row = non_digit_rows(j);
    placeholder_value = patterns.expectedMatches{current_row};
    
    %find the row with matching name
    matching_row = find(strcmpi(placeholder_value, patterns.name));
    if length(matching_row) < 1
        error('Row with name corresponding to expectedMatches placeholder value for %s not in directory_patterns.csv',strjoin(patterns{current_row,'name'}));
    elseif length(matching_row) > 1
        error('More than 1 row has a corresponding name to expectedMatches placeholder value for %s not in directory_patterns.csv',strjoin(patterns{current_row,'name'}));
    end
    
    patterns.expectedMatches{current_row} = num2str(patterns.present(matching_row));
end
patterns.expectedMatches = cellfun(@str2num, patterns.expectedMatches);



%- identify the subject number
if ~isempty(strfind(subj,'NIH')); %#ok<STREMP>
    subject_ID_split = strsplit(subj,'NIH');
elseif ~isempty(strfind(subj,'TRE')); %#ok<STREMP>
    subject_ID_split = strsplit(subj,'TRE');
end
subject_number = str2num(subject_ID_split{2});



%---------------------------------------------------------------------------------------------------------
% loop over all items in the pattern table and see if they are present/absent
%---------------------------------------------------------------------------------------------------------
any_warnings = 0;
for row_index = 1:size(patterns,1)
    
    current_expected = patterns.expectedMatches(row_index);
    current_present = patterns.present(row_index);
    current_name = patterns.name{row_index};
    current_target = patterns.target{row_index};
    current_matches = patterns.matches{row_index};
    current_minimum_subj_thresh = patterns.minimumPatientNumber(row_index);
    
    
    if current_expected ~= -1%is this an optional pattern
        
        if subject_number >= current_minimum_subj_thresh | strcmp(subj(1:3),'TRE')%does this pattern apply to this subject
            
            if current_expected ~= current_present
                
                any_warnings = 1;
                
                logM = sprintf('WARNING: Expected %d %s (%s), found %d instead.', current_expected, current_name, current_target, current_present);
                fprintfWrapper(sprintf('%s\n', logM));
                log_table = updateLog(log_table, current_name, logM, 'warning');
                
                %print a message if there is one
                if length(current_matches) > 0 && ~isequal(patterns.messageIfPresent(row_index),{''})
                    
                    fprintfWrapper(sprintf('%s\n',patterns.messageIfPresent{row_index}));
                    
                end
                
                %print out the matches if printMatchesIfWarningPresent is
                %set to 1
                if length(current_matches) > 0 && patterns.printMatchesIfWarningPresent(row_index) == 1
                    
                    fprintfWrapper(sprintf('-------%s\n', all_files{current_matches}));
                    
                end
                
            else%means we found what we expected
                
                logM = sprintf('GOOD: Expected %d %s (%s), found %d.', current_expected, current_name, current_target, current_present);
                fprintfWrapper(sprintf('%s\n', logM));
                log_table = updateLog(log_table, current_name, logM, 'good');
            end
            
        else%this pattern does not apply this subject
            
            logM = sprintf('SKIP: %s only applies to subjects above %d.', current_name, current_minimum_subj_thresh);
            fprintfWrapper(sprintf('%s\n', logM));
            log_table = updateLog(log_table, current_name, logM, 'skip');
            
        end
        
    else% -1, so it's an optional match, but maybe there is still a message
        
        if length(current_matches) > 0 && ~isequal(patterns.messageIfPresent(row_index),{''})
            fprintfWrapper(sprintf('%s\n',patterns.messageIfPresent{row_index}));
        end
    end
end


%---------------------------------------------------------------------------------------------------------
% check for warnings in the alignmentStats.txt
%---------------------------------------------------------------------------------------------------------
align_stats_row = find(strcmp('alignStats', patterns.name));
align_stats_present = patterns{align_stats_row, 'present'};

stats_warning_count = 0;

if align_stats_present == 1
    
    fid_stats = fopen(fullfileEEG(rootEEGdir,subj,'behavioral/alignmentStats.txt'));
    
    tline = fgets(fid_stats);
    while ischar(tline)
        
        matches = regexp(tline,'WARNING: \d+ of \d+ events');
        if length(matches) > 0
            any_warnings = 1;
            stats_warning_count = stats_warning_count + 1;
        end
        
        tline = fgets(fid_stats);
    end
    
    fclose(fid_stats);
end

if stats_warning_count > 0
    logM =sprintf('WARNING: %d unaligned event warnings present in alignmentStats.txt', stats_warning_count);
    fprintfWrapper(sprintf('%s\n', logM));
    log_table = updateLog(log_table, 'alignmentStats_warnings', logM, 'warning');
else
    logM = sprintf('GOOD: no alignment warnings present in alignmentStats.txt');
    fprintfWrapper(sprintf('%s\n', logM));
    log_table = updateLog(log_table, 'alignmentStats_warnings', logM, 'good');
end


%---------------------------------------------------------------------------------------------------------
% check alignment summary good alignment pairings
%---------------------------------------------------------------------------------------------------------
align_summary_row = find(strcmp('alignSummary', patterns.name));
align_summary_present = patterns{align_summary_row, 'present'};

if align_summary_present == 1
    
    summary_text = fileread(fullfileEEG(rootEEGdir,subj,'behavioral/alignmentSummary.txt'));
    
    tasks_aligned = regexp(summary_text,'\d+ of \d+ task-sessions aligned', 'match');
    tasks_aligned_nums = regexp(tasks_aligned, '\d+', 'match');
    if isequal(tasks_aligned_nums{1}{1},tasks_aligned_nums{1}{2})
        
        logM = sprintf('GOOD: %s in alignmentSummary.txt', strjoin(tasks_aligned));
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table, 'alignmentSummary_task_sessions', logM, 'good');
    else
        logM = sprintf('WARNING: %s in alignmentSummary.txt. The numbers should match.', strjoin(tasks_aligned));
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table, 'alignmentSummary_task_sessions', logM, 'warning');
        any_warnings = 1;
    end
    
    align_summary_patterns = {'\d+ missing events\.mat' '\d+ missing eeg RAW data' '\d+ missing eeg sync' '\d+ alignment possible' '\d+ unpaired eeg RAW data'};
    align_summary_fields = {'alignmentSummary_missing_events.mat' 'alignmentSummary_missing_eeg_raw' 'alignmentSummary_missing_eeg_sync' 'alignmentSummary_alignment_possible' 'alignmentSummary_unpaired_eeg'};
    
    for pi = 1:length(align_summary_patterns)
        align_matches = regexp(summary_text, align_summary_patterns{pi}, 'match');
        align_matches_num = regexp(align_matches, '\d+', 'match');
        if isequal(align_matches_num{1}{1},'00')
            logM = sprintf('GOOD: %s in alignmentSummary.txt', strjoin(align_matches));
            fprintfWrapper(sprintf('%s\n', logM));
            log_table = updateLog(log_table, align_summary_fields{pi}, logM, 'good');
        else
            logM = sprintf('WARNING: %s in alignmentSummary.txt. There should be 0.', strjoin(align_matches));
            fprintfWrapper(sprintf('%s\n', logM));
            log_table = updateLog(log_table, align_summary_fields{pi}, logM, 'warning');
            any_warnings = 1;
        end
    end
    
    
end


%---------------------------------------------------------------------------------------------------------
% check eventEegFilePaths.txt
% eventEegFilePaths.txt for 2 last lines... should read "1 unique" and path set to "Volumes/Shares/FRNU/data/eeg/NIHXXX/eeg.noreref"
%---------------------------------------------------------------------------------------------------------
eeg_filepaths_row = find(strcmp('eventEegFilePaths', patterns.name));
eeg_filepaths_present = patterns{eeg_filepaths_row, 'present'};

if eeg_filepaths_present == 1
    
    fid_eeg_filepaths = fopen(fullfileEEG(rootEEGdir,subj,'behavioral/eventEegFilePaths.txt'));
    
    tline = fgets(fid_eeg_filepaths);
    second_to_last_line = '';
    last_line = '';
    
    while ischar(tline)
        second_to_last_line = last_line;
        last_line = tline;
        
        tline = fgets(fid_eeg_filepaths);
    end
    
    fclose(fid_eeg_filepaths);
    
    second_to_last_match = regexp(second_to_last_line, 'SUMMARY: 1 unique eegfile paths in events\.mat files across all experiments', 'match');
    last_match = regexp(last_line, strcat('\/Volumes\/Shares\/FRNU\/data\/eeg\/',subj,'\/eeg.noreref\/'), 'match');
    
    if isempty(second_to_last_match)
        logM = sprintf('WARNING: 2nd to last line of eventEegFilePaths.txt does not specify "1 unique" eegfile path.');
        fprintfWrapper(sprintf('%s\n', logM));
        %remove newline
        logM = sprintf('WARNING: 2nd to last line of eventEegFilePaths.txt does not specify "1 unique" eegfile path.');
        log_table = updateLog(log_table,'eventEEGFilePaths_2ndToLastLine' , logM, 'warning');
        any_warnings = 1;
    else
        logM = sprintf('GOOD: 2nd to last line of eventEegFilePaths.txt specifies "1 unique" eegfile path.');
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table,'eventEEGFilePaths_2ndToLastLine' , logM, 'good');
    end
    
    if isempty(last_match)
        logM = sprintf(strcat('WARNING: last line of eventEegFilePaths.txt does not specify subject eeg.noreref file path.'));
        fprintfWrapper(sprintf('%s\n', logM));
        %remove newline
        logM = sprintf(strcat('WARNING: last line of eventEegFilePaths.txt does not specify subject eeg.noreref file path.'));
        log_table = updateLog(log_table,'eventEEGFilePaths_lastLine' , logM, 'warning');
        any_warnings = 1;
    else
        logM = sprintf('GOOD: last line of eventEegFilePaths.txt specifies subject eeg.noreref file path.');
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table,'eventEEGFilePaths_lastLine' , logM, 'good');
    end
end


%---------------------------------------------------------------------------------------------------------
% check for elementInfo warnings -- isResected, isIctal, isInterictal should usually NOT be empty (these are added by Kareem after patient discharge)
%---------------------------------------------------------------------------------------------------------
elementInfo_exists_row = find(strcmp('element_info', patterns.name));
elementInfo_exists_present = patterns{elementInfo_exists_row, 'present'};
%keyboard
if elementInfo_exists_present == 1
    
    fieldList = {'resected', 'ictal', 'interictal'};
    
    %- parse patient info to see if resected electrodes are expected.
    patient_info = parsePatientInfo(subj, rootEEGdir);
    
    %- first confirm patient info has a resection field filled out (whether None or other)
    warnStr = 'patientInfo_containsResection';
    if isempty(patient_info.resection),
        logM = sprintf('WARNING: patient_info resection field is blank');
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table, warnStr , logM, 'warning');
        any_warnings = 1;
    else
        logM = sprintf('GOOD: patient_info resection field is not blank');
        fprintfWrapper(sprintf('%s\n', logM));
        log_table = updateLog(log_table, warnStr , logM, 'good');
    end
    
    %- now confirm that resected/ictal/interictal channels are indicated in element_info
    for iF = 1:length(fieldList),
        
        %- if patient info says "none" for resected, then dont expect any channels
        expectNoChan = 0;
        if strcmp(fieldList{iF},'resected')         & contains(lower(patient_info.resection),'none'),
            expectNoChan = 1;
        elseif strcmp(fieldList{iF},'ictal')        & patient_info.Ictal_isNone==1,
            expectNoChan = 1;
        elseif strcmp(fieldList{iF},'interictal')   & patient_info.Interictal_isNone==1,
            expectNoChan = 1;
        end
        
        warnStr   = sprintf('elementInfo_emptyField_%s', fieldList{iF});
        chanNames = getLeads(subj, rootEEGdir, 'markedAs', fieldList{iF});
        
        if isempty(chanNames) & expectNoChan==0
            logM = sprintf('WARNING: no electrodes specified as %s in elementInfo.csv', fieldList{iF});
            fprintfWrapper(sprintf('%s\n', logM));
            log_table = updateLog(log_table, warnStr , logM, 'warning');
            any_warnings = 1;
        else
            if expectNoChan==1
                logM = sprintf('GOOD: 0 %s electrodes in elementInfo is noted in patientInfo', fieldList{iF});
            else
                logM = sprintf('GOOD: 1 or more electrodes specified as %s in elementInfo.csv', fieldList{iF});
            end
            fprintfWrapper(sprintf('%s\n', logM));
            log_table = updateLog(log_table, warnStr , logM, 'good');
        end
    end
end


%replace blank notes with '--'
for i = 1:size(log_table,1)
    row = log_table{i,log_table.Properties.VariableNames};
    
    if isequal(row{4}, '')
        row{4} = '--';
        log_table{i,log_table.Properties.VariableNames} = row;
    end
    
end


%error count and reminders
fprintfWrapper(sprintf('-----------------------------------------------------------\n'));

if any_warnings == 0
    fprintfWrapper(sprintf('Everything seems in order. 0 warnings.\n'));
else
    fprintfWrapper(sprintf('At least 1 warning to fix or annotate.\n'));
end

fprintfWrapper(sprintf('REMINDER: check that clinical_summary.pdf contains no patient IDs.\n'));
fprintfWrapper(sprintf('REMINDER: check that clinical_montage.xlsx contains no patient IDs and that electrode list are all on a single sheet.\n'));

if p.Results.writeOut
    closeLogfile(log_table,checkListFile);
end

if p.Results.returnTable
    returnLogTable = log_table;
else
    returnLogTable = 'returnTable parameter set to 0, returning this string';
end

end




%-------------------------------------------------------------------------------------------------------%
function log_table=openLogFile(checkListFile, fieldsInPattern, fieldsNotInPattern)

global writeOut;

fprintfWrapper(sprintf('checking for existing finalChecklist.csv ...\n'));
fields = {fieldsInPattern{:} fieldsNotInPattern{:}};

if exist(checkListFile,'file')
    
    log_in = fopen(checkListFile);
    
    log_table = textscan(log_in, '%q%q%q%q', 'Delimiter',',', 'HeaderLines',1);
    log_table = table(log_table{1},log_table{2},log_table{3},log_table{4});
    log_table.Properties.VariableNames = {'name' 'outputMessage' 'status' 'notes'};
    
    fclose(log_in);
    
    if ~isequal(log_table(:,1), fields)
        
        fprintfWrapper(sprintf('existing log file does not contain all fields present in finalChecklist_fields.txt\n'));
        
        new_log_table = newLog(fields);
        
        for i = 1:height(new_log_table)
            current_name = find(strcmp(new_log_table.name{i}, log_table.name));
            if ~isempty(current_name)
                new_log_table.notes{i} = log_table.notes{current_name};
            end
        end
        log_table = new_log_table;
        
    end
else
    
    log_table = newLog(fields);
end


end




%-------------------------------------------------------------------------------------------------------%
function log_table=newLog(fields)

fprintfWrapper(sprintf('making log from scratch.\n'));
log_table = table(fields');


log_table = [log_table cell(length(fields), 3)];
log_table.Properties.VariableNames = {'name' 'outputMessage' 'status' 'notes'};
end




%-------------------------------------------------------------------------------------------------------%
function closeLogfile(log_table, checkListFile)

fprintfWrapper(sprintf('\n\nwriting log to %s\n',checkListFile));
%log_out = fopen(fullfileEEG(rootEEGdir,subj,'finalChecklist.csv'),'w');
%fprintfWrapper(sprintf(log_out,'"%s","%s","%s","%s"\n', log_table.Properties.VariableNames{1:4});

writetable(log_table, checkListFile, 'QuoteStrings', 1);

%making a copy on FRNU, deleting the original prevents a weird 'file not accessible error when opening with excel'
system(sprintf('cp %s %s', checkListFile,[checkListFile '.copy'] ));
system(sprintf('rm %s', checkListFile ));
system(sprintf('mv %s %s', [checkListFile '.copy'],checkListFile));

%fclose(log_out);
end




%-------------------------------------------------------------------------------------------------------%
function log_table=updateLog(log_table, name, message, status)

current_row = find(strcmpi(log_table.name, name));

if length(current_row) ~= 1
    error('attempting to update field in log_table that either does not exist or is not unique.\n');
end

log_table{current_row, 'outputMessage'} = {message};
log_table{current_row, 'status'} = {status};

end



%-------------------------------------------------------------------------------------------------------%
function fprintfWrapper(s)
global printOut;

if printOut
    fprintf(s);
end
end

