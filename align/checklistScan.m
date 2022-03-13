function checklistScan(rootPath, varargin)
%
% checklistScan
% rootPath - path where NIHXXX subject folders are present
% writeScanSteps - binary flag (default false) for writing list of each subject's changes (e.g., 47 items the same, 2 good -> warning, 3 warning --> good)
% writeWarningSummary - binary flag (default TRUE) for writing out master warning list
% verifyChecklists - binary flag (default FALSE) that will trigger a re-run of subjectDirCheck on each subject
% deleteInvalid - default TRUE.  if TRUE and verifyingChecklist, delete any finalChecklist.invalid.csv and dont create new ones (notes will be copied from invalid to new checklist before deleting)
% debug - binary flag (default FALSE) all written files are output to debug directory instead of rootPath
%
% this function loops over the subjects in rootPath and creates a master warning list
%
% 10/1 CZ - created
% 12/13 CZ - heavy edits, no resemblance to previous version
%
% 4/4/18 -- JW tweaked to shape published summary doc within eeg
%   4/12/18 -- JW removes requirement for finalChecklistFields.txt... field list is now dynamically creatd from patterns.xls
% 12/16/2019 -- SJ change subjBuild.xls to .xlsx
% 2/12/2021 --SJ: exclude NIH085

p = inputParser;
p.addParameter('writeScanSteps', 0, @isint);  %- JW added flag to keep/cut the subject-level update of changes.
p.addParameter('writeWarningSummary', 1, @isint);
p.addParameter('verifyChecklists',0, @isint);
p.addParameter('deleteInvalid',1, @isint);
p.addParameter('debug',0, @isint);
parse(p, varargin{:});

global dbg;
global dbg_dir;
dbg = p.Results.debug;
dbg_dir = '~/checklistScan_debug/';


fprintf('\n\n');

if dbg
    
    stars = '***********************************************************************************';
    msg1 = 'Running in debug mode';
    msg2 = sprintf('all written files are saved to dbg_dir (%s)', dbg_dir);
    
    fprintf('%s\n', stars);
    fprintf('%s\n', stars);
    fprintf('%s\n', msg1);
    fprintf('%s\n', msg2);
    fprintf('%s\n', stars);
    fprintf('%s\n', stars);
    
    if ~exist(dbg_dir, 'dir')
        mkdir(dbg_dir);
    end
    
end

global dataPath;
dataPath = rootPath;

global outFID;

if dbg
    if p.Results.writeWarningSummary
        outFID = fopen([dbg_dir '/subjProcessingWarnings.txt'], 'w');
    end
    logFileStr = [dbg_dir '/_checklistScan_verify.log'];
else
    if p.Results.writeWarningSummary
        outFID = fopen([rootPath '/subjProcessingWarnings.txt'], 'w'); %- JW tweaked this name
    end
    logFileStr = [rootPath '/_checklistScan_verify.log'];
end
logFID = fopen(logFileStr, 'a');


%- generate list of fields using finalChecklist
global fields;
%fields_file = which('finalChecklist_fields.txt'); %- old way required a separate text file
%fields_readin = fileread(fields_file);
%fields = strsplit(fields_readin,'\n');

patterns_filename = which('directory_patterns.csv');
if isempty(patterns_filename)
    error('%s directory_patterns.csv file does not exist or cannot be accessed.\n');
end
patterns = readtable(patterns_filename,'ReadVariableNames',1, 'Delimiter', ','); %- read table now to genearte field list for open function JW 4/12/18
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

fields = {fieldsInPattern{:} fieldsNotInPattern{:}};


if ~exist(rootPath,'dir'),
    fprintf('\n error, rootPath not found: %s',rootPath);
    keyboard
    return;
end
dir_res = dir(rootPath);
tableRows = [];
for i = 1:length(dir_res)
    if ( length(dir_res(i).name) > 3 ) && (isequal( dir_res(i).name(1:3), 'NIH' ) | isequal( dir_res(i).name(1:3), 'TRE' )) && (~strcmpi(dir_res(i).name,'NIH085')) %- JW added TRE, SJ added NIH085
        tableRows = [tableRows ; dir_res(i).name];
    end
end


%%% if verifyChecklists flag is passed!!
if p.Results.verifyChecklists
    
    existingTables = gatherExistingWarnings(rootPath, tableRows, 0, 0, 0);
    
    statusChanges = table(fields');
    statusChanges = [statusChanges cell(length(fields), size(tableRows,1))];
    statusChanges.Properties.VariableNames = [cellstr('field');cellstr(tableRows)]';
    
    good2warning_count = zeros(1,size(tableRows,1));
    warning2good_count = zeros(1,size(tableRows,1));
    absent2warning_count = zeros(1,size(tableRows,1));
    absent2good_count = zeros(1,size(tableRows,1));
    same_count = zeros(1,size(tableRows,1));
    
    mismatchTables = [];
    
    fprintf('\n\nverifying NIH subjects in this directory... hold tight\n\n')
    
    for z = 1:size(tableRows,1)
        
        s = tableRows(z,:);
        
        outstr = '\n...reading existing and recomputing finalChecklist.csv for %s';
        if z == 1
            fprintf(outstr, s);
        else
            %- backspace trick is slick when it works, but warnings and messages for some subjects screw it up... so back to old school
            %for v = 1:length(outstr)+4
            %    fprintf('\b');
            %end
            fprintf(outstr, s);
        end
        
        
        if ~exist([dataPath '/' s], 'dir')
            continue;
        end
        
        

        
        %create new checklist from scratch for each subject
        %if deleteInvalid, delete any finalChecklist.invalid.csv WITHIN subjDirCheck, because then the notes can be copied before delete
        newLogTable = subjectDirCheck(s, dataPath, 'writeOut', 0, 'returnTable', 1, 'printOut', 0, 'deleteInvalid', p.Results.deleteInvalid);
        
        existingLogTable = existingTables{1,z};
        
        statusChanges{:,s} = compareTables(newLogTable, existingLogTable);
        
        %in this case invalidate and prompt in the next loop
        if any(cellfun(@(x) isequal(x,'good2warning'),  statusChanges{:,s} ) | cellfun(@(x) isequal(x,'absent2warning'),  statusChanges{:,s} ))
            
            if p.Results.deleteInvalid ~= 1
                invalidateChecklist(s);
            end
            
            invalid_field_row = find(cellfun(@(x) isequal(x, 'invalid_finalChecklist'), newLogTable{:,'name'}));
            newLogTable{invalid_field_row, 'outputMessage'} = {'WARNING: Expected 0 invalid_finalChecklist (docs/finalChecklist.invalid.csv), found 1 instead.'};
            newLogTable{invalid_field_row, 'status'} = {'warning'};
            %recompute with new warning
            statusChanges{:,s} = compareTables(newLogTable, existingLogTable);
            
            writeChecklist(newLogTable, s);
            
            mismatchTables = [mismatchTables;s];
            
            %otherwise just overwrite the existing file, things could only have gotten better
        else
            writeChecklist(newLogTable, s);
        end
        
        good2warning_count(1,z) = sum(cellfun(@(x) isequal(x,'good2warning'),  statusChanges{:,s} ));
        warning2good_count(1,z) = sum(cellfun(@(x) isequal(x,'warning2good'),  statusChanges{:,s} ));
        absent2warning_count(1,z) = sum(cellfun(@(x) isequal(x,'absent2warning'),  statusChanges{:,s} ));
        absent2good_count(1,z) = sum(cellfun(@(x) isequal(x,'absent2good'),  statusChanges{:,s} ));
        same_count(1,z) = sum(cellfun(@(x) isequal(x,'same'),  statusChanges{:,s} ));
        
    end
    
    fprintf('\n\nupdated subject finalChecklists.csv files -- %s\n', datestr(datetime('now','Format','yyyy-MM-dd HH:mm')));
    fprintf(logFID, '\n\n\nupdated subject finalChecklists.csv files using checklistScan.m -- %s\n', datestr(datetime('now','Format','yyyy-MM-dd HH:mm')));
    
    star_str = '******';
    star_empty = '';
    
    for h = 1:size(tableRows,1)
        fprintf('%s updated finalChecklist contains status changes:\n\tsame status %d\n\twarning -> good %d\n\tgood -> warning %d  %s\n\tabsent -> good %d\n\tabsent -> warning %d  %s\n\n',...
            tableRows(h,:), same_count(1,h), warning2good_count(1,h), good2warning_count(1,h), ifelse(good2warning_count(1,h) > 0, star_str, star_empty), absent2good_count(1,h), absent2warning_count(1,h), ifelse(absent2warning_count(1,h) > 0, star_str, star_empty));
        fprintf(logFID,'%s updated finalChecklist contains status changes:\n\tsame status %d\n\twarning -> good %d\n\tgood -> warning %d  %s\n\tabsent -> good %d\n\tabsent -> warning %d  %s\n\n',...
            tableRows(h,:), same_count(1,h), warning2good_count(1,h), good2warning_count(1,h), ifelse(good2warning_count(1,h) > 0, star_str, star_empty), absent2good_count(1,h), absent2warning_count(1,h), ifelse(absent2warning_count(1,h) > 0, star_str, star_empty));
    end
    
    fprintf('\n\n\n');
    
    if ~isempty(mismatchTables)
        
        done_looking = 0;
        
        while ~done_looking
            
            fprintf('\n\nSubjects with newly noted warnings: [%s]\n', strjoin(cellstr(mismatchTables)));
            fprintf('You may want to inspect the checklist files for those subjects to see what has changed\n');
            fprintf('Type ''NIHXXX'' to open finalChecklists for that subject\nOR type ''quit'' to continue to the master warning list\n');
            inp = input('.... ','s');
            
            if isequal(inp, 'quit')
                done_looking = 1;
                
            else
                
                subj_search_idx = 'NIHXXX';
                for i = 1:size(tableRows,1)
                    if isequal(tableRows(i,:), inp)
                        subj_search_idx = tableRows(i,:);
                    end
                end
                
                if isequal(subj_search_idx, 'NIHXXX')
                    
                    fprintf('%s is not a valid subject in %s\n', inp, rootPath);
                    continue;
                    
                else
                    
                    isMismatch = 0;
                    for i = 1:size(mismatchTables,1)
                        if isequal(subj_search_idx, mismatchTables(i,:))
                            isMismatch = 1;
                        end
                    end
                    
                    fprintf('opening finalChecklist/s for %s\n',subj_search_idx);
                    
                    if dbg
                        system(sprintf('open %s', [dbg_dir '/' subj_search_idx '_finalChecklist.csv']));
                        if isMismatch
                            system(sprintf('open %s', [dbg_dir '/' subj_search_idx '_finalChecklist.invalid.csv']));
                        end
                    else
                        system(sprintf('open %s/%s/docs/finalChecklist.csv', rootPath, subj_search_idx));
                        if isMismatch
                            system(sprintf('open %s/%s/docs/finalChecklist.invalid.csv', rootPath, subj_search_idx));
                        end
                    end
                end
            end
        end%end while
        
    end
end


writeXLSsummary = 1;
existingTables = gatherExistingWarnings(rootPath, tableRows, p.Results.writeWarningSummary, 1, writeXLSsummary);



if p.Results.writeWarningSummary
    fclose(outFID);
end

fclose(logFID);


%- delete the "scanSteps" file
if p.Results.writeScanSteps == 0,
    delete(logFileStr);
end

fprintf('\n\n\n');

end




%-------------------------------------------------------------------------------------------------------%
function changes = compareTables(newLogTable, existingLogTable)

changes = cell(1,height(newLogTable));

for i = 1:height(newLogTable)
    
    newField = newLogTable{i,'name'};
    newField = newField{1,1};
    
    newFieldStatus = newLogTable{i,'status'};
    newFieldStatus = newFieldStatus{1,1};
    
    existingFieldStatus = 'good';
    
    if isempty(existingLogTable)
        existingFieldStatus = 'warning';
    else
        
        existingFieldMatch = find(cellfun( @(x) isequal(x, newField) , existingLogTable{:,'name'}));
        if isempty(existingFieldMatch)
            existingFieldStatus = 'absent';
        else
            existingFieldStatus = existingLogTable{existingFieldMatch, 'status'};
            existingFieldStatus = existingFieldStatus{1,1};
        end
        
    end
    
    stat_switch = 'same';
    if ~isequal(newFieldStatus,existingFieldStatus)
        if isequal(newFieldStatus, 'good') && isequal(existingFieldStatus, 'warning')
            stat_switch = 'warning2good';
        elseif isequal(newFieldStatus, 'warning') && isequal(existingFieldStatus, 'good')
            stat_switch = 'good2warning';
        elseif isequal(newFieldStatus, 'good') && isequal(existingFieldStatus, 'absent')
            stat_switch = 'absent2good';
        elseif isequal(newFieldStatus, 'warning') && isequal(existingFieldStatus, 'absent')
            stat_switch = 'absent2warning';
        end
    end
    
    changes{1,i} = stat_switch;
    
end

changes = changes';
end




%-------------------------------------------------------------------------------------------------------%
function saveToLog(d, type)

global outFID;

if type == 'table'
    
    extra_col_space = 4;
    max_col_space = zeros(1,length(d.Properties.VariableNames));
    %find the max length for each column
    
    %write the header
    for i = 1:length(d.Properties.VariableNames)
        max_col_space(1,i) = max(cellfun(@length, d{:,i}));
        
        if length(d.Properties.VariableNames{i}) > max_col_space(1,i)
            max_col_space(1,i) = length(d.Properties.VariableNames{i});
        end
        
        fprintf(outFID, '%s', d.Properties.VariableNames{i});
        
        for j = 1:( max_col_space(1,i) - length(d.Properties.VariableNames{i}) + extra_col_space )
            fprintf(outFID, ' ');
        end
    end
    
    fprintf(outFID, '\n');
    
    for g = 1:length(max_col_space)
        for h = 1:max_col_space(1,g)
            fprintf(outFID, '-');
        end
        for h = 1:extra_col_space
            fprintf(outFID, ' ');
        end
    end
    
    fprintf(outFID, '\n');
    
    %write the rest of the table
    for i = 1:height(d)
        for j = 1:length(d.Properties.VariableNames)
            
            max_space = max_col_space(1,j);
            val = d{i,j};
            
            fprintf(outFID, '%s', val{1} );
            
            for k = 1:( max_space - length(val{1}) + extra_col_space)
                fprintf(outFID, ' ');
            end
            
        end
        fprintf(outFID, '\n');
    end

    
elseif type == 'text'
    fprintf(outFID, '%s', d);
end

end




%-------------------------------------------------------------------------------------------------------%
function invalidateChecklist(s)

global dataPath;
global dbg;
global dbg_dir;

thisFileSrc = fullfileEEG(dataPath, s,'docs/finalChecklist.csv');

if dbg == 0
    
    thisFileDst = fullfileEEG(dataPath, s,'docs/finalChecklist.invalid.csv');
    if exist(thisFileSrc,'file')
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(thisFileSrc,thisFileDst,'f');
        whtSpc   = ''; if length(thisFileSrc)<80, whtSpc(80-length(thisFileSrc))=' '; end
        if SUCCESS, fprintf('\n %s %s --> renamed to docs/finalChecklist.invalid.csv', thisFileSrc, whtSpc);
        else        fprintf('\n %s %s --> error when moving', thisFileSrc, whtSpc); keyboard; end;
    end
    
else
    
    thisFileDst = [dbg_dir '/' s '_finalChecklist.invalid.csv'];
    copyfile(thisFileSrc, thisFileDst);
    
end

end




%-------------------------------------------------------------------------------------------------------%
function writeChecklist(log_table, s)

global dataPath;
global dbg;
global dbg_dir;

if dbg == 0
    
    checkListFile = [dataPath '/' s '/docs/finalChecklist.csv'];
    writetable(log_table, checkListFile, 'QuoteStrings', 1);
    
    %making a copy on FRNU, deleting the original prevents a weird 'file
    %not accessible error when opening with excel'
    system(sprintf('cp %s %s', checkListFile,[checkListFile '.copy'] ));
    system(sprintf('rm %s', checkListFile ));
    system(sprintf('mv %s %s', [checkListFile '.copy'],checkListFile));
    
else
    
    checkListFile = [dbg_dir '/' s '_finalChecklist.csv'];
    writetable(log_table, checkListFile, 'QuoteStrings', 1);
    
end

end




%-------------------------------------------------------------------------------------------------------%
function existingTables = gatherExistingWarnings(rootPath, tableRows, write_flag, disp_flag, writeXLS_flag)

global outFID;
global fields;
global dbg;
global dbg_dir;

existingWarningTable = cell2table({'x' 'x' 'x' 'x'});
existingWarningTable.Properties.VariableNames = {'subject' 'warning' 'notes' 'existingExpired'};

existingTables = cell(1, size(tableRows,1));

if disp_flag == 1
    fprintf('\n\n-- looking at existing finalChecklist.csv files ... hold tight\n\n\n');
end

count = 1;
expired_present = 0;

for h = 1:size(tableRows,1)
    
    s = char(tableRows(h,:));
    
    if ~exist([rootPath '/' s], 'dir')
        existingWarningTable = [existingWarningTable ; {s '-- not a valid subject --' '-- not a valid subject --' '-- not a valid subject --'}];
        count = count + 1;
        continue;
    end
    
    if dbg == 1
        checkListFile = [dbg_dir, '/', s, '_finalChecklist.csv'];
        if ~exist(checkListFile, 'file')
            checkListFile = strcat(rootPath, '/',s, '/docs/finalChecklist.csv');
        end
    else
        checkListFile = strcat(rootPath, '/',s, '/docs/finalChecklist.csv');
    end
    
    
    if exist(checkListFile,'file')
        
        log_in = fopen(checkListFile);
        
        existingLogTable = textscan(log_in, '%q%q%q%q', 'Delimiter',',', 'HeaderLines',1);
        existingLogTable = table(existingLogTable{1},existingLogTable{2},existingLogTable{3},existingLogTable{4});
        existingLogTable.Properties.VariableNames = {'name' 'outputMessage' 'status' 'notes'};
        
        existingTables{1,count} = existingLogTable;
        
        fclose(log_in);

        if ~isequal(existingLogTable.name, fields')
            expired_present = 1;
            expirationStatus = 'existing finalChecklist.csv does not contain all the correct fields/rows, rerun subjectDirCheck';
        else
            expirationStatus = 'all fields/rows are present';
        end
        
        warningsExistingTable = strcmpi('warning', existingLogTable.status);
        
        if size(warningsExistingTable,1) ~= 0
            for i = 1:size(warningsExistingTable,1)
                if warningsExistingTable(i,1) == 1
                    existingWarningTable = [existingWarningTable ; {s  existingLogTable{i,2}  existingLogTable{i,4} expirationStatus}];
                end
            end
        end
        
    else
        existingWarningTable = [existingWarningTable ; {s '-- no finalChecklist.csv found --' '-- no finalChecklist.csv found --' '-- no finalChecklist.csv found --'}];
        existingTables{1,count} = table();
    end
    
    count = count + 1;
    
end

existingWarningTable(1,:) = [];
missing_notes = sum(cellfun(@(x) isequal(x, '--'), existingWarningTable{:,'notes'}));

invalid_warning = 'WARNING: Expected 0 invalid_finalChecklist (docs/finalChecklist.invalid.csv), found 1 instead.';
invalid_warning_rows = cellfun(@(x) isequal(x, invalid_warning), existingWarningTable{:,'warning'});
invalid_warning_subjs = unique(existingWarningTable{invalid_warning_rows, 'subject'});


%- write a summary of the warnings to a text file
if write_flag ==1
    
    fprintf(outFID, 'Warnings from existing finalChecklist.csv files --- %s\n\n',datestr(now));   %- JW added a datestring
    saveToLog(existingWarningTable, 'table');
    fprintf(outFID,'\n\n');
    
    if expired_present == 1
        fprintf(outFID, '     !!!!!!!!!! At least 1 subject has a finalChecklist that is out of date,\n     !!!!!!!!!! re-run subjectDirCheck on those subject or run checklistScan(''%s'', ''verifyChecklists'', 1)\n', rootPath);
    end
    
    fprintf(outFID, '\n     Subject folders with invalidated finalChecklist files: [ %s ]\n\n', strjoin(invalid_warning_subjs));
    fprintf(outFID, '     Use the invalidated versions \n\t-> to indentify newly noted errors \n\t-> delete the invalid versions\n\t-> rerun subjectDirCheck or checklistScan(''%s'', ''verifyChecklists'', 1)\n', rootPath);

    fprintf(outFID,'\n\n     Of %d warnings, %d have empty notes (''--'') that need to be filled\n', height(existingWarningTable), missing_notes);
    if missing_notes>0,
        fprintf(outFID,'\n------------------------------------------------------------------------------------------------------------------------------------\n');
        fprintf(outFID,'------------------------------------------ JUST THE WARNINGS WITH EMPTY NOTES: -----------------------------------------------------\n\n');
        iMissingNotes = find(strcmp(existingWarningTable.notes,'--'));
        saveToLog(existingWarningTable(iMissingNotes,:), 'table');
        fprintf(outFID,'\n\n');
    end
    
end


%- display the warning list on the command line
if disp_flag ==1
    
    disp(existingWarningTable);
    
    fprintf('\n\n');
    
    if expired_present == 1
        fprintf('     !!!!!!!!!! At least 1 subject has a finalChecklist that is out of date,\n     !!!!!!!!!! re-run subjectDirCheck on those subject or run checklistScan(''%s'', ''verifyChecklists'', 1)\n', rootPath);
    end
    
    fprintf('\n     Of %d warnings, %d have empty notes (''--'') that need to be filled\n', height(existingWarningTable), missing_notes);
    
    fprintf('\n     Subject folders with invalidated finalChecklist files: [ %s ]\n\n', strjoin(invalid_warning_subjs));
    fprintf('     Use the invalidated versions \n\t-> to indentify newly noted errors \n\t-> delete the invalid versions\n\t-> rerun subjectDirCheck or checklistScan(''%s'', ''verifyChecklists'', 1)\n', rootPath);
    
    if dbg
        fprintf('\n     Debugging: Output also written to %s\n', [dbg_dir '/_master_warning_list.txt']);
    else
        fprintf('\n     Output also written to %s\n', [rootPath '/_master_warning_list.txt']);
    end
end



%- write a summary of all good/warning to a single xls spreadsheet
if writeXLS_flag == 1,
    
    %-
    nameXLfile = fullfileEEG(rootPath,sprintf('subjBuild_%s.xlsx',datestr(now,'yymmdd')));
    XLsheet    = 'SUBJ_DIR_CHECK';
    
    
    %- generate the data:  rows=errors x subjects=columns
    numSubj  = length(existingTables);
    numCheck = length(existingTables{1}.name);
    
    xlOUT = {'check' '' '' existingTables{1}.name{:}}'; %- create the rows
    for iSubj = 1:numSubj,
        xlOUT{iSubj+3,1} = tableRows(iSubj,:);
        for iRow = 1:numCheck,
            xlOUT{      1,iRow+3} = existingTables{1}.name{iRow};
            thisStatus = existingTables{iSubj}.status{iRow};
            thisNote   = existingTables{iSubj}.notes{iRow};
            if strcmp(thisStatus,'warning') & ~isempty(thisNote) & ~strcmp(thisNote,'--'),
                thisStatus = 'noted warning';
            end
            xlOUT{iSubj+3,iRow+3} = thisStatus;
        end
    end
    %- add summaries per row and column
    warnOrError = strcmp(xlOUT,'warning');
    xlOUT{        2,1} = 'warnPerCheck';
    xlOUT{numSubj+5,1} = 'warnPerCheck';
    for iRow = 1:numCheck,
        xlOUT{        2,iRow+3} = sum(warnOrError(:,iRow+3)); %- show sum at top and bottom of data
        xlOUT{numSubj+5,iRow+3} = sum(warnOrError(:,iRow+3));
    end
    xlOUT{1,2} = 'warnPerSubj';
    xlOUT{1,numCheck+5} = 'warnPerSubj';
    for iSubj = 1:numSubj,
        xlOUT{iSubj+3,2}          = sum(warnOrError(iSubj+3,:));
        xlOUT{iSubj+3,numCheck+5} = sum(warnOrError(iSubj+3,:));
    end
  

    %- always call "createSubjTaskList" first to create a fresh xls file.  If existing sheet found throw and error
    if exist(nameXLfile,'file'),
        [STATUS,SHEETS] = xlsfinfo(nameXLfile);
        if contains(SHEETS,XLsheet),
            fprintf('\n WARNING: summary excel file %s already contains SHEET %s; \n  best to delete xls file and create fresh',nameXLfile,XLsheet);
            keyboard;
        end
        %delete(nameXLfile);
    end
    
    %- write the data
    warning('OFF', 'xlwrite:AddSheet')
    status = xlwrite(nameXLfile, xlOUT, XLsheet, 'B2');
    if status==0,
        fprintf('\n ERROR WRITING');
        keyboard;
    end
    status = xlwrite(nameXLfile, xlOUT', [XLsheet '2'], 'B2');
    
    %- procesing date
    dateTxt = {sprintf('Created %s',datestr(now))};
    status = xlwrite(nameXLfile, dateTxt, XLsheet, 'A1:A1');
    
    fprintf('\n %s created with %d subjects and %d tasks\n',nameXLfile, length(existingTables), length(existingTables{iSubj}.name));
    
end



end %- end main



%-------------------------------------------------------------------------------------------------------%
function res = ifelse(cond, out1, out2)

if all(cond)
    res = out1;
else
    res = out2;
end
end
