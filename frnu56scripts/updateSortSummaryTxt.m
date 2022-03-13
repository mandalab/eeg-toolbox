function [sortNoteSummary, tSortNotes, flag_sortsIncomplete] =  updateSortSummaryTxt(subj,sortedByPath,sessID,sortTxtFilesCln)

% subj = 'NIH069';
% sortedByPath = '/Volumes/Seagate Backup Plus Drive/local56/56PROC/micro_behavioral/micro_pristine/NIH069/sorts_manual/190117_1336/reref(stimMask5)_sortedByDY';
% sessID = '190117_1336';

% sortRootName = fullfile(sortedByPath, sprintf('sortNotes(%s)_sortedBy',sessID));
% dSortNotes = dir([sortRootName '*']);
% tSortNotes = {};
% Anything you change in here must be changed in extractSpikeInfo_v3b.m!!
% Eventually use this function in there.

fprintf('\n%s','Updating sort summary text file ... ');

sortRootName = fullfile(sortedByPath, sprintf('sortNotes(%s)_sortedBy',sessID));
dSortNotes = dir([sortRootName '*']);
tSortNotes = {};

if length(dSortNotes)>=1,
    for iNote=1:length(dSortNotes),
        tSortNotes{iNote,1} = dSortNotes(iNote).name;
        tSortNotes{iNote,2} = readtable( fullfile(sortedByPath,dSortNotes(iNote).name));
    end   
    if length(dSortNotes)>1,
        fprintf('\n kinda weird to have two or more sort notes here... should just be one.  What is going on?');
        keyboard;
    end
else
    fprintf('\n Uh Oh, no sort notes found, expected %s', sortRootName);
    keyboard
    tSortNotes{1,1} = 'no sort notes found';
    tSortNotes{1,2} = table; %- initailize to empty table
end


%- crate a "sortNoteSummary"... how many A's Bs Cs, etc... then see how many of those have sort files.  Want to detect if sorts are missing

%- first check for and delete the old summaries
old_sort_summary = dir(fullfile(sortedByPath,'sorts_*.txt'));
if ~isempty(old_sort_summary), for iOld = 1:length(old_sort_summary), delete(fullfile(sortedByPath,old_sort_summary(iOld).name)); end; end

if istable( tSortNotes{1,2} ),

    %- assume we should be looking at the first sort note if more than one was found
    iNote = 1;

    sortNoteSummary = sprintf('Summary of grades found in %s/%s, %s',subj,tSortNotes{iNote,1},datetime);
    tSorts = tSortNotes{iNote,2};
    cGrade = find(contains(tSorts.Properties.VariableNames,'Grade'));
    grades = tSorts{:,cGrade};
    
    if ~iscell(grades), %- grades should always be a cell array, but if blank might be read as "nan".  convert to cell if so
        if sum(isnan(grades))==length(grades),
            grades = {};  for ii=1:height(tSorts), grades{ii,1} = 'nan'; end
        else
            fprintf('\n wasnt cell, but values are not all nan... did somebody type in numerical grades?');
            keyboard;
        end 
    end

    %- define a grade as the highest grade of multiple measurements
    grades = upper(grades);
    grades(cellfun('isempty',grades)) = {'Z'};                

    if size(grades,2)>1 % Normally should be
        grades = sort(grades')'; %- force all to upper case so unique sorting doesn't separate them below
        grades = grades(:,1); % This is the higher of the 2 column grades
    else % In older cases, like NIH069, there was only 1 grade
    % Do nothing, because you already have your 1 column!
    end


    % Check for blanks - strncmp below only checks that the
    % first character in the string matches (so ignores -/+)
    blanks = find( ~strncmp(grades,'A',1) & ~strncmp(grades,'B',1) & ~strncmp(grades,'C',1) & ~strncmp(grades,'D',1) & ~strncmp(grades,'F',1)) ;
    if length(blanks)>0,
        for ii=1:length(blanks),
            grades{blanks(ii),1} = 'F(was blank)'; %
        end
    end
    tSorts.maxGrade = grades; %- create a new column
    %%%%%%%

    %- now create a column in tSorts that indicates which channels have a sort.txt
    hasSortTxt = zeros(height(tSorts),1);
    for iChan=1:height(tSorts)
        hasSortTxt(iChan) = sum(strcmp(sortTxtFilesCln,sprintf('%s.txt',tSorts.SortChanName{iChan}))); %- use the "Clean" version of sortTxtFiles, which is scrubbed of "(grabbedXXYYZZ)" 
    end
    if any(hasSortTxt>1), fprintf('\n how did this happen? should be zero or one.'); keyboard; end
    tSorts.hasSortTxt = hasSortTxt; %- creates new column

    %- put the modified version of the table (two new columns) back into the tSortNotes cell array so it is appended to spikeInfo
    tSortNotes{iNote,2} = tSorts; 

    %- grade list
    flag_sortsIncomplete = 0;
    gradeList = unique({tSorts.maxGrade{:},'A','B','C','D','F'}); %- add A,B,C,D,F here in case none were reported, but still use unique incase other things were written there
    sortsComplete = NaN(length(gradeList),1); 
    sortsRedo = NaN(length(gradeList),1); 
    for iG=1:length(gradeList),
        thisGrade = gradeList{iG};

        %- search for sort files... were all the "a's" sorted? and b's, etc???
        iChan = find(strcmp(tSorts.maxGrade,thisGrade));
        gradeCount = length(iChan);
        sortCount  = sum(tSorts.hasSortTxt(iChan));

        %- check whether each grade is sorted or not.  If A or B, then flag as incomplete sorts overall
        if   sortCount==gradeCount,  
            if sortCount>0,
                summaryStr = '(completely sorted)';
                sortsComplete(iG) = 3;
            else
                summaryStr = '(none to sort)';
                sortsComplete(iG) = 0;
            end
            sortsRedo(iG) = 0;

        else
            if contains(thisGrade(1),'A') || contains(thisGrade(1),'B'),  
                flag_sortsIncomplete = 1; 
                sortsRedo(iG) = 1;
                summaryStr = '<<< incomplete "A" or "B" sort. Need to add those!!!!';
            else
                summaryStr = '<<< incomplete sorts?!?';
            end
            if sortCount == 0
                sortsComplete(iG) = 1;
            else
                sortsComplete(iG) = 2;
            end
        end
        sortNoteSummary = sprintf('%s\n   Max Grade from nPlay %s: %d identified, %d sorted %s ',sortNoteSummary, thisGrade, gradeCount,sortCount,summaryStr);
    end

    for gg = 1:numel(sortsRedo) %sortsComplete)
        if isnan(sortsRedo(gg))
            comp_below = sortsComplete(gg+1:end);
            if ismember(3,comp_below) || ismember(2,comp_below)
                sortsRedo(gg) = 1;
            else
                sortsRedo(gg) = 0;
            end
        end
    end
    if length(blanks)>0
        if length(blanks)==length(grades),
            sortNoteSummary = sprintf('%s\n   PROBLEM: no grades marked in sortNotes, so cant tell if all units were sorted or not',sortNoteSummary);
        else
            sortNoteSummary = sprintf('%s\n   PROBLEM: %s grades blank in sortNotes, so cant tell if all units were sorted or not',sortNoteSummary,num2str(length(blanks)));
        end
        flag_sortsIncomplete = 1;
    end    
    redoChanList = {};
    if sum(sortsRedo)>0
        sortNoteSummary = sprintf('%s\n   PROBLEM: Sorts for higher grades not completed before lower grades (and/or A/B not completely sorted).',sortNoteSummary);
        flag_sortsIncomplete = 1;
        % Suggest which channels to sort
        reSortSummary = sprintf('\n%s\n','Incomplete channels to sort: ');
        redoChanList = {};
        for gg2 = 1:numel(sortsRedo)
            if sortsRedo(gg2) == 1
                thisGrade2 = gradeList{gg2};
                reSortSummary = sprintf('%s\n%s\n',reSortSummary,['Grade: ' thisGrade2]);
                redoChanList = tSorts.SortChanName(strcmp(tSorts.maxGrade,thisGrade2) & tSorts.hasSortTxt == 0);
                %sprintf('\t%s\n',redoChanList{:});
                %reSortSummary = sprintf('%s%s\n',reSortSummary,redoChanList{:});
                reSortSummary = sprintf('%s%s',reSortSummary,sprintf('\t%s\n',redoChanList{:}));
                %redoChanList = [redoChanList; tSorts.SortChanName(strcmp(tSorts.maxGrade,thisGrade2) & tSorts.hasSortTxt == 0)];
            end
        end
        sortNoteSummary = sprintf('%s\n%s\n',sortNoteSummary,reSortSummary);
    end
                        %- output a text file for easy determination of whether sorts are complete or not
    if flag_sortsIncomplete, strComplete = 'INCOMPLETE'; else strComplete = '(complete)'; end
    fid_sortSum = fopen(fullfile(sortedByPath,sprintf('sorts_%s.txt',strComplete)),'w+');
    fprintf(fid_sortSum,'%s',sortNoteSummary);
    fclose(fid_sortSum);

else
    sortNoteSummary = sprintf('%s\n table not read?',sortNoteSummary);
    flag_sortsIncomplete = 1; %- unknown if complete or not without the table... assume incomplete
    keyboard;

    fid_sortSum = fopen(fullfile(sortedByPath,'sorts_MISSING_SORT_NOTES.txt'),'w+');
    fclose(fid_sortSum);
end

fprintf('%s\n','Done. Sort note summary: ');

disp(sortNoteSummary);
end

