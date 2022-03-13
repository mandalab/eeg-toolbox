%
%
%  Loop over previously sorted .txt files and rename them according to the new name scheme
%
%  1st step:  organize all sorts into session folders.  Should have same organization/hierarchy as if these were created by getSortFiles
%
%  rootSessDir/ session_folder / sorted / XXreref_sortedBYY / sort_txt / sortTextFiles.txt
%         (e.g. /190117_1336/sorted/noreref_sortedByDY/sort_txt/utah_AMTG02.txt)
%
%
%  2nd step: create a renameOldSorts.csv.   
%         Start with a new copy of micro_RenamedChanInfo.csv for that subject
%             cut out the first two rows, and the columns listing the sessions 
%             add a new column: "SortNameOld".   This should have the list of channel names that were originally split out. 
%


%-
rootSessDir   = '/Volumes/JW24TB/data24TB/localFRNU56/NIH069/first passes/copy of sorts OG names rename'; 
renameKeyPath = fullfile(rootSessDir,'renameOldSorts_NIH069.csv');


%- set these parameters
OVERWRITE_DEST_FILE = 0;

%- definitions of subfoldernames within sorts
SRC_SORT_TXT_DIR = 'sort_txt_oldSortNames';  %- SOURCE: put the old sort files (that need to be renamed) in here.  This is 
DST_SORT_TXT_DIR = 'sort_txt';               %- DESTINATION: this is where freshly generated sorts should go... after renaming the old ones will be in here


%- this is it.  The key for going from old sort name to new sort name.  Confirm that none of the SortChanNames match any SortNameOlds (very unlikely)
renameKey = readtable(renameKeyPath);

if ~isempty(intersect(renameKey.SortChanName, renameKey.SortNameOld)),
   fprintf('\n No good... rename key uses one or more of the same name in "old" and "new"... this will cause a naming conflict');
   keyboard;
   return;
end


sessList = dir(rootSessDir);  
sessList = {sessList([sessList.isdir] & ~strncmp({sessList.name},'.',1)).name}; 


%- loop over sessions
for iSess=1:length(sessList),
    
    fprintf('\n  In Session %s:', sessList{iSess});
    
    sortSubDir = dir(fullfile(rootSessDir,sessList{iSess},'sorted'));
    sortSubDir = {sortSubDir([sortSubDir.isdir] & ~strncmp({sortSubDir.name},'.',1)).name}; 
   
    %- loop over sorts within that session (can be more than one, e.g. different sorters or different referencing used
    for iSortDir=1:length(sortSubDir),
        
        fprintf('\n     In sortDir %s:', sortSubDir{iSortDir});
        sortSubDirPath = fullfile(rootSessDir,sessList{iSess},'sorted',sortSubDir{iSortDir});
    
        %- SOURCE:      path to OLD sorts that need to be renamed
        sortTxtPathOld  = fullfile(sortSubDirPath,SRC_SORT_TXT_DIR); %- this is the source directory
        sortTxtFilesOld = dir([sortTxtPathOld '/*.txt']);   sortTxtFilesOld = {sortTxtFilesOld.name};
        if length(sortTxtFilesOld)==0,
            fprintf('\n     No Sort files found. Make sure they are in the correct subdirectory: \n%s',sortTxtPathOld);
            keyboard;
            continue;
        end
        
        
        %
        %-  around here is where you could check for .bin files from old sort... those would be compared to bin files from new sort to see if timing shoudl be shifted
        %
        
        
        %- DESTINATION: path to correctly named sorts. 
        sortTxtPathNew  = fullfile(sortSubDirPath,DST_SORT_TXT_DIR);  %- this is the target output directory of the rename
        sortTxtFilesNew = dir([sortTxtPathNew '/*.txt']);   sortTxtFilesNew = {sortTxtFilesNew.name};
        if length(sortTxtFilesNew)==0,
            if ~exist(sortTxtPathNew,'dir'),
                fprintf('\n     Created an output directory to contain the corrected sorts: %s',sortTxtPathNew);
                mkdir(sortTxtPathNew);
            end
        end
        
        
        %- loop over the actual sort text files
        iKeyCopied = [];
        for iSortFile = 1:length(sortTxtFilesOld),
           
            oldSortName = sortTxtFilesOld{iSortFile}(1:end-4); %- remove .txt from end so it matches the jacksheet entry
            oldSortPath = fullfile(sortSubDirPath,SRC_SORT_TXT_DIR,sprintf('%s.txt',oldSortName));
            
            iKey = find(strcmp(renameKey.SortNameOld,oldSortName));
            
            if length(iKey)==1,
                newSortName = renameKey.SortChanName{iKey};
                newSortPath = fullfile(sortSubDirPath,DST_SORT_TXT_DIR,sprintf('%s.txt',newSortName));
                
                if     exist(newSortPath,'file') & OVERWRITE_DEST_FILE==0,
                    fprintf('\n         file %s --> %s [new file already exists, so no action taken]', oldSortName,newSortName);
                else
                    if exist(newSortPath,'file'),
                        fprintf('\n         file %s --> %s [new file already exists, but was overwritten]', oldSortName,newSortName);
                    else
                        fprintf('\n         file %s --> %s', oldSortName,newSortName);
                    end
                    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(oldSortPath, newSortPath,'f');
                    if ~SUCCESS,
                        fprintf('  --- Uh oh... copy error');
                        keyboard;
                    else
                        fprintf('  --- copied');
                        iKeyCopied = [iKeyCopied iKey];
                    end
                end
            elseif length(iKey)==0
                fprintf('\n         file %s --> No match found in key', oldSortName);
            else
                fprintf('\n         file %s --> ERROR: multiple matches found in key', oldSortName);
            end
        end
        if length(iKeyCopied)>0,
            %- save a table that lists the channels that were copied
            renameKeyCopied = renameKey(iKeyCopied,:);
            fileNameSave = fullfile(sortSubDirPath,SRC_SORT_TXT_DIR,sprintf('renameLookupTable[files copied %s].csv',datestr(now,'yymmdd_HHMM')));
            writetable(renameKeyCopied,fileNameSave);
            fprintf('\n  %d files copied, renameLookupTable.csv saved to %s', height(renameKeyCopied), fileNameSave);
        end
        
    end
    
end

fprintf('\n all done \n');



