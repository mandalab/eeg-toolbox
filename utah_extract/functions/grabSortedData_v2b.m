function grabSortedData(subj,rootMicroDir,sorted_micro_dir,sessFolders)
%
% Go through a processed folder, look for any sorts (.txt files) in the spikes_XXX folders, and move them to a "sorted" folder for subsequent processing
%
% utahDir should point to root directory of subjects for micro processing.
%
%  NOTE: this function has added functionality intended to convert or deal with odl
%
%  This function looks for _ready sessions, which the sorter has marked.  
%
%
% JW 3/2019
% MK 4/2019... started layering in functionality of grabbing legacy sorts... possibly should get split to a different function?
% JW 9/2019 tweak all path pointers to prepare for new plan of 56PROC for processing, 56PUB for sorting
% SJ 3/13/2020 - added 'end' at the bottom, just in case this solves the freezing of 56 problem


RENAME_GRABBED_SORTS = 1;  % default (1):  set to 0 if you DONT want to remane the grabbed text files to sort_grabbed.txt


manualReadyString = '_ready';  %- this string is manually appended to a session once manual sorting is complete... flag that it is ready to grab


grabDateStr = datestr(now,'YYmmdd');



fprintf('processing subject %s (current datestring = %s)',subj,grabDateStr)


%- if source is 56PUB/micro_forSorting, the session directories are flatter (no "processed_SPK" first)
%- however, we want this code to be able to work for sorts that are stored on pristine as well... how to we make that branch happen?

%- specify path to target (processing folder) and source (spike sort folder)... roots *can* be the same, but intent is to move to 56PUB for sorts and 56PROC for process
sortSrcRoot = fullfile(sorted_micro_dir,subj,'processed_SPK');
if ~exist(sortSrcRoot,'dir'),
    %- this is intended to catch the flattened version of processed_SPK data on 56PUB, where sessions are one directory up
    sortSrcRoot = fullfile(sorted_micro_dir,subj);
end
sortDstRoot = fullfile(rootMicroDir,subj,'sorts_manual');  %- manual sorts all go here... kilosort/mountainsort/etc will go in equivalent root
fprintf('\n grab sorts from %s \n         move to %s',sortSrcRoot,sortDstRoot);


%- if user doesn't pass in session folder list, then check all of the folders
if nargin<4 || length(sessFolders)==0,
    sessFolders = dir(sortSrcRoot);
    sessFolders = {sessFolders([sessFolders.isdir]).name};
    sessFolders = sessFolders(~strcmp(sessFolders,'.') & ~strcmp(sessFolders,'..') & ~strcmp(sessFolders,'_extraction_notes'));
    numAllSess  = length(sessFolders);
    sessFolders = sessFolders(contains(sessFolders,manualReadyString));
    fprintf('\n      >> Found %d potential session folders; %d designated as "_ready" to be grabbed and processed', numAllSess, length(sessFolders));
end

numSamples = 20; % samples to use for finding matches with legacy .bins



%- loop over data_processedSKP/session folders, find the spike folders within (e.g., spikes_noreref; spikes_reref; etc)
for s = 1:length(sessFolders)
    
    %- look inside this session folder
    sessID       = sessFolders{s};
    sessIDclean  = strrep(sessID,manualReadyString,''); %- remove "_ready"... ready will get replaced with "grabbed" once this processs is done
    spikeDirList = dir(fullfile(sortSrcRoot,sessID,'sort_*'));
    spikeDirList = {spikeDirList([spikeDirList.isdir]).name};
    
    %- little check to make sure old naming convention isn't being used by mistake
    spikeDirList_old = dir(fullfile(sortSrcRoot,sessID,'spikes_*'));
    if ~isempty(spikeDirList_old),
        fprintf('\n Hold up... looks like this has old naming conventions.  Physio splits for sorts should be saved in sort_reref, sort_noreref, etc)');
        keyboard;
        error('\n fix this... code below assumes 4 character name and will cut things off');
    end

    if isempty(spikeDirList),
        fprintf('\n No spike directory found: %s',sessID);
        keyboard;
        continue;
    end

    
    %- is there a sortNotes.xls in this session folder?  By default there should be one with "??" at the end.
    %   we want to take the one with user initials and then use those initials to name the subsequent sorts folder
    %   i.e.,  if sortNotes_sortedByJW.xlsx,  move the sorts to sorts_manual/session/reref_sortedByJW
    sortNotesList = dir(fullfile(sortSrcRoot,sessID,'sortNotes(*)_sortedBy*.xlsx'));
    if length(sortNotesList)==1,
        sortNotesUse = fullfile(sortSrcRoot,sessID,sortNotesList(1).name);
        %- pull the sorter initials from the sortNotes name
        sortInitialsUse = sortNotesList(1).name(strfind(sortNotesList(1).name,'sortedBy')+8:strfind(sortNotesList(1).name,'.xlsx')-1);
    elseif length(sortNotesList)>1,
        %- if more than one, try cutting out "??" one if it exists
        sortNotesList = sortNotesList(~contains({sortNotesList.name},'sortedBy??')); %- get rid of unnamed version
        if length(sortNotesList)>1,
            fprintf('\n\n ERROR: more than one sort notes with sorter initials in %s\n decide which one is valid and delete the other', fullfile(sortSrcRoot,sessID));
            keyboard
            error('\n two sortNotes found in %s... remove the invalid one', fullfile(sortSrcRoot,sessID));
        end
        sortNotesUse = fullfile(sortSrcRoot,sessID,sortNotesList(1).name);
        %- pull the sorter initials from the sortNotes name
        sortInitialsUse = sortNotesList(1).name(strfind(sortNotesList(1).name,'sortedBy')+8:strfind(sortNotesList(1).name,'.xlsx')-1);
    else
        fprintf('\n no sort notes found in %s. Should at least be the "??" one.  Kinda weird.  Maybe ok for legacy?', fullfile(sortSrcRoot,sessID));
        keyboard;
        sortNotesUse    = '';
        sortInitialsUse = '??';
    end
    
    
    %- this function will move sorts into a parallel path... instead of data_processedSPK/sess it will be sorts_manual/sess.
    %    check to see if that has already been done, and if so potentially add new sorts
    %    Note: the sort folders *should* have been modifed with sorter initials, so use dir instead of just specifying full target
    sortDirSess = fullfile(sortDstRoot,sessIDclean);
    if exist(sortDirSess,'dir'),
        sortDirList = dir(fullfile(sortDirSess));
        sortDirList = {sortDirList([sortDirList.isdir]).name};
        sortDirList = sortDirList(~strcmp(sortDirList,'.') & ~strcmp(sortDirList,'..'));
    else
        sortDirList = {};
    end
    
    
    
    %- loop over spike direcotries within this session (e.g., spikes_reref/spikes_noreref/spikes_reref(mask5ms))
    %     and look for sort.txt files.  If found, rename them (sort_grabbed.txt) and copy to sorted/session
    for iDir=1:length(spikeDirList)
        
        fprintf('\ngoing into %s, %s\n',spikeDirList{iDir},sessID)
        
        spikeDir    = spikeDirList{iDir}; %- e.g., spikes_reref/spikes_noreref/spikes_reref(mask5ms)
        spikeDirCln = spikeDir(6:end); %- cut off "sort_" from the string
        
        
        %- find the matching sort directory or create one
        iSortDir    = find(strncmp(sortDirList,sprintf('%s_',spikeDirCln),length(spikeDirCln)+1));
        if length(iSortDir)>1,
            fprintf('\n uh oh, multiple sort directories with same name. not sure what to do here');
            keyboard;
        elseif length(iSortDir)==1,
            sortDir = fullfile(sortDirSess,sortDirList{iSortDir},'sort_txt');
        else
            %- sort directory not created yet... get ready to create it in case we find sorts
            sortDir = fullfile(sortDirSess,sprintf('%s_sortedBy%s',spikeDirCln,sortInitialsUse),'sort_txt');
        end
        
        
        %- if legacy, make sure we are pointing to the correct legacy version (some sessions have multiple legacy sorts)
        if strcmpi(spikeDirCln,'legacy')
            spikeDir = sprintf('%s/%s',spikeDir,'v_use');
            fprintf('\n code for legacy processing commented out below... need to make sure everything still makes sense after tweaking file locations');
            keyboard;  %- JW says Mo was working on this... not clear what organization is supposed to be for this to work
        end
        
        
        %- look for text files in the spikes folder
        binFiles  = dir(fullfile(sortSrcRoot,sessID,spikeDir,'*.bin'));
        binFiles  = {binFiles.name};
        sortFiles = dir(fullfile(sortSrcRoot,sessID,spikeDir,'*.txt'));
        sortFiles = {sortFiles.name};
        sortFiles = sortFiles(~contains(sortFiles,'readme') & ~contains(sortFiles,'params')); %- should just be sort text files left
        
          
        iSortGrab = find(contains(sortFiles,'(grabbed'));
        iSortNeed = find(~contains(sortFiles,'(grabbed'));
        
        
        
        
        %- At this point we know which files there are to move... so move then!
        fprintf('\n  ******   Session %s, %s:  %d bins, %d text files total, %d already marked as grabbed   *****', sessID, spikeDir, length(binFiles), length(sortFiles), length(iSortGrab));
        
        if length(sortFiles)==0,  %  |  length(iSortNeed)==0,
            fprintf('\n  no sort.txt files found, skipping this session'); 
            
        else
            %- getting here means there are text files that haven't been grabbed
            fprintf('\n %d sort text files to copy and rename;  %d sort text files that were previously grabbed but we will confirm we have a copy of in PROC', length(iSortNeed), length(iSortGrab));
            
            iSortNeed = 1:length(sortFiles); %- forces the code below to treat every text file as a potential "need"... checking that the preivously grabbed ones are in fact good to go.
            
            
            %- create the sort folder if not created yet: "sorts_manual/session/reref_sortedByXX"
            if ~exist(sortDir,'dir'),
                
                %- sorter initials are derived from sortNotes, else are "LE" if a legacy sort.  If sortNotes were not renamed (i.e., still called sortNotes_sortedBy??.xlsx), then prompt here.
                if strcmp(sortInitialsUse,'??')
                    sortInitialsUse = input(sprintf('\n\n About to create a new "sorted" foder for %s, but sortNotes were not labeled with sorter... \n\n Please enter the initials of the sorter: (?? default)  ',tempDir),'s');
                    if isempty(reply),
                        sortInitialsUse = '??';
                    else
                        sortDir = fullfile(sortDirSess,sprintf('%s_sortedBy%s',spikeDirCln,sortInitialsUse),'sort_txt');
                    end
                end
                
                %- make the directory
                mkdir(sortDir);
                fprintf(' --> created directorty %s\n', sortDir);
                
                %- copy over the extrationInfo folder and all its contents
                srcDir = fullfile(sortSrcRoot,sessID,spikeDir,'_extractionInfo');
                dstDir = fullfile(sortDir(1:strfind(sortDir,'sort_txt')-1),'_extractionInfo');
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(srcDir,dstDir,'f');
                if ~SUCCESS,
                    fprintf('\n uh oh... extrationInfo folder didnt copy to sorts folder');
                end
                
            end
            
            if ~isempty(sortNotesUse),
                %- always overwrite... if doing new sorts get the most recent version of the notes file
                srcDir = sortNotesUse;                                       %- specify source file
                dstDir = fullfile(sortDir(1:strfind(sortDir,'sort_txt')-1)); %- specify destination directory
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(srcDir,dstDir,'f');   %- copy file to new dir
            end
            
            
            %- loop over the sort.txt files, copy them, and then rename the ones left behind so we know we got them
            previouslyGrabbedAreOK = 0;
            
            for iF=1:length(iSortNeed),
                
                %- source files will be renamed "(grabbedXXYYZZ)"... once renamed they will not be grabbed a second time
                %  destination files will be named "(grabbedXXYYZZ)".  
                %    For each newly grabbed source, we need to identify the "clean" (non-grabbed name) of the destination files to see if there is a conflict
                
                srcFile = fullfile(sortSrcRoot,sessID,spikeDir,sortFiles{iSortNeed(iF)});
                if contains(sortFiles{iSortNeed(iF)},'(grabbed'),
                    renameFile = sortFiles{iSortNeed(iF)}; %- previously grabbed... dont relabel an old grab
                    dstFileChecked = fullfile(sortDir,renameFile);           %- look to see if this old grab is present
                else
                    renameFile = sprintf('%s_(grabbed%s).txt',sortFiles{iSortNeed(iF)}(1:end-4),grabDateStr);
                    dstFileChecked = fullfile(sortDir,sprintf('%s_(grabbed*).txt',sortFiles{iSortNeed(iF)}(1:end-4)));
                end
                srcFileRenamed = fullfile(sortSrcRoot,sessID,spikeDir,renameFile);
                dstFileRenamed = fullfile(sortDir,renameFile);
                
                
                
                %- check for existance.  if conflict exists, move the old version to "old sorts renamed"
                dstFileCheck = dir( dstFileChecked );
                if ~isempty(dstFileCheck),
                    if length(dstFileCheck)>1,
                        fprintf('\n not sure how to deal with 2 collisions... below assumes one');
                        keyboard;
                    end
                    dstInfo = dstFileCheck(1); %- assume just one collision
                    srcInfo = dir(srcFile);
                    
                    if dstInfo.bytes ~= srcInfo.bytes | ~strcmp(dstInfo.date,srcInfo.date),
                        %- source and destination files are different.  Update the destination
                        fprintf('\n HEADS UP: detected an updated sort that is colliding with an old sort in destination folder. Old sort will be saved in _old_sorts_redone',srcFile);
                        
                        redone_sort_folder = '_old_sorts_redone';
                        redone_sort_path = fullfile(sortDir,redone_sort_folder);
                        if ~exist(redone_sort_path,'dir'), mkdir(redone_sort_path); end
                        SUCCESS = movefile( fullfile(sortDir,dstInfo.name), fullfile(sortDir,redone_sort_folder,dstInfo.name) );
                        if SUCCESS==0,
                            fprintf('\n move to old failed... do it manually?');
                            keyboard
                        end
                    elseif  contains(sortFiles{iSortNeed(iF)},'(grabbed') && strcmp(sortFiles{iSortNeed(iF)},dstInfo.name),
                        %- source and destination files are identical, and they were grabbed previously
                        previouslyGrabbedAreOK = previouslyGrabbedAreOK+1;
                        %fprintf('\n %d of %d) %s previously grabbed file matches (no action)', iF, length(iSortNeed), sortFiles{iSortNeed(iF)});
                        continue;
                    else
                        %- collision is actually the same file, but with different files names or not grabbed previously
                        fprintf('\n kinda a weird situation... there is a potential collision (same filename root in PUB and PROC)\n File contents are the same, yet name is different');
                        fprintf('\n this can happen if an older sort is copied into a folder with the identical "grabbed" file, but does that actually happen?');
                        fprintf('\n not sure how the code below will deal with this or if we even have to think about it... possible course of action would be to skip this channel');
                        keyboard;
                        error('\n decide what to do then fix it');
                    end
                end
                
                %- do the copy and check for errors
                [SUCCESS1,MESSAGE,MESSAGEID] = copyfile(srcFile,dstFileRenamed,'f');
                if SUCCESS1==0, fprintf('\n uh oh. copy error'); keyboard; end
                
                %- for some reason move failed below.  make an empty rename file then delete the real file
                %[SUCCESS2,MESSAGE,MESSAGEID] = copyfile(srcFile,renameFile,'f');
                %if SUCCESS2==0, fprintf('\n uh oh. copy error'); keyboard; end
                
                %- less desirable version: make an empty rename file then delete the real file
                %[FID, MESSAGE] = fopen(renameFile,'a+');
                %if FID==0, fprintf('\n cant fopen, probably a permissions issue'); keyboard;
                %else       fclose(FID); end
                
                %- hang up on this one during debugging...
                if RENAME_GRABBED_SORTS && ~strcmp(srcFile,srcFileRenamed),
                    [SUCCESS2,MESSAGE,MESSAGEID] = movefile(srcFile,srcFileRenamed,'f');
                    if SUCCESS2==0,
                        fprintf('\n first move attempt didnt work... create the file then move (does this still happen?)');
                        keyboard
                        [FID, MESSAGE] = fopen(srcFileRenamed,'a+');
                        if FID==0,
                            fprintf('\n cant fopen, probably a permissions issue');
                            keyboard;
                        else
                            fclose(FID);
                            delete(srcFileRenamed)
                            [SUCCESS2,MESSAGE,MESSAGEID] = copyfile(srcFile,srcFileRenamed,'f');
                            if SUCCESS2==0,
                                fprintf('\n second attempt didnt work...');
                                keyboard
                            end
                        end
                    end
                end %- if RENAME_GRABBED_SORTS,
                
                %- double check new files were made, if so, delete the original
                if exist(srcFileRenamed,'file') & exist(dstFileRenamed,'file'),
                    fprintf('\n %d of %d) %s copied and renamed', iF, length(iSortNeed), sortFiles{iSortNeed(iF)});
                    
                    if exist(srcFile) && ~strcmp(srcFile,srcFileRenamed),
                        delete(srcFile);
                        fprintf('\n Uh oh, cant delete source file: %s', srcFile);
                        keyboard;
                        
                    else
                        fprintf(' (and source deleted)');
                    end
                    
                end %if exist rename file
                
            end %- for sortNeed
             
            if length(iSortGrab)>0 & previouslyGrabbedAreOK>0,
                fprintf('\n ** %d of %d previously grabbed files were confirmed to be the same as in the destination (no action taken)', length(iSortGrab),previouslyGrabbedAreOK);
            end
                      
            
        end %- if there are files to sort
        
    end %- for iDir=spkDirList (look over reref/ reref(mask5ms) / noreref, etc)
    
    
    %- now rename the session folder from _ready to _(grabbed)
    if contains(sessID,manualReadyString),
        scrPath = fullfile(sortSrcRoot,sessID);
        dstPath = fullfile(sortSrcRoot,sprintf('%s_(grabbed)',sessIDclean));
        [SUCCESS4,MESSAGE,MESSAGEID] = movefile(scrPath,dstPath,'f');
    end
    
end %- for 1:length(sessFolders)


fprintf('\n  <done grabbing %d sessions>\n', length(sessFolders));

end %SJ
