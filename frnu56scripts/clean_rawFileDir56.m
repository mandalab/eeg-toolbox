function SUCCESS = clean_rawFileDir56(path2data_raw, SUBJ)
%
%%  function organizeDataRawCervello
%
%  this is the function that clean up the messy outputs from Cervello by putting sessions into session folders 
%
%
%  OLD VERSION... used to rename files
%
%  NEW VERSION (this version), just lumps the files together into the same session folder
%
%
%% Code takes files that are saved by cervello system on the research
% computer, and renames them so that they are have folders and titles
% corresponding to proper date-time.
%
% Code should be carefully run and output inspected with REAL_MOVE set to
% 0, as Cervello file names can change, causing incorrect and irreversible
% file name changes.
%
%  from Cervello-Controlled Saves, example file name:  yJEidTiJzzz-20170310-115925-INST1.ccf  --->  170310_1159
%%
%
%  JW 4/2019 -- slightly tweaked working version and commited to SVN  (created about a year before that but maintained locally)
%

%clear all;


%- Copy from full local eeg folder to eegTrim
%SUBJ = 'NIH069';
%path2data_raw = sprintf('/Volumes/56E/UTAH_E/%s/data_raw',SUBJ);     %- office-local (could change this to server)
%path2data_raw = sprintf('/Volumes/72E/UTAH_E/%s/data_raw',SUBJ);     %- office-local (could change this to server)

SUCCESS = 0;


%
GET_EEG_ROOT = path2data_raw;

%- First confirm that NSP1 and NSP2 data are combined and in "data_raw"
if contains(GET_EEG_ROOT,'dataNSP1') | contains(GET_EEG_ROOT,'dataNSP1'),
    fprintf('\n ERROR: EEG_ROOT (%s) contains dataNSP1 or dataNSP2 folder.  Aborting rename.\n Combine all loose blackrock files from *both* NSP1 and NSP2 in NIHXXX/data_raw then rerun. \n',GET_EEG_ROOT);
    keyboard;
    return;
end
parts = strsplit(GET_EEG_ROOT,'/');
if ~strcmp(parts{end},'data_raw'),
    fprintf('\n ERROR: EEG_ROOT (%s) does not end in "data_raw".  Aborting rename.\n Combine all loose blackrock files from *both* NSP1 and NSP2 in NIHXXX/data_raw then rerun. \n',GET_EEG_ROOT);
    keyboard;
    return;
end




%REAL_MOVE = 1;  %- Now just always use 1   %%% old way... set to 0 initally... only set to 1 when you know everything is working ok.
MULTI_SUBJ_NAME_OK = 0;

str2repeat = '';
for CHECK_FILES_IN_FOLDERS = [0 1],
    
    %- for CHECK_FILES_IN_FOLDERS, first do a quick check that all the files have been moved to folders
    if CHECK_FILES_IN_FOLDERS,
        %- quick check that there are not loose files in the folder first
        fList  = dir(GET_EEG_ROOT);
        fDirs  = {fList(~[fList.isdir]).folder};    %- track the actual directory... this will be used later for the CHECK_FILES_IN_FOLDER option
        fNames = {fList(~[fList.isdir]).name};
        
        iF = find(strcmp(fDirs,GET_EEG_ROOT));
        foundBad = 0;
        for iiF=iF,
            if ~(strcmp(fNames{iiF}(1),'.') | strcmp(fNames{iiF}(1),'_'))
                fprintf('\n UNEXPECTED FILE IN ROOT: %s ... need to put all files into folders before running "CHECK_FILES_IN_FOLDERS"\n',fNames{iiF});
                foundBad=1;
            end
        end
        if foundBad,
            fprintf('\n breaking out\n');
            keyboard;
            return;
        end
    end
    
    
    
    %- get the files
    fprintf('\n populating list of folders or files now...');
    if CHECK_FILES_IN_FOLDERS,
        fList       = dir([GET_EEG_ROOT '/*/*']); %- get the files in the folders
        fileNameMod = '_post';
        REAL_MOVE   = 0;
    else
        fList       = dir(GET_EEG_ROOT);
        fileNameMod = '';%
        REAL_MOVE   = 1;
    end
    fDirs  = {fList(~[fList.isdir]).folder};    %- track the actual directory... this will be used later for the CHECK_FILES_IN_FOLDER option
    fNames = {fList(~[fList.isdir]).name};
    [~,iSort] = sort(fNames);  iSort = iSort(end:-1:1);
    fNames = fNames(iSort); fDirs=fDirs(iSort);  %- order in reverse... this trick will correctly identify the short session when two sessions start in the same minute
    
    
    %- check to see if anything there
    if ~exist(GET_EEG_ROOT,'dir') | length(fList)==0,
        fprintf('\n Uh oh.... nothing in directory: %s\n', GET_EEG_ROOT);
        keyboard;
        return;
    end
    
    
    %- Force a not-real-move first
    %     if REAL_MOVE == 1
    %         masterFileNamesFile = sprintf('%s/_fileNameCounts%s.txt',GET_EEG_ROOT,fileNameMod);
    %
    %         if ~exist(masterFileNamesFile,'file'),
    %             fprintf('\n Looks like this is the first pass. Setting "REAL_MOVE" to zero so you can check first');
    %             REAL_MOVE = 0;
    %         else
    %             fprintf('\n About to do a real move... continue if youre sure you want to risk this fuck-up (risk is pretty low, but present)\n')
    %             keyboard
    %         end
    %     end
    
    
    %- create a master text file that tracks all the changes
    if CHECK_FILES_IN_FOLDERS,
        fTrackChangesPath = sprintf('%s/_fileNameChanges%s.txt',GET_EEG_ROOT,fileNameMod);
        fTrackChanges     = fopen(fTrackChangesPath,'w+');
        fprintf(fTrackChanges,'\n\n******************%s; real-move=%d ***********************\n*******************  generated: %s  *****************',GET_EEG_ROOT,REAL_MOVE,datestr(now));
        fTrackCountsPath  = sprintf('%s/_fileNameCounts%s.txt',GET_EEG_ROOT,fileNameMod);
        fTrackCounts      = fopen(fTrackCountsPath,'w+');
        fprintf(fTrackCounts,'\n\n******************%s; real-move=%d ***********************\n*******************  generated: %s  *****************',GET_EEG_ROOT,REAL_MOVE,datestr(now));
    end
    
    newFileNameList = {};  thisSubjRoot = {};  newSessList={}; numNewNames=0; numCopied=0;
    suffixListCV = {'INST0.nev','INST0.ns2','INST0.ns3','INST0.ns5','INST0.ns6',... %SJ added .ns6
        'INST1.nev','INST1.ns2','INST1.ns3','INST1.ns5','INST1.ns6',...%SJ added .ns6
        'INST0.ccf','INST0.csr','INST0.toc','INST0.sif','INST1.ccf','INST1.csr','INST1.toc','INST1.sif','catch'};
    % suffixes that we dont use but could in the future: {'INST0.ns1','INST0.ns4','INST1.ns1','INST1.ns4'}
    suffixListMAT1 = {'ieeg.nev','ieeg.ns2','ieeg.ns3','ieeg.ns5','ieeg.ns6',...%SJ added .ns6
        'utah.nev','utah.ns2','utah.ns3','utah.ns5','utah.ns6',...%SJ added .ns6
        'ieeg.ccf','ieeg.csr','ieeg.toc','ieeg.sif','utah.ccf','utah.csr','utah.toc','utah.sif','catch'};
    suffixListMAT2 = {'ieeg1.nev','ieeg1.ns2','ieeg1.ns3','ieeg1.ns5','ieeg1.ns6',...%SJ added .ns6
        'ieeg2.nev','ieeg2.ns2','ieeg2.ns3','ieeg2.ns5','ieeg2.ns6',...%SJ added .ns6
        'ieeg1.ccf','ieeg1.csr','ieeg1.toc','ieeg1.sif','ieeg2.ccf','ieeg2.csr','ieeg2.toc','ieeg2.sif','catch'};
    suffixListMAT3 = {'micro.nev','micro.ns2','micro.ns3','micro.ns5','micro.ns6',...%SJ added .ns6
        'utah.nev','utah.ns2','utah.ns3','utah.ns5','utah.ns6',...%SJ added .ns6
        'micro.ccf','micro.csr','micro.toc','micro.sif','utah.ccf','utah.csr','utah.toc','utah.sif','catch'};      %- NIH050
    suffixListMAT4 = {'utah1.nev','utah1.ns2','utah1.ns3','utah1.ns5','utah1.ns6',...%SJ added .ns6
        'utah2.nev','utah2.ns2','utah2.ns3','utah2.ns5','utah2.ns6',...%SJ added .ns6
        'utah1.ccf','utah1.csr','utah1.toc','utah1.sif','utah2.ccf','utah2.csr','utah2.toc','utah2.sif','catch'};  %- NIH034
    suffixListMAT5 = {'ieeg.nev','ieeg.ns2','ieeg.ns3','ieeg.ns5','ieeg.ns6',...%SJ added .ns6
        'utah1.nev','utah1.ns2','utah1.ns3','utah1.ns5','utah1.ns6',...%SJ added .ns6
        'ieeg.ccf','ieeg.csr','ieeg.toc','ieeg.sif','utah1.ccf','utah1.csr','utah1.toc','utah1.sif','catch'};      %- NIH036
    suffixListMAT6 = {'utah_u.nev','utah_u.ns2','utah_u.ns3','utah_u.ns5','utah_u.ns6',...%SJ added .ns6
        'utah_m.nev','utah_m.ns2','utah_m.ns3','utah_m.ns5','utah_m.ns6',...%SJ added .ns6
        'utah_u.ccf','utah_u.csr','utah_u.toc','utah_u.sif','utah_m.ccf','utah_m.csr','utah_m.toc','utah_m.sif','catch'};  %- NIH037
    
    %
    suffixList = suffixListCV; %- assume cervello based suffix, with INST0 and INST1

    if strcmp(SUBJ,'NIH029'), suffixList = suffixListMAT1; end %-
    if strcmp(SUBJ,'NIH030'), suffixList = suffixListMAT1; end %-
    if strcmp(SUBJ,'NIH034'), suffixList = suffixListMAT4; end %-
    if strcmp(SUBJ,'NIH036'), suffixList = suffixListMAT5; end %-
    if strcmp(SUBJ,'NIH037'), suffixList = suffixListMAT6; end %-
    
%     ind_nev_ns = cell2mat(regexp(suffixList,'.*\.((ns\d)|(nev))'));
%     ind_nev_ns2 = regexp(suffixList,'.*\.((ns\d)|(nev))');
%     ind_nev_ns2(2) = {[]}
%     ind_nev_ns3 = cell2mat(ind_nev_ns2)
    
    foundAnyUnexpected = 0;
    suffixCounts = zeros(length(fNames),length(suffixList));
    for iF=1:length(fNames)
        fileName = fNames{iF};
        if fileName(1)=='.' | contains(fileName,'LogFile') | contains(fileName,'rename') | contains(fileName,'fileName') | contains(fileName,'jacksheetBR') | strcmpi(fileName,'notes.txt'),
            continue;
        end
        fileDirActual = fDirs{iF};
        
        
        %%- EXTRACT OUT THE SESSION FOLDER FROM THE NAME, and THE SUFFIX to CONFIRM EXPECTED SUFFIX
        parts         = strsplit(fileName,'-');
        if length(parts)==1 & (contains(fileName,'ieeg') | contains(fileName,'utah') | contains(fileName,'micro')),
            if CHECK_FILES_IN_FOLDERS,
                subjRoot = '';
                newSessStrOG = fileName(1:11); %- YYMMDD_HHMM
                fileSuffixOG = fileName(13:end);
                fileSuffixOGnoBEH = regexprep(fileSuffixOG,'_beh','');
                fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'_premie',''); %- make ns5 and ns5 fall into the same slot <<- subj<NIH050 used ns6 for micro
                %fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'.ns6','.ns5'); %- make ns5 and ns5 fall into the same slot <<- subj<NIH050 used ns6 for micro
                fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'.ns4','.ns3'); %- make ns4 and ns3 fall into the same slot <<- subj<NIH050 used ns4 for pulses
                if sum(contains(suffixList,fileSuffixOGnoBEH))>0,
                    suffixList = suffixList; %- for unusual subjects, like NIH036, override happens above
                elseif sum(contains(suffixListMAT1,fileSuffixOGnoBEH))>0,
                    suffixList = suffixListMAT1;
                elseif sum(contains(suffixListMAT2,fileSuffixOGnoBEH))>0,
                    suffixList = suffixListMAT2;
                elseif sum(contains(suffixListMAT3,fileSuffixOGnoBEH))>0,
                    suffixList = suffixListMAT3;
                elseif sum(contains(suffixListMAT4,fileSuffixOGnoBEH))>0,
                    suffixList = suffixListMAT4;
                else
                    fprintf('\n uh oh... no matches: %s \n',fileName);
                    keyboard;
                end
                
            else
                fprintf('\n shoudlnt happen. break out %s\n', fileName);
                keyboard;
                return;%
            end
            %- matlab based file or a filename that was already converted
            %subjRoot = '';
            
        else
            %- standard approach... parse the file name to see what session folder it SHOULD go into
            subjRoot      = parts{1};%
            newSessStrSEC = sprintf('%s_%s',parts{end-2},parts{end-1}); %- Y%YYYMMDD_HHMMSS
            newSessStrOG  = newSessStrSEC(3 : end-2); % YYMMDD_HHMM
            fileSuffixOG  = parts{end};%
            suffixList    = suffixListCV; %- cervello list
            
        end
        
        
        %- confirm just one subject in this folder...
        if iF==1,
            thisSubjRoot = parts{1};
            sessStrCnt = 1;
        elseif ~strcmp(thisSubjRoot, subjRoot) & length(subjRoot)>0 & MULTI_SUBJ_NAME_OK==0,
            fprintf('\n\n Uh oh... looks like this folder contains multiple subjects: %s, %s\n found the second subject name in %s\n', thisSubjRoot, parts{1},newSessStrOG);
            reply = input('Should we allow two subject names for one subject (almost never true)','s');
            if isempty(reply) | strcmpi(reply,'n'),
                reply = 'N';
                fprintf('\n remove 2nd subject folder (or burry within a subfolder) and rerun.\n aborting\n');
                return;
            elseif strcmpi(reply,'y'),
                MULTI_SUBJ_NAME_OK = 1;
            end
        end
        
        
        %- keep a tally of all the files found for each session... are some missing in specific files?
        unexpectedSuffix = '';
        fileSuffixOGnoBEH = regexprep(fileSuffixOG,'_beh','');
        fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'_premie',''); %- make ns5 and ns5 fall into the same slot <<- subj<NIH050 used ns6 for micro
        %fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'.ns6','.ns5'); %- make ns5 and ns5 fall into the same slot <<- subj<NIH050 used ns6 for micro
        fileSuffixOGnoBEH = regexprep(fileSuffixOGnoBEH,'.ns4','.ns3'); %- make ns4 and ns3 fall into the same slot <<- subj<NIH050 used ns4 for pulses
        if sum(strcmp(suffixList,fileSuffixOGnoBEH))==0,
            if foundAnyUnexpected==0,
                fprintf('\n unexpected suffix found: %s/%s..."%s"\n shouldnt happen... take a look and then continue if seems OK\n', fileDirActual,fileName,fileSuffixOG);
                %keyboard;
                foundAnyUnexpected = 1;
            end
            iSuffList = find(strcmp(suffixList,'catch'));
            unexpectedSuffix = ' << unexpected file suffix';
            keyboard;
        else
            iSuffList = find(strcmp(suffixList,fileSuffixOGnoBEH));
        end
        
        
        %- carry over suffixes to the new file name
        if CHECK_FILES_IN_FOLDERS,
            newSessDir = fileDirActual;
        else
            if contains(fileName,'_beh'),    addBeh='_beh';       else addBeh=''; end
            if contains(fileName,'_premie'), addPremie='_premie'; else addPremie=''; end
            newSessDir = sprintf('%s/%s%s%s',GET_EEG_ROOT,newSessStrOG,addBeh,addPremie);
        end
        newSessStr = sprintf('%s/%s',newSessDir,fileName);
        
        if sum(strcmp(newSessList,newSessDir))==0,
            newSessList{end+1} = newSessDir;
        end
        iSL = find(strcmp(newSessList,newSessDir));
        suffixCounts(iSL,iSuffList)=suffixCounts(iSL,iSuffList)+1;
        
        
        %- create a loop to add "x's" when two or more files were created within the same minute...
        % NEW VERSION that preserves Cervello filename wyith seconds SHOULDNT HAVE THIS PROBLEM... session folder will just hold both without collision
        extraSessPerMin = '';
        while sum(strcmp(newFileNameList,newSessStr))>0,
            newSessDir = sprintf('%s_shortSess',newSessDir);  %- new way is to call the earlier session "shortSess"
            newSessStr = sprintf('%s/%s%s',newSessDir,newSessStrOG,fileSuffix);
            
            extraSessPerMin = ' << extra session in same minute';
        end
        newFileNameList{iF} = newSessStr;
        numNewNames = numNewNames+1;
        
        
        %- actual copy/rename
        if REAL_MOVE,
            if ~exist(newSessDir,'dir'), 
                mkdir(newSessDir); 
            elseif exist(fullfile(newSessDir,'jacksheetBR_complete.xls'),'file'), 
                fprintf('\n moving a file into a session folder that already has a jacksheetBR.  delete the jacksheet so it can be made fresh');
                delete(fullfile(newSessDir,'jacksheetBR_complete.xls')); 
            end;
            if exist(newSessStr),
                fprintf('\n uh oh... possibly two sets of files in same minute that wasnt caught? Dont continue until you know whats up!!!');
                keyboard;
            end
            [success,message,messageid] = movefile([GET_EEG_ROOT '/' fileName],newSessStr);
            if success==0,
                fprintf('\n uh oh... copy error');
                keyboard;
            else
                numCopied = numCopied+1;
                %- no need to save a renameLog anymore... not actually renaming, just moving into a folder. but lets keep this for the time being just in case
                %fid=fopen([newSessDir '/renameLog2.txt'],'a+');   %- named renameLog2 because this log is just moving the file into a folder, not actually renaming
                %fprintf(fid,'\n\n %s --> %s %s %s', fileName, newSessStr,extraSessPerMin,unexpectedSuffix);
                %fclose(fid);
            end
        end
        
        
        fprintf('\n %s --> %s %s %s', fileName, newSessStr,extraSessPerMin,unexpectedSuffix);
        
    end
    
    suffixCounts = suffixCounts(1:length(newSessList),:); %- clear out spacers
    strConclusion = sprintf('\n\n %d files in directory, %d files with potential rename, %d files actually moved \n\n',length(fNames),numNewNames,numCopied);
    fprintf(strConclusion)
    
    if CHECK_FILES_IN_FOLDERS,
        delete(fTrackChangesPath);
        
        %-sanity check.. session list should be in order because files strings are in order
        [~,iSess]=sort(newSessList); iDiff = iSess(diff(iSess)~=-1);
        if length(iDiff)>0, fprintf('\n WARNING: session list is not in alphabetical order... might mean a file name is off. check these\n'); disp(newSessList([iDiff(end:-1:1)])'); end
        
        fprintf(             '\n\n Session Dir:           sum key file counts NSP0, NSP1:   nev/ns2/ns3/ns5/ns6;  nev/ns2/ns3/ns5/ns6;  [catch] '); %SJ added ns6
        fprintf(fTrackCounts,'\n\n Session Dir:           sum key file counts NSP0, NSP1:   nev/ns2/ns3/ns5/ns6;  nev/ns2/ns3/ns5/ns6;  [catch] '); %SJ added ns6
        iFixUnexpected = [];
        iSplitMultiple = [];
        for iSC=length(newSessList):-1:1,
            catchStr = '';
            iComp = [1:8];
            if suffixCounts(iSC,end)>0,                                                  catchStr = '<<< UNEXPECTED FILE TYPES';            iFixUnexpected = [iFixUnexpected iSC]; end
            if (iSC>1) && (sum(suffixCounts(iSC,iComp) - suffixCounts(iSC-1,iComp))~=0), catchStr = [catchStr ' ** change in file counts']; end
            if max(suffixCounts(iSC,iComp))>1,                                           catchStr = '<<< MULTIPLE COPIES, MUST SPLIT DIR';  iSplitMultiple = [iSplitMultiple iSC];  end
            %fprintf(             '\n %s: %d INST0, %d INST1:   (%d/%d/%d/%d/%d,  %d/%d/%d/%d/%d)  [%d] %s', newSessList{iSC}, sum(suffixCounts(iSC,1:4)), sum(suffixCounts(iSC,5:8)), suffixCounts(iSC,[1:8 end]),catchStr);
            fprintf(             '\n %s: %d INST0, %d INST1:   (%d/%d/%d/%d/%d,  %d/%d/%d/%d/%d)  [%d] %s', newSessList{iSC}, sum(suffixCounts(iSC,1:5)), sum(suffixCounts(iSC,6:10)), suffixCounts(iSC,[1:10 end]),catchStr); %SJ changed to account for .ns6
            fprintf(fTrackCounts,'\n %s: %d INST0, %d INST1:   (%d/%d/%d/%d/%d,  %d/%d/%d/%d/%d)  [%d] %s', newSessList{iSC}, sum(suffixCounts(iSC,1:5)), sum(suffixCounts(iSC,6:10)), suffixCounts(iSC,[1:10 end]),catchStr); %SJ changed to account for .ns6
        end
        fclose(fTrackCounts);
        
        %- now force a folder-name change on the violators
        iFix = unique([iSplitMultiple iFixUnexpected]);
        if length(iFix)>0,
            
            %- make this step automatic
            fprintf('\n\n HEADS UP, looks like %d session directories need to be dealt \n with because of file splits (%d) and unexpected files (%d). See list above. \n', length(iFix), length(iSplitMultiple), length(iFixUnexpected));
            ADD_DIR_SUFFIX = 1;
            %reply = input('If you agree type "Y" to have the directories automatically renamed to make it easier to manually fix, else hit return: y or [n]  ','s');
            %if isempty(reply) | strcmpi(reply,'n'),
            %    ADD_DIR_SUFFIX = 0;
            %elseif strcmpi(reply,'y'),
            %    ADD_DIR_SUFFIX = 1;
            %end
            
            if ADD_DIR_SUFFIX,
                fprintf('\n\n RENAMING %d session directories that need to be dealt with:',length(iFix));
                for iSC = 1:length(iFix),
                    if exist(newSessList{iFix(iSC)},'dir'),
                        sessSuff = '';
                        if sum(ismember(iSplitMultiple,iFix(iSC))), sessSuff = [sessSuff '_FixSPLIT'];      mkdir([newSessList{iFix(iSC)} '_shortSess']); end
                        if sum(ismember(iFixUnexpected,iFix(iSC))), sessSuff = [sessSuff '_FixUNEXPECTED']; mkdir([newSessList{iFix(iSC)} '/bad_file_names']); end
                        
                        [success,message,messageid] = movefile(newSessList{iFix(iSC)},[newSessList{iFix(iSC)} sessSuff]);
                        if success==0,
                            fprintf('\n uh oh... copy error');
                            keyboard;
                        else
                            fprintf('\n     %s --> add %s', newSessList{iFix(iSC)}, sessSuff);
                        end
                    else
                        fprintf('\n shouldnt get here');
                        keyboard
                    end
                end
                
                %-
                
                if length(iSplitMultiple)>0,
                    fprintf('\n\n For any "_FixSPLIT" folder, move the short segment to the newly created _shortSess folder then delete the string "_FixSPLIT" from the dir');
                end
                if length(iFixUnexpected)>0,
                    fprintf('\n\n For any "_FixUNEXPECTED" folder, correct the name error if possible, else move them to a sub-sub-folder called "bad_file_names" then delete the string "_FixUNEXPECTED" from the dir');
                end
                fprintf('\n\n ^^^ FIX the tagged directories and then RERUN  ^^^ \n');
                
                %keyboard
                return; %- error return
            end
        end
        
        
    elseif numCopied>0,
        str2repeat = strConclusion;
        fprintf('\n about to double check the file contents after that move');
        pause(2);
    end
    
    
end %- CHECK_FILES_IN_FOLDERS

if length(str2repeat)>0,
    fprintf(str2repeat);
else
    fprintf('\n\n Completed with %d files in %d sesssions; %d copies executed this run. ',length(fNames),length(newSessList),numCopied);
end
fprintf('\n all done \n');
SUCCESS = 1;

