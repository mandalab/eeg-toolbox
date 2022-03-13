%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   JW's calls to eeg_toolbox alignment functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------------------%
%----------------------  Choose the subject of interest and root EEG directory -----------------------%




%rootEEGdir = '/Volumes/56PROC/private/.eeg_pristine/update/tempUpdate';  %- new path from frnu56v2 (partial update)
%rootEEGdir = '/Volumes/56PROC/eeg_processing/.eeg_pristine/update/eeg';         %- new path from frnu56v2 (complete update)
rootEEGdir = '/Volumes/56PROC/eeg_processing/.update/eeg';  %- new path SJ 10/8/19
rootPUBdir = '/Volumes/56PUB/readWrite/micro_forSorting';   %- path for generating sortBuild file SJ 10/17/19

%% ------------------------------------------------------------------------------------------------------------------------%%
%% SELECT WHICH PROCESS TO RUN... ONLY THE FIRST ONE PICKED WILL GET DONE
%% ------------------------------------------------------------------------------------------------------------------------%%

%- INDIVIDUAL SUBJECT CHECK and/or PREPandALIGN (not required, just here for convenience if having issues with ALL_SUBJ CHECKS)
INDIV_SUBJ_CHECK        = 0;
INDIV_SUBJ_PREPandALIGN = 0;
INDIV_SUBJ_STR          = 'NIH077';


%- FULL DIRECTORY CHECK and TASK_LIST
ALL_SUBJ_FINAL_CHECK    = 0;  %- creates a clean subject task list then runs check list to populate the xlsx file
ALL_SUBJ_SORT_BUILD     = 1;  %- creates the sortBuild.xlsx file, which summarizes which sorts have been completed; put in PUB and .update/eeg

FIND_LOCKED_FILES       = 0;  %- only need to run this is chronosync is having problems with transfering due to locked files





%% ------------------------------------------------------------------------------------------------------------------------%%
%% EXECUTION OF BLOCKS SELECTED ABOVE... ONLY ONE WILL RUN (each has a return statement)
%% ------------------------------------------------------------------------------------------------------------------------%%



if INDIV_SUBJ_CHECK,
    tabOut = subjectDirCheck(INDIV_SUBJ_STR,rootEEGdir,'writeout',1,'printout',1,'returnTable',1,'copyInvalidNotes',1);
    % function returnLogTable=subjectDirCheck(subj, rootEEGdir, varargin)
    %   'writeOut' (1 default) ---  saves to subj/
    %   'printOut' (1 default) --- does what?
    %   'returnTable' (0 default) ---
    %   'copyInvalidNotes' (1 default)
     return;
end


if INDIV_SUBJ_PREPandALIGN,
    fprintf('\n fresh jack split or just run with existing split... continue or run a particular line depending on what you want');
    keyboard;
    eegPrepAndAlign(INDIV_SUBJ_STR, rootEEGdir,  'freshJackandSplit',1, 'freshExtractAlign',1, 'batch',0, 'skipProcessReref', 1); %- first run through and try to find issues
    eegPrepAndAlign(INDIV_SUBJ_STR, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',0, 'batch',0, 'skipProcessReref', 0);%- then the slow step of making processed
    return;
end




if ALL_SUBJ_FINAL_CHECK,
    %-  destroys then creates a new subjBuild... run this before checklistScan adds a sheet to the xlsx file
    createSubjTaskList(rootEEGdir);   % just pass in the root directory and this does the rest
    
    %- final run of checkListScan should come immediately after createSubjTaskList
    checklistScan(rootEEGdir,'writeScanSteps',0,'writeWarningSummary',1,'verifyChecklists',1,'deleteInvalid',1,'debug',0);
    % rootPath - path where NIHXXX subject folders are present
    % writeScanSteps - binary flag (default false) for writing list of each subject's changes (e.g., 47 items the same, 2 good -> warning, 3 warning --> good)
    % writeWarningSummary - binary flag (default TRUE) for writing out master warning list
    % verifyChecklists - binary flag (default FALSE) that will trigger a re-run of subjectDirCheck on each subject
    % deleteInvalid - default TRUE.  if TRUE and verifyingChecklist, delete any finalChecklist.invalid.csv and dont create new ones
    % debug - binary flag (default FALSE) all written files are output to debug directory instead of rootPath
   
    fprintf('\n now if moving to pristine, do a full directory check on permissions\n');
    %return;
    
end

if ALL_SUBJ_SORT_BUILD
    % Creates the sortBuild file, which summarizes per subject and session,
    % which sorts are complete or not
    % This will generate the file in the PUB directory and then copy it
    % over to .update/eeg
    rootUPDATEdir = rootEEGdir;
    status = createSubjSortList(rootPUBdir,rootUPDATEdir);
    fprintf('\n now if moving to pristine, do a full directory check on permissions\n');
end


if FIND_LOCKED_FILES,
    % shouldn't have to run this often... mostly important if chronosync is having issues making a clean copy
    searchAndDestroyLockedFiles(rootEEGdir)
    return;
end

