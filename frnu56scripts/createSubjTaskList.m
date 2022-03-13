function createSubjTaskList(rootEEGdir)
%
%%  function createSubjTaskList
%
%  This is a helper function for updating the eeg database on frnu.  last step before pushing an update is to update an xls that says which subjects did which tasks
%  Loop over subjects and create table of completed tasks, saved in an xls folder at the root
%
%
%  JW 4/2019 -- commited to SVN (created about a year before that but maintained locally)
%
%  SJ 12/16/2019 -- Change subjBuild from .xls to .xlsx
% 

%clear all;

%rootEEGdir = '/Volumes/Shares/FRNU/data/eeg';                        %office-local
%rootEEGdir = '/Volumes/JW24TB/data24TB/eeg_new';                        %office-local

nameXLfile = fullfileEEG(rootEEGdir,sprintf('subjBuild_%s.xlsx',datestr(now,'yymmdd'))); %SJ change to .xlsx
%nameXLfile = fullfileEEG(rootEEGdir,'subjTaskSummary.xls');
XLsheet    = 'NUM_SESS_FOLDERS';

allTasksFound   = {};
sessCntAllTasks = [];




%- loop over subjects

subjList = dir(fullfileEEG(rootEEGdir,'*')); %- catch TRE and NIH
subjListCln = {}; for iT=1:length(subjList), if ~strcmp(subjList(iT).name(1),'.') & subjList(iT).isdir, subjListCln{end+1} = subjList(iT).name; end; end; subjList = subjListCln;
xlOUT = {};
for subj = subjList,
    
    subj = subj{1};
    iSubj=find(strcmp(subjListCln,subj));
    
    
    xlOUT{1+iSubj,1} = subj;
    
    
    %- loop over task folders and check events for each task
    prePostStr = '';
    taskListDir = dir(fullfileEEG(rootEEGdir,subj,sprintf('behavioral%s',prePostStr)));
    taskListCln = {}; for iT=1:length(taskListDir), if ~strcmp(taskListDir(iT).name(1),'.') & ~strcmp(taskListDir(iT).name(1),'_') & taskListDir(iT).isdir, taskListCln{end+1} = sprintf('%s%s',taskListDir(iT).name,prePostStr); end; end
    taskListFound = taskListCln;
    
    %- check for preOp/postOp version
    prePostStr = '_preOp';
    taskListDir = dir(fullfileEEG(rootEEGdir,subj,sprintf('behavioral%s',prePostStr)));
    taskListCln = {}; for iT=1:length(taskListDir), if ~strcmp(taskListDir(iT).name(1),'.') & ~strcmp(taskListDir(iT).name(1),'_') & taskListDir(iT).isdir, taskListCln{end+1} = sprintf('%s%s',taskListDir(iT).name,prePostStr); end; end
    taskListFoundPreOp = taskListCln;
   
    prePostStr = '_postOp';
    taskListDir = dir(fullfileEEG(rootEEGdir,subj,sprintf('behavioral%s',prePostStr)));
    taskListDir = dir(fullfileEEG(rootEEGdir,subj,'behavioral_postOp'));
    taskListCln = {}; for iT=1:length(taskListDir), if ~strcmp(taskListDir(iT).name(1),'.') & ~strcmp(taskListDir(iT).name(1),'_') & taskListDir(iT).isdir, taskListCln{end+1} = sprintf('%s%s',taskListDir(iT).name,prePostStr); end; end
    taskListFoundPostOp = taskListCln;
   
    subjTasksFound = unique({taskListFound{:} taskListFoundPreOp{:} taskListFoundPostOp{:}},'stable');
    %if length(allTasksFound)==0, allTasksFound = taskListFound; end
    allTasksFound = unique({allTasksFound{:} subjTasksFound{:}},'stable');
    %disp(allTasksFound);
    
    
    %- loop over found "names", which includes preOp and postOp suffix, then modify so the correct directory is searched
    numSessFolder = zeros(1,length(allTasksFound));
    for iTask = 1:length(subjTasksFound),
        
        thisTaskName = subjTasksFound{iTask};
        if    contains(thisTaskName,'_preOp'),
            thisBehDir = 'behavioral_preOp';
            thisTaskDir = thisTaskName(1:strfind(thisTaskName,'_preOp')-1);
        elseif contains(thisTaskName,'_postOp'),
            thisBehDir = 'behavioral_postOp';
            thisTaskDir = thisTaskName(1:strfind(thisTaskName,'_postOp')-1);
        else
            %- usual situation
            thisBehDir = 'behavioral';
            thisTaskDir = thisTaskName;
        end
        sessFolder = dir(fullfileEEG(rootEEGdir,subj,thisBehDir,thisTaskDir,'session_*'));
        %keyboard;
        
        %- save to the meta structure
        iAllTask = find(strcmp(allTasksFound,thisTaskName));
        sessCntAllTasks(iSubj,iAllTask) = length(sessFolder);
        xlOUT{1,1+iAllTask}             = allTasksFound{iAllTask};
        xlOUT{1+iSubj,1+iAllTask}       = length(sessFolder);
    end
    
    
    %- use this code to get different information about each task... e.g., num trials, numEEG files, etc
    if 0,
        numSess       = zeros(1,length(allTasksFound));
        for iTask = 1:length(taskListFound),
            evFile = fullfileEEG(rootEEGdir,subj,'behavioral',taskListFound(iTask),'events.mat');
            if exist(evFile,'file'),
                events = [];
                load(evFile);
                
                disp(evFile(1));
                keyboard;
                
                %grab stats about each task
                % numSessions (non-integer to represent number of completed sessions)
                % numTrials
                % numEEG files
                % numEvents
            end
        end% for iTask
    end% if 0
    
end


%- sort based on task names
[Y,I] = sort(upper(allTasksFound));
allTasksFound   = allTasksFound(I);
sessCntAllTasks = sessCntAllTasks(:,I);
subjPerTask     = sum(sessCntAllTasks>0);

xlOUT2 = xlOUT;
for iTask=1:length(allTasksFound),
    for row=1:length(subjList)+1,
        xlOUT{row,1+iTask} = xlOUT2{row,1+I(iTask)};
        if isempty(xlOUT{row,1+iTask}), 
            xlOUT{row,1+iTask} = ' '; 
        end
    end
    
    %- add summary at the bottom
    
    xlOUT{length(subjList)+3,1} = 'SUBJ>0';
    xlOUT{length(subjList)+3,1+iTask} = subjPerTask(iTask);
end


%
if exist(nameXLfile,'file'),
    [STATUS,SHEETS] = xlsfinfo(nameXLfile);
    if contains(SHEETS,XLsheet) & 0,
        fprintf('\n WARNING: summary excel file %s already contains SHEET %s; \n  best to delete xls file and create fresh',nameXLfile,XLsheet);
        keyboard;
    end
   delete(nameXLfile);
   fprintf('\n deleting existing summary: %s',nameXLfile);
end

status = xlwrite(nameXLfile, xlOUT, XLsheet, 'B2');
if status==0,
    fprintf('\n ERROR WRITING');
    keyboard;
end

%- procesing date
dateTxt = {sprintf('Created %s',datestr(now))};
status = xlwrite(nameXLfile, dateTxt, XLsheet, 'A1:A1'); 

fprintf('\n %s created with %d subjects and %d tasks\n',nameXLfile, length(subjList), length(allTasksFound));

