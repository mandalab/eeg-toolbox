
function [events] = ltm_Strupt_Extract_Event(behDir,sessionName)
% analyze the data
% Zane 2014.2.11

% clear all;
% close all;
% clc;

% set the data path as the current data directory
datapath = [behDir filesep sessionName];

if ~strcmp(sessionName,'session_0')

% cd(datapath)
ExcludeSub = [];
% subject 1 2 3 study 150 objects per condition, and therefore guessing rat
% is very high

% subject 9 and 15 have a guessing rate greater than 80%
%ExcludeSub = [1 3 4 5 6  10 11 13];
%ExcludeSub = [1:40];
dataname = ('Record');
datafilename = sprintf('*%s_*.mat',dataname);
groupfile = dir(fullfile(datapath, datafilename));
nfile = length(groupfile);
% fprintf('%s\n',datapath);
% fprintf('There are %d files in the current directory\n\n', nfile);

%% search from the ls
if strcmp(sessionName,'session_1')
    for ifile = 1:nfile
        analysis_file = groupfile(ifile).name;
        infname = analysis_file;
        readingfile = fullfile(datapath, infname);
        scanfile = sscanf(infname,'Record_S%d_ID_sn%d_b%d.mat');
        subjectID(ifile)= scanfile(1,1);
        Session(ifile) = floor(scanfile(2,1));
        BlockID(ifile) = scanfile(3,1);
    end
else  % for later test
        for ifile = 1:nfile
        analysis_file = groupfile(ifile).name;
        infname = analysis_file;
        readingfile = fullfile(datapath, infname);
        scanfile = sscanf(infname,'DelayTest_Record_S%d_ID_sn%d_b%d.mat');
        subjectID(ifile)= scanfile(1,1);
        Session(ifile) = floor(scanfile(2,1));
        BlockID(ifile) = scanfile(3,1);
        end
end

subID = unique(subjectID(~ismember(subjectID,ExcludeSub)));
BlockID = unique(BlockID);
SessionID = unique(Session(~ismember(Session,ExcludeSub)));

%update nsub, nblock
nsub = length(subID);
nblock = length(BlockID);
% fprintf('we have %d subjects in this analysis\nEach has %d blocks\n',nsub, nblock);

%% We are going to load the file, subject by subject with all conditions
datawhole={};
DataAll = {};

%subID = unique(subjectID(~ismember(subjectID,ExcludeSub)));
for isub = 1:nsub
    
    AllStudyStruct=[];
    AllPractStruct=[];
    AllTestStruct=[];
    AllDelayTestStruct=[];
    
    if strcmp(sessionName,'session_1')
        for iblk = 1:nblock
            clear Datastruct_Study  Datastruct_Practice  Datastruct_Test
            loadingfile = [datapath sprintf('/Record_S%d_ID_sn%d_b%d.mat',subID(isub),SessionID(isub),BlockID(iblk))];
            load(loadingfile);
            %             fprintf('loading:%s\n',loadingfile);
            AllStudyStruct = [AllStudyStruct Datastruct_Study];
            AllPractStruct = [AllPractStruct Datastruct_Practice];
            AllTestStruct = [AllTestStruct Datastruct_Test];
            %             fprintf('loading:session %d block%d\n',ise,iblk);
        end
    else
        for iblk = 1:nblock
            clear Datastruct_Study  Datastruct_Practice  Datastruct_Test
            loadingfile = [datapath sprintf('/DelayTest_Record_S%d_ID_sn%d_b%d.mat',subID(isub),SessionID(isub),BlockID(iblk))];
            load(loadingfile);
            %             fprintf('loading:%s\n',loadingfile);
            AllDelayTestStruct = [AllDelayTestStruct Datastruct_DelayTest];
            %             fprintf('loading:session %d block%d\n',ise,iblk);
        end
    end
    
    if isempty(AllStudyStruct) StudyField=[];else StudyField= fieldnames(AllStudyStruct); end
    if isempty(AllPractStruct) PractField=[];else PractField= fieldnames(AllPractStruct); end
    if isempty(AllTestStruct) TestField=[];else TestField= fieldnames(AllTestStruct); end
    if isempty(AllDelayTestStruct) DelaytestField=[];else DelaytestField= fieldnames(AllDelayTestStruct); end
    
    if strcmp(sessionName,'session_1')
        overallapingfields=StudyField(ismember(StudyField,PractField));
        overallapingfields=overallapingfields(ismember(overallapingfields,TestField));   
        for itr=1:size(AllPractStruct,2)
            for i = 1:length(overallapingfields)
                AllStudyStruct1(itr).(overallapingfields{i})  = getfield(AllStudyStruct(itr),overallapingfields{i});
                AllPractStruct1(itr).(overallapingfields{i})  = getfield(AllPractStruct(itr),overallapingfields{i});
                AllTestStruct1(itr).(overallapingfields{i})  = getfield(AllTestStruct(itr),overallapingfields{i});
            end
        end
        AllstructToEvent=[AllStudyStruct1 AllPractStruct1 AllTestStruct1 ];
    else
        overallapingfields=DelaytestField;
        AllstructToEvent=[AllStudyStruct AllPractStruct AllTestStruct AllDelayTestStruct]; 
    end

    
    Events={AllstructToEvent.Events};
    
    clear events
        curEvent = 1;
        for itr = 1:length(Events)
            for ie = 1:size(Events{itr},2)
                
            clear thisEvent
            thisEvent.block =AllstructToEvent(itr).block;
            thisEvent.trial = AllstructToEvent(itr).trial;
            thisEvent.experiment = 'ltm_Strupt';
            thisEvent.subject = subID(isub);        
            thisEvent.session = SessionID;
            thisEvent.eventType = AllstructToEvent(itr).Events(ie).name;
            thisEvent.mstime = AllstructToEvent(itr).Events(ie).unixStartTime;
            
            %variables
            thisEvent.PractCond = AllstructToEvent(itr).PractCond;
            thisEvent.IsSame = AllstructToEvent(itr).IsSame;
  
            % memory
            thisEvent.Resp = AllstructToEvent(itr).Resp;
            thisEvent.ConfidenceRating = AllstructToEvent(itr).ConfidenceRating;           
            thisEvent.Corr_Memory = AllstructToEvent(itr).Corr_Memory;
            thisEvent.RT_Memory = AllstructToEvent(itr).RT_Mem;
            
            % categorization task
            thisEvent.Corr_Category1 = AllstructToEvent(itr).Corr_Category1;
            thisEvent.RT_Category1 = AllstructToEvent(itr).RT_Category1;
            thisEvent.Corr_Category2 = AllstructToEvent(itr).Corr_Category2;
            thisEvent.RT_Category2 = AllstructToEvent(itr).RT_Category2;
          
            thisEvent.CorrectMemoItem = AllstructToEvent(itr).CorrectMemoItem;
            thisEvent.MemNonTestedItem = AllstructToEvent(itr).MemNonTestedItem;
            thisEvent.PracticeItemPresent = AllstructToEvent(itr).PracticeItemPresent;

            events(curEvent) = thisEvent;
            curEvent = curEvent+1;
            end
        end   
    
  
        save([datapath '/events.mat'], 'events');
        
end
    
 fprintf('\nloading behavioral task finished\n\n');
   
else
    events=[]; % do not process the practice block
end

end
