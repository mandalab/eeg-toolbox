
function [events] = stm_color_stimulation_eventExtract(behDir,sessionName)
% analyze the data
% Zane 2014.2.11

% clear all;
% close all;
% clc;

% set the data path as the current data directory
datapath = [behDir filesep sessionName];


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
    for ifile = 1:nfile
        analysis_file = groupfile(ifile).name;
        infname = analysis_file;
        readingfile = fullfile(datapath, infname);
        scanfile = sscanf(infname,'Record_S%d_ID_sn%d_b%d.mat');
        subjectID(ifile)= scanfile(1,1);
        Session(ifile) = floor(scanfile(2,1));
        BlockID(ifile) = scanfile(3,1);
    end


subID = unique(subjectID(~ismember(subjectID,ExcludeSub)));
BlockID = unique(BlockID);
SessionID = unique(Session(~ismember(Session,ExcludeSub)));

%update nsub, nblock
nsub = length(subID);
nblock = length(BlockID);
% fprintf('we have %d subjects in this analysis\nEach has %d blocks\n',nsub, nblock);

%% We are going to load the file, subject by subject with all conditions


%subID = unique(subjectID(~ismember(subjectID,ExcludeSub)));
for isub = 1:nsub
    

    AllstructToEvent=[];
    
    for iblk = 1:nblock
        clear Datastruct_Study  Datastruct_Practice  Datastruct_Test
        loadingfile = [datapath sprintf('/Record_S%d_ID_sn%d_b%d.mat',subID(isub),SessionID(isub),BlockID(iblk))];
        load(loadingfile);
        %             fprintf('loading:%s\n',loadingfile);
        AllstructToEvent = [AllstructToEvent Datastruct];
        %             fprintf('loading:session %d block%d\n',ise,iblk);
    end
    

    Events={AllstructToEvent.Events};
    
    clear events
     
    curEvent = 1;
        for itr = 1:length(Events)
            for ie = 1:size(Events{itr},2)
                
            clear thisEvent
            thisEvent.block =AllstructToEvent(itr).block;
            thisEvent.trial = AllstructToEvent(itr).trial;
            thisEvent.experiment = 'stm_color_stimulation';
            thisEvent.subject = subID(isub);        
            thisEvent.session = SessionID;
            thisEvent.eventType = AllstructToEvent(itr).Events(ie).name;
            thisEvent.mstime = AllstructToEvent(itr).Events(ie).unixStartTime;
            
            %variables
            thisEvent.Stimulation = AllstructToEvent(itr).Stimulation;
            thisEvent.M1Size = AllstructToEvent(itr).M1Size;
            thisEvent.MemRot = AllstructToEvent(itr).MemRot;
            thisEvent.estimate = AllstructToEvent(itr).estimate;
            thisEvent.Error = AllstructToEvent(itr).Error;
            thisEvent.RT = AllstructToEvent(itr).RT;           
            thisEvent.Mem = AllstructToEvent(itr).Mem;
            thisEvent.Test = AllstructToEvent(itr).Test;
            
            events(curEvent) = thisEvent;
            curEvent = curEvent+1;
            end
        end   
    
        save([datapath '/events.mat'], 'events');
        
end
    
 fprintf('\nloading behavioral task finished\n\n');
   

end
