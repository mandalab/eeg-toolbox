
function [events] = Color_STM_Extract_Event(behDir,sessionName)
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
datafilename = sprintf('%s_*.mat',dataname);
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
datawhole={};
DataAll = {};

%subID = unique(subjectID(~ismember(subjectID,ExcludeSub)));
for isub = 1:nsub
    
    for ise = 1:length(SessionID)
        
        for iblk = 1:nblock
            
            loadingfile = [datapath sprintf('/Record_S%d_ID_sn%d_b%d.mat',subID(isub),SessionID(ise),BlockID(iblk))];
            load(loadingfile);
%             fprintf('loading:%s\n',loadingfile);
            ntrail = size(Datastruct);
            ntrail = ntrail(1,2);
            
            for i = 1:ntrail
                datawhole{i+(iblk-1)*ntrail,isub} = Datastruct(i);
            end
            
%             fprintf('loading:session %d block%d\n',ise,iblk);
        end
    end
    
    
    for ise = 1:length(SessionID)
        
        
        for iblock = 1:nblock
            for itr = 1:ntrail
                %Get general info
                block(itr+(iblock-1)*ntrail,isub) = datawhole{itr+(iblock-1)*ntrail,isub}.block; % change no change, 1 no change, 0 for change
                trial(itr+(iblock-1)*ntrail,isub) = datawhole{itr+(iblock-1)*ntrail,isub}.trial; % change no change, 1 no change, 0 for change
       
                %Get the independent variable
                TempMemPresent(:,itr+(iblock-1)*ntrail)= datawhole{itr+(iblock-1)*ntrail,isub}.Mem; %  1 First Gen, 2 recent Gen                
                ToRem(itr+(iblock-1)*ntrail,isub) = TempMemPresent(1,itr+(iblock-1)*ntrail); % 
                ToForget(itr+(iblock-1)*ntrail,isub) = TempMemPresent(2,itr+(iblock-1)*ntrail,isub); %  
                Radius(itr+(iblock-1)*ntrail,isub) = datawhole{itr+(iblock-1)*ntrail,isub}.Radius; %  1 First Gen, 2 recent Gen
                %Get the dependent variable
                TestItem(itr+(iblock-1)*ntrail,isub)= datawhole{itr+(iblock-1)*ntrail,isub}.TestItem;                
                ClutRotateDeg(itr+(iblock-1)*ntrail,isub)= datawhole{itr+(iblock-1)*ntrail,isub}.ClutRotateDeg;
                Error(itr+(iblock-1)*ntrail,isub) = datawhole{itr+(iblock-1)*ntrail,isub}.Error; %  1 First Gen, 2 recent Gen
                RT(itr+(iblock-1)*ntrail,isub) = datawhole{itr+(iblock-1)*ntrail,isub}.RT;
                Study(itr+(iblock-1)*ntrail,isub)= datawhole{itr+(iblock-1)*ntrail,isub}.Study;               
                MemRot(itr+(iblock-1)*ntrail,isub)= datawhole{itr+(iblock-1)*ntrail,isub}.MemRot(1);               
                estimate(itr+(iblock-1)*ntrail,isub)= datawhole{itr+(iblock-1)*ntrail,isub}.estimate;
                Events{itr+(iblock-1)*ntrail,isub}= datawhole{itr+(iblock-1)*ntrail,isub}.Events;
            end
        end
        
        DataAll{isub}.Events = Events;
        DataAll{isub}.block = block(1:ntrail*nblock,isub);
        DataAll{isub}.trial = trial(1:ntrail*nblock,isub);
        DataAll{isub}.ToRem = ToRem(1:ntrail*nblock,isub); % this is the hiden true color
        DataAll{isub}.Study = Study(1:ntrail*nblock,isub); % this color is rotated by a random jitter
        DataAll{isub}.ToForget = ToForget(1:ntrail*nblock,isub);
        DataAll{isub}.TestItem = TestItem(1:ntrail*nblock,isub);
        DataAll{isub}.TempMemPresent=TempMemPresent(:,1:ntrail*nblock); 
        DataAll{isub}.MemRot = MemRot(1:ntrail*nblock,isub); % this color is rotated true color, based on ClutRotateDeg
        DataAll{isub}.ClutRotateDeg = ClutRotateDeg(1:ntrail*nblock,isub);
        DataAll{isub}.Error = Error(1:ntrail*nblock,isub);
        DataAll{isub}.RT = RT(1:ntrail*nblock,isub);
        DataAll{isub}.estimate = estimate(1:ntrail*nblock,isub); % error = estimate - MemRot 
        DataAll{isub}.Radius = Radius(1:ntrail*nblock,isub);
        DataAll{isub}.subID(isub) = subID(isub);
        DataAll{isub}.SessionID(ise) = SessionID(ise);

        DataSub = DataAll(isub);
        
        save([datapath '/DataAll_' num2str(subID(isub)) '_Sn' num2str(SessionID(ise)) '.mat'], 'DataSub');
        
        %create events.mat file 
        clear events;
        
        curEvent = 1;
        for itr = 1:size(Events,1)
            for ie = 1:size(Events{1},2)
                
            clear thisEvent
            thisEvent.block = DataAll{isub}.block(itr);
            thisEvent.trial = DataAll{isub}.trial(itr);
            thisEvent.experiment = 'Color_STM';
            thisEvent.subject = DataAll{isub}.subID(isub);        
            thisEvent.session = DataAll{isub}.SessionID(ise);
            thisEvent.eventType = DataAll{isub}.Events{itr}(ie).name;
            thisEvent.TestItem = DataAll{isub}.TestItem(itr);
            thisEvent.Radius = DataAll{isub}.Radius(itr);
            thisEvent.TempMemPresent = DataAll{isub}.TempMemPresent(:,itr);
            thisEvent.ClutRotateDeg = DataAll{isub}.ClutRotateDeg(itr);
            thisEvent.ToRem =  DataAll{isub}.ToRem(itr);
            thisEvent.ToForget = DataAll{isub}.ToForget(itr);
            thisEvent.Study = DataAll{isub}.Study(itr);
            thisEvent.estimate = DataAll{isub}.estimate(itr);
            thisEvent.MemRot =  DataAll{isub}.MemRot(itr);
            thisEvent.RT = DataAll{isub}.RT(itr);
            thisEvent.Error = DataAll{isub}.Error(itr);            
            thisEvent.mstime = DataAll{isub}.Events{itr}(ie).unixStartTime;

            events(curEvent) = thisEvent;
            curEvent = curEvent+1;
            end
        end
        
        save([datapath '/events.mat'], 'events');
        
    end
    
end

fprintf('\nloading behavioral task finished\n\n');

end
