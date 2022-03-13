%makeNewestEvStructoldSubj
% converts newest event structure subjects 37 to 43? to better format:
 % New file format
%  eventscrs.subject = events.subject; *
%  eventscrs.expFolder = events.sessionName;*
%  eventscrs.amplitude = events.amplitude;*
%  eventscrs.frequency = events.frequency;*
%  eventscrs.isWN_REP = events.isWN_REP;*
%  eventscrs.eegfile = events.eegfile;*
%  eventscrs.cumulative_pulses = events.cumulative_pulses;*
%  eventscrs.dc_Pulses = events.num_dc_pulses;*
%  eventscrs.numPulses = events.numPulses;*
% eventscrs.session % 0, 1,2,3,4 (now changes if type or stimpair changes)*
% eventscrs.stimType % CC0,CC1, SM, WN,WNStack %
% eventscrs.monopolar %0 1 yes no. for old type instead 1 may be the value
% of difference... *
 % In old version. values are nans if these are bipolar sessions 0 if cc0; 1 if cc1...
 % should i just change to all nans... I think yes
 %eventscrs.stimLocTag % make for all! *
 %eventscrs.WNRepNum % nan. 0, 1,2,3, *
 %eventscrs.NKreference % as TAG *
 %eventscrs.annotation
 % for NIH040, NIH041, maybe NIH042 can go to their stim log 30-May-2016
 % and edit from there.
oldrootDir = '/Volumes/Shares/FRNU/people/tim/wn_events/';
funDir = '/Volumes/Shares/FRNU/people/steinhardt/';
labfunDir = '/Users/steinhardtcr/Documents/eeg_toolbox_trunk';
addpath(genpath(funDir));
addpath(genpath(labfunDir));
rootEEGdir =   '/Users/steinhardtcr/Documents/OrganizedStimAnalysis/newEvent3sGood';
%vs (['/Volumes/Shares/FRNU/dataWorking/eeg/' subj '/behavioral/stimMapping'])
for subjend = 40%%[37:39 41 43]
subj = ['NIH0' num2str(subjend)];
cd([rootEEGdir '/' subj '/behavioral/stimMapping'])
allSess= dir('sess*');
%cd([oldrootDir subj]);
%load('wn_events.mat'); % not complete and referencein all the events strutures i have but those aren't complete enough to od this...
% 
% eegfileDivs = (unique({events(:).eegfile}));
% 
%  q = cellfun(@(x) strsplit(x,'/'),eegfileDivs,'UniformOutput',0);
% timeSes= cellfun(@(x) str2num(x{end}(end-3:end)),q,'UniformOutput',0);
% [srtedtimes idxsrted]=sort([timeSes{:}]);

sesscnt = 0 ;
for sessSplitFiles = 1:length(allSess)
    cd([rootEEGdir '/' subj '/behavioral/stimMapping'])

cd(allSess(sessSplitFiles).name)
clear events eventscrs
%cd([oldrootDir subj]);
load('events.mat');
%% go "session" folder by "session" folder - more like blocks- and sort by time at which occurred.

%redo 49 events in the way i did today then do this and try the repetitions

 % New file format
 [eventscrs(1:length(events)).subject] =   deal(events(1).subject);
 [eventscrs.expFolder] =  deal(events.session);%events(1).expFolder);
 [eventscrs.stimLocTag] = deal(events.stimLocTag);

 %%
 [eventscrs.amplitude] =  deal(events.amplitude);
 [eventscrs.frequency] =  deal(events.frequency);
 [eventscrs.isWN_REP] =  deal(-1); %is all nanas in old subjects
 [eventscrs.wnRepNum] = deal(-1);
 [eventscrs.polarity] = deal(events.polarity); % polarity super confusin
% see 0,1,-1. 10
%%
% for old subjects stim pair or big offsets are going to be the way that we separate
% sessions
%events([events.eegoffset] == 0) = [];
offsetSizes=abs(diff([events(:).eegoffset]));
stckpulses = find(abs(diff([events.eegoffset])) < 100 & abs(diff([events.eegoffset])) >1 );
%%Will mark all stckulse pulses and the following as stack because they'd be stack corrupted
[events(stckpulses).stimType] = deal('WN_Stack');
[events(stckpulses+1).stimType] = deal('WN_Stack');
[eventscrs.stimType] = deal(events.stimType);
%sessSeparators= find(offsetSizes > mean(offsetSizes) & offsetSizes >2000); %start of next trial
 % for now use session breaks instead of offsets because of weird timing
 % glitches...
 %load('events')
 
UsessNames = unique({events(:).stimType});
rmdbrks =cellfun(@(x) strncmp(x,'WN_REP_D',8),UsessNames);
UsessNames([rmdbrks]) = [];
UsessNames=UsessNames([cellfun(@(x) length(x)>20,UsessNames)]) ;
    sepTimes= cellfun(@(x) max(find(strcmp({events(:).stimType},x))),UsessNames,'UniformOutput',0);
   sessSeparators = sort([sepTimes{:}]);
   sessSeparators(sessSeparators == 1)=[];
 
 clear sessIdxsfin
paircombs = unique({events.stimLocTag});
 if size(paircombs,1) ==1,
     cnt =1;
     [eventscrs([1:(sessSeparators(1)-1)]).session] = deal(cnt); 
     sessIdxsfin{cnt} = [1:sessSeparators(1)-1];  
     cnt = cnt +1;
    for sessDivs = 1:length(sessSeparators)-1
        for m = sessSeparators(sessDivs):(sessSeparators(sessDivs+1)-1)
             eventscrs(m).session = cnt;
        end
    sessIdxsfin{cnt} = [sessSeparators(sessDivs):(sessSeparators(sessDivs+1)-1)];
        cnt = cnt +1;
    end
    for m = sessSeparators(sessDivs+1):length(eventscrs)
    eventscrs(m).session = cnt;
    end
     sessIdxsfin{cnt} = [sessSeparators(sessDivs+1):length(events)];  
 else % if stim locations vary and we can only use offsets to tell things apart
     cnt =1;
     curevNum = 1;
     while curevNum <=length(events)
         if isempty(sessSeparators(min(find(sessSeparators > curevNum))))
             curSessEvs =[curevNum:length(events)];
         else
     curSessEvs =[curevNum:sessSeparators(min(find(sessSeparators > curevNum)))];
         end
     curStimSessEvs = find(strcmp({eventscrs(curSessEvs).stimLocTag},eventscrs(curevNum).stimLocTag)); 
     [eventscrs([curSessEvs(curStimSessEvs)]).session] = deal(cnt); 
     sessIdxsfin{cnt} = [curSessEvs(curStimSessEvs)];
     cnt =cnt +1;
     curevNum = curSessEvs(curStimSessEvs(end)) +1;
     end
 end
 %eventscrs(find([events.eegoffset] ==0))=[];
 [ eventscrs.eegoffset] = deal(events.eegoffset);
 eventscrs([events.eegoffset] == 0) = []
%  if min(abs(unique(offsetSizes))) < 100
%     % diffOffsets = unique(offsetSizes);
%      for n = unique([eventscrs(find(abs(diff([eventscrs.eegoffset])) < 100 & abs(diff([eventscrs.eegoffset])) > 1 )).session])
%          for stacksess= find([eventscrs.session] == n)
%         eventscrs(stacksess).stimType ='WN_Stack';   
%          end
%      end
%           sprintf(['we have stack pulses in ' subj ' session' num2str(n)]) 
%  end
 sessions = unique([eventscrs.session]);
if  sessions(1) ~=1
    for n = sessions-sessions(1) +1
        for m = find([eventscrs.session] == sessions(n))
            eventscrs(m).session = n;
        end
    end
end
%eventscrs([events.eegoffset] ==0) = [];
%
%evtest = eventscrs;
%eventscrs([events.eegoffset] == 0) = []
clear sessLength
 repcnt = 1;
 for n = 1:max([eventscrs(:).session])
  %  [eventscrs(find([eventscrs(:).session] == n)).session] =  deal(n-1);
    sessLength(n) = length(find(([eventscrs(:).session] == n)));
 end
 allSes= 1:max([eventscrs(:).session]);
 %[events([eventscrs.session] == possiblereps).amplitude]
 
 %eventscrs([events.eegoffset] ==0) = [];
 for testReps = 1:max([eventscrs(:).session])
     if unique([eventscrs([eventscrs.session] == testReps).wnRepNum] ) == -1
     choosePool = allSes;
     choosePool(testReps)= [];
    possiblereps = choosePool(find(sessLength(choosePool) == sessLength(testReps)));
    if sum(possiblereps)~=0 & ~isempty(possiblereps)
        repcnt= 2;
        for curTestSes =possiblereps 
            
          if sum( [eventscrs([eventscrs.session] == curTestSes).amplitude] -[eventscrs([eventscrs.session] == testReps).amplitude]) == 0
            for strtRep =[find([eventscrs(:).session] == testReps)]
              eventscrs(strtRep).wnRepNum =  1;
              eventscrs(strtRep).isWN_REP = 1;
            end
            for matchReps = find([eventscrs(:).session] == curTestSes)
            eventscrs(matchReps).wnRepNum =  repcnt;
            eventscrs(strtRep).isWN_REP = 1;
            
            end
            repcnt =  repcnt +1;
          end
        
        end
    else
        'did nothing'; % no matches but same length
    end
     else
         'nothing'; % nothing even the same length
     end 
 end
 %% get the eegfiel path correct

 events([events.eegoffset] == 0) = [];
%   tmpdr = strsplit(events(1).eegfile,'/');
%   if strcmp(tmpdr(5),'dataWorking') & ~isempty(find(strncmp('data',tmpdr,4)))
%   else
%  [events] = ev_part_switch(events,'data','dataWorking');
%   end
%  %  if strcmp(tmpdr{8}(5),'r')
%  % [events] = ev_part_switch(events,'eeg.reref','eeg.processed');
%  %  end
%  
% 
%      if strncmp(tmpdr(end),'NIH',3)
%          [events] = ev_part_switch(events,tmpdr(end),tmpdr{end}(length(subj)+2:end));
%        
%      else
        [eventscrs.eegfile] = deal(events.eegfile);
%     end

 %% less important pulse info
 [eventscrs.cumulative_pulses] = deal(events.cumulative_pulses);
 [eventscrs.dc_Pulses] = deal(events.num_dc_pulses);
 [eventscrs.numPulses] = deal(events.numPulses);
 [eventscrs.annotation] = deal(events.annotation);
 % these two are super necessary!
 [eventscrs.eegfile] = deal(events(1).eegfile);
 [eventscrs.eegoffset] = deal(events.eegoffset);
 
 cd( '/Volumes/Shares/FRNU/people/steinhardt/');
if ~exist('newEventStructures', 'dir')
    mkdir('newEventStructures')
end

clear events
events = eventscrs;
cd([funDir  '/newEventStructures'])
 save([subj '_B_' num2str(sesscnt) '_events'],'events')
 sesscnt =sesscnt +1;
end
end
 %% testing if the electrode names were referencing the right places... seems like yes...
%  curEvSet= find(ismember(vertcat(events(:).stimLocNum),allstimLocCombos(end-1,:),'rows'));
%      wnevents = find(strcmp({events(:).stimType},'WN'));
%      cursess = intersect(wnevents,curEvSet);
% 
% 
% % stename = {chan2name(subj,'/Volumes/Shares/FRNU/dataWorking/eeg',10)}; -
% % same as current cha2name
% 
% 
% curEvSet= find(ismember(vertcat(events(:).stimLocNum),allstimLocCombos(end-1,:),'rows'));
% wnevents = find(strcmp({events(:).stimType},'WN'));
% cursess = intersect(wnevents,curEvSet);
% curEvSet = cursess;
%  if isnan(events([curEvSet(1)]).eegoffset), curEvSet(1) = [] ;end
%  if isnan(events([curEvSet(end)]).eegoffset),curEvSet(end)=[];end
% 
% 
% Sesduration = events([curEvSet(end)]).eegoffset-events([curEvSet(1)]).eegoffset;
% %events = ev_part_switch(events,'data','dataWorking');
% %tmp = strsplit(events(1).eegfile,'/');
% [events.eegfile] = deal(strjoin(tmp(1:9),'/'));
% [eegdis]=gete(cell2mat(AllTAGS(30)),events([curEvSet(1)]),Sesduration);
% [eegdis2]=gete(cell2mat(AllTAGS(31)),events([curEvSet(1)]),Sesduration);
% [eegdis3]=gete(cell2mat(AllTAGS(allstimLocCombos(3,2))),events([curEvSet(1)]),Sesduration);
% figure; plot(eegdis,'r'); hold on;
% plot(eegdis2,'b');
% plot(eegdis3,'g');
% xlim([3000 6000])
% 
% for i = 1:length(paircombos)
%     i
% stimEnames{i} =[chan2name(subj,currootpath,num2str(allstimEs(i)))];
% end
% 
% if str2num(subj(5:6)) < 44
% stchn1 = name2chan(subj,currootpath,StimLocCombs{curstimSes}{1});
% stchn2= name2chan(subj,currootpath,StimLocCombs{curstimSes}{2});
% rndchk= name2chan(subj,currootpath,randChkElec);
% 
% else
%     if str2num(subj(5:6)) == 39 % might be for all of them..
%         curPairLocs{1}= StimLocCombs{curstimSes}{1};
%          curPairLocs{2}=StimLocCombs{curstimSes}{2}; 
%     end
%     
%  stchn1 =curPairLocs{1};
%  stchn2 =curPairLocs{2};
%  rndchk=chan2name(subj,currootpath,34);%randChkElec;
% end

%% make plotable and plot! (Reps to stim)
% 
% if str2num(subj(5:6)) > 43
%     stimprnums = curPairLocs;
% else
%     stimprnums = StimLocCombs{1};
% end
% if isnan(cursessEv(1).eegoffset) | isnan(cursessEv(end).eegoffset)
% pulllength1 =  cursessEv(end-1).eegoffset -  cursessEv(2).eegoffset;
% else
% pulllength1 =  cursessEv(end).eegoffset -  cursessEv(1).eegoffset;  
% end
% [eegstim1, ~] = gete_ms(stchn1, cursessEv(2),pulllength1);
% [eegstim2, ~] = gete_ms(stchn2, cursessEv(2),pulllength1);
% [eego1, ~] = gete_ms(rndchk, cursessEv(2),pulllength1);
% 
% plotEEGOut(eegstim1,eegstim2, eego1,subj,stimprnums,rndchk,[2000 12000],': WN MonoPolar Session')


