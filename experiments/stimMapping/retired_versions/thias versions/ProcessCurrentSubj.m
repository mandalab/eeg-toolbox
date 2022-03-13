%% Initialize to run functions
addpath(genpath('/Users/steinhardtcr/Documents/eeg_toolbox_trunk'));
addpath(genpath('/Users/steinhardtcr/Documents/dungeon_toolbox_17a'));
%rootEEGdir  = '/Volumes/Macintosh HD 2/';                        %office-local
%rootEEGdir = '/Users/steinhardtcr/Documents/offlinefiles/';%;'/Volumes/Shares/FRNU/dataWorking/eeg';
%%rootEEGdir ='/Users/steinhardtcr/Documents/StimProcessingTesting/ExtractionPackage';
rootEEGdir = '/Users/steinhardtcr/Documents/OrganizedStimAnalysis/newEvent3sGood/';%'/Volumes/Shares/FRNU/dataWorking/eeg';                %office-server
addpath(genpath('Volumes/Shares/EEG'));
subj       = 'NIH045';
%prefixMachine = 'CA211' ;  % 'HA866' 'CA211'
prefixMachine = 'CA211';

strDateStart = '7/10/16'; 
strDateEnd = '7/12/16';  


%% ------- LIST THE CONTENTS OF INATI's SERVER
nktdir = '/Volumes/Shares/EEG/LTVEEG_DATA/Archive3/NKT/EEG2100';
eegRawServerDir(subj, rootEEGdir, prefixMachine, strDateStart, strDateEnd,nktdir);
%% ------- GRAB SPECIFIC FILES FROM INATI's SERVER
grabSessions = {'2NE'};  %212','216','21I','21S','21Y'};
%nih050 - CA2114S7.EEG 
%nih051 -  CA2114RL.EEG , CA2114RN.EEG, CA2114RP.EEG , 
eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, grabSessions, {'21E','EEG','LOG'},nktdir,1,0);% 1 means stim folder version 
%% PLOT THE DC09 CHANNEL FROM EACH
rawPath = fullfileEEG(rootEEGdir,subj,'raw');  eegPulseVisualize(rawPath, {'DC10'});   % STIM {'21E','EEG','LOG'}

%% ------- RAM files
% grabSessions = {'1MJ',};
% eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, grabSessions, {'21E','EEG'},nktdir);  
% %% plot DC09
% rawPath = fullfileEEG(rootEEGdir,subj,'raw/_RAM');  eegPulseVisualize(rawPath, {'DC09'});   % STIM {'21E','EEG','LOG'}
% 
% 
% %% grab stim Sessions
% grabSessions = {'1T7','1T8'}; %3DZ 3E1 3A7
% eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, grabSessions, {'21E','EEG','Log'});  
% %% Just Log
% grabSessions = {'160'}; %3DZ 3E1 3A7
% eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, grabSessions, {'LOG'});  
% %% visualize
% rawPath = fullfileEEG(rootEEGdir,subj,'raw/STIM');
% eegPulseVisualize(rawPath, {'DC10'}); 
% 
% %% Behavioral Processing
% taskList = {'attentionTask', 'palRam', 'palRamStim', 'paRemap',  'stimMapping','stimMapping2',  'SequenceMem',  'playPass', 'auditoryLexicalDecision', 'auditoryVowels', 'goAntiGo', 'pa3', 'languageTask'};
% for iTask = [5],
%     behavioralProcessing(subj, rootEEGdir, taskList{iTask});  % attentionTask languageTask moveTask pa3 paRepeat playPass
% end   
% 

%% EEG prep & align
subj= 'NIH045';

tmpeegPrepAndAlign(subj, rootEEGdir,'useStim',0);

eegPrepAndAlign(subj, rootEEGdir,'useStim',1);

%% Create leads.txt
% fid = fopen('leads.txt','w');
% for k = 1:95,
%     fprintf(fid,'%d\n',k);
% end




