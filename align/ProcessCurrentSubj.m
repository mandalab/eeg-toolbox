%% Initialize to run functions
addpath(genpath('/Users/steinhardtcr/Documents/eeg_toolbox/tags'));

%rootEEGdir  = '/Volumes/Macintosh HD 2/';                        %office-local
rootEEGdir = '/Volumes/Shares/FRNU/dataWorking/eeg';
%%rootEEGdir ='/Users/steinhardtcr/Documents/StimProcessingTesting/ExtractionPackage';
%rootEEGdir = '/Volumes/Shares/FRNU/data/eeg';                %office-server
addpath(genpath('Volumes/Shares/EEG'));
addpath(genpath('/Users/steinhardtcr/Documents/eeg_toolbox_trunk'));
subj       = 'NIH049';
%prefixMachine = 'CA211' ;  % 'HA866' 'CA211'
prefixMachine = 'HA866' ;

strDateStart = '5/24/17'; 
strDateEnd = '5/27/17';  


%% ------- LIST THE CONTENTS OF INATI's SERVER
nktdir = '/Volumes/Shares/EEG/LTVEEG_DATA/Archive5/NKT/EEG2100';
eegRawServerDir(subj, rootEEGdir, prefixMachine, strDateStart, strDateEnd,nktdir);
%% ------- GRAB SPECIFIC FILES FROM INATI's SERVER
grabSessions = {'4A5'}'%{'49X','4A5','49Z','4A7'};

eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, grabSessions, {'21E','EEG','LOG'},nktdir,1);% 1 means stim folder version 
%% PLOT THE DC09 CHANNEL FROM EACH
rawPath = fullfileEEG(rootEEGdir,subj,'raw/STIM');  eegPulseVisualize(rawPath, {'DC10'});   % STIM {'21E','EEG','LOG'}

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
eegPrepAndAlign(subj, rootEEGdir,'useStim',1);

%% Create leads.txt
% fid = fopen('leads.txt','w');
% for k = 1:95,
%     fprintf(fid,'%d\n',k);
% end




