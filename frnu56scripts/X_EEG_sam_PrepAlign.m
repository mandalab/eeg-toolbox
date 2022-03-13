%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to run PrepAndAlign for multiple (or 1) patients as a check after processing, or when implementing a change.

subject_list = {'NIH088'};

rootEEGdir = '/Volumes/56PROC/eeg_processing/.update/eeg';
%rootEEGdir = '/Volumes/EEG_56STAGE-1/stage/working';
%rootEEGdir = '/Volumes/EEG_56STAGE-1/stage/ready';
%rootEEGdir = '/Volumes/EEG_56STAGE-1/stage/grabbed';

just_behavioralProcessing   = 0; taskType = 'stm_face'; % 

freshJackandSplit           = 1; %YES
freshExtractAlign           = 1; %YES
batch                       = 0;
skipProcessReref            = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick checks to alert the user about rootEEGDir
dir_correct = input(['rootEEGdir is set to: ' rootEEGdir '. Is this correct? (Press any key to Continue) ']);

if contains(rootEEGdir,'56PROC')
    fprintf('%s\n',['Caution: rootEEGdir is set to ' rootEEGdir ', which is on 56PROC. Are you sure you mean to do this?']);
    keyboard
elseif contains(rootEEGdir,'stage')
    if contains(rootEEGdir, 'working')
        fprintf('%s\n',['Caution: rootEEGdir is set to ' rootEEGdir ', which is in the working folder! Make sure you have another copy or switch to ready!']);
        keyboard
    end
else
    fprintf('%s\n','Caution: rootEEGdir is not set to 56PROC or STAGE. Please confirm you are in the right spot before continuing.');
    keyboard
end

if ~exist(rootEEGdir,'dir')
    fprintf('%s\n',['ERROR! rootEEGdir is set to ' rootEEGdir ', which does not exist. Are you sure you are connected to the Volume?']);
end

% Looping through subjects
for sub = 1:numel(subject_list)
    subject = subject_list{sub};
    
    if ~exist(fullfile(rootEEGdir,subject),'dir')
        fprintf('%s\n',['ERROR! ' fullfile(rootEEGdir,subject) ' does not exist!!!']);
        keyboard
    end
    
    if just_behavioralProcessing == 1
        input('User selected JUST behavioral processing. Is this correct? (Press any key to Continue) ');
        behavioralProcessing(subject,rootEEGdir,taskType);
        fprintf('\n%s\n\n',['Behavioral Processing complete for ' subject ': ' taskType]);
    else
        eegPrepAndAlign(subject,rootEEGdir,'freshJackandSplit',freshJackandSplit,'freshExtractAlign',freshExtractAlign,'batch',batch,'skipProcessReref',skipProcessReref);
    end
          
    fprintf('%s\n','----------------------------------------------------------------------------------------------------');
    fprintf('%s\n',['-------------------------------------- ' subject ' Complete [' num2str(sub) '/' num2str(numel(subject_list)) '] ---------------------------------------']);
    fprintf('%s\n','----------------------------------------------------------------------------------------------------');
end
    
fprintf('\n\n%s\n','----------------------------------------------------------------------------------------------------');
fprintf('\t\t%s\n\n','RUN COMPLETE');


