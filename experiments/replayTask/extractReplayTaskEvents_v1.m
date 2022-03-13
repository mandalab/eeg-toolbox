function [events] = extractReplayTaskEvents_v1(rootEEGdir, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Function for extracting behavioral data for replayTask %%%%%%%%%%%%% 
%   
% Created by Jiali Zhang 12/2020
% Inputs:
%   rootEEGdir: eeg directory path 
%   subject: e.g. 'NIHXXX'
%   sessionName: e.g. 'session_0'

% % COMMENT OUT BEFORE PUTTING IN BEHAVIORALPROCESSING:
% rootEEGdir  = '/Users/zhangj34/Documents/experiments/replayTask/trunk/data';
% % sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/attentionFeedback', sessionName);
% sessionDir = '/Users/zhangj34/Desktop/replayTask/NIH088/session_1';
% subject     = 'NIH088';   % NIHXXX
% sessionName = 'session_1'; 

% Version 1: NIH084(alternative copy), NIH086 

%% set up pathways

%%% define %%%
sessionDir = fullfileEEG(rootEEGdir, subject, 'behavioral','replayTask',sessionName);
sessLogFile = fullfileEEG(sessionDir,'session.log');

%%% load CSV files containing image labels and paths %%%
% extract by category 
animalsCSV = readtableSafe(fullfile(sessionDir,'animals.csv'),'ReadVariableNames',false);
facesCSV   = readtableSafe(fullfile(sessionDir,'famousfaces.csv'),'ReadVariableNames',false);
foodsCSV   = readtableSafe(fullfile(sessionDir,'foods.csv'),'ReadVariableNames',false);
objectsCSV = readtableSafe(fullfile(sessionDir,'objects.csv'),'ReadVariableNames',false);
% combine 
fullCSV = [animalsCSV ; facesCSV ; foodsCSV; objectsCSV];

%% open and read files 

%%% open and read session log %%% 

% open 
fid = fopen(sessLogFile,'r'); % if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); % - used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

% convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find(sessionName >= '0' & sessionName <= '9'); % gets index/indices of numbers in sessionName 
if max(diff(strNumeric))>1, fprintf('\n ERROR: problem converting session name into a numeric.. %s --> %s', sessionName, sessionName(strNumeric)); keyboard; end
sessNum = str2double(sessionName(strNumeric)); % store indices of sessionName             
if isempty(sessNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end % shouldn't need this catch...

% read session.log line-by-line and convert to events structure
events = [];
index  = 0;
type   = '';
while true
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        %if ~contains(type,'BLANK') && ~strcmp(type,'E')
        %    fprintf('\n%s\n','Warning: the last line of sessionlog is not BLANK. Last line BLANK is needed to update msDuration in events. Break out, check sessionlog, and ask Courtney or Sam for instructions.')
        %    keyboard;
        %end
        break;
    end
    
    % generic text scan to get time, offset, and type
    xTOT          = textscan(thisLine,'%f%d%s');
    mstime        = xTOT{1}(1);   % must add (1) because numbers after the string in the above line cause overflow to first %f
    msoffset      = xTOT{2}(1);
    type          = xTOT{3}{1};
    
    % default Parameters (details will be filled out/altered based on type)
    experiment    = 'replayTask';
    subject       = subject;
    sessionName   = sessionName;
    sessionNum    = sessNum;  % store in state var so all events are assigned a sessionNum
    
    mstime        = mstime;  
    msoffset      = msoffset;   %#ok<*NASGU,*ASGSL>
    
    % parameters 
    
    taskStr       = ' '; % (char) specifying start of task section  
%     taskSec       = ' '; % 'FUNCLOC' for functional localizer or 'RECALL' for main replay/recall task 
    taskStartms   = nan; % start time of task 
    blockStr      = ' '; % (char) specifying task block 
    blockType     = ' '; % (char) specifying block type within task section
    runStr        = ' '; % (char) specifying repetition number for stimuli presentation 
    runNum        = nan; % (double) repetition number for presentation of stimuli 
    
    trialStr      = ' '; % (char) specifying trial number 
    trialNum      = nan; % (double) specifying trial number 
    
    stimType      = ' '; % (char) specifying stimuli type: funcloc, replay sequence, answer choices, target
    stimStr       = ' '; % (char) sepcifying stimuli number 
    stimNum       = nan; % (double) specifying stimuli number 
    stimlabel     = ' '; % 
    imageFile     = ' '; % (char) of image file 
    imageID       = ' '; % (char) of image ID,'object1.jpg','food2.jpg' etc
    imageCat      = ' '; % (char) specifying image category: ANIMAL, FACE, FOOD, or OBJECT
   
    stimStartms   = nan; % start time of stimuli presentation 
    stimEndms     = nan; % (double) time stimuli was removed from screen 
    stimDurms     = nan; % (double) stimuli duration
%     stimDurms     = nan; 
    ISIDurStr     = ' '; % (char) specifying 'DELTA' or the duration of ISI +/- jitter 
    ISIDurms      = nan; % (double) specifying duration of ISI (+/- jitter)
    
    blockStartms  = nan; % (double) specifying start time of block periods: rest, distractor, recall, full recall etc. 
    blockEndms    = nan; % (double) specifying end time of block periods 
    blockDurms    = nan; % (double) specifying duration of block 
    
    questionsNum  = nan; % (double) specifying recall question, 0 or 1 ('earliest' or 'latest,' respectively)
    stimPos       = nan; % (double) specifying display position of image in recall question task block
    stimInd       = nan; % (double) specifying index of image from original sequence of same trial: 0, 1, 2, 3
    subjRespStr   = ' '; % (char) specifying subject's keyboard response on session.log 
    subjResp      = ' '; % (char) specifying subject's keyboard response: 'LEFT 'RIGHT' or 'UP' 
    correctness   = ' '; % (char) specifying subject's response correctness 
    corr          = nan; % (logic) specifying whether subject's response is correct 
    corrRecallStr = ' '; % (char) specifying number of recall questions correct
    corrNum       = nan; % (double) specifying number of recall questions correct 
    rxnTimeStr    = ' '; % (char) specifying reaction time 
    rxnTime       = nan; % (double) specifying reaction time
    
    vocalLabels   = {};  % (empty cell) to store vocalized image labels 
    vocalTimes    = [];  % (empty array) to store vocalization times
    
    
    TRIAL = sprintf('TRIAL_%d,', (0:20)); % create var TRIAL to store all trial names 'TRIAL_0' etc...
    TRIAL = split(TRIAL, ','); 
    DEC_Q = sprintf('DECISION_QUESTION_%d,', (0:10));
    DEC_Q = split(DEC_Q, ','); 
    
    switch type
        
        case 'SESSION_START' 
            
            xTOT = textscan(thisLine, '%f%d%s%s%s%s%s'); % grab the session number 
            newSession = sessionNum;
            sessStartms = mstime; 
        
        case {'FUNCLOC_START', 'REPLAY_START'}
            
            xTOT = textscan(thisLine, '%f%d%s'); % grab the tast section 
            taskStr = xTOT{3}; 
            taskSec = extractBefore(taskStr, '_START'); % assign task section
            taskStartms = mstime;
            index = index + 1;
                         
        case {'RUN_0', 'RUN_1'} 
            
            xTOT = textscan(thisLine, '%f%d%s%s%s%s%s');
            runStr = xTOT{3};
            stimStr = xTOT{4}; 
            imageFile = xTOT{6};
            
            blockType = 'STIM_SHOW'; % functional localizer stim 
            
            runNum = extractAfter(runStr, 'RUN_'); 
            runNum = str2double(runNum{1}); 
            
            stimType = 'FL'; 
            stimNum = extractAfter(stimStr, 'FUNCLOC_STM_');
            stimNum = str2double(stimNum{1}); 
            
            stimStartms = mstime; 
            
            imageID = extractAfter(imageFile, 'stimuli/'); 
            if contains(imageID, 'animal') == 1, imageCat = 'ANIMAL'; end 
            if contains(imageID, 'face')   == 1, imageCat = 'FACE'  ; end
            if contains(imageID, 'food')   == 1, imageCat = 'FOOD'  ; end
            if contains(imageID, 'object') == 1, imageCat = 'OBJECT'; end
            
            stimlabel = upper(char(fullCSV.Var1(strcmp(imageFile,fullCSV.Var2)))); % extract image name from CSV
            
            index = index + 1;
            
        case {'ISI_START', 'ISI_END'}
            
            xTOT = textscan(thisLine, '%f%d%s%s');
            
            if contains(type, 'START')
                
                ISIStartms = mstime; 
                if ~isempty(events)
                    if isnan(events(end).stimEnd), events(end).stimEnd = mstime; end % should always be true, image ends when ISI starts
                    if isnan(events(end).stimDuration), events(end).stimDuration = events(end).stimEnd - events(end).stimStart; end
                end
                
            elseif contains(type, 'END')
                ISIDurStr = xTOT{4}; 
                ISIDurms = extractAfter(ISIDurStr, 'DELTA_'); 
                ISIDurms = str2double(ISIDurms{1}); 
                if isnan(events(end).ISIDuration), events(end).ISIDuration = ISIDurms; end 
            end 
            
        case {'FUNCLOC_REST_START', 'FUNCLOC_REST_END'}

            blockType = 'REST'; 
            if contains(type, 'FUNCLOC_REST_START')
                blockStartms = mstime; 
                index = index + 1;    
            elseif contains(type, 'FUNCLOC_REST_END')
                events(end).blockEnd = mstime;
                events(end).blockDuration = events(end).blockEnd - events(end).blockStart; 
            end 
         
        case {'FREE_RECALL_PROMPT', 'FREE_RECALL_END'} % start at 'prompt' because subjects often begin recalling immediately 
            
            blockType = 'FREE_RECALL'; 
            if contains(type, 'FREE_RECALL_PROMPT')
                blockStartms = mstime;  
                index = index + 1;
            elseif contains(type, 'FREE_RECALL_END') 
                events(end).blockEnd = mstime;  
                events(end).blockDuration = events(end).blockEnd - events(end).blockStart;
            end 
           
            
        case {TRIAL{:}} % containing any 'TRIAL_X' X is between 0 & 20 
            
            xTOT = textscan(thisLine, '%f%d%s%s%s%s%s%s');
            trialStr = xTOT{3};
            
            trialNum = extractAfter(trialStr, 'TRIAL_');
            trialNum = str2double(trialNum{1}); 
            
            block = xTOT{4}{1}; 
            
            switch block 
                
                case {'REPLAY_STM_0', 'REPLAY_STM_1', 'REPLAY_STM_2', 'REPLAY_STM_3'}

                    stimStr = xTOT{4};
                    imageFile = xTOT{6};

                    blockType = 'STIM_SHOW'; % replay stim 

                    stimType = 'RP'; 
                    stimNum = extractAfter(stimStr, 'REPLAY_STM_');
                    stimNum = str2double(stimNum{1});
                    stimStartms = mstime; 

                    imageID = extractAfter(imageFile, 'stimuli/'); 
                    if contains(imageID, 'animal') == 1, imageCat = 'ANIMAL'; end 
                    if contains(imageID, 'face')   == 1, imageCat = 'FACE'  ; end
                    if contains(imageID, 'food')   == 1, imageCat = 'FOOD'  ; end
                    if contains(imageID, 'object') == 1, imageCat = 'OBJECT'; end
                    
                    stimlabel = upper(char(fullCSV.Var1(strcmp(imageFile,fullCSV.Var2))));
                    
                    index = index + 1;
                
                case {'REST_START', 'REST_END'}
                    
%                     blockType = 'REST'; 
%                     if contains(block, 'REST_START'), blockStartms = mstime;    
%                     elseif contains(block, 'REST_END') blockEndms = mstime; blockDurms = blockEndms - events(index).blockStart; end 
%                     index = index + 1;

                    blockType = 'REST'; 
                    if contains(block, 'REST_START')
                        blockStartms = mstime; 
                        index = index + 1;    
                    elseif contains(block, 'REST_END')
                        events(end).blockEnd = mstime;
                        events(end).blockDuration = events(end).blockEnd - events(end).blockStart; 
                    end 
                    
                case {'DISTRACTOR_START', 'DISTRACTOR_END'} 
                    
%                     blockType = 'DISTRACTOR'; 
%                     if contains(block, 'DISTRACTOR_START'), blockStartms = mstime;    
%                     elseif contains(block, 'DISTRACTOR_END') blockEndms = mstime; blockDurms = blockEndms - events(index).blockStart; end 
%                     index = index + 1;
                    
                    blockType = 'DISTRACTOR'; 
                    if contains(block, 'DISTRACTOR_START')
                        blockStartms = mstime; 
                        index = index + 1;    
                    elseif contains(block, 'DISTRACTOR_END')
                        events(end).blockEnd = mstime;
                        events(end).blockDuration = events(end).blockEnd - events(end).blockStart; 
                    end 
                
                case {'RECALL_QUESTION_0', 'RECALL_QUESTION_1'} 
                    
                    if block == 'RECALL_QUESTION_0', questionsNum = '0'; end %#ok<BDSCA> % which image appeared earliest in sequence?
                    if block == 'RECALL_QUESTION_1', questionsNum = '1'; end %#ok<BDSCA> % which image appeared latest in sequence? 
                    stimPosStr = xTOT{7};
                    stimIndStr = xTOT{8};
%                     recallQStart 

                    if contains(xTOT{5}, 'IMAGE_ON_SCREEN') 
                        
                        blockType = 'RECALL_Q'; 
                        
                        imageFile = xTOT{6};
                        
                        imageID = extractAfter(imageFile, 'stimuli/'); 
                        if contains(imageID, 'animal') == 1, imageCat = 'ANIMAL'; end 
                        if contains(imageID, 'face')   == 1, imageCat = 'FACE'  ; end
                        if contains(imageID, 'food')   == 1, imageCat = 'FOOD'  ; end
                        if contains(imageID, 'object') == 1, imageCat = 'OBJECT'; end
                        
                        stimlabel = upper(char(fullCSV.Var1(strcmp(imageFile,fullCSV.Var2))));

                        stimType = 'CHOICE'; 
                        stimPos = extractAfter(stimPosStr, 'POSITION_'); 
                        stimInd = extractAfter(stimIndStr, 'INDEX_'); 
                        stimInd = str2double(stimInd{1}); 
                        
                        stimStartms = mstime;
                        index = index + 1; 
                    
                    end 
                    
                    if contains(xTOT{5}, 'SUBJECT_RESPONSE')
                        
                        blockType = 'RECALL_R'; 
                        
                        subjRespStr = xTOT{5}; 
                        subjResp = extractAfter(subjRespStr, 'SUBJECT_RESPONSE_'); 
                        index = index + 1;
                        
                    end                             

                    if contains(xTOT{5}, 'CORRECT')
                        
                        correctness = xTOT{5}(1); 
                        if contains(correctness{1}, 'IN'), events(end).correct = false; 
                        else, events(end).correct = true; end
                        rxnTimeStr = xTOT{6}; 
                        rxnTime = extractAfter(rxnTimeStr{1}, 'REACTION_TIME_'); 
                        events(end).reactionTime = str2double(rxnTime);
                            
                    end
        

%                 case {'CORRECT_RECALL_0', 'CORRECT_RECALL_1', 'CORRECT_RECALL_2'} 
%                     
%                     blockType = 'RECALL_Q_SUM'; % recall questions summary for trial 
%                     corrRecallStr = xTOT{4}; 
%                     corrNum = extractAfter(corrRecallStr, 'CORRECT_RECALL_'); 
%                     corrNum = str2double(corrNum{1}); 
%                     index = index + 1; 
                    
                case {'FULLRECALL_PROMPT', 'FULLRECALL_END'}
                    
%                     blockType = 'FULL_RECALL';
%                     if contains(block, 'FULLRECALL_START'), blockStartms = mstime;    
%                     elseif contains(block, 'FULLRECALL_END') blockEndms = mstime; blockDurms = blockEndms - events(index).blockStart; end 
%                     index = index + 1; 
                    
                    blockType = 'FULL_RECALL'; 
                    if contains(block, 'FULLRECALL_PROMPT')
                        blockStartms = mstime; 
                        index = index + 1;    
                    elseif contains(block, 'FULLRECALL_END')
                        events(end).blockEnd = mstime;
                        events(end).blockDuration = events(end).blockEnd - events(end).blockStart; 
                    end 

                                 
                case {DEC_Q{:}}  %#ok<*CCAT1>
                    
                    xTOT = textscan(thisLine, '%f%d%s%s%s%s%s');
                    questionStr = xTOT{4}; 
                    questionsNum = extractAfter(questionStr{1}, 'DECISION_QUESTION_'); 
                    questionsNum = str2double(questionsNum); 
                    
                    if contains(xTOT{5}, 'IMAGE_ON_SCREEN') 
                        
                        blockType = 'DECISION_Q';
                        
                        imageFile = xTOT{6};
                        stimIndStr = xTOT{7};
                        
                        imageID = extractAfter(imageFile, 'stimuli/'); 
                        if contains(imageID, 'animal') == 1, imageCat = 'ANIMAL'; end 
                        if contains(imageID, 'face')   == 1, imageCat = 'FACE'  ; end
                        if contains(imageID, 'food')   == 1, imageCat = 'FOOD'  ; end
                        if contains(imageID, 'object') == 1, imageCat = 'OBJECT'; end
                        
                        stimlabel = upper(char(fullCSV.Var1(strcmp(imageFile,fullCSV.Var2))));
                        
                        stimType = 'CHOICE'; 
                        stimInd = extractAfter(stimIndStr{1}, 'INDEX_'); 
                        stimInd = str2double(stimInd); 
                        
                        stimStartms = mstime;
                        
                        if contains(xTOT{5}, 'TARGET_IMAGE_ON_SCREEN'), stimType = 'TARGET'; end
                        
                        index = index + 1; 
                    end
                    
                    if contains(xTOT{5}, 'SUBJECT_RESPONSE')
                        
                        blockType = 'DECISION_R';
                        
                        subjRespStr = xTOT{5}(1); 
                        subjResp = extractAfter(subjRespStr, 'SUBJECT_RESPONSE_'); 
                        
                        index = index + 1; 
                        
                    end
                    
                    if contains(xTOT{5}, 'CORRECT')
                        
                        correctness = xTOT{5}(1); 
                        if contains(correctness{1}, 'IN'), events(end).correct = false; 
                        else, events(end).correct = true; end
                        rxnTimeStr = xTOT{6}; 
                        rxnTime = extractAfter(rxnTimeStr{1}, 'REACTION_TIME_'); 
                        events(end).reactionTime = str2double(rxnTime);
                    end
                    
                    
            end
        
        case {'RECORDING_START'}
            
            MISSING_ANN = 0;
            
            xTOT = textscan(thisLine, '%f%d%s%s');
            annFileRoot = xTOT{4}{1}; 
            annFileName = sprintf('%s.ann', annFileRoot);
            annFile = fullfile(sessionDir, annFileName);
            
            if ~exist(annFile,'file')
                if MISSING_ANN == 0
                    fprintf('\n >>> %s (and possibly others) were not found in %s',annFileName,sessionDir);
                end
                MISSING_ANN = MISSING_ANN + 1;
                resultStr = sprintf('ANN: no ann file found'); 
            
            else
                
                % ann file present... process it
                fid2 = fopen(annFile,'r');
                
                if fseek(fid2,1,'bof')==-1 %annotation file is empty
                    fprintf('\n%s is empty',annFile); keyboard;
                else
                    
                    fseek(fid2,0,'bof');
                    while true
                        tmpAnnLine = fgetl(fid2);
                        if ~ischar(tmpAnnLine)       
                            break
                        end
                        if numel(tmpAnnLine)==0      
                            continue 
                        end
                        if strcmp(tmpAnnLine(1),'#') 
                            continue
                        end % advance past comments and empty lines
                        
                        xANN = textscan(tmpAnnLine,'%f%f%s');

                        if contains(annFileRoot, 'FUNCLOC')
                
                            blockType = 'VOCALIZED';
                            
                            wordTime = round(xANN{1}); 
                            wordNum = xANN{2}; 
                            imgWord = xANN{3}{1}; 
                            
                            vocalLabels{end+1} = imgWord;   %#ok<*AGROW>
                            vocalTimes(end+1) = wordTime + mstime; 

                        end 
            
                        if contains(annFileRoot, 'FREE_RECALL')
                            
                            blockType = 'FREE_RECALL';

                            wordTime = round(xANN{1}); 
                            wordNum = xANN{2}; 
                            imgWord = xANN{3}{1}; 
                            
                            vocalLabels{end+1} = imgWord; 
                            vocalTimes(end+1) = wordTime + mstime; 
                            
                            events(end).vocalizedImages = vocalLabels; 
                            events(end).vocalizedTimes = vocalTimes; 
                            
                        end
            
                        if contains(annFileRoot, 'FULLRECALL')
                            
                            blockType = 'FULL_RECALL';
                            wordTime = round(xANN{1}); 
                            wordNum = xANN{2}; 
                            imgWord = xANN{3}{1}; 
                            
                            vocalLabels{end+1} = imgWord; 
                            vocalTimes(end+1) = wordTime + mstime; 
                            
                            events(end).vocalizedImages = vocalLabels; 
                            events(end).vocalizedTimes = vocalTimes; 
                            
                        end 
                        
                    end
                    
                    fclose(fid2);
                    
                end
                
            end
            
            if contains(annFileRoot,'FUNCLOC')
                 index = index + 1; 
            end

    end 

    % assign values to events array 
    if index > length(events)
        
        % create dummy event structure that is upddated below based on type
        
        clear thisEvent
        
        % experiment & section, task, block info: name, start/end times, durations
        thisEvent.experiment        = experiment;
        thisEvent.subject           = subject;
        thisEvent.sessionName       = sessionName;
        thisEvent.sessionNum        = sessionNum;   % store in state var so all events are assigned a sessionNum
        thisEvent.sessionStart      = sessStartms; 
        thisEvent.mstime            = mstime;
        thisEvent.taskSection       = taskSec{1};
        thisEvent.taskStart         = taskStartms; 
        thisEvent.blockType         = blockType; 
        thisEvent.blockStart        = blockStartms; 
        thisEvent.blockEnd          = blockEndms; 
        thisEvent.blockDuration     = blockDurms; 
        
        % run, trial, or stim info: number, file, category, start/end times, duration 
        thisEvent.runNumber         = runNum; 
        thisEvent.trialNumber       = trialNum; 
        thisEvent.stimType          = stimType; 
        thisEvent.stimNumber        = stimNum; 
        thisEvent.stimFile          = imageFile;
        thisEvent.stimCategory      = imageCat; 
        thisEvent.stimlabel         = stimlabel; 
        thisEvent.stimStart         = stimStartms; 
        thisEvent.stimEnd           = stimEndms; 
        thisEvent.stimDuration      = stimDurms; 
        thisEvent.ISIDuration       = ISIDurms; 
        
        % question info: number, choices, response, reaction time, correctness
        thisEvent.questionNumber    = questionsNum;  
        thisEvent.choicePosition    = stimPos; 
        thisEvent.choiceIndex       = stimInd; 
        thisEvent.subjResponse      = subjResp; 
        thisEvent.correct           = corr; 
        thisEvent.reactionTime      = rxnTime; 
        thisEvent.numberCorrect     = corrNum; 
   
        % annotation file info, voclized images & vocalization times 
        thisEvent.vocalizedImages   = vocalLabels; 
        thisEvent.vocalizedTimes    = vocalTimes; 
        
        if (index==1)
            events = thisEvent; % before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
        
    end 
    
end

fclose(fid);  % close session.log

end



