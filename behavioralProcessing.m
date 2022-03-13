function behavioralProcessing(subject,rootEEGdir,taskType,altBehDir)

% BehavioralProcessing.m€™ extracts events for each task.
%
% Input:
%      subject = 'NIHXXX',
%      rootEEGdir = ~/data/eeg',
%      taskType = 'languageTask'
%      altBehDir = 'behavioral_preOp', 'behavioral_postOp', 'behavioral'  < optional parameter, just pass first 3 to use default 'behaioral' directory
% Output:
%      saves out event.mat to behavioral session directories
%


% optional input arguement to all for preOp and postOp behavioral directory analysis

if nargin==3,  behFolder = 'behavioral';
else           behFolder = altBehDir;
    if ~strcmp(behFolder(1:10),'behavioral'),  %- full string should be: 'behavioral', 'behavioral_preOp', 'behavioral_postOp'
        fprintf('\n Uh oh, alternate behavioral folder (%s) may not be specified correctly', altBehDir);
    end
end


subjDir  = fullfileEEG(rootEEGdir,subject);
behDir   = fullfileEEG(subjDir,behFolder,taskType);  % 'behavioral', but could also have 'behavioral_preOp' or 'behavioral_postOp' here
sessions = dir(fullfileEEG(behDir,'session_*'));

subjNum  = str2num(subject(4:6)); %- NIHXXX or BEHXXX


fprintf('\n\n\n------------------------------------------------------------------------------\nAttempting to create event.mat files in %d sessions of: %s\n', length(sessions), behDir);



%%- quick loop through sessions to correctly order double digit session numbers (else 10 comes after 1 and before 2)
if length(sessions)==0,  fprintf('\n experiment folder found, but no valid session folders within; exiting');   return; end
for iSess = 1:length(sessions),
    sessName    = sessions(iSess).name;
    strNumeric  = find( sessName >= '0' & sessName <= '9');  % following conditional statement pulls out the first chunk of numbers and assumes that is the session number
    if max(diff(strNumeric))>1, 
        iKeep=[1:find(diff(strNumeric)>1,1,'first')]; 
        %fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessName, sessName(strNumeric), sessName(strNumeric(iKeep))); 
        strNumeric=strNumeric(iKeep); 
    end;
    sessNum     = str2num( sessName(strNumeric) );  if isempty(sessNum), sessNum=iSess; end; %shouldn't need this catch...
    sessNumAr(iSess) = sessNum;
end
if length(sessions)>0,
    [sortVal sortInd] = sort(sessNumAr);
    sessions = sessions(sortInd);
    sessNumAr = sortVal;
end

if length(sessNumAr)>length(unique(sessNumAr)),
    fprintf('\n Warning: at least two session folders with same numeric value found in behavioral processing;\n');
    %fprintf('       "session numbers" will be modified to reflect index into session array\n');  sessNumAr=[1:length(sessions)];
end
%%- loop through sessions, create events file for each session, then master events.mat in root dir
allEvents = [];
allUItables = []; % JIC

for iSess = 1:length(sessions)
    
    %- construct strings or session numbers from "sessions" to pass to individual extraction functions
    sessName    = sessions(iSess).name;
    sessNum     = sessNumAr(iSess);
    %sessNum     = str2num( sessName(find(sessName=='_')+1:end) );
    %if isempty(sessNum), sessNum=iSess-1; fprintf('\n Warning: %s getting assigned session number %d \n',sessName, sessNum); end; %shouldn't need this catch... but in case session_1 directory renamed session_1trim
    sessionDir  = fullfileEEG(behDir,sessions(iSess).name);
    sessFileStr = fullfileEEG(behDir,sessions(iSess).name,'session.log');
    eventfile   = fullfileEEG(sessionDir,'events.mat');
    events = [];
    switch taskType
        case {'attentionFeedback'}
            if subjNum>=89
                events = extractAttnFeedbackEvents_v2(sessFileStr, subject, sessName, allEvents);  %added 3/2021 by CB
            else
                events = extractAttnFeedbackEvents_v1(sessFileStr, subject, sessName, allEvents);  %added 9/2020 by CB
            end
        case {'attentionTask' 'attentionZap'}
            if (strcmp(subject(1:3),'NIH') & str2num(subject(4:6))>=14) | strcmp(subject(1:3),'BEH'),
                events = jwAttnTaskEvents_v3b(sessFileStr, subject, sessName, allEvents);  % version 3a extracts free recall info
            else
                events = extractAttentionEvents_v2(subject, behDir, sessNum);  % this *may* not work for all subj < NIH016...
            end

        case 'auditoryVowels'
            events = extractAudVowEvents(subject,behDir, sessNum);
            
        case 'baseline'
            [events, event_table] = create_BL_eventsPrepAndAlign(subject,sessionDir,rootEEGdir,sessNum);            
            % JIC: if the bit below causes any issues in the pipeline feel free to remove
            hgexport(event_table, fullfile(sessionDir,'events_table'), hgexport('factorystyle'), 'Format', 'png'); % this is an undocumented hack so it may break in future releases
            delete(event_table)
            
%             % save master table
%             if iSess == length(sessions)
%                 allUItables(1).Data = vertcat(allUItables.Data);
%                 allUItables = allUItables(1);
%             end
                        
        case 'delayMatching'
            events = extractdelayMatching_v1(sessFileStr, subject, sessName,[]);
            
        case {'freeRecall','freeRecallStim'}
            %if strcmp(subject(1:3),'NIH') & str2num(subject(4:6))<43,  %- jw think's NIH043 is when freeRecall session.log changed (moved to system2.0)
            %  events = RAM_FR_CreateAllEvents(subject, behDir, sessNum);
            %elseif strcmp(subject(1:3),'NIH') & str2num(subject(4:6))<48,  %- jw think's NIH048 is when freeRecall session.log changed (moved to system3.0)
            if subjNum<48,  %- jw think's NIH048 is when freeRecall session.log changed (moved to system3.0)
                events = RAM_FR_sys2_CreateTASKEvents(subject, behDir, sessNum, sessName);  %- not sure if this works for pre-system 2.0... if not, uncomment lines above
            else
                events = RAM_FR_sys3_CreateTASKEvents(subject, behDir, sessNum, sessName);
            end
            
        case 'goAntiGo'
            events    = extractGoAntiGoEvents(subject,behDir,num2str(sessNum)); % extract the events from the session.log file
            
        case 'languageTask'
            if strcmp(subject(1:3),'NIH') & str2num(subject(4:6))<20,
                events    = extractLangEvents_v1(subject,behDir,sessNum);
            else
                events    = extractLangEvents(subject,behDir,sessNum);   %- seems to have problems with NIH025 and NIH026...
            end
            
        case 'ltmStrupt'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/ltmStrupt'];
            [events] = ltm_Strupt_Extract_Event(behDir,sessName);
            
        case 'memoryDecision'
            events    = extractMemoryDecisionEvents_v1(sessFileStr, subject, sessName, []);
            
        case 'manipMemExp'
            if strcmp(subject,'NIH053')
                events    = extractManipMemEventsExp_v01(subject, behDir, sessNum);
            else
                events    = extractManipMemEventsExp_v02(subject, behDir, sessNum);
            end
            
        case 'manipMemImp'
            events = extractManipMemEventsImp_v01(subject,behDir,sessNum);
        case 'moveTask'
            events    = extractMoveTask(subject,behDir,sessNum);
            
        case 'motivation'
            events    = extractMotivationEvents(sessFileStr,subject,sessName); %(sessLogFile, subject, sessionName)
            
        case 'mSeqGUI'
            % This will be the research stim behavioral folder for NIH077
            % onwards. will use same extraction code as stimMapGUI
            events = extractStimEventsFromAll(subject, sessName, behDir, rootEEGdir);
            
        case {'multiModal' 'multiModalAV'}
            % For audiovisual version
            if strcmp(subject,'NIH045')|| strcmp(subject,'NIH046')
                events = extractAVMultimodalEvents(rootEEGdir, subject, sessName, allEvents);
                % Old version
            else
                stimKeyDir = fullfile(behDir,'stimKey.txt');
                events = extractMultiModal(stimKeyDir,sessFileStr, subject, sessName);
            end
            
        case 'namingTask'
            events = extractNamingTask(sessFileStr, subject, sessName);
            
        case 'pa3'
            numBlocks = 15;
            numTrials = 4;
            %touchAnnFiles(sessionDir,numBlocks,numTrials);  %- creates fake (empty) annotation files if ann not found... dont do this! better to wait for annotation!
            try
                events    = extractPA3events(subject,behDir,sessNum,sessionDir);
            catch err
                disp(getReport(err,'extended'));
                events    = [];
            end
            
        case {'palRam', 'palRamStim', 'PAL'}
            if subjNum < 48 || subjNum > 58
                %- palRAM, system 1.0 and system 2.0:            NIH026 to 47
                %- "RAM", which is the same as system 1.0/2.0:   NIH062 to
                events = RAM_PAL_CallCreateTaskFunctions(subject,behDir, sessNum, sessName);
            else
                %- palRAM, system 3.0 (comma value pairs in session.log)
                %- NIH048 to 058
                events = RAM_PAL_CallCreateTaskFunctions_sys3(subject,behDir, sessNum, sessName);
            end
            
        case  {'PALrepeat'}
            events = RAM_PALrepeat_CallCreateTaskFunctions(subject,behDir, sessNum, sessName);
            
            
        case 'paRemap'
            events    = paRemap_ExtractEvents(sessFileStr, subject, sessName);
            
        case 'paRepeat'
            events    = paRepeat_ExtractEvents(sessFileStr, subject, sessionDir);
            
        case 'playPass'
            events    = extractPlayPassEvents(subject,behDir, sessNum);
            
        case 'prob_sel'
            events    = [];
            fprintf('\n WARNING: prob_sel does not have extraction code yet, so no events.mat being created');
            
        case 'restingState'
            events = extractRestingState(sessFileStr, subject, sessName);
            
        case 'seizure'
            events = load_SZ_info_PrepAndAlign(subject,rootEEGdir);
            
        case 'semanticSpan'
            events = semanticSpan_extractEvents(sessFileStr, sessName);
            
        case 'SequenceMem'
            events    = extractSequenceMemEvents(subject,behDir,num2str(sessNum)); % extract the events from the session.log file
            
        case 'SerialRT'
            events    = extractSerialRTEvents(subject,behDir,num2str(sessNum));
            
        case 'stimMapAnn'
            %- this should be called for NIH020 to NIH045.... NIH035 to NIH045 have a mix of stimMapAnn and stimMapGUI
            events = extractStimEventsFromAnn_v10(sessFileStr, subject, sessName, sessNum, rootEEGdir); %- used to convert NK annotations into a session.log and eeg.eeglog (requires raw/STIM_MAP NK files are initailly processed by createStimMapSessFromAnn_v01(rootrootEEGdir,subj);
            
        case 'stimMapGUI'
            %- this should be called for >=NIH035.... NIH035 to NIH045 have a mix of stimMapAnn and stimMapGUI, NIH045 and up just has stimMapGUI
            if subjNum >= 62 || contains(subject,'TRE'),
                %- Dave updated the stimLogs and event extraction code starting at NIH062, this code will deal with small differences between subjects above 62
                events = extractStimEventsFromAll(subject, sessName, behDir, rootEEGdir);
            elseif subjNum < 62,
                %- Tim and JW created the original version of stimMapGUI and event extraction code.
                %- for subjects 35 to 61, excluding 40
                if subjNum==40,
                    events = extractStimEventsFromOGstimLog_v10(sessFileStr, subject, sessName, sessNum, rootEEGdir);  %- stim session using Tim's matlab GUI to control the cerestim
                else
                    %- this should be called for >=NIH035.... NIH035 to NIH045 have a mix of stimMapAnn and stimMapGUI, NIH045 and up just has stimMapGUI
                    events = extractStimEventsFromStimLog_v10(subject, sessName, behDir, rootEEGdir);  %- stim session using Tim's matlab GUI to control the cerestim
                end
            end
            
        case 'stimulusBlast'
            events = stimBlast_ExtractEvents_v2(sessFileStr, subject, sessName);
            
        case 'stimulusGuess'
            events = stimGuess_ExtractEvents_v3(sessFileStr, subject, sessName);
            
        case 'theDoors'
            if strcmp(subject,'NIH042')|| strcmp(subject,'NIH043')
                events    = extractDoorsEvents(sessFileStr, subject, sessName, allEvents);
            elseif strcmp(subject,'NIH041')
                events    = extractDoorsEvents_NIH041(sessFileStr, subject, sessName, allEvents);
            end
            
        case 'zapList'
            events = extractzapList_v1(sessFileStr, subject, sessName,[]);
            
        case 'stm_Face'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/stm_Face'];
            events = Face_STM_Extract_Event(behDir, sessName);
            
        case 'stm_Color'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/stm_Color'];
            events = Color_STM_Extract_Event(behDir, sessName);
            
        case 'stm_Color_v2'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/stm_Color_v2'];
            events = Color_STM_Extract_Event(behDir, sessName);     
            
        case 'stm_Color_stimulation'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/stm_Color_stimulation'];
            [events] = stm_color_stimulation_eventExtract(behDir,sessName);
            
        case 'stm_Replay'
            behDir = [rootEEGdir '/' subject '/behavioral/' '/stm_Replay'];
            [events] = stm_Replay_eventExtract(behDir,sessName);  
            
        case 'replayTask'
            [events] = extractReplayTaskEvents(rootEEGdir, subject, sessName); 
    end
    
    
    fprintf('\n%d) extracted %d events from %s', iSess, length(events), sessFileStr);
    if length(events)>0
        save(eventfile,'events', '-v7');  %- version 7 file format is compact... only need version 7.3 if >2GB
    else
        fprintf(' -- NO events, so not creating events.mat');  
    end
    
    if length(allEvents)>0 && length(fieldnames(allEvents))~=length(fieldnames(events))
        fprintf('\n field mismatch... fix it');
        keyboard;
    end
    allEvents=[allEvents, events];
end


%%- Confirm success and save new MASTER matfile at root
if isempty(allEvents)
    fprintf('\nNo events for %s. Exiting without creating master events.mat.', taskType);
else
    % Change name and save to root events.mat
    events=allEvents;
    rootevntFileStr = sprintf('%s/events.mat',behDir);
    save(rootevntFileStr, 'events', '-v7');  %- version 7 file format is compact... only need version 7.3 if >2GB
    fprintf('\n --- %d events saved in %s ---\n', length(events), rootevntFileStr);
end

