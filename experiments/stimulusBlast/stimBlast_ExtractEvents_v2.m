function [events] = stimBlast_ExtractEvents(sessLogFile, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from stimulusBlast %%%%
%
%   extraction designed for stimulusBlast
%
%   create an event for every presented word and every response (no words from the training sections should be included)
%
% JHW 3/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% uncomment following lines to directly run script
% clear all
%
% %rootEEGdir  = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg';
% %rootEEGdir  = '/Volumes/Shares/FRNU/dataWorking/eeg';
% %subject     = 'NIH031';   % EEG002  NIH016
% %sessionName = 'session_1';
% %sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/stimulusBlast',sessionName);
%
% subject     = 'fullTest3b';   % EEG002  NIH016
% sessionName = 'session_0';
% sessionDir  = fullfileEEG('../data/',subject,sessionName);
%
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function needs to be edited to account for a field called BLOCK_BLUR',
% which will contain 5 numerical values separated by a space. 1st value =
% block, and other 4 values = the blur levels corresponding to that target
% category, as updated after that block's performance



%fprintf('\nOn session: %s ', sessionName)


VERBOSE = 0;

%---------------------------------------------------------------------------------------%
%- Confirm Session File Exists & Create cleaned session number string
fid    = fopen(sessLogFile,'r');
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
    sessLogWithKeyFile = fullfileEEG(sessionDir,'sessionWithKeyLog.log');
    sessLogOriginal    = fullfileEEG(sessionDir,'sessionOriginalNoKeyLog.log');
    
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1,
    iKeep=[1:find(diff(strNumeric)>1,1,'first')];
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep)));
    strNumeric=strNumeric(iKeep);
end;
sessionNum = str2num( sessionName(strNumeric) );
if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...
%---------------------------------------------------------------------------------------%


%- make a new session log that integrates keylog entries (more accurate key-press timing)
MAKE_SESS_KEY_LOG = 1;
if exist('sessionDir','var') && ~exist(sessLogWithKeyFile,'file') & MAKE_SESS_KEY_LOG == 1,
    
    keyLogFile = fullfileEEG(sessionDir,'keyboard.keylog');
    if exist(keyLogFile),
        
        fidKL = fopen(keyLogFile,'r');
        spacePressMS   = [];
        spacePressLine = {};
        while 1
            thisLine = fgetl(fidKL);
            if ~ischar(thisLine); break; end
            
            %- convert string into cell array for each element in string
            xTOT = textscan(thisLine,'%s',20);
            thisLineStr = {};
            for k = 1:length(xTOT),
                if ~isempty(xTOT{k}),
                    thisLineStr = cat(2,thisLineStr,xTOT{k});
                end
            end
            
            %- Generic text scan to get time, offset, and type
            mstime  = str2num(thisLineStr{1}); %- dont add the offset
            mstime2 = str2num(thisLineStr{1}) + str2num(thisLineStr{2}); %- add the error term
            
            %- fair way to incorporate error.  Actual event happens sometime between first and second number.  So lets just take 1/2 of second number
            mstime3 = str2num(thisLineStr{1}) + str2num(thisLineStr{2})/2.0;
            
            
            %-only include SPACE bar "P"resses
            if strcmp(thisLineStr{3},'P') &  strcmp(thisLineStr{4},'SPACE')
                %if strcmp(thisLineStr{4},'SPACE')
                spacePressMS(end+1)   = mstime3;
                thisLineOut = sprintf('%s\t%s\t%s\t%s',thisLineStr{1:2},'KEY_LOG_PRESS',thisLineStr{4});
                spacePressLine{end+1} = thisLineOut;
            end
        end
        fclose(fidKL);
        
        
        %- now loop over session.log and get strings + mstimes
        if exist(sessLogOriginal,'file'), fclose(fid); fid = fopen(sessLogOriginal,'r'); end %- use the "original" if it exists
        sessLogMS   = [];
        sessLogLine = {};
        while 1
            thisLine = fgetl(fid);
            if ~ischar(thisLine); break; end
            
            %- convert string into cell array for each element in string
            xTOT = textscan(thisLine,'%s',20);
            thisLineStr = {};
            for k = 1:length(xTOT),
                if ~isempty(xTOT{k}),
                    thisLineStr = cat(2,thisLineStr,xTOT{k});
                end
            end
            
            %- Generic text scan to get time, offset, and type
            mstime       = str2num(thisLineStr{1});
            sessLogMS(end+1)   = mstime;
            sessLogLine{end+1} = thisLine;
        end
        fclose(fid);
        
        %- now combine to create "sesion log with keyboard log"
        fidNewLog = fopen(sessLogWithKeyFile,'w+');
        
        allMS     = [spacePressMS sessLogMS];
        allLines  = {spacePressLine{:} sessLogLine{:}};
        [~,iSort] = sort(allMS);
        allLines  = allLines(iSort);
        %disp(allLines');
        fprintf(fidNewLog,'%s\n',allLines{:});
        fclose(fidNewLog);
        
        successM = movefile(sessLogFile,sessLogOriginal);
        successC = copyfile(sessLogWithKeyFile,sessLogFile);
        fid      = fopen(sessLogFile,'r'); %- this is the updated session log, which is a copy of the combined log
        if successM==0 | successC==0 | fid==-1,
            fprintf('\n Move or copy or open failed... what up?');
            keyboard;
            return;
        end
        
    end
end



%- lists of possible event tags (3rd column of session.log)
sessEvents  = {'B','SESSION_START','STIMULUS_LIST','SESSION_TYPE','SESS_DURATION_TOTAL','SESS_END','E'};
blockEvents = {'BLOCK_PROMPT','BLOCK_START','BLOCK_TARGET','BLOCK_SUMMARY','BLOCK_DURATION','BLOCK_END'};
namingEvent = {'REC_START','NAMING_QUERY','COUNTDOWN_DONE','REC_STOP'};
stimEvents  = {'TXT_ON_SCREEN','IMG_ON_SCREEN','BLANK_SCREEN','KEYBOARD_PRESS'};
stimCats    = {'NOISE','ANIMAL','OBJECT','PERSON','PLACE'};
targetTypes = {'CORE_FIXED','CORE_RAND','UNIQUE'};
nonTargets  = {'CORE_IMG','CORE_TXT','CORE_BOTH','UNIQUE'};
respTypes   = {'CORRECT(hit)','INCORRECT(false_alarm)','INCORRECT(MISS)'};
keyLogEvent = {'KEY_LOG_PRESS'};

%- key for conversion from string to number
%catStr2Num  = {'ANIMAL', 'OBJECT', 'PERSON', 'PLACE','','','','','','NOISE'}; %- make 'NOISE' a 10
catStr2Num  = {'ANIMAL', 'OBJECT', 'PERSON', 'PLACE','BUSH','CLINTON','KENNEDY','TRUMP','','NOISE'}; %- make 'NOISE' a 10
            
%- real responses should never be less than 200-300 ms (visual-motor latency)
MIN_RXN_TIME = 300; %- if rxn time is <300 ms, then assume it is a late response from the previous stimulus


%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 1;
useKeyPress = 0; %- only try to append keypress to stimulus if in stimulus-presentation period
moreToRead  = true;
while 1
    thisLine = fgetl(fid);
    if ~ischar(thisLine); break; end
    
    %- convert string into cell array for each element in string
    xTOT = textscan(thisLine,'%s',20);
    thisLineStr = {};
    for k = 1:length(xTOT),
        if ~isempty(xTOT{k}),
            thisLineStr = cat(2,thisLineStr,xTOT{k});
        end
    end
    
    %- Generic text scan to get time, offset, and type
    mstime       = str2num(thisLineStr{1});
    eventType    = thisLineStr{3};
    
    
    %- default event field values
    makeNewEvent  = 0;  %- set to 1 to create new event, otherwise updating previous event
    isBlockTarget = 0;
    thisRespMade  = nan; %- was the space bar pressed or not
    thisRespCorr  = nan; %- was the press correct
    %-
    switch eventType
        
        
        case {'SESSION_TYPE'}
            thisSessType   = thisLineStr{4}; % TRAIN_SLOW, TRAIN_FAST, TEST
            
            
        case {'BLOCK_TARGET'}
            targetCategory    = thisLineStr{5}; % ANIMAL, OBJECT, PERSON, PLACE
            targetCategoryNum = find(strcmp(targetCategory,catStr2Num));
            
        case {'TXT_ON_SCREEN' 'IMG_ON_SCREEN'}
            %- update duration of pervious stimulus
            if length(events)>0 && isnan(events(end).duration),
                events(end).duration = mstime - events(end).mstime;
            end
            
            makeNewEvent   = 1; %-increment event counter
            thisBlockStr   = thisLineStr{4}; % 0_22, 30_0, etc
            thisCategory   = thisLineStr{5}; % NOISE, ANIMAL, OBJECT, PERSON, PLACE
            thisCategoryNum = find(strcmp(thisCategory,catStr2Num));
            thisTargetStr  = thisLineStr{6}; % TARGET, NONTARGET
            thisRepeatType = thisLineStr{7}; % CORE_IMG, CORE_TXT, CORE_FIXED, CORE_RAND, UNIQUE
            thisStimName   = thisLineStr{8}; % Dog, White House, etc
            stimIDstr      = thisLineStr{9};
            thisStimIDnum  = thisCategoryNum*100 + str2num(stimIDstr(find( stimIDstr >= '0' & stimIDstr <= '9')));
            
            thisIsText     = strcmp(eventType(1:3),'TXT');
            thisIsImage    = strcmp(eventType(1:3),'IMG');
            thisIsStimulus = ~strcmp(thisCategory, 'NOISE');
            thisIsNoise    = strcmp(thisCategory, 'NOISE');
            %if thisIsNoise, makeNewEvent=0; end %- dont bother outputing the noise events for now
            
            thisIsTarget   = strcmp(thisTargetStr,'TARGET');
            thisIsRecall   = 0;
            
            iHyphen        = strfind(thisBlockStr,'_');
            thisBlockNum   = str2num(thisBlockStr(1:iHyphen-1));
            thisStimNum    = str2num(thisBlockStr(iHyphen+1:end));
            
            thisDuration   = nan;
            thisRxnTime    = nan;
            thisRecallTime = nan; %- requires analyzing annotation
            thisIsRecalled = 0;   %- assume not recalled
            %- assume no key-press... that will be overwritten later if one happens
            if thisIsTarget,
                thisResponse = 'MISS'; %- assume no response, which is an unlabeled correct rejection for non-targets'
                thisRespMade = 0; %-
                thisRespCorr = 0; %-
            else
                thisResponse = 'CORRECT_REJECTION'; %- assume no response, which is an unlabeled correct rejection for non-targets'
                thisRespMade = 0; %-
                thisRespCorr = 1; %-
            end
            
            
            
        case 'KEYBOARD_PRESS'
            %- ignore these... switch to using the "key log" press entries... those seem more accurate and dont get lost
            %events(end).response = thisLineStr{6};
            %if ~strcmp(thisLineStr{6},'INCORRECT(MISS)'),
            %    events(end).rxnTime  = str2num(thisLineStr{7});
            %end
            
            
        case 'KEY_LOG_PRESS'
            if index>1 && events(end).isStimulus, 
                pressTime = mstime + 0.5*str2num(thisLineStr{2});  %- keypress events can have big (>100ms) error term, assume in the middle
                rxnTime   = pressTime - events(end).mstime;
                
                %- if response was <200 - 300 ms, then assume it is a late response to the previous stimulus
                if rxnTime < MIN_RXN_TIME && index>2 && events(end-1).isStimulus,
                    iEvResp = length(events)-1;
                else
                    iEvResp = length(events);
                end
                events(iEvResp).rxnTime = pressTime - events(iEvResp).mstime;
                
                if events(iEvResp).isTarget,
                    events(iEvResp).responseStr = 'HIT';
                    events(iEvResp).isRespMade  = 1;
                    events(iEvResp).isRespCorr  = 1;
                else
                    events(iEvResp).responseStr = 'FALSE_ALARM';
                    events(iEvResp).isRespMade  = 1;
                    events(iEvResp).isRespCorr  = 0;
                end
            end

            
        case 'BLANK_SCREEN'
            if length(events)>0 && isnan(events(end).duration),
                events(end).duration = mstime - events(end).mstime;
            end
            
            
        case 'REC_START'
            
            annFileRoot  = thisLineStr{5};
            annFileName  = sprintf('%s.ann',annFileRoot);
            
            
            %- prep for "recall" event
            thisStimName   = annFileName;
            thisStimIDnum  = nan;
            thisCategory   = '';
            thisCategoryNum = nan;
            thisTargetStr  = '';
            thisRepeatType = ''; %- was this stimulus target or not?
            thisBlockNum   = str2num(thisLineStr{4});
            thisStimNum    = nan;
            
            thisDuration   = nan;
            thisResponse   = ''; %- assume no response, which is an unlabeled correct rejection for non-targets'
            thisRxnTime    = nan;
            thisRecallTime = nan; %- requires analyzing annotation
            thisIsRecalled = 0;   %- assume not recalled
            thisRespMade   = 0; %-
            thisRespCorr   = 0; %-
            
            %- booleans for selecting events
            thisIsText     = 0;
            thisIsImage    = 0;
            thisIsStimulus = 0;
            thisIsNoise    = 0;
            thisIsTarget   = 0; %- only applies to stimuli, not statement of what the target is
            thisIsRecall   = 1;

            
            %- initialize varibles for storing free recall info
            frWords  = {};
            frTimes  = [];
            frStimEv = {};
            frTarStr = {};
            
            % go through annotation file and get all the recalls here.
            annFile      = fullfile(sessionDir,annFileName);
            if exist(annFile,'file'),
                
                %- ann file present... process it
                fidKL = fopen(annFile,'r');
                if fseek(fidKL,1,'bof')==-1 %annotation file is empty
                    fprintf('\n%s is empty',annFile); keyboard;
                else
                    fseek(fidKL,0,'bof');
                    while true
                        tmpAnnLine=fgetl(fidKL);
                        if ~ischar(tmpAnnLine);      break;    end
                        if numel(tmpAnnLine)==0;     continue; end
                        if strcmp(tmpAnnLine(1),'#');continue; end %- advance past comments and empty lines
                        
                        x2=textscan(tmpAnnLine,'%f%f%s');
                        thisRT = round(x2{1});
                        thisWordNum = x2{2};
                        thisRecWord = x2{3}{1};
                        
                        if thisRT>20000,
                            fprintf('\n crazy RT... probably a bug in total recall file save; reannoate this one NOW');
                            keyboard;
                        end
                        
                        %- find the corresponding stimulus
                        iStimEv = find([events.blockNum]==thisBlockNum & strcmp(upper(thisRecWord),upper({events.stimName})));  %- total recall changes capitolization
                        
                        %- if stimulus was found in this block mark it
                        if     length(iStimEv)==1,
                            if isnan(events(iStimEv).recallTime),
                                events(iStimEv).recallTime = thisRT+mstime;  %- go back and label the encoding event (only with the first instance of recall)
                                events(iStimEv).isRecalled = 1;  %- it was recalled (different from "isRecall", which means this is a recall event
                            end
                            if events(iStimEv).isTarget,
                                thisTargetStr = 'TargetFromThisBlock';
                            else
                                thisTargetStr = 'NonTargetFromThisBlock';
                            end
                            
                        elseif length(iStimEv)==0,
                            iStimEvPrior = find([events.blockNum]<thisBlockNum & strcmp(upper(thisRecWord),upper({events.stimName})));  %- total recall changes capitolization
                            if length(iStimEvPrior)>0,
                                thisTargetStr = 'IntrusionFromPreviousBlock';
                            else
                                thisTargetStr = 'IntrusionNotSeenBefore';
                            end
                            %fprintf('\n didnt find match');
                            iStimEv = [];
                            
                        else
                            fprintf('\n shouldnt happen...'); keyboard;
                        end
                        
                        %- save cell array of info so recall events can get saved as well
                        frWords{end+1}  = thisRecWord;
                        frTimes(end+1)  = thisRT+mstime;
                        frStimEv{end+1} = iStimEv;
                        frTarStr{end+1} = thisTargetStr;
                    end
                end
                fclose(fidKL);
                
                makeNewEvent = length(frWords); %- new event for each recalled word
             
            else  %- ann file not found
                thisIsRecall    = nan;          %- only condition where this is nan
                frWords{end+1}  = annFileName;  %-
                frTimes(end+1)  = mstime;
                frStimEv{end+1} = [];
                frTarStr{end+1} = 'AnnotationNotFound';
                makeNewEvent = length(frWords); %- empty annotation event created for each list with a missing annotation
                
            end % if ~exist(annFile,'file'),
            
            
        otherwise
            %- nothing to process on this line
            
    end
    
    
    
    %- if this is a new stimulus, create an event...
    for ii=1:makeNewEvent,
        
        if thisIsRecall==1 | isnan(thisIsRecall),
            mstime        = frTimes(ii);
            thisStimName  = frWords{ii};
            thisTargetStr = frTarStr{ii};
            if ~isempty(frStimEv{ii}),
                thisStimNum  = events(frStimEv{ii}).stimNum;
            else
                thisStimNum  = nan;
            end
        end
        
        clear thisEvent
        thisEvent.experiment        = 'stimulusBlast';
        thisEvent.subject           = subject;
        thisEvent.sessionName       = sessionName;
        thisEvent.sessionNum        = sessionNum;
        thisEvent.sessionType       = thisSessType;
        thisEvent.mstime            = mstime;
        thisEvent.eventType         = eventType; %
        thisEvent.logStr            = thisLine;
        thisEvent.stimName          = thisStimName;
        thisEvent.stimIDnum         = thisStimIDnum;     %- 1 through 60 * 2 * 4 (numStim * img vs text * num categories)
        thisEvent.targetCategory    = targetCategory;    %- target for this block
        thisEvent.targetCategoryNum = targetCategoryNum; %- 
        thisEvent.stimCategory      = thisCategory;     %- category of this stimulus
        thisEvent.stimCategoryNum   = thisCategoryNum;  %- [nan, 0,1,2,3,4 -->  nonStim, noise,animal,object,person,place]
        thisEvent.targetStr         = thisTargetStr;    %- was this stimulus target or not?
        thisEvent.repeatType        = thisRepeatType;   %- CORE_BOTH CORE_IMG CORE_TXT UNIQUE
        thisEvent.blockNum          = thisBlockNum;
        thisEvent.stimNum           = thisStimNum;
        thisEvent.duration          = thisDuration;
        thisEvent.responseStr       = thisResponse;     %- string describing response (hit/miss/correct rejection/false alarm)
        thisEvent.isRespMade        = thisRespMade;     %- was the space bar pressed or not
        thisEvent.isRespCorr        = thisRespCorr;     %- was the response correct or not
        thisEvent.isRecalled        = thisIsRecalled;   %- was the response correct or not
        thisEvent.rxnTime           = thisRxnTime;
        thisEvent.recallTime        = thisRecallTime;
        %- booleans for selecting events
        thisEvent.isText            = thisIsText;
        thisEvent.isImage           = thisIsImage;
        thisEvent.isStimulus        = thisIsStimulus;   %- image or text that is NOT noise
        thisEvent.isNoise           = thisIsNoise;
        thisEvent.isTarget          = thisIsTarget;     %- stimulus that is in the target category
        thisEvent.isRecall          = thisIsRecall;     %- is this a recall event (i.e., ms time = vocalization time); nan per blast if annotation not found
        
        %- add the event to the event list
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end %- if (index==1)
        index = index+1;
        
    end % if isNewEvent
    
    
end %- while 1
fclose(fid);  % close session.log


if VERBOSE==0, return; end

%- enumerate possible categories: PEOPLE, PLACES, OBJECTS, ANIMALS
ev = events; %- all stimuli
catList  = unique({ev.stimCategory});

%- organize the stimuli base on categories
stimList = {};  stimCat = {};  stimCatNum = [];
for thisCat=catList,
    thisCatStim = unique({events(strcmp({ev.stimCategory},thisCat)).stimName});
    
    tmpCnt = [];
    for iStim=1:length(thisCatStim),
        tmpCnt(iStim) = sum([ev(strcmp({ev.stimName},thisCatStim{iStim})).isText] + [ev(strcmp({ev.stimName},thisCatStim{iStim})).isImage]);
    end
    [Y,I] = sort(tmpCnt,2,'descend');
    thisCatStim = thisCatStim(I);
    
    stimList(end+1:end+length(thisCatStim))   = thisCatStim; %- master list of all stimuli... organize this by category
    stimCat(end+1:end+length(thisCatStim))    = thisCat;
    stimCatNum(end+1:end+length(thisCatStim)) = find(strcmp(catList,thisCat));
end


%- create a structure that contains counts for each stimulus, txt and img
stimCnt.numStim     = length(stimList);
stimCnt.numCat      = length(catList);
stimCnt.stimList    = stimList;
stimCnt.stimListNum = [1:stimCnt.numStim];
stimCnt.stimCat     = stimCat;
stimCnt.stimCatNum  = stimCatNum;
for iStim=1:stimCnt.numStim,
    %- all events
    ev = events;
    stimCnt.txt_all(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isText]);
    stimCnt.img_all(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isImage]);
    
    %- just non-targets
    ev = events([events.isTarget]==0);
    stimCnt.txt_nonT(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isText]);
    stimCnt.img_nonT(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isImage]);
    
    %- just targets
    ev = events([events.isTarget]==1);
    stimCnt.txt_tar(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isText]);
    stimCnt.img_tar(iStim) = sum([ev(strcmp({ev.stimName},stimList{iStim})).isImage]);
    
end
fprintf('\n STIMULUS COUNTS...\n      INFO   | ALL Trials | NonTarget | Target\n    cat;  ID;  txt,  img;  txt,  img;  txt,  img\n');
counts = [stimCnt.stimCatNum' stimCnt.stimListNum' stimCnt.txt_all' stimCnt.img_all' stimCnt.txt_nonT' stimCnt.img_nonT' stimCnt.txt_tar' stimCnt.img_tar'];
disp(counts)
disp(sum(counts))


%- All Text vs All Images
ev = events; %- all events
numTxt = sum([ev.isText]);
numImg = sum([ev.isImage]);
fprintf('\n All stim: %d txt, %d img',numTxt,numImg);

%- control for category and attention state
evNonTar = events([events.isTarget]==0);
for thisCat=catList,
    ev = evNonTar(strcmp({evNonTar.stimCategory},thisCat{1}));;
    numTxt = sum([ev.isText]);
    numTxtUnique = length(unique({ev([ev.isText]).stimName}));
    numImg = sum([ev.isImage]);
    numImgUnique = length(unique({ev([ev.isImage]).stimName}));
    fprintf('\n Non-Target %s: %d txt (%d unique); %d img (%d unique)', thisCat{1},numTxt,numTxtUnique,numImg,numImgUnique);
end
%- control for category and attention state
evNonTar = events([events.isTarget]==1);
for thisCat=catList,
    ev = evNonTar(strcmp({evNonTar.stimCategory},thisCat{1}));;
    numTxt = sum([ev.isText]);
    numTxtUnique = length(unique({ev([ev.isText]).stimName}));
    numImg = sum([ev.isImage]);
    numImgUnique = length(unique({ev([ev.isImage]).stimName}));
    fprintf('\n     Target %s: %d txt (%d unique); %d img (%d unique)', thisCat{1},numTxt,numTxtUnique,numImg,numImgUnique);
end



%- Distribution of resopnse times and stimulus durations
figure(1); clf
subplot(211)
hist([events.duration],100)
xlabel('Stim Duration (ms)')

subplot(212)
hist([events.rxnTime],100)
xlabel('Response Time (ms)')



% %- now some sanity checks... make sure expected trial counts are there
% iCore = find(strcmp({events.targetType},'CORE_IMG'));
% events(iCore(1))
% iOwl = find(strcmp({events.stimName},'Owl'));
% length(iOwl)
% iOwl = find(strcmp({events.stimName},'Owl') & [events.isTarget]==1);
% length(iOwl)
%
% iCore = find(strcmp({events.targetType},'CORE_BOTH'));
% length(iCore)
% events(iCore(1))
% iOwl = find(strcmp({events.stimName},'Pigeon') & [events.isTarget]==1);
% sum(strcmp({events.stimName},'Pigeon') & [events.isText]==1 & [events.isTarget]==0)
% sum(strcmp({events.stimName},'Pigeon') & [events.isImage]==1 & [events.isTarget]==0)
% events(iOwl(1))



% %%% CUT FROM SESSION LOG
% 1488912049711	0	B	Logging Begins
% 1488912051109	0	SESSION_START	SESSION_0	data/jwTest2/session_0 stimulusBlast v4.0b
% 1488912051109	0	STIMULUS_LIST	stim_list/stimList_50stim_40blocks_6targets.csv
% 1488912051109	0	SESSION_TYPE	TEST
% 1488912058795	1	BLOCK_PROMPT	0
% 1488912062818	0	BLOCK_START	0
% 1488912063295	1	BLOCK_TARGET 	0 	PERSON
% 1488912065290	1	BLANK_SCREEN
% 1488912066040	1	TXT_ON_SCREEN 	0_0 	NOISE 	MASK 	TEXT 	XXXXXX
% 1488912066792	1	BLANK_SCREEN
% 1488912066792	1	IMG_ON_SCREEN 	0_1 	NOISE 	MASK 	IMAGE 	noise 	images/noise1.jpg
% 1488912067567	1	BLANK_SCREEN
% 1488912067567	1	TXT_ON_SCREEN 	0_2 	ANIMAL 	NONTARGET 	UNIQUE 	Kangaroo
% 1488912068320	1	BLANK_SCREEN
% 1488912068320	1	IMG_ON_SCREEN 	0_3 	PERSON 	TARGET 	UNIQUE 	Matt_Damon 	images/FamousFace/face33.jpg
% 1488912068874	0	KEYBOARD_PRESS 	0_3 	SPACE 	CORRECT(hit) 	507
% 1488912069118	1	BLANK_SCREEN
% 1488912069118	1	TXT_ON_SCREEN 	0_4 	ANIMAL 	NONTARGET 	CORE_TXT 	Frog
% 1488912069869	1	BLANK_SCREEN
% 1488912069869	1	IMG_ON_SCREEN 	0_5 	ANIMAL 	NONTARGET 	CORE_IMG 	Owl 	images/Animal/animal5.jpg
% 1488912070674	1	BLANK_SCREEN
% 1488912070674	1	TXT_ON_SCREEN 	0_6 	OBJECT 	NONTARGET 	UNIQUE 	Laptop
% 1488912071427	1	BLANK_SCREEN
% 1488912071427	1	IMG_ON_SCREEN 	0_7 	ANIMAL 	NONTARGET 	CORE_BOTH 	Pigeon 	images/Animal/animal10.jpg
% 1488912072212	1	BLANK_SCREEN
% ...
% 1488912096905	1	TXT_ON_SCREEN 	0_40 	PLACE 	NONTARGET 	CORE_TXT 	Great_Wall_of_China
% 1488912097659	1	BLANK_SCREEN
% 1488912097659	1	IMG_ON_SCREEN 	0_41 	OBJECT 	NONTARGET 	CORE_IMG 	Cup 	images/Object/object6.jpg
% 1488912098444	1	BLANK_SCREEN
% 1488912098444	1	TXT_ON_SCREEN 	0_42 	PERSON 	TARGET 	CORE_RAND 	Emma_Watson
% 1488912099068	0	KEYBOARD_PRESS 	0_42 	SPACE 	CORRECT(hit) 	619
% 1488912099198	1	BLANK_SCREEN
% 1488912099198	1	IMG_ON_SCREEN 	0_43 	OBJECT 	NONTARGET 	CORE_IMG 	Scissors 	images/Object/object8.jpg
% 1488912099982	1	BLANK_SCREEN
% 1488912099982	1	TXT_ON_SCREEN 	0_44 	NOISE 	MASK 	TEXT 	XXXXXX
% 1488912100734	1	BLANK_SCREEN
% 1488912100734	1	IMG_ON_SCREEN 	0_45 	NOISE 	MASK 	IMAGE 	noise 	images/noise1.jpg
% 1488912101507	1	BLANK_SCREEN
% 1488912103448	1	REC_START 	0 	BLOCK_0_148891210
% 1488912103451	1	NAMING_QUERY 	0 	PERSON
% 1488912114649	1	COUNTDOWN_DONE 	0
% 1488912115648	1	REC_STOP 	0 	BLOCK_0_148891210
% 1488912115654	1	BLOCK_SUMMARY 	0 	hits 6	 mis 0	 fa 2	 cr 38
% 1488912118648	0	BLOCK_DURATION 	0 	55.8 sec
% 1488912118648	0	BLOCK_END	0
% 1488912118653	1	BLOCK_PROMPT	1
% 1488912120211	0	BLOCK_START	1
% 1488912120710	1	BLOCK_TARGET 	1 	PLACE
% 1488912122706	1	BLANK_SCREEN
% 1488912123456	1	TXT_ON_SCREEN 	1_0 	NOISE 	MASK 	TEXT 	XXXXXX
% 1488912124208	1	BLANK_SCREEN
% 1488912124208	1	IMG_ON_SCREEN 	1_1 	NOISE 	MASK 	IMAGE 	noise 	images/noise1.jpg
% 1488912124981	1	BLANK_SCREEN
% 1488912124981	1	TXT_ON_SCREEN 	1_2 	PERSON 	NONTARGET 	CORE_TXT 	Elvis_Presley
% 1488912125735	1	BLANK_SCREEN
% 1488912125735	1	IMG_ON_SCREEN 	1_3 	ANIMAL 	NONTARGET 	CORE_IMG 	Cat 	images/Animal/animal6.jpg
% 1488912126520	1	BLANK_SCREEN
% 1488912126520	1	TXT_ON_SCREEN 	1_4 	PLACE 	TARGET 	UNIQUE 	Lincoln_Memorial
% 1488912127092	0	KEYBOARD_PRESS 	1_4 	SPACE 	CORRECT(hit) 	567
% 1488912127274	1	BLANK_SCREEN
% ...
% 1488912156725	1	IMG_ON_SCREEN 	1_43 	ANIMAL 	NONTARGET 	CORE_IMG 	Gorilla 	images/Animal/animal8.jpg
% 1488912157511	1	BLANK_SCREEN
% 1488912157511	1	TXT_ON_SCREEN 	1_44 	NOISE 	MASK 	TEXT 	XXXXXX
% 1488912158263	1	BLANK_SCREEN
% 1488912158263	1	IMG_ON_SCREEN 	1_45 	NOISE 	MASK 	IMAGE 	noise 	images/noise1.jpg
% 1488912159036	1	BLANK_SCREEN
% 1488912160779	1	REC_START 	1 	BLOCK_1_148891216
% 1488912160782	1	NAMING_QUERY 	1 	PLACE
% 1488912171980	1	COUNTDOWN_DONE 	1
% 1488912172979	1	REC_STOP 	1 	BLOCK_1_148891216
% 1488912172985	1	BLOCK_SUMMARY 	1 	hits 6	 mis 0	 fa 0	 cr 40
% 1488912175979	0	BLOCK_DURATION 	1 	55.8 sec
% 1488912175979	0	BLOCK_END	1
% 1488912175979	0	TASK_DONE
% 1488912178979	0	SESS_DURATION_TOTAL 2.1 min
% 1488912178979	0	SESS_END
% 1488912179009	0	E	Logging Ends
%
