function [events] = stimGuess_ExtractEvents_v1(sessLogFile, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from stimulusGuess %%%%
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
% subject     = 'jwReal';   % EEG002  NIH016
% %sessionName = 'session_0_ninetyBlur';
% sessionName = 'session_2_sixteyBlur';
% %sessionName = 'session_4_sixteyBlur+2';
% %sessionName = 'session_5_sixteyBlur-2 newLog';
% sessionDir  = fullfileEEG('../data/',subject,sessionName);
%
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sessLogFile, subject, sessionName

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
%if exist('sessionDir','var') && MAKE_SESS_KEY_LOG == 1,
    
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
            mstime2  = str2num(thisLineStr{1}) + str2num(thisLineStr{2}); %- add the error term
            
            %- fair way to incorporate error.  Actual event happens sometime between first and second number.  So lets just take 1/2 of second number
            mstime3 = floor(str2num(thisLineStr{1}) + str2num(thisLineStr{2})/2.0);
            
            
            %-only include SPACE bar "P"resses
            if strcmp(thisLineStr{3},'P') &  (strcmp(thisLineStr{4},'LEFT') | strcmp(thisLineStr{4},'RIGHT') | strcmp(thisLineStr{4},'UP') | strcmp(thisLineStr{4},'DOWN'))
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
sessEvents  = {'B','SESSION_START','STIMULUS_LIST','BLUR_DIR','SESSION_TYPE','SESSION_TIMING','SESS_DURATION_TOTAL','SESS_END','E'};
blockEvents = {'BLOCK_PROMPT','BLOCK_START','BLOCK_FEEDBACK','BLOCK_SUMMARY','CATEG_SUMMARY','CATEG_SUMMARY2','BLOCK_DURATION','BLOCK_END'};

stimEvents  = {'IMG_ON_SCREEN','IMG_OFF_SCREEN','BLANK_SCREEN','KEYBOARD_PRESS'};
stimCats    = {'NOISE','ANIMAL','OBJECT','PERSON','PLACE','BUSH','CLINTON','KENNEDY','TRUMP'};
respTypes   = {'CORRECT(hit)','INCORRECT(fa)','INCORRECT(skip)'};

%- key for conversion from string to number
catStr2Num  = {'ANIMAL', 'OBJECT', 'PERSON', 'PLACE','BUSH','CLINTON','KENNEDY','TRUMP','','NOISE'}; %- make 'NOISE' a 10
catStr2NumB = {'ANIMAL', 'OBJECT', 'PERSON', 'PLACE','BUSH','CLINT','JFK','TRUMP'}; %- make 'NOISE' a 10

%- presummed mapping between response key and response category... confirmed below, and used if key_log response differs from keyBoard
respCatList   = {'DOWN(ANIMAL)', 'RIGHT(OBJECT)', 'UP(PERSON)', 'LEFT(PLACE)','DOWN(BUSH)','UP(CLINTON)','LEFT(KENNEDY)','RIGHT(TRUMP)'};
respCatListV5 = {'DOWN(ANIMAL)', 'RIGHT(OBJECT)', 'UP(PERSON)', 'LEFT(PLACE)','DOWN(BUSH)','LEFT(CLINTON)','UP(KENNEDY)','RIGHT(TRUMP)'}; %- old version V5


%- predefine some globals that will get overwritten below if newer version of session.log
TIME_BEFORE_BUTTON = 200; %- assume that task rejects button presses before 200ms
thisSessIsZap = 0; %- assume non-stim session

%- Read session.log line-by-line and convert to events structure
moreToRead    = true;
events        = [];
index         = 1;
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
    %isBlockTarget = 0;
    thisRespMade  = nan; %- was the space bar pressed or not
    thisRespCorr  = nan; %- was the press correct
    %-
    switch eventType
        
        
        case {'BLUR_DIR'}
            thisSessBlurDir = thisLineStr{4};
            if contains(thisSessBlurDir,'ninetyBlur'),
                isNinteyBlurSess = 1;  %- blur levels are set at 16 (17 for exemplars)
                thisSessBlurLevel = 'ninetyBlur';
            else
                isNinteyBlurSess = 0;  %- expect to read the blur level from each image
                if contains(thisSessBlurDir,'->'), %- V5 session hack
                    thisSessBlurDir = 'multiBlur'; thisSessBlurLevel=5; %- call it 5... but really it is a mix of 0 and 5
                else
                    thisSessBlurDir = sprintf('%s,%s,%s',thisLineStr{[4:6]});
                    thisSessBlurLevel = thisLineStr{5};
                end
            end
            
            
        case {'SESSION_TYPE'}
            thisSessType   = thisLineStr{4}; % TRAIN_SLOW, TRAIN_FAST, TEST
            
            
        case {'SESSION_TIMING'}
            iStr = find(contains(thisLineStr,'TIME_BEFORE_BUTTON'));
            if length(iStr==1),
                TIME_BEFORE_BUTTON = str2num(thisLineStr{iStr}(20:end));
                if VERBOSE, fprintf('\n got TIME_BEFORE_BUTTON from session log = %d',TIME_BEFORE_BUTTON'); end
            else
                fprintf('\n SESSION_TIMING row does not contain "TIME_BEFORE_BUTTON".  weird.');
                keyboard;
            end
            
            
        case {'SESSION_STIM=True'}
            thisSessIsZap = 1;
        
        case {'SESSION_STIM'}
            thisSessIsZap  = strcmp(thisLineStr{4},'IS_STIM_SESS=True');
            
            
        case {'IMG_ON_SCREEN'}
            %- catch from old version... shouldnt happen much
            if strcmp(thisLineStr{5},'REMOVE_IMAGE')
                %- this happend for a few sessions for NIH064 and NIH066 preOp before I updated the task script...
                %- update duration of previous stimulus here because "IMAGE_OFF" doesn't exist in this log
                if length(events)>0 && isnan(events(end).durationImg),
                    events(end).durationImg = mstime - events(end).mstime;
                end
                continue;
            end
            
            %- update duration of pervious stimulus
            if length(events)>0,
                if isnan(events(end).durationTrial) events(end).durationTrial = mstime - events(end).mstime; end %- should always be true
                if isnan(events(end).durationImg)     events(end).durationImg = mstime - events(end).mstime; end %- if flashed stimulus this will not be true
            end
            
            
            makeNewEvent   = 1; %-increment event counter
            thisBlockStr   = thisLineStr{4}; % 0_22, 30_0, etc
            thisCategory   = thisLineStr{5}; % NOISE, ANIMAL, OBJECT, PERSON, PLACE
            thisCategoryNum = find(strcmp(thisCategory,catStr2Num));
            
            
            %thisTargetStr  = thisLineStr{6}; % TARGET, NONTARGET
            %thisRepeatType = thisLineStr{7}; % CORE_IMG, CORE_TXT, CORE_FIXED, CORE_RAND, UNIQUE
            thisStimName   = thisLineStr{8}; % Dog, White House, etc
            stimIDstr      = thisLineStr{9};
            if isNinteyBlurSess | strcmp(thisCategory, 'NOISE')
                thisStimIDnum  = thisCategoryNum*100 + str2num(stimIDstr(find( stimIDstr >= '0' & stimIDstr <= '9')));
                thisStimBlur = 15;  %- for NIH062 initial
                
                %% TRUTH FOR NIH062 initial eCog sessions and NIH063 preop... after that it was 15.
                %if strcmp(thisCategory,'ANIMAL') | strcmp(thisCategory,'OBJECT') | strcmp(thisCategory,'PERSON') | strcmp(thisCategory,'PLACE'),
                %    thisStimBlur = 16;
                %else
                %    thisStimBlur = 17; %- exemplars: presidends
                %end
                
            else
                stimIDstrPreBlur = stimIDstr(1:end-6);
                thisStimIDnum  = thisCategoryNum*100 + str2num(stimIDstrPreBlur(find( stimIDstrPreBlur >= '0' & stimIDstrPreBlur <= '9')));
                thisStimBlur = str2num(stimIDstr(end-5:end-4));
            end
            
            thisIsImage    = strcmp(eventType(1:3),'IMG');
            thisIsStimulus = ~strcmp(thisCategory, 'NOISE');
            thisIsNoise    = strcmp(thisCategory, 'NOISE');
            
            iHyphen        = strfind(thisBlockStr,'_');
            thisBlockNum   = str2num(thisBlockStr(1:iHyphen-1));
            thisStimNum    = str2num(thisBlockStr(iHyphen+1:end));
            
            thisDurationImg   = nan;
            thisDurationTrial = nan;
            thisRxnTime       = nan;
            
            %- assume no key-press... that will be overwritten later if one happens
            thisRespBut    = 'NONE';
            thisRespCat    = 'NONE'; %- assume no response, which is an unlabeled correct rejection for non-targets'
            thisRespCatNum = nan;
            thisRespMade   = 0; %-
            thisRespCorr   = 0; %-
            
            thisKeyLogBut  = 'NONE'; %- placeholders for keylog values
            thisKeyLogRT   = nan;
            
            %- zap (electrical stimulation) stuff
            if thisSessIsZap,
                thisZapType    = thisLineStr{10};
                thisZapTarg    = thisLineStr{11};
                thisZapTargNum = str2num(thisLineStr{12}(2:end-1));
            else
                thisZapType    = 'n/a';
                thisZapTarg    = 'n/a';
                thisZapTargNum = [];
            end
            
            
        case {'IMG_OFF_SCREEN'}
            %- update duration of pervious stimulus
            if length(events)>0 && isnan(events(end).durationImg),
                events(end).durationImg = mstime - events(end).mstime;
            end
            
            
        case 'KEY_LOG_PRESS'
            %- in stimBlast we needed to use key_log to catch existance and timing of all presses... maybe not with stimGuess
            if length(events)>0 && isnan(events(end).keyLogRxnTime),  %- dont use keystrock from block prompt in events.mat, and dont overwrite previous entry
                
                %- check RT against "TIME_BEFORE_BUTTON"... dont log it here if too quick
                keyLogRT = str2num(thisLineStr{1})-events(end).mstime;
                if keyLogRT>TIME_BEFORE_BUTTON,
                    events(end).keyLogButton    = thisLineStr{4};  %- save a copy of the whole line
                    events(end).keyLogRxnTime   = str2num(thisLineStr{1})-events(end).mstime;
                    
                    %- sometimes keyboard_press comes before key_log... if so do correction now
                    if events(end).isRespMade,
                        if VERBOSE, fprintf('\n key log after keyboard press: event %d, \n%s',length(events),events(end).logStrResp); end
                        
                        if ~strcmp(events(end).respButton, events(end).keyLogButton),
                            %- need to update: respCategory, respCatNum, and isRespCorr
                            if VERBOSE, fprintf('\n updating response from keylog value, event %d: %s --> %s  << from keylog AFTER keypress\n%s',length(events),events(end).respButton,events(end).keyLogButton,thisLine); end
                            events(end).respButton      = events(end).keyLogButton;
                            thisRespStr = respCatList(find(contains(respCatList,events(end).respButton)));
                            if events(end).stimCategoryNum <= 4,
                                %-categories
                                thisRespStr = thisRespStr(find(contains(thisRespStr,catStr2Num{1})|contains(thisRespStr,catStr2Num{2})|contains(thisRespStr,catStr2Num{3})|contains(thisRespStr,catStr2Num{4})));
                            else
                                %-presidents
                                thisRespStr = thisRespStr(find(contains(thisRespStr,catStr2Num{5})|contains(thisRespStr,catStr2Num{6})|contains(thisRespStr,catStr2Num{7})|contains(thisRespStr,catStr2Num{8})));
                            end
                            if length(thisRespStr)==1, thisRespStr = thisRespStr{1};
                            else fprintf('\n uh oh... didnt return exactly 1 match'); keyboard; end
                            events(end).respCategory    = thisRespStr(strfind( thisRespStr,'(')+1:end-1);
                            events(end).respCategoryNum = find(strcmp(events(end).respCategory,catStr2Num));
                            events(end).isRespCorr      = strcmp(events(end).respCategory,events(end).stimCategory);
                        end
                        if events(end).rxnTime ~= events(end).keyLogRxnTime,
                            if VERBOSE, fprintf('\n updating rxn time from keylog value, event %d:  %d --> %d  << from keylog AFTER keypress', length(events), events(end).rxnTime, events(end).keyLogRxnTime); end
                            events(end).rxnTime = events(end).keyLogRxnTime;
                        end
                    end
                end %- if >TIME_BEFORE_BUTTON
            end
            
            
        case 'KEYBOARD_PRESS'
            %- ignore these... switch to using the "key log" press entries... those seem more accurate and dont get lost
            %events(end).response = thisLineStr{6};
            
            if sum(strcmp(respCatList,thisLineStr{5})) ~= 1,
                respCatList = respCatListV5;
                if sum(strcmp(respCatList,thisLineStr{5})) == 1,
                    fprintf('\n respCatList didnt match new version but good with V5');
                else
                    fprintf('\n respCatList not matching session.log (%s)... clear the cat list and re-populate',thisLineStr{5});
                    keyboard;
                    if length(respCatList)==8,respCatList={}; end
                    respCatList{end+1} = thisLineStr{5};
                end
            end
            
            events(end).logStrResp      = thisLine;  %- save a copy of the whole line
            events(end).respButton      = thisLineStr{5}(1:strfind( thisLineStr{5},'(')-1);
            events(end).respCategory    = thisLineStr{5}(strfind( thisLineStr{5},'(')+1:end-1);
            events(end).respCategoryNum = find(strcmp(events(end).respCategory,catStr2Num));
            events(end).isRespMade      = 1;
            events(end).isRespCorr      = strcmp(events(end).respCategory,events(end).stimCategory);
            events(end).rxnTime         = str2num(thisLineStr{7});
            
            if strcmp(thisLineStr{6},'CORRECT(hit)') ~= events(end).isRespCorr,
                fprintf('\n log "CORRECT" does not match computed correct... look into this JW');
                keyboard;
            end
            
            
            %- confirm that key_log_press VALUE and TIMING match with keyboard press... if not update values recorded above
            if ~isnan(events(end).keyLogRxnTime),
                if ~strcmp(events(end).respButton, events(end).keyLogButton),
                    %- need to update: respCategory, respCatNum, and isRespCorr
                    if VERBOSE, fprintf('\n updating response from keylog value, event %d: %s --> %s\n%s',length(events),events(end).respButton,events(end).keyLogButton,thisLine); end
                    events(end).respButton      = events(end).keyLogButton;
                    thisRespStr = respCatList(find(contains(respCatList,events(end).respButton)));
                    if events(end).stimCategoryNum <= 4,
                        %-categories
                        thisRespStr = thisRespStr(find(contains(thisRespStr,catStr2Num{1})|contains(thisRespStr,catStr2Num{2})|contains(thisRespStr,catStr2Num{3})|contains(thisRespStr,catStr2Num{4})));
                    else
                        %-presidents
                        thisRespStr = thisRespStr(find(contains(thisRespStr,catStr2Num{5})|contains(thisRespStr,catStr2Num{6})|contains(thisRespStr,catStr2Num{7})|contains(thisRespStr,catStr2Num{8})));
                    end
                    if length(thisRespStr)==1, thisRespStr = thisRespStr{1};
                    else fprintf('\n uh oh... didnt return exactly 1 match'); keyboard; end
                    events(end).respCategory    = thisRespStr(strfind( thisRespStr,'(')+1:end-1);
                    events(end).respCategoryNum = find(strcmp(events(end).respCategory,catStr2Num));
                    events(end).isRespCorr      = strcmp(events(end).respCategory,events(end).stimCategory);
                end
                if events(end).rxnTime ~= events(end).keyLogRxnTime,
                    if VERBOSE, fprintf('\n updating rxn time from keylog value, event %d:  %d --> %d', length(events), events(end).rxnTime, events(end).keyLogRxnTime); end
                    events(end).rxnTime = events(end).keyLogRxnTime;
                end
            end
            
            
        case 'BLANK_SCREEN'
            if length(events)>0,
                if isnan(events(end).durationTrial) events(end).durationTrial = mstime - events(end).mstime; end %- should always be true
                if isnan(events(end).durationImg)     events(end).durationImg = mstime - events(end).mstime; end %- if flashed stimulus this will not be true
            end
            
            
            
        otherwise
            %- nothing to process on this line
            %fprintf('\n
    end
    
    
    
    %- if this is a new stimulus, create an event...
    for ii=1:makeNewEvent,
        
        clear thisEvent
        thisEvent.experiment        = 'stimulusGuess';
        thisEvent.subject           = subject;
        thisEvent.sessionName       = sessionName;
        thisEvent.sessionBlurDir    = thisSessBlurDir;
        thisEvent.sessionBlurVal    = thisSessBlurLevel;
        thisEvent.sessionNum        = sessionNum;
        thisEvent.sessionType       = thisSessType;
        thisEvent.sessionIsZap      = thisSessIsZap;
        thisEvent.mstime            = mstime;
        thisEvent.eventType         = eventType; %
        thisEvent.logStr            = thisLine;
        thisEvent.stimFileName      = stimIDstr;
        thisEvent.stimNickName      = thisStimName;
        thisEvent.stimIDnum         = thisStimIDnum;     %- 1 through 60 * 2 * 4 (numStim * img vs text * num categories)
        thisEvent.stimBlurVal       = thisStimBlur;
        %thisEvent.targetCategory    = targetCategory;    %- target for this block
        %thisEvent.targetCategoryNum = targetCategoryNum; %-
        thisEvent.stimCategory      = thisCategory;     %- category of this stimulus
        thisEvent.stimCategoryNum   = thisCategoryNum;  %- [nan, 0,1,2,3,4 -->  nonStim, noise,animal,object,person,place]
        %thisEvent.targetStr         = thisTargetStr;    %- was this stimulus target or not?
        %thisEvent.repeatType        = thisRepeatType;   %- CORE_BOTH CORE_IMG CORE_TXT UNIQUE
        thisEvent.blockNum          = thisBlockNum;
        thisEvent.stimNum           = thisStimNum;    %- presentation order within the list
        thisEvent.durationImg       = thisDurationImg;
        thisEvent.durationTrial     = thisDurationTrial;
        thisEvent.logStrResp        = '';
        thisEvent.respButton        = thisRespBut;     %- string describing response button ('UP','DOWN', etc)
        thisEvent.respCategory      = thisRespCat;     %- string describing response category ('ANIMAL', etc)
        thisEvent.respCategoryNum   = thisRespCatNum;     %- string describing response category ('ANIMAL', etc)
        thisEvent.isRespMade        = thisRespMade;     %- was the space bar pressed or not
        thisEvent.isRespCorr        = thisRespCorr;     %- was the response correct or not
        thisEvent.rxnTime           = thisRxnTime;
        thisEvent.keyLogButton      = thisKeyLogBut;   %- placeholder for keylog data...
        thisEvent.keyLogRxnTime     = thisKeyLogRT;
        thisEvent.zapTypeStr        = thisZapType;
        thisEvent.zapTargStr        = thisZapTarg;
        thisEvent.zapTargNum        = thisZapTargNum;
        
        %- booleans for selecting events
        thisEvent.isImage           = thisIsImage;
        thisEvent.isStimulus        = thisIsStimulus;   %- image or text that is NOT noise
        thisEvent.isNoise           = thisIsNoise;
        
        
        
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

keyboard

%- quick plots of some sanity checks.  RXN Times, Trial counts per stim and per category

%- Distribution of resopnse times and stimulus durations
figure(1); clf
subplot(212)
hist([events.durationImg],100)
xlabel('Stim Duration (ms)')

subplot(211)
hist([events.rxnTime],100)
xlabel('Response Time (ms)')


%- distribution of trial counts vs image and vs category
figure(2); clf
subplot(211)
hist([events.stimCategoryNum],100)
xlabel('Category Num')

subplot(212)
hist([events.stimIDnum],5000)
xlabel('Stimulus ID Num');


%- image-by-image measure of response AND rxnTime (two separate color plots)
blocks  = unique( [events.blockNum]  );  %- create "metaBlockNum" when event is multiple sessions combined
stimIDs = unique( [events.stimIDnum] );
respByStim = nan( length(blocks), length(stimIDs)); %- max stimIDs shoud be 1000 (noise)... compress out empty columns later
rxntByStim = respByStim;
blurByStim = respByStim;
for iB=1:length(blocks),
    evList = events([events.blockNum]==blocks(iB));
    
    for iEv=1:length(evList),
        iStim = find(stimIDs == evList(iEv).stimIDnum);
        respByStim(iB,iStim) = evList(iEv).respCategoryNum;
        rxntByStim(iB,iStim) = evList(iEv).rxnTime;
        blurByStim(iB,iStim) = evList(iEv).stimBlurVal;
    end
end


fS = 20; strTitle = sprintf('%s(%s)',subject,sessionName); strTitle(find(strTitle=='_')) = ' ';
xLine = 0.5+[0:60:240 240+[15:15:60]];


figure(3); clf; set(gcf,'color','w');
%subplot(211)
imagesc(1:length(stimIDs),blocks,respByStim,[0 8])
set(gca,'xtick',[30 90 150 210 240+[7:15:60]],'xticklabel',catStr2NumB(1:8),'tickdir','out','box','off','fontsize',fS); hold on;
colormap([1 1 1; lines(8)]);
colorbar
title(strTitle);
for iL=1:length(xLine), plot([1 1]*xLine(iL),get(gca,'ylim'),'k--','linewidth',2); end

figure(4); clf; set(gcf,'color','w');
%subplot(212)
imagesc(1:length(stimIDs),blocks,rxntByStim,[0 3000])
set(gca,'xtick',[30 90 150 210 240+[7:15:60]],'xticklabel',catStr2NumB(1:8),'tickdir','out','box','off','fontsize',fS); hold on;
colorbar
colormap([1 1 1; parula(100)]);
for iL=1:length(xLine), plot([1 1]*xLine(iL),get(gca,'ylim'),'k--','linewidth',2); end

figure(5); clf; set(gcf,'color','w');
%subplot(212)
imagesc(1:length(stimIDs),blocks,blurByStim,[0 10])
set(gca,'xtick',[30 90 150 210 240+[7:15:60]],'xticklabel',catStr2NumB(1:8),'tickdir','out','box','off','fontsize',fS); hold on;
colorbar
colormap([1 1 1; parula(11)]);
for iL=1:length(xLine), plot([1 1]*xLine(iL),get(gca,'ylim'),'k--','linewidth',2); end

%
%




return;

%- enumerate possible categories: PEOPLE, PLACES, OBJECTS, ANIMALS
ev = events; %- all stimuli
%catList  = unique({ev.stimCategory}); %- mucking up order set above...
catList = catStr2Num;

%- organize the stimuli base on categories
stimList = {};  stimCat = {};  stimCatNum = [];
for thisCat=catList,
    thisCatStim = unique({events(strcmp({ev.stimCategory},thisCat)).stimFileName});
    
    tmpCnt = [];
    for iStim=1:length(thisCatStim),
        tmpCnt(iStim) = sum([ev(strcmp({ev.stimFileName},thisCatStim{iStim})).isImage]);
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
    stimCnt.img_all(iStim) = sum([ev(strcmp({ev.stimFileName},stimList{iStim})).isImage]);
end

%- this is confusing because the stimCatNum is no
fprintf('\n STIMULUS COUNTS...\n      INFO   | ALL Trials | NonTarget | Target\n    cat;  ID;  txt,  img;  txt,  img;  txt,  img\n');
counts = [stimCnt.stimCatNum' stimCnt.stimListNum' stimCnt.img_all' ];
disp(counts)
disp(sum(counts))


%- All Text vs All Images
ev = events; %- all events
numImg = sum([ev.isImage]);
fprintf('\n All stim: %d img',numImg);


evNonTar = events;
for thisCat=catList,
    ev = evNonTar(strcmp({evNonTar.stimCategory},thisCat{1}));;
    numImg = sum([ev.isImage]);
    numImgUnique = length(unique({ev([ev.isImage]).stimFileName}));
    fprintf('\n     Target %s: %d img (%d unique)', thisCat{1},numImg,numImgUnique);
end


