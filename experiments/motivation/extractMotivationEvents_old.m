function events = extractMotivationEvents(sessLogFile, subject, sessionName)
% EXTRACTMOTIVATIONEVENTS extracts one event per trial of a Motivation session
%
% events = extractMotivationEvents(sessionLogFile, subject, sessionName)
%
% This old version creates 1 event per trial. The new version does 1 event
% per line
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% uncomment following lines to directly run script
% clear all
%
% rootEEGdir  = '/Users/wittigj/DataJW/AnalysisStuff/dataLocal/eeg';
% rootEEGdir  = '/Volumes/Shares/FRNU/dataWorking/eeg';
% subject     = 'NIH031';   % EEG002  NIH016
% sessionName = 'session_1';
%
% sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/paRemap',sessionName);
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'motivation';

fid = fopen(sessLogFile, 'r'); %if loops check to see that sessFile exists
if (fid == -1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric)) > 1
    iKeep=1:find(diff(strNumeric)>1,1,'first'); 
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', ...
        sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); 
    strNumeric=strNumeric(iKeep); 
end
sessionNum = str2double(sessionName(strNumeric));               
if isempty(sessionNum)
    fprintf('\n ERROR: problem converting session name into a numeric'); 
    keyboard;  
end; %shouldn't need this catch...


%- Read session.log line-by-line and convert to events structure
index = 0;
events = [];
word = '';
type = '';

while true
    thisLine            = fgetl(fid); % get new line
    if ~ischar(thisLine); break; end  % reached EOF
    
    pieces = strsplit(thisLine);
    timestamp = str2double(pieces{1});
    clockType = str2double(pieces{2});
    eventType = pieces{3};
    eventInfo = pieces(4:end);
    
    switch eventType
        case 'SESS_START'
            % INFO: session number (0-based)
            assert(str2double(eventInfo{1}) == sessionNum,...
                'session number given (%s) does not match log file (%d)', eventInfo{1}, sessionNum);
            
        case 'BLOCK_START'
            % INFO: block number (0-based)
            block = str2double(eventInfo{1});
            
        case 'BLANK'
            % INFO: trial
            blankBeginTime = timestamp;
    
            % signal end of study event
            type = 'STUDY';
            if ~isempty(word)
                index = index + 1;
            end

        case 'MASK_PRESENT'
            % INFO: trial
            maskTime = timestamp;
            
        case 'VAL_PRESENT'
            % INFO: trial, value
            value = str2double(eventInfo{2});
            valueTime = timestamp;
            
        case 'WORD_PRESENT'
            % INFO: trial, word
            trial = str2double(eventInfo{1});
            word = eventInfo{2};
            wordTime = timestamp;

        case 'WORD_TEST'
            % INFO: [SHOWN/FOIL], word
            isFoil = strcmpi(eventInfo{1}, 'FOIL');
            word = eventInfo{2};
            wordTime = timestamp;
            
        case 'RESPONSE'
            % INFO: [SEEN/UNSEEN/PASS], points
            response = eventInfo{1};
            pointsGiven = str2double(eventInfo{2});
            responseTime = timestamp;
            
            % signal end of test event
            index = index + 1;
            type = 'TEST';
            
    end
    
    
    if index > length(events)        
       
        % initialize blank structure
        fields = {...
            'experiment';
            'subject';
            'sessionName';
            'sessionNum';
            'isClock';
            'mstime';
            'blockNum';
            'type';
            'word';
            'wordTime';
            'responseTime';
            'pointsGiven';
            'response';
            'isFoil';
            'value';
            'valueTime';
            'maskTime';
            'blankBeginTime';
            'trialNum'};
        thisEvent = struct();
        for i = 1 : length(fields)
            thisEvent.(fields{i}) = [];
        end
            
        % set the structure
        thisEvent.experiment        = experiment  ;
        thisEvent.subject           = subject     ;
        thisEvent.sessionName       = sessionName ;
        thisEvent.sessionNum        = sessionNum  ;
        thisEvent.isClock           = clockType     ;
        thisEvent.blockNum          = block       ;
        thisEvent.type              = type        ;
        
        if strcmpi(type, 'TEST')
            thisEvent.word = word;
            thisEvent.wordTime = wordTime;
            thisEvent.responseTime = responseTime;
            thisEvent.pointsGiven = pointsGiven;
            thisEvent.response = response;
            thisEvent.isFoil = isFoil;
            thisEvent.mstime = wordTime;
            
            
        elseif strcmpi(type, 'STUDY')
            thisEvent.word = word;
            thisEvent.maskTime = maskTime;              % first
            thisEvent.valueTime = valueTime;            % second
            thisEvent.wordTime = wordTime;              % third
            thisEvent.blankBeginTime = blankBeginTime;  % fourth
            thisEvent.value = value;
            thisEvent.trialNum = trial;
            thisEvent.mstime = valueTime;
            assert(valueTime < wordTime && wordTime < blankBeginTime,...
                'Word %s has timing:\n%d - value\n%d - word\n%d - blank',...
                word, valueTime, wordTime, blankBeginTime);
        else
            error('Index incremented unexpectedly')
        end
            
        
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent; %#ok<*AGROW>
        end
    end
    
end
fclose(fid);  % close session.log