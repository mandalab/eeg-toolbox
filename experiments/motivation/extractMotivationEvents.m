function events = extractMotivationEvents(sessLogFile, subject, sessionName)
% EXTRACTMOTIVATIONEVENTS extracts one event per line of a Motivation session
%
% events = extractMotivationEvents(sessionLogFile, subject, sessionName)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% uncomment following lines to directly run script
% clear all
%
% rootEEGdir  = '/Users/trottams/NIH/dataWorking/AnalysisStuff/dataLocal/eeg';
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
events = [];
word = '';
phase = '';
trial = -1;

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
            block = 0;
            assert(str2double(eventInfo{1}) == sessionNum,...
                'session number given (%s) does not match log file (%d)', eventInfo{1}, sessionNum);
            continue;
            
        case 'BLOCK_START'
            % INFO: block number (0-based)
            block = str2double(eventInfo{1});
            continue;
            
        case 'BLANK' % end trial
            % INFO: trial
            phase = 'STUDY';
            

        case 'MASK_PRESENT' % start trial
            % INFO: trial
            phase = 'STUDY';
            word = '';
            value = [];
            trial = trial + 1;
            
        case 'VAL_PRESENT'
            % INFO: trial, value
            value = str2double(eventInfo{2});
            phase = 'STUDY';
            
        case 'WORD_PRESENT' 
            % INFO: trial, word
            
            trial = str2double(eventInfo{1});
            word = eventInfo{2};
            phase = 'STUDY';
            
        case 'WORD_TEST' % start of test trial
            % INFO: [SHOWN/FOIL], word
            
            isFoil = strcmpi(eventInfo{1}, 'FOIL');
            word = eventInfo{2};
            phase = 'TEST';
            pointsGiven = [];
            response = '';
            soundType = '';
            
        case 'RESPONSE'
            % INFO: [SEEN/UNSEEN/PASS], points
            
            response = eventInfo{1};
            pointsGiven = str2double(eventInfo{2});
            phase = 'TEST';
            
        case 'SOUND'
            % INFO: [CORRECT/INCORRECT]
            phase = 'TEST'; 
            soundType = eventInfo{1};
            
        case 'SCORE_PRESENT'
            % INFO: +/- points given and shown
            phase = 'TEST';
            pointsGiven = str2double(eventInfo{1});
            
        otherwise
            %disp(thisLine);
            continue;
    end
    
    

    % initialize blank structure
    fields = {...
        'experiment';
        'subject';
        'sessionName';
        'sessionNum';
        'clockType';
        'mstime';
        'blockNum';
        'type';
        'phase';
        'word';
        'pointsGiven';
        'response';
        'isFoil';
        'value';
        'trialNum';
        'soundType'};
    thisEvent = struct();
    for i = 1 : length(fields)
        thisEvent.(fields{i}) = [];
    end

    % set the structure
    thisEvent.experiment        = experiment  ;
    thisEvent.subject           = subject     ;
    thisEvent.sessionName       = sessionName ;
    thisEvent.sessionNum        = sessionNum  ;
    thisEvent.clockType         = clockType   ;
    thisEvent.blockNum          = block       ;
    thisEvent.mstime            = timestamp   ;
    thisEvent.type              = eventType   ;
    thisEvent.phase             = phase       ;
    
    if strcmpi(phase, 'TEST')
        thisEvent.word = word;
        thisEvent.pointsGiven = pointsGiven;
        thisEvent.response = response;
        thisEvent.isFoil = isFoil;
        thisEvent.soundType = soundType;

    elseif strcmpi(phase, 'STUDY')
        thisEvent.word = word;
        thisEvent.value = value;
        thisEvent.trialNum = trial;
        
    end


    events = [events, thisEvent]; %#ok<*AGROW>

    
end
fclose(fid);  % close session.log