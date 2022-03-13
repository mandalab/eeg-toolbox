function events = extractAVMultimodalEvents(rootEEGdir, subject, sessionName, priorEvents)

%rootEEGdir  = '/Volumes/shares/FRNU/data/eeg';
%rootEEGdir='/Users/sreekumarv/Documents/MyCode/theDoors';
% rootEEGdir = '/Users/sreekumarv/Desktop/localDataWorking';
% subject     = 'NIH046';   % 
% sessionName = 'session_0';

sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/multiModalAV',sessionName);
sessLogFile = fullfileEEG(sessionDir,'session.log');
eventFile   = fullfileEEG(sessionDir,'events.mat');
priorEvents = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');   if max(diff(strNumeric))>1, fprintf('\n ERROR: problem converting session name into a numeric.. %s --> %s', sessionName, sessionName(strNumeric)); keyboard; end;
sessNum    = str2num( sessionName(strNumeric) );               if isempty(sessNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...


if (~isempty(priorEvents))
    priorBlockCount  = max([priorEvents.blockCount])+1;
else
    priorBlockCount  = 0;
end

%- Init some state variables
newBlock        = nan;   %- used to track block increments
block           = nan;   %- set using TRIAL X_X (for each stim type, and given oddball class, iterate through the list "block" number of times)
msBlockStart    = -1 ;

%newSession      = -1;      %- previously used to isolate the first SES_START in the file (should only be 1, but could be multiple if session incomplete) JW cut on 2/2015
sessionNum      = sessNum; %- defined from session folder name above
msSessStart     = -1;      %- should be set to mstime of the last "SESSION_START" event in the sesion.log file
blockCount      = priorBlockCount ;

%stimulusDur         = -999;
%correct             = nan       ; 	%- 0-error, 1=correct
%RT                  = -999      ;
%measuredTime        = -999      ;   %- used for ISI and ITI (see logging in pyepl code, some discrepancy (1-20ms?) between timestamps and recorded ISI and ITI within the pyepl code)
trial               = -999      ;   %- trial number within each block of a given stim category (and chosen oddball categ)
iterationNum        = -999      ;   %- iteration through list (typically two times, with oddball classes reversed the second time)
soundvideo          = ''        ;   %- sound or video
stimCateg           = ''        ;
oddballCateg        = ''        ;
nonoddballCateg     = ''        ;
stimulusPath = '';
%%%%
%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 0;
while true
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); break; end
    
    
    %- Generic text scan to get time, offset, and type
    xTOT                = textscan(thisLine,'%f%d%s','delimiter', '\t');
    mstime              = xTOT{1}(1);   %- must add (1) because numbers after the string in the above line cause overflow to first %f
    msoffset            = xTOT{2}(1);
    type                = xTOT{3}{1};
    
    %- default Parameters (details will be filled out/altered based on type)
    experiment          = 'multiModal';
    subject             = subject   ;
    sessionName         = sessionName ;
    sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
    msDuration          = nan       ;
    mstime              = mstime    ;
    msoffset            = msoffset  ;
    type                = type      ;
    subtype             = ''        ;
    block               = block     ;
    stimCateg           = stimCateg ;
    oddballCateg        = oddballCateg;
    nonoddballCateg     = nonoddballCateg ;
    iterationNum        = iterationNum; 
    soundvideo          = soundvideo;
    stimulusDur         = -999;
    correct = nan;
    RT = nan;
    respKey = '';

    switch type
        case 'SESSION_START'
            xTOT=textscan(thisLine,'%f%d%s%s','delimiter', '\t'); % grab the session number
            textAbrev           = sprintf('SESSION %s',sessionName);
            textStr             = sprintf('session %s',sessionName);
            resultStr           = sprintf('%s',sessionName);
            msSessStart         = mstime;
   
%         case 'ISI_END'
%             xTOT=textscan(thisLine,'%f%d%s%d');            
%             measuredTime = xTOT{4}(1);

%         case 'BLOCK'
%             xTOT=textscan(thisLine,'%f%d%s%d%s');
%             block        = xTOT{4}(1)+1;
%             subtype      = xTOT{5}{1};
            
        case 'STIMULUS TYPE'
            xTOT=textscan(thisLine,'%f%d%s%s%s%s%s%s%s%d','delimiter', '\t');
            stimCateg = xTOT{4}{1};
            % Spanish to English
            if strcmp(stimCateg,'los siguientes videos conteniendo acciones ')
                stimCateg = 'action videos ';
            elseif strcmp(stimCateg,'los siguientes acordes de piano ')
                stimCateg = 'piano chords ';
            elseif strcmp(stimCateg,'los siguientes sonidos de instrumentos ')
                stimCateg = 'instrument sounds '; 
            elseif strcmp(stimCateg,'las siguientes melodias ')
                stimCateg = 'melodies ';                 
            end
            
            oddballCateg = xTOT{6}{1};                
            if strcmp(oddballCateg,'una melodia ascendente ')
                oddballCateg = 'upward melody ';
            elseif strcmp(oddballCateg,'una melodia descendente ')
                oddballCateg = 'downward melody ';
            elseif strcmp(oddballCateg,'con una mano ')
                oddballCateg = 'one hand ';    
            elseif strcmp(oddballCateg,'con dos manos ')
                oddballCateg = 'two hands '; 
            elseif strcmp(oddballCateg,'una trompeta ')
                oddballCateg = 'trumpet '; 
            elseif strcmp(oddballCateg,'una flauta ')
                oddballCateg = 'flute ';  
            elseif strcmp(oddballCateg,'triste ')
                oddballCateg = 'mellow '; 
            elseif strcmp(oddballCateg,'alegre ')
                oddballCateg = 'bright ';                 
            end

            iterationNum = xTOT{8}(1);
            
            if strcmp(oddballCateg,'upward melody ')
                nonoddballCateg = 'downward melody ';
                soundvideo = 'sound';
            elseif strcmp(oddballCateg,'downward melody ')
                nonoddballCateg = 'upward melody ';
                soundvideo = 'sound';
            elseif strcmp(oddballCateg,'one hand ')
                nonoddballCateg = 'two hands ';  
                soundvideo = 'video';
            elseif strcmp(oddballCateg,'two hands ')
                nonoddballCateg = 'one hand ';    
                soundvideo = 'video';
            elseif strcmp(oddballCateg,'trumpet ')
                nonoddballCateg = 'flute ';
                soundvideo = 'sound';
            elseif strcmp(oddballCateg,'flute ')
                nonoddballCateg = 'trumpet '; 
                soundvideo = 'sound';
            elseif strcmp(oddballCateg,'mellow ')
                nonoddballCateg = 'bright ';
                soundvideo = 'sound';
            elseif strcmp(oddballCateg,'bright ')
                nonoddballCateg = 'mellow '; 
                soundvideo = 'sound';
            end
            
        case 'TRIAL'
            xTOT=textscan(thisLine,'%f%d%s%s%s%d%s','delimiter', '\t');
            BLOCK_TRIAL = xTOT{4}{1};
            blocktrial = strsplit(BLOCK_TRIAL,'_');
            block=str2double(blocktrial(1));
            trial = str2double(blocktrial(2));
            subtype = xTOT{7}{1};
            
            switch subtype                
                case {'SOUND_PLAYING','VIDEO_PLAYING'}
                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%d%s%s%s%s%s%d');
                    stimulusPath = xTOT{8}{1};
                    stimulusDur  = xTOT{10}(1);
                    
                case 'RESPONSE'
                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%d%s%d%s%s%s%s%s%d');
                    respKey = xTOT{8}{1};
                    RT = xTOT{10}(1);
                    correct = xTOT{12}(1);                    
                case {'BLANK_SCREEN','QUERY_TEXT_ON'}
                    % do nothing; modify these to be TYPE in the pyepl code                    
            end           

            
        case {'B','ISI_START ','ISI_END ','STIMULUS_LIST'...
                'TASK_DONE','RETURN PRESS','RETURN'...
                'SESS_END', 'SESS_END_REWARD','SESS_DURATION',...
                'SESS_DURATION_ACTIVE', 'SESS_DURATION_TOTAL','BLOCK_END','E'}
            % experimenter feedback. do nothing
        otherwise
                fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sessionName)

    end
    index=index+1;
%    if index==1001; break; end

 
    
    
    
    %- asign values to events array
    if index>length(events),        
%         %keyboard
%         if newBlock ~= block,
%             msBlockStart    = mstime ;
%             block           = newBlock ;
%             if ~isnan(newBlock),
%                 blockCount      = blockCount+1;
%             end
%         end
        
        %- create dummy event structure that is upddated below based on type
        clear thisEvent
        thisEvent.experiment        = experiment  ;
        thisEvent.subject           = subject     ;
        thisEvent.sessionName       = sessionName ;
        thisEvent.sessionNum        = sessionNum  ;   % store in state var so all events are assigned a sessionNum  %% JW updated 2/2015
        thisEvent.block             = double(block) ; % store in state var so all events are assigned a block
        thisEvent.type              = type        ;
        thisEvent.subtype           = subtype     ;
        thisEvent.msoffset          = msoffset    ;
        thisEvent.mstime            = mstime      ;
%         thisEvent.mstimeEnd         = mstime + msDuration ; %msDuration isn't really being calculated for everything, so infer duration during analysis using the mstime for the next event
%         thisEvent.msDuration        = msDuration  ;   %- this ms time to following blank time
%         

        thisEvent.trial             = trial ;   %- Cycle number within a block 
        thisEvent.responseCorrect   = double(correct) ;  %- 0-error, 1=correct, -1=skipped, -999=n/a    % force double so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
        thisEvent.RT                = RT          ;
        thisEvent.stimulusDur           = stimulusDur;

        %- strings
        thisEvent.soundorvideo     = soundvideo    ;
        thisEvent.stimCateg        = stimCateg     ; % melodies, instrument sound, videos, piano chords
        thisEvent.respKey          = respKey       ;
        thisEvent.stimulusPath     = stimulusPath  ; %- filepath to stimulus
        thisEvent.oddballCateg     = oddballCateg  ;
        thisEvent.nonoddballCateg  = nonoddballCateg;
        
        %- timing relative to local time points
        thisEvent.msSessStart       = msSessStart  ;
        thisEvent.msBlockStart      = msBlockStart ;
        


            %iSrc = find([events.sampleListIndex]==sampleListIndex & [events.block]==block & [events.isSample]==1);
            %events(iSrc).responseCorrect = double(responseCorrect) ;  % need to force double here so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
            %events(iSrc).resultStr       = resultStr ;
            %events(iSrc).RT              = RT ; %propogate the reaction time back to the samples so sample events can easily be sorted by reaction time

         
            
            
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
    end    
end
fclose(fid);  % close session.log

save(eventFile,'events');

end
