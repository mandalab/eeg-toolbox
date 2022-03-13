function events = extractdelayMatching_v1(sessLogFile, subject, sessionName, priorEvents)


%eventFile = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncomment following lines to directly run script
% clear all
% addpath(genpath('/Users/vishnusreekumar/Dropbox/eeg_toolbox/trunk'))
% 
% %rootEEGdir  = '/Volumes/shares/FRNU/data/eeg';
% %rootEEGdir='/Users/sreekumarv/Documents/MyCode/theDoors';
% rootEEGdir = '/Users/sreekumarv/Desktop/localDataWorking/delayMatching';
% subject     = 'NIH047';%'DBS049_left';   % 
% sessionName = 'session_1';
% 
% sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral',sessionName);
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];

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
block           = nan;   %- set using BLOCK_CONTEXT and/or BLOCK_TRIAL
msBlockStart    = -1 ;

%newSession      = -1;      %- previously used to isolate the first SES_START in the file (should only be 1, but could be multiple if session incomplete) JW cut on 2/2015
sessionNum      = sessNum; %- defined from session folder name above
msSessStart     = -1;      %- should be set to mstime of the last "SESSION_START" event in the sesion.log file
blockCount      = priorBlockCount ;

%textAbrev          = ''        ;   %- abbreviated string describing condition for plots


%- default Parameters (details will be filled out/altered based on type)
experiment          = 'delayMatching';
subject             = subject   ;
sessionName         = sessionName ;
sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
msDuration          = nan       ;

% block               = block     ; 	%- store in state var so all events are assigne a block
correct             = nan       ; 	%- 0-error, 1=correct
RT                  = -999      ;
delay               = -999      ;   %- the delay on a given trial after stim goes off the screen.
measuredTime        = -999      ;   %- used for ISI and ITI (see logging in pyepl code, some discrepancy (1-20ms?) between timestamps and recorded ISI and ITI within the pyepl code)
trial               = -999      ;   %- trial number within each context at study, and within test
leftChoice          = ''        ;
rightChoice         = ''        ;
stimIm              = ''        ;
respKey             = ''        ;
corrKey             = ''        ;


%%%%
%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 0;
currentBlock = nan;
while true
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); break; end   
    
    %- Generic text scan to get time, offset, and type
    xTOT                = textscan(thisLine,'%f%d%s');
    mstime              = xTOT{1}(1);   %- must add (1) because numbers after the string in the above line cause overflow to first %f
    msoffset            = xTOT{2}(1);
    type                = xTOT{3}{1};
    subtype             = '';
    
  

    switch type
        case 'SESSION_START'
            xTOT=textscan(thisLine,'%f%d%s%s'); % grab the session number
            textAbrev           = sprintf('SESSION %s',sessionName);
            textStr             = sprintf('session %s',sessionName);
            resultStr           = sprintf('%s',sessionName);
            msSessStart         = mstime;
            
        case 'BLOCK_PROMPT'
            xTOT=textscan(thisLine,'%f%d%s%d'); % grab the session number
            block = xTOT{4}(1)+1;     
            
        case 'BLOCK_START'
            xTOT=textscan(thisLine,'%f%d%s%d'); % grab the session number
            block = xTOT{4}(1)+1;
            
        case 'ITI_END'
            xTOT=textscan(thisLine,'%f%d%s%d');           
            measuredTime = xTOT{4}(1);

        case 'BLOCK'
            xTOT=textscan(thisLine,'%f%d%s%d%s%s');
            block        = xTOT{4}(1)+1;
            subtype      = xTOT{6}{1};

            switch subtype
                case {'STIM_PRESENT_ONSET','DELAY_BEGIN','DELAY_END','CHOICE_END'}
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%s%d%s%s%s%d%s%s%s%s%s%s%s%s%s%d%s%d');
                    trial   = xTOT{8}(1)+1;
                    stimIm  = xTOT{10}{1};
                    delay   = xTOT{12}(1);
                    leftChoice = xTOT{14}{1};
                    rightChoice = xTOT{16}{1};
                    respKey     = xTOT{18}{1};
                    corrKey     = xTOT{20}{1};
                    correct     = xTOT{22}(1);
                    RT          = xTOT{24}(1);

               

                    
                    
                otherwise %word
                    fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sessionName)
 
            end
            
        case {'B','ITI_START',...
                'TASK_DONE','RETURN PRESS','RETURN','REST_PROMPT','CONTINUE_BLOCK',...
                'SESS_END', 'SESS_END_REWARD','SESS_DURATION',...
                'SESS_DURATION_ACTIVE', 'SESS_DURATION_TOTAL','','BLOCKS_DURATION_ACTIVE','BLOCK_END','E'}
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
        thisEvent.mstimeEnd         = mstime + msDuration ; %msDuration isn't really being calculated for everything, so infer duration during analysis using the mstime for the next event
        thisEvent.msDuration        = msDuration  ;   %- this ms time to following blank time
        

        thisEvent.blockCount        = blockCount  ;
        thisEvent.trial             = trial ;   %- Cycle number within a context
        thisEvent.stimImage         = stimIm;
        thisEvent.delay             = delay;
        thisEvent.leftChoice        = leftChoice;
        thisEvent.rightChoice       = rightChoice;
        thisEvent.subResp           = respKey;
        thisEvent.corrResp          = corrKey;
        %- test/response characteristics
        thisEvent.responseCorrect   = double(correct) ;  %- 0-error, 1=correct, -1=skipped, -999=n/a    % force double so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
        thisEvent.RT                = RT          ;

        thisEvent.measuredTime      = measuredTime; %- measured in exp code as opposed to calculating it from log times (probably gets rid of lag between log command and actual logging)
        

        
        %- timing relative to local time points
        thisEvent.msSessStart       = msSessStart  ;
        thisEvent.msBlockStart      = msBlockStart ;
        
            
            
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
    end    
end
fclose(fid);  % close session.log

% save(eventFile,'events');

end
