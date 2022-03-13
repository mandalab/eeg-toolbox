function events = extractzapList_v1(sessLogFile, subject, sessionName, priorEvents)


%eventFile = '';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% uncomment following lines to directly run script
% clear all
% addpath(genpath('/Users/vishnusreekumar/Dropbox/eeg_toolbox/trunk'))
% 
% %rootEEGdir  = '/Volumes/shares/FRNU/data/eeg';
% %rootEEGdir='/Users/sreekumarv/Documents/MyCode/theDoors';
% rootEEGdir = '/Users/sreekumarv/Desktop/localDataWorking/zapList';
% subject     = 'NIH047';%'DBS049_left';   % 
% sessionName = 'session_2';
% 
% sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral',sessionName);
% sessLogFile = fullfileEEG(sessionDir,'session.log');
% eventFile   = fullfileEEG(sessionDir,'events.mat');
% priorEvents = [];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%------------------------------------------------------------------------------------------
%- fix the eeg.log here... cut out the "s" at the end of it
eegLogName = [sessLogFile(1:end-11) 'eeg.eeglog'];
eegLogWith  = [eegLogName '.originalWithS'];
eegLogTemp  = [eegLogName '.temp'];
if ~exist(eegLogWith,'file'),
    
    fprintf('\n EEG.EEGLOG has not been copied to ".withS" version... do that now');
    fidEEGold = fopen(eegLogName,'r');
    if fidEEGold==-1, error('\n couldnt move and open eeglog with S at the end'); end;
    fidEEGnew = fopen(eegLogTemp,'w');
    if fidEEGnew==-1, error('\n couldnt move new eeg log'); end;
    
    %- copy everything except for the final "s"
    while true
        thisLine = fgetl(fidEEGold);
        if ~ischar(thisLine),    break; end
        if ~strcmp(strtrim(thisLine),'s'),
            fprintf(fidEEGnew, '%s\n', thisLine);
        else
            %keyboard
        end
    end
    fclose(fidEEGold);
    fclose(fidEEGnew);
    
    %- only copy after its done
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(eegLogName,eegLogWith,'f');
    if SUCCESS==0, error('\n copy failed'); end
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(eegLogTemp,eegLogName,'f');
    if SUCCESS==0, error('\n copy failed'); end
    delete(eegLogTemp); % remove temp
end
%------------------------------------------------------------------------------------------





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
experiment          = 'zapList';
subject             = subject   ;
sessionName         = sessionName ;
sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
msDuration          = nan       ;

% block               = block     ; 	%- store in state var so all events are assigne a block

block_context1_FileName = ''    ;   %- context 1 study trial info for this block
block_context2_FileName = ''    ;   %- context 2 study trial info for this block
block_test_FileName = ''        ;   %- test trial info for this block
contextIm           = ''        ;   %- context filename
context          = -999      ;   %- 1-1, 2-2

correct             = nan       ; 	%- 0-error, 1=correct
RT                  = -999      ;
conf                = -999      ;   % 0- not confident, 1- confident.
measuredTime        = -999      ;   %- used for ISI and ITI (see logging in pyepl code, some discrepancy (1-20ms?) between timestamps and recorded ISI and ITI within the pyepl code)
trial               = -999      ;   %- trial number within each context at study, and within test
contextDependent    = -999      ;   %- 1 if the higher prob choice is context dependent, 0 if correct choice is context independent
condition           = -999      ;
context1_Im         = ''        ;
context2_Im         = ''        ;
word                = ''        ;
studyStim           = ''        ;
testStim            = ''        ;




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
            
        case 'BLOCK_CONTEXT'
            xTOT=textscan(thisLine,'%f%d%s%s%s%s'); % grab the session number
            block_context = xTOT{4}{1};
            strNumeric = strsplit(block_context,'_'); 
            block = str2num(strNumeric{1})+1;
            contextNum = str2num(strNumeric{2})+1; %1 or 2
            contextIm = xTOT{6}{1};
            if contextNum==1;context1_Im = contextIm;
            elseif contextNum==2;context2_Im = contextIm;
            end
            
          
        case 'ITI_END'
            xTOT=textscan(thisLine,'%f%d%s%d');           
            measuredTime = xTOT{4}(1);

        case 'BLOCK'
            xTOT=textscan(thisLine,'%f%d%s%d%s');
            block        = xTOT{4}(1)+1;
            subtype      = xTOT{5}{1};

            switch subtype
                case 'STUDY'
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%s%s%d%s%s%s%s%s%s%s%s%s%s%s%s%s%d');
                    %[]...]
                    %xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%d%s%s%s%s%f%s%f%s%s%s%d%s%s%s%f','delimiter',' ');
                    context = xTOT{7}(1); % 1 or 2
                    contextIm = xTOT{9}{1};
                    trial   = xTOT{11}(1)+1;
                    word = xTOT{13}{1};
                    s1=xTOT{15}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{16}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{17}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{18}{1};s4(regexp(s4,'[[,],,]'))=[];
                    
                    studyStim = strcat(s1,s2,s3,s4);
                    
                    s1=xTOT{20}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{21}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{22}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{23}{1};s4(regexp(s4,'[[,],,]'))=[];                    
                    
                    testStim = strcat(s1,s2,s3,s4);
                    condition = xTOT{25}(1);
               

                case {'TEST_WORD_ONSET','TEST_CHOICE_ONSET','TEST_RESPONSE','_EARLY_TEST_RESPONSE'}   %%- JW put _EARLY TEST RESPONSE in this list... is that OK?
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%s%s%d%s%s%s%s%s%s%s%s%s%s%s%s%s%d%s%d%s%d');
                    context = xTOT{7}(1); % 1 or 2
                    contextIm = xTOT{9}{1};
                    trial   = xTOT{11}(1)+1;
                    word = xTOT{13}{1};
                    s1=xTOT{15}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{16}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{17}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{18}{1};s4(regexp(s4,'[[,],,]'))=[];
                    
                    studyStim = strcat(s1,s2,s3,s4);
                    
                    s1=xTOT{20}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{21}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{22}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{23}{1};s4(regexp(s4,'[[,],,]'))=[];                    
                    
                    testStim = strcat(s1,s2,s3,s4);
                    condition = xTOT{25}(1);
                    correct = xTOT{27}(1); % 1 if correct, 0 if not
                    RT = xTOT{29}(1);      

                    
                case 'TEST_RESPONSE_CONF'
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%s%s%d%s%s%s%s%s%s%s%s%s%s%s%s%s%d%s%d%s%d%s%d');
                    context = xTOT{7}(1); % 1 or 2
                    contextIm = xTOT{9}{1};
                    trial   = xTOT{11}(1)+1;
                    word = xTOT{13}{1};
                    s1=xTOT{15}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{16}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{17}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{18}{1};s4(regexp(s4,'[[,],,]'))=[];
                    
                    studyStim = strcat(s1,s2,s3,s4);
                    
                    s1=xTOT{20}{1};s1(regexp(s1,'[[,],,]'))=[];
                    s2=xTOT{21}{1};s2(regexp(s2,'[[,],,]'))=[];
                    s3=xTOT{22}{1};s3(regexp(s3,'[[,],,]'))=[];
                    s4=xTOT{23}{1};s4(regexp(s4,'[[,],,]'))=[];                    
                    
                    testStim = strcat(s1,s2,s3,s4);
                    condition = xTOT{25}(1);
                    correct = xTOT{27}(1); % 1 if correct, 0 if not
                    RT = xTOT{29}(1);   
                    conf = xTOT{31}(1); 
                    
                    
                otherwise %word
                    fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s\n', type,index,sessionName)
 
            end
            
        case {'B','ITI_START',...
                'TASK_DONE','RETURN PRESS','RETURN'...
                'SESS_END', 'SESS_END_REWARD','SESS_DURATION',...
                'SESS_DURATION_ACTIVE', 'SESS_DURATION_TOTAL','','BLOCKS_DURATION_ACTIVE','BLOCK_END','E'...
                'START', 'PROB', 'STOP'...  %- math stuff
                'You'}
            % experimenter feedback. do nothing
        otherwise
                fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sessionName)

    end
    index=index+1;
%    if index==1001; break; end

 
    
    
    
    %- asign values to events array
    if index>length(events),        
        
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
        thisEvent.context           = context;
        thisEvent.word              = word;
        thisEvent.condition         = condition;
        %- test/response characteristics
        thisEvent.responseCorrect   = double(correct) ;  %- 0-error, 1=correct, -1=skipped, -999=n/a    % force double so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
        thisEvent.RT                = RT          ;
        thisEvent.conf              = conf        ;

        thisEvent.measuredTime      = measuredTime; %- measured in exp code as opposed to calculating it from log times (probably gets rid of lag between log command and actual logging)
        
        %- strings
        thisEvent.contextIm         = contextIm   ; %- filename of background context
        thisEvent.studyStim         = studyStim   ;
        thisEvent.testStim          = testStim    ;
        
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

%save(eventFile,'events');

end
