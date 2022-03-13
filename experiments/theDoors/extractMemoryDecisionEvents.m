function events = extractMemoryDecisionEvents(sessLogFile, subject, sessionName, priorEvents)


%eventFile = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncomment following lines to directly run script
clear all
addpath(genpath('/Users/sreekumarv/Dropbox/eeg_toolbox/trunk'))

%rootEEGdir  = '/Volumes/shares/FRNU/data/eeg';
%rootEEGdir='/Users/sreekumarv/Documents/MyCode/theDoors';
rootEEGdir = '/Users/sreekumarv/Desktop/localDataWorking';
subject     = 'Wittig_Pilot';%'DBS049_left';   % 
sessionName = 'session_0';

sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/theDoors',sessionName);
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
block           = nan;   %- set using BLOCK_CONTEXT and/or BLOCK_TRIAL
msBlockStart    = -1 ;

%newSession      = -1;      %- previously used to isolate the first SES_START in the file (should only be 1, but could be multiple if session incomplete) JW cut on 2/2015
sessionNum      = sessNum; %- defined from session folder name above
msSessStart     = -1;      %- should be set to mstime of the last "SESSION_START" event in the sesion.log file
blockCount      = priorBlockCount ;

%textAbrev          = ''        ;   %- abbreviated string describing condition for plots


%- default Parameters (details will be filled out/altered based on type)
% experiment          = 'theDoors';
% subject             = subject   ;
% sessionName         = sessionName ;
% sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
% msDuration          = nan       ;
% samediff            = -999      ;   %- 1 for paths associated with a target preserved across pair of contexts
% mstime              = mstime    ;
% msoffset            = msoffset  ;
% type                = type      ;
% block               = block     ; 	%- store in state var so all events are assigne a block

corrResp            = ''        ;   %- the higher prob stimulus (LEFT OR RIGHT)
rewResp             = ''        ;   %- the side that was rewarded (LEFT OR RIGHT)
block_context1_FileName = ''    ;   %- context 1 study trial info for this block
block_context2_FileName = ''    ;   %- context 2 study trial info for this block
block_test_FileName = ''        ;   %- test trial info for this block
%subtype             = ''        ;
contextIm           = ''        ;   %- context filename
context          = -999      ;   %- 1-1, 2-2

correct             = nan       ; 	%- 0-error, 1=correct
RT                  = -999      ;
leftProb            = -999      ;   %- prob associated with stim appearing on the left
rightProb           = -999      ;   %- prob associated with stim appearing on the right 
measuredTime        = -999      ;   %- used for ISI and ITI (see logging in pyepl code, some discrepancy (1-20ms?) between timestamps and recorded ISI and ITI within the pyepl code)
trial               = -999      ;   %- trial number within each context at study, and within test
contextDependent    = -999      ;   %- 1 if the higher prob choice is context dependent, 0 if correct choice is context independent
ACateg              = ''        ;
BCateg              = ''        ;
CCateg              = ''        ;
DCateg              = ''        ;
context1_Im         = ''        ;
context2_Im         = ''        ;
changeCategs        = {}  ;
constCategs         = {} ;
leftIm = '';
rightIm = '';


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
    
    %- default Parameters (details will be filled out/altered based on type)
    experiment          = 'theDoors';
    subject             = subject   ;
    sessionName         = sessionName ;
    sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
    msDuration          = nan       ;
    samediff            = -999      ;   %- 1 for paths associated with a target preserved across pair of contexts
    mstime              = mstime    ;
    msoffset            = msoffset  ;
    type                = type      ;
    subtype             = '';
    block               = block     ;
    changeCategs        = changeCategs  ;
    constCategs         = constCategs  ;

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
            % Read the csv file that contains category information for this block
            % This is used to determine which choices are context dependent and
            % which aren't
            categFileName = sprintf('block%d_CategoryInfo_v3.3.csv',block-1);
            blockCategFullFile = fullfileEEG(sessionDir,categFileName);
            Categs = textread(blockCategFullFile,'%s','delimiter',',');
            ACateg = Categs{1}; BCateg = Categs{2}; CCateg = Categs{3}; DCateg = Categs{4};
            % The correct A or B choice depends on context whereas the
            % correct C or D choice remains the same across contexts.
            changeCategs = {ACateg,BCateg};
            constCategs = {CCateg,DCateg}; 
            
            block_context1_FileName = sprintf('block%d_0v3.3.csv',block);
            block_context2_FileName = sprintf('block%d_1v3.3.csv',block);
            block_test_FileName = sprintf('block%d_test_v3.3.csv',block);
            
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
                case {'STUDY','STUDY_CHOICE','STUDY_FEEDBACK_START','STUDY_FEEDBACK_END'}
                    
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%d%s%s%s%s%f%s%f%s%s%s%d%s%s%s%d');
                    %xTOT=textscan(thisLine,'%f%d%s%d%s%s%d%s%d%s%s%s%s%f%s%f%s%s%s%d%s%s%s%f','delimiter',' ');
                    context = xTOT{7}(1)+1; % 1 or 2
                    if context==1; contextIm = context1_Im;
                    elseif context==2; contextIm = context2_Im;
                    end
                    trial   = xTOT{9}(1)+1;
                    leftIm  = xTOT{11}{1};
                    rightIm = xTOT{12}{1};
                    leftProb = xTOT{14}(1);
                    rightProb = xTOT{16}(1);
                    corrResp = xTOT{18}{1}; % LEFT or RIGHT
                    correct = xTOT{20}(1); % 1 if correct, 0 if not, 999 default/non-response event logging
                                           % for 'STUDY', read from nextline to get info
                                           % about whether this choice was
                                           % correct - fix this in the code
                                           % later to directly log this
                                           % info on this 'STUDY' line
                    rewResp = xTOT{22}{1}; % which side was rewarded on this trial.
                                           % If corrResp and rewResp
                                           % mismatch, that is a prediction
                                           % error trial when probs have
                                           % been learned well.
                    RT = xTOT{24}(1);      % default = 2000000  non-response event logging
                                           % no response = 1000000
                                           
                    % get categories of left and right images, compare them
                    % with ACateg,BCateg and CCateg,DCateg.
                    % contextDependent = 1 if match with {ACateg,BCateg}.
                    % Else 0.
                    lSplit = strsplit(leftIm,'/');
                    leftCateg = lSplit{2};
                    rSplit = strsplit(rightIm,'/');
                    rightCateg = rSplit{2};  
                    currCategs = {leftCateg,rightCateg};
                    if isempty(setxor(currCategs,changeCategs)); contextDependent = 1;
                    elseif isempty(setxor(currCategs,constCategs)); contextDependent = 0;
                    else contextDependent = 12345;
                    end
                    
                         

                case {'TEST','TEST_CHOICE'}
                    xTOT=textscan(thisLine,'%f%d%s%d%s%s%s%s%d%s%s%s%s%s%s%d%s%d');
                    contextIm = xTOT{7}{1}; 
                    if strcmp(contextIm,context1_Im);context=1;
                    elseif strcmp(contextIm,context2_Im); context=2;
                    end
                    trial   = xTOT{9}(1)+1;
                    leftIm  = xTOT{11}{1};
                    rightIm = xTOT{12}{1};
                    corrResp = xTOT{14}{1}; % LEFT or RIGHT
                    correct = xTOT{16}(1); % 1 if correct, 0 if not, 999 default/non-response event logging
                                           % for 'STUDY', read from nextline to get info
                                           % about whether this choice was
                                           % correct - fix this in the code
                                           % later to directly log this
                                           % info on this 'STUDY' line
                    RT = xTOT{18}(1);      % default = 2000000  non-response event logging
                                           % no response = 1000000
                    lSplit = strsplit(leftIm,'/');
                    leftCateg = lSplit{2};
                    rSplit = strsplit(rightIm,'/');
                    rightCateg = rSplit{2};  
                    currCategs = {leftCateg,rightCateg};
                    if isempty(setxor(currCategs,changeCategs)); contextDependent = 1;
                    elseif isempty(setxor(currCategs,constCategs)); contextDependent = 0;
                    else contextDependent = 12345;
                    end
                    %%% Add leftProb and rightProb info based on
                    %%% contextDependent and context identity (or context 1
                    %%% vs 2) - higher diff always in context 1?
                    
                    
                    
                otherwise %word
                    fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sessionName)
 
            end
            
        case {'B','ITI_START',...
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
        thisEvent.mstimeEnd         = mstime + msDuration ; %msDuration isn't really being calculated for everything, so infer duration during analysis using the mstime for the next event
        thisEvent.msDuration        = msDuration  ;   %- this ms time to following blank time
        

        thisEvent.blockCount        = blockCount  ;
        thisEvent.trial             = trial ;   %- Cycle number within a context
        thisEvent.context           = context;
        thisEvent.contextIm         = contextIm;
        thisEvent.contextDependent  = contextDependent;
        %- test/response characteristics
        thisEvent.responseCorrect   = double(correct) ;  %- 0-error, 1=correct, -1=skipped, -999=n/a    % force double so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
        thisEvent.RT                = RT          ;
        thisEvent.leftProb          = leftProb   ;
        thisEvent.rightProb         = rightProb   ;
        thisEvent.corrResp          = corrResp    ;
        thisEvent.rewResp           = rewResp ;
        thisEvent.measuredTime      = measuredTime; %- measured in exp code as opposed to calculating it from log times (probably gets rid of lag between log command and actual logging)
        
        %- strings
        thisEvent.contextIm         = contextIm   ; %- filename of background context
        thisEvent.blockcontext1FName = block_context1_FileName;   %- context 1 study trial info for this block
        thisEvent.block_context2FName = block_context2_FileName;   %- context 2 study trial info for this block
        thisEvent.block_testFName    = block_test_FileName;   %- test trial info for this block
        thisEvent.leftIm = leftIm;
        thisEvent.rightIm = rightIm;
        
        %- timing relative to local time points
        thisEvent.msSessStart       = msSessStart  ;
        thisEvent.msBlockStart      = msBlockStart ;
        
        %structs
        thisEvent.changeCategs    =     changeCategs;
        thisEvent.constCategs     =     constCategs;

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
