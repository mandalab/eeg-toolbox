function events = extractDoorsEvents_NIH043(sessLogFile, subject, sessionName, priorEvents)


addpath(genpath('/Users/sreekumarv/Dropbox/eeg_toolbox/trunk'))

%eventFile = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rootEEGdir = '/Users/sreekumarv/Desktop/localDataWorking';
%rootEEGdir  = '/Volumes/shares/FRNU/data/eeg';
%% uncomment following lines to directly run script
% clear all
% rootEEGdir='/Users/sreekumarv/Documents/MyCode/theDoors';
% subject     = 'NIH041';   % 
% sessionName = 'session_5';
% 
% sessionDir  = fullfileEEG(rootEEGdir,subject,sessionName);
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
msDuration      = nan;
blockCount      = priorBlockCount ;

%textAbrev          = ''        ;   %- abbreviated string describing condition for plots
contextIm           = ''        ;   %- context filename
targetIm            = ''        ;
contextNum          = -999      ;   %- 1-1, 2-2
presentationNum     = -999      ;   %- All steps (0-3) get assigned a cycle number (number of times it has been presented in the past + 1). Only type=MSG will get the previous presentationNum since this is propagated. If not, log presentationNum for 'MSG'

%%%%
%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 0;
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
    block               = block     ; 	%- store in state var so all events are assigne a block

    corrResp            = ''        ;   %- 'green' or 'brown' at step 0 or ['green','girl'], etc at step 1
    subResp             = ''        ;   %- ""
    subRespKey          = ''        ;   %- LEFT, RIGHT
    foundTarget         = ''        ;   %- target filename which is obtained
    blockFileName       = ''        ;   %- need to save these files too in the individual session folders (modify pyepl code)
    onScreen            = ''        ;   %- what is being displayed on the screen
    subtype             = ''        ;

    correct             = nan       ; 	%- 0-error, 1=correct, -1=skipped, -999=n/a, gets assigned for steps 1,2, and 3 (overall path correct/wrong at reward presentation - basically, correct at steps 1 AND 2).
    RT                  = -999      ;
    respLoc             = -999      ;   %- 1-Left, 2-Right, 3 rest (shouldn't happen)    
    rewProb             = -999      ;   %- prob of getting a quarter. 0 for incorrect
    measuredTime        = -999      ;   %- used for ISI and ITI (see logging in pyepl code, some discrepancy (1-20ms?) between timestamps and recorded ISI and ITI within the pyepl code)
    stepNum             = -999      ;   %- 1-1, 2-2, 0=target presentation, 3=reward presentation, -999 otherwise
    trialNum            = -999      ;   %- TRIAL in BLOCK_TRIAL gives you overall serial pos within the context presentation


    contextStimuli      = struct();
    %targetIm           = ''        ;   %- target filename, Not logging this
    %currently!! Do it. Though it can be inferred from CORR_RESP at step 1
    %using the contextStimuli struct.  

%     index=index+1
%     if index==6;break;end
    
%     disp('----')
%     disp(mstime)
%     disp(type)
%     disp('----')
    

    
    switch type
        case 'SESSION_START'
            xTOT=textscan(thisLine,'%f%d%s%s'); % grab the session number

            textAbrev           = sprintf('SESSION %s',sessionName);
            textStr             = sprintf('session %s',sessionName);
            resultStr           = sprintf('%s',sessionName);
            msSessStart         = mstime;

%         case {'BLOCK_START','BLOCK_PROMPT'}
%             xTOT=textscan(thisLine,'%f%d%s%s');
        
        case {'STIMULUS_LIST_CONTEXT1','STIMULUS_LIST_CONTEXT2'}
            xTOT=textscan(thisLine,'%f%d%s%s');
            blockFileName = xTOT{4}{1};
            %disp(blockFileName)
      
        case {'STIMULI_CONTEXT1','STIMULI_CONTEXT2'}
            xTOT=textscan(thisLine,'%f%d%s%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q','delimiter',' []','MultipleDelimsAsOne',true);
            contextStimuli = struct();
            for i=1:4
                StimuliName = struct();
                for j= i*4:i*4+2
                    k=j-4*(i-1)-3;
                    if k ==1
                        stimName = 'context';
                        if i ==1,StimuliName.(stimName)=char(xTOT{j}{1}(4:end-2));
                        else StimuliName.(stimName) = char(strsplit(xTOT{j}{1},','));end
                    elseif k==2
                        stimName = 'corrSeq';
                        first = char(xTOT{j}{1}(3:7));
                        second = char(xTOT{j}{1}(12:15));
                        StimuliName.(stimName) = strcat(strrep(first,'''',''),',', strrep(second,'''',''));
                    elseif k==3
                        stimName = 'targetPath';
                        StimuliName.(stimName) = char(strsplit(xTOT{j}{1},','));   
                    end
                    %StimuliName{k,1} = char(strsplit(xTOT{j}{1},',')); 
                end
                field=sprintf('target%d',i);
                contextStimuli.(field) = StimuliName;
            end
                
            %disp(blockFileName)*

        case 'ISI_END'
            xTOT=textscan(thisLine,'%f%d%s%d');            
            measuredTime = xTOT{4}(1);
            
        case 'ITI_END'
            xTOT=textscan(thisLine,'%f%d%s%d');            
            measuredTime = xTOT{4}(1);
            
        case 'BLOCK_CONTEXT'
            xTOT=textscan(thisLine,'%f%d%s%s%s%s');
            block_context       = xTOT{4}{1};
            subtype             = xTOT{5}{1};
            contextIm           = xTOT{6}{1};
            strNumeric          = strsplit(block_context,'_'); 
            block               = str2num(strNumeric{1});
            contextNum          = str2num(strNumeric{2});  
       
            
        case 'BLOCK_TRIAL'
            xTOT=textscan(thisLine,'%f%d%s%s%s');
            block_trial         = xTOT{4}{1};
            strNumeric          = strsplit(block_trial,'_');   
            block               = str2num(strNumeric{1});
            trialNum            = str2num(strNumeric{2}); 
            subtype             = xTOT{5}{1};


            switch subtype
                case 'MSG'
                    xTOT=textscan(thisLine,'%f%d%s%s%s%s%s%s%s%s%s%s%s%d%s%.6f','delimiter',' ');
                    %targetIm = xTOT{6}{1}; %log the image filename!
                    onScreen = 'cue target';
                    %samediff = xTOT{12}(1); %Don't bother. Get this at the
                    %first step below
                    %rewProb  = xTOT{10}(1);
                    
                case 'TARGET_IMAGE'
                    xTOT=textscan(thisLine,'%f%d%s%s%s%s');
                    targetIm           =  xTOT{6}{1};
                    
                case 'STEP'
                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s');
                    stepNum  = xTOT{6}(1);
                    onScreen = xTOT{7}{1}; % What's on screen
                    %disp(onScreen)
                    
                    switch stepNum                        
                        case 0
                            switch onScreen
                                case 'IMAGES_ON_SCREEN'
                                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%s');
                                    presNumOrSubResp = xTOT{10}{1};
                                    
                                    switch presNumOrSubResp
                                        case 'PRESENTATION_NUM'
                                            xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%s%d');
                                            presentationNum = xTOT{11}(1);
                                        case 'SUB_RESP'
                                            xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%s%q%q%d');
                                            subRespParseStr = xTOT{11}{1};
                                            subResp = subRespParseStr(isstrprop(subRespParseStr,'alpha'));
                                    end
                            end
 
                        case 1                    
                            switch onScreen
                                case 'IMAGES_ON_SCREEN'
                                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%s%d');
                                    presentationNum = xTOT{11}(1);
                                case 'TARGET_ON_SCREEN'
                                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s%s%s%s%s%s%s%s%d%s%d');    
                                    corrRespParseStr1 = xTOT{10}{1};
                                    corrRespParseStr2 = xTOT{11}{1};
                                    corrResp = [corrRespParseStr1(isstrprop(corrRespParseStr1,'alpha')), corrRespParseStr2(isstrprop(corrRespParseStr2,'alpha'))];
                                    subRespParseStr1  = xTOT{13}{1};
                                    subRespParseStr2  = xTOT{14}{1};
                                    subResp = [subRespParseStr1(isstrprop(subRespParseStr1,'alpha')),subRespParseStr2(isstrprop(subRespParseStr2,'alpha'))];
                                    correct  = xTOT{16}(1);
                                    RT       = xTOT{18}(1);
                                    msDuration = RT;
                                case 'REWARD'
                                    xTOT=textscan(thisLine,'%f%d%s%s%s%d%s%s');
                                    foundTarget = xTOT{8}{1};
                            end
                    end
                  

                otherwise %word
                    fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sesionName)
 
            end
            
        case {'B','BLOCK_PROMPT','BLOCK_START','ISI_START','ITI_START',...
                'TASK_DONE','RETURN PRESS','RETURN'...
                'SESS_END', 'SESS_END_REWARD','SESS_DURATION',...
                'SESS_DURATION_ACTIVE', 'SESS_DURATION_TOTAL','BLOCK_END','E'}
            % experimenter feedback. do nothing
        otherwise
                fprintf('warning: session.log entry not parsed: type:%s\n index:%d\n session:%s/n', type,index,sesionName)

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
        thisEvent.targetPath        = samediff    ;   %- 1 if same, 0 if diff
        
        %- event counters
        thisEvent.blockCount        = blockCount  ;
        thisEvent.presentationNum   = presentationNum  ;   %- Cycle number within a context
        thisEvent.trialCount        = trialNum ;   %- trial number from beginning of block
        thisEvent.stepNum           = stepNum;
        thisEvent.contextNum        = contextNum;
        
        %- test/response characteristics
        thisEvent.responseCorrect   = double(correct) ;  %- 0-error, 1=correct, -1=skipped, -999=n/a    % force double so NaN's can be found using bracket notations (e.g., isnan([events.responseCorrect])
        thisEvent.RT                = RT          ;
        thisEvent.rewProb           = rewProb     ;
        thisEvent.respLoc           = respLoc     ; %- 1=left, 2=right
        thisEvent.subResp           = subResp     ;
        thisEvent.corrResp          = corrResp    ;
        thisEvent.measuredTime      = measuredTime; %- measured in exp code as opposed to calculating it from log times (probably gets rid of lag between log command and actual logging)
        
        %- strings
        thisEvent.foundTarget       = foundTarget ; %- filename of target found
        thisEvent.contextIm         = contextIm   ; %- filename of background context
        thisEvent.blockFileName     = blockFileName;
        thisEvent.onScreen          = onScreen    ;
        %- timing relative to local time points
        thisEvent.msSessStart       = msSessStart  ;
        thisEvent.msBlockStart      = msBlockStart ;
        thisEvent.target            = targetIm     ;
        
        %structs
        thisEvent.contextStimuli    = contextStimuli;

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

end
