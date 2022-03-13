function events=RAM_PAL_CreateTASKEvents_v3(subject,expDir,session,sessionDir,forceSESSION, startingElec)
%
% FUNCTION:
%  extractPA3events.m
%
% DESCRIPTION:
%  extracts the events associated with pa3.
%
% INPUTS:
%  subject.......... 'UP020'
%  expDir........... '/data/eeg/UP020/behavioral/pa3/'
%  session.......... 0
%  subject= 'NIH049';
%  expDir = '/Users/sreekumarv/Desktop/localDataWorking/NIH049/behavioral/palRAM/';
%  session=1;
%  sessionDir='session_0';
%  forceSESSIONS = 1;
%  startingElec=0;
%  forceSESSIONS.... [optional] 1 = forces session to this number
%  startingElec..... [optional] for systems with recording systems that
%                    start at 0, enter 0. Defaults to 1.
% OUTPUTS:
%  events=the events structure
%
% CORRECT FIELD:
%  the study orientation, study pair, test orientation,
%  and test pair all contain information whther the item
%  is successfully recalled or not.  If the item was
%  successfully recalled, then the item recalled as well
%  as the reaction time also will appear in all four of
%  these fields
%
% INTRUSION FIELD:
%  0 if the word was correct, was a vocalization, or was a 'pass'.
%  -1 if it was an XLI.  PLI have a number.  The XLI/PLI
%  information is taken from the .ann file.
%
% ISTRIGGER FIELD: Depends on feedback condition
%   'none': always -999
%   'STIM': a trial during which stim was delivered.  NOTE that stim is
%           delivered for a set duration.  The recalled events
%           during a stim trial will all have isTrigger=1, but
%           if took a long time for the participant to respond,
%           they vocalozations may have occured after the stim
%           was turned off.
%   'REAL_TIME': a trial in which word presentation was triggered
%           of a measured oscilation.
%
%
% OTHER NOTES:
%  (1) Written by jfburke 4/11 (john.fred.burke@gmail.com)
%  (2) 5/11 (jfb): added all the fields (correct, resp_word, etc) to the
%  orient events so a user can easily filter the orient period for
%  correct vs. incorrect recalls
%  (3) 5/11: added all the fields (correct, resp_word, etc) to stim
%  events
%  (4) changed the way intrusions are scored (5/12/11)
%
%  (5) updated by VS 10/2018 to fix teh following issues:
%    1. Pythonic indexing used for serial position entries for test events whereas matlab indices used for the same info for study events.
%    2. Probe positions were not correct
%    3. study events had no details filled in (probe pos, resp word, rt, etc which are useful if people just want to do [studyevents.rt] for instance) which had to be done posthoc after processing the test events
%    4. REC start end, probe end, etc weren't being processed right. I populated them with the same test event info just in case people want to extract these events relating to a test word.
%
%

clear global
global SUBJECT SESSION STIMTRIAL STIMTYPE STIMLOC events versionNum stimParams currParamSet
SUBJECT = subject;
SESSION = session;
STIMTYPE='';
STIMTRIAL=nan;
currParamSet = 0;
versionNum = '';
%thisSessDir = sprintf('session_%d',SESSION);  %
thisSessDir = sessionDir;  % NINDS allows for text suffix to session folder (e.g. session_1a, for a broken up session)
sessFile    = fullfile(expDir,thisSessDir,'session.log');

fid = fopen(sessFile,'r');
if fid==-1
    fprintf('session %d..no session.log file found.\n',SESSION);
    fprintf('EXITING\n\n');
    return
end

% you can change the session
if exist('forceSESSION','var') && ~isempty(forceSESSION)
    SESSION=forceSESSION;
end
if ~exist('startingElec','var') || isempty(startingElec)
    startingElec = 1;
end

% get experimental variables
% NOUNPOOL  = getNounPool_local(expDir);

evCounter    = 0;
trialCounter = 0;
events       = [];
ListNumber   = [];
while true
    %try
    thisLine = fgetl(fid);
    if ~ischar(thisLine); fclose(fid); return;end
    
    % get the third string before the underscore
    xTOT=textscan(thisLine,'%f\t%f\t%s');
    if isempty(xTOT{1})
        xTOT= textscan(thisLine, '(%fL, %f)\t%*s\t%s','delimiter','\t');
    end
    thisTYPE = xTOT{3}{1};
    thisMSTIME  = xTOT{1}(1);
    thisMSOFF   = xTOT{2};
    
    %sprintf('%s %f',thisTYPE,thisMSTIME)
    % based on the type write different fields for this event
    switch upper(thisTYPE)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case {'E','B','STIM_ON','FORCED_BREAK','INSTRUCT_END','COUNTDOWN_START','COUNTDOWN_END'}
            % make the event
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF);
            appendNewEvent_local(evCounter,'type',thisTYPE);
            if strcmp(thisTYPE,'STIM_ON')
                appendNewEvent_local(evCounter,'serialpos',pairCounterThisList);
            end
            if strcmp(thisTYPE,'E') && trialCounter>0
                trialCounter = trialCounter -1 ;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'STIM_PARAMS'
            %keyboard
            x=textscan(thisLine,'%f%f%s%d%s%s%s%s%s%f');
            ParamSet = x{4}+1;
            AnodeStr = x{6}{1};
            CathodeStr = x{8}{1};
            
            [AnodeNum, AnodeTag, CathodeNum, CathodeTag] = getStimParams(expDir, subject, AnodeStr, CathodeStr, startingElec);
            
            stimParams(ParamSet).Anode = AnodeNum;
            stimParams(ParamSet).AnodeTag = AnodeTag;
            stimParams(ParamSet).Cathode = CathodeNum;
            stimParams(ParamSet).CathodeTag = CathodeTag;
            stimParams(ParamSet).Amp = x{10};
            
            
            
        case 'INSTRUCT_START'
            x=textscan(thisLine,'%f%f%s%d%s%s');
            thisSess  = x{4};
            if (thisSess-1)~=session
                %error('sessions dont match');
                fprintf('\nWARNING: SESS_START doesnt match session folder'); keyboard;
            end
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,x{1},x{2});
            appendNewEvent_local(evCounter,'type',x{3}{1});
            trialCounter = false;
            
        case {'ENCODING_START','RETRIEVAL_START'}
            % start/increment counters
            %sprintf('thisType %s, trialCounter %d',thisTYPE,trialCounter)
            pairCounterThisList  = 0;
            probeCounterThisList = 0;
            probePos             = 0;
            %             if strcmp(thisTYPE,'ENCODING_START') && trialCounter~=false
            %                 trialCounter = trialCounter + 1;
            %             end
            % make the event
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF)
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'list',trialCounter)
            
            
        case {'TRIAL'}
            evCounter = evCounter+1;
            
            x=textscan(thisLine,'%f%f%s%d%s');
            % Reading in practice, non-stim, or stim
            [subchunk] = regexp(thisLine,'\".*?\"','match','ignorecase');
            temp = subchunk{2}; temp(1) = []; temp(end)=[];
            isStim = strcmp(temp,'STIM');
            
            nums=regexp(thisLine,['\d'],'match');
            ListNumber = str2double(cell2mat([nums(end-1) nums(end)]));
            %thisTrial = x{4}(1);
            %             thisTrial = thisLine(end-1);
            %             thisTrial = str2double(thisTrial);
            %              ListNumber = thisTrial;
            %             if thisTrial~=trialCounter && trialCounter~=false
            %                 %sprintf('mstime %f, thisTrial %d, trialCounter %d',x{1}(1),thisTrial,trialCounter)
            %                 error('trials out of order?')
            %             elseif trialCounter==false
            %                 trialCounter=thisTrial;
            %             end
            
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF);
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'list',trialCounter)
            appendNewEvent_local(evCounter,'isStim',isStim);
            
            % melkalliny - removed rec_start from this list, it's processed
            % below
            % vs - additionally removed rec_end, probe_end. Need to populate
            % those with event features as well.
        case {'STUDY_ORIENT','STUDY_ORIENT_OFF','RETRIEVAL_START','RETRIEVAL_ORIENT_OFF',...
                'PROBE_OFF','DISTRACT_START','DISTRACT_END',...
                'ENCODING_END','PAIR_OFF','PRACTICE_ORIENT',...
                'PRACTICE_ORIENT_OFF','PRACTICE_PAIR_OFF','WAITING_START','WAITING_END',...
                'PRACTICE_RETRIEVAL_ORIENT','ORIENT_START','ORIENT_END','RECALL_END'}; %'REC_START','RECALL_END','PROBE_END','REC_END'
            evCounter = evCounter+ 1;
            mkNewEvent_local(evCounter, thisMSTIME, thisMSOFF);
            %probeCounterThisList = probeCounterThisList + 1;
            %sprintf('type %s ,mstime %f, trialNum %d',thisTYPE,thisMSTIME);
            %sprintf('mstime%f, trialNum %d,trialCounter %d',thisMSTIME,trialNum,trialCounter)
            
            %Wtf is this? What is trialNum?
            %             if trialNum~=trialCounter
            %                 error('trials out of order?')
            %             end
            
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'list',trialCounter)
            
            
            if strcmp(thisTYPE,'RECALL_END')
                trialNum = trialNum + 1;
                if trialCounter~=false,trialCounter=trialCounter+1; end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case {'STUDY_PAIR_START','PRACTICE_PAIR'}
            % make the event
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF)
            pairCounterThisList = pairCounterThisList + 1;
            % read the data again to extract the trial number and words
            
            if strcmp(thisTYPE,'STUDY_PAIR_START')
                x=textscan(thisLine,'%f%f%s%d%s%s%s%s'); %Obsolete here
                [subchunk] = regexp(thisLine,'\".*?\"','match','ignorecase');
                %Finding if stim or not
                temp = subchunk{end}; temp(1)=[]; temp(end)=[];
                isStim = strcmp(temp,'STIM');
                %Finding list number
                temp = isspace(thisLine); temp = find(temp==1);
                temp = thisLine(temp(4):(temp(5))); temp(1)=[]; temp(end-1:end)=[];
                trialNum = str2double(temp);
                listNo = trialNum;
                %trialNum = str2double(lastPart(x{5}{1}))+1;
                temp = isspace(thisLine); temp = find(temp==1);
                temp = thisLine(temp(10):(temp(11))); temp(1:2)=[]; temp(end-2:end)=[];
                word1 = temp;
                %word1 = lastPart(x{6}{1});
                temp = isspace(thisLine); temp = find(temp==1);
                temp = thisLine(temp(12):(temp(13))); temp(1:2)=[]; temp(end-2:end)=[];
                word2 = temp;
                %word2 = lastPart(x{7}{1});
            else
                x=textscan(thisLine,'%f%f%s%d%s%s');
                isStim = -999; %practice
                trialNum=0; %practice
                word1 = lastPart(x{5}{1});
                word2 = lastPart(x{6}{1});
            end
            %Find serial position
            temp = isspace(thisLine); temp = find(temp==1);
            temp = thisLine(temp(8):(temp(9))); temp(1)=[]; temp(end-1:end)=[];
            serialPos = str2double(temp);
            %serialPos = x{4}+1;
            
            
            % check if the trial number is what we expect
            %sprintf('study pair x{1}=%f,trialNum=%d,trialCounter=%d',x{1}(1),trialNum,trialCounter)
            %             if trialNum~=trialCounter
            %                 error('trials out of order?')
            %             end
            %
            %             if serialPos~=pairCounterThisList
            %                 error('Pairs out of order?')
            %             end
            
            % now append
            appendNewEvent_local(evCounter,'serialpos',pairCounterThisList);
            appendNewEvent_local(evCounter,'study_1',word1);
            appendNewEvent_local(evCounter,'study_2',word2);
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'list',trialNum);
            appendNewEvent_local(evCounter,'isStim',isStim);
            
            % go back and fill in this events oreient
            thisPair = [events.list]==trialNum & ...
                [events.serialpos]==serialPos;
            thisOrient = thisPair & strcmp({events.type},'STUDY_ORIENT');
            appendNewEvent_local(thisOrient,'study_1',word1);
            appendNewEvent_local(thisOrient,'study_2',word2);
            
            % if there was a stim, fill that in too
            thisStim = thisPair & strcmp({events.type},'STIM_ON');
            if any(thisStim)
                appendNewEvent_local(thisStim,'study_1',word_1);
                appendNewEvent_local(thisStim,'study_2',word_2);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'PRACTICE_TRIAL' %No longer exists
            trialNum=0;
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF);
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'list',trialCounter)
            
        case {'REC_START','PROBE_START','PROBE_END','REC_END'}
            % make the event
            evCounter = evCounter + 1;
            mkNewEvent_local(evCounter,thisMSTIME,thisMSOFF);
            
            % VS- increment only for REC_START since REC starts before
            % every probe
            if strcmp(thisTYPE,'REC_START')
                probeCounterThisList = probeCounterThisList + 1;
            end
            probePos = probeCounterThisList;
            % read the data again to extract the trial number and words
            x=textscan(thisLine,'%f%f%s%d%s%s%s%s%s%s');
            [subchunk] = regexp(thisLine,'\".*?\"','match','ignorecase');
            
            %Looks like stim info is not present in this type of line
            %Will need to use encoding info
            try
                isStim = strcmp(x{10}{1},'STIM');
            catch
                isStim = 0;
                %x=textscan(thisLine,'%f%f%s%d%s%s%s%s%s');
                %sprintf('set isStim to %d',isStim)
                %sprintf('x{7} = %s',x{7}{1})
            end
            spaces = isspace(thisLine); spaces = find(spaces==1);
            temp = thisLine(spaces(8):(spaces(9))); temp(1)=[]; temp(end-1:end)=[];
            serialPos = temp;
            if ischar(serialPos)  %- JW adding conversion to integer to avoid warning below
                serialPos = str2num(serialPos);
            end
            
            %probePos = x{4}+1; No probePos here
            %serialPos = str2double(lastPart(x{5}{1}))+1;
            %             trialNumStr = lastPart(x{6}{1});
            %             if isstrprop(trialNumStr,'alpha')
            %                 trialNum = 0; %practice trial
            %             else
            %                 trialNum = str2double(lastPart(x{6}{1}))+1;
            %             end
            %sprintf('%f,%d,%d',x{1}(1),trialNum,trialCounter)
            temp = thisLine(spaces(4):(spaces(5))); temp(1:2)=[]; temp(end-2:end)=[];
            probe = temp;
            %probe = lastPart(x{7}{1});
            temp = thisLine(spaces(10):end); temp(1:2)=[]; temp(end-1:end)=[];
            expecting = temp;
            % expecting = lastPart(x{8}{1});
            temp = thisLine(spaces(6):(spaces(7))); temp(1)=[]; temp(end-1:end)=[];
            cue_direction = temp;
            %cue_direction = str2double(lastPart(x{9}{1}));
            
            trialNum = ListNumber;

            % lets go ahead and assign the info into this word. we will
            % populate it later 
            appendNewEvent_local(evCounter,'type',thisTYPE);
            appendNewEvent_local(evCounter,'probepos',probePos);
            appendNewEvent_local(evCounter,'serialpos',serialPos+1);
            appendNewEvent_local(evCounter,'list',trialNum);
            appendNewEvent_local(evCounter,'probe_word',probe);
            appendNewEvent_local(evCounter,'expecting_word',expecting);
            appendNewEvent_local(evCounter,'cue_direction',cue_direction);
            appendNewEvent_local(evCounter,'isStim',isStim);
            
            
            if ~strcmp(thisTYPE,'PROBE_START')
                continue
            end
            
            
            %Correct/incorrect
            annFileName  = sprintf('%d_%d.ann',ListNumber,probeCounterThisList-1);
            annFile      = fullfile(expDir,thisSessDir,annFileName);
            if exist(annFile);
                fid2=fopen(annFile,'r');
                thisRT = [];
                lastRT = -999;
                lastRecWord = -999;
                lastCorrect =  false;
                while true
                    tmpAnnLine=fgetl(fid2);
                    if ~ischar(tmpAnnLine);break;end
                    if numel(tmpAnnLine)==0;continue;end
                    if strcmp(tmpAnnLine(1),'#');continue;end
                    x2=textscan(tmpAnnLine,'%f%f%s');
                    thisRT = round(x2{1});
                    thisWordNum = x2{2};
                    thisRecWord = x2{3}{1};
                    isVOC    =  ~strcmpi(thisRecWord,'<>'); %0 if no voc
                    thisRecWord = x2{3}{1};
                    correct  = strcmp(thisRecWord,expecting);
                    
                    if isVOC
                        lastRT = thisRT;
                        lastRecWord = thisRecWord;
                        lastCorrect = correct;
                    end
                    
                    if ~correct && isVOC
                        intrusion=thisWordNum;
                    else
                        intrusion=0;
                    end
                    
                    % M.E. edits - make new event fields for any
                    % events from annotation
                    
                    evCounter = evCounter + 1;
                    
                    mkNewEvent_local(evCounter,thisMSTIME+thisRT,20) % why 20?
                    
                    appendNewEvent_local(evCounter,'type','REC_EVENT');
                    appendNewEvent_local(evCounter,'probepos',probePos);
                    appendNewEvent_local(evCounter,'serialpos',serialPos+1);
                    appendNewEvent_local(evCounter,'list',trialNum);
                    appendNewEvent_local(evCounter,'probe_word',probe);
                    appendNewEvent_local(evCounter,'expecting_word',expecting);
                    appendNewEvent_local(evCounter,'cue_direction',cue_direction);
                    appendNewEvent_local(evCounter,'isStim',isStim);
                    appendNewEvent_local(evCounter,'RT',thisRT);
                    appendNewEvent_local(evCounter,'correct',correct);
                    appendNewEvent_local(evCounter,'vocalization',isVOC);
                    appendNewEvent_local(evCounter,'intrusion',intrusion);
                    appendNewEvent_local(evCounter,'resp_word',thisRecWord);
                    
%                     %Apply corr/incorr label to encoding event as well
%                     temp = find(strcmp({events.study_1},expecting));
%                     if (size(temp,2) ~=0);
%                         temp = find(strcmp({events.study_1},expecting));
%                         if (size(temp,2) > 1);
%                             [events(temp).correct] = deal(correct);
%                         else
%                             events(temp).correct = correct;
%                         end
%                     else
%                         temp = find(strcmp({events.study_2},expecting));
%                         
%                         if (size(temp,2) > 1);
%                             [events(temp).correct] = deal(correct);
%                         else
%                             events(temp).correct = correct;
%                         end
%                     end
                end
                fclose(fid2);
            else
                thisRecWord = 'NoAnn';
                thisRT = -999;
                fprintf('\n >>> warning: missing annotation file %s <<<',annFile);
            end
            
            
            % now, lets assign RT/corr/resp_word into encoding and retrieval events
            thisPair = [events.list]==trialNum & ...
                ([events.serialpos]==serialPos+1 | ...
                [events.probepos]==probePos);
            theseEvents = thisPair & ...
                (strcmp({events.type},'STUDY_PAIR_START') | ...
                strcmp({events.type},'PROBE_START'));   
            appendNewEvent_local(theseEvents,'correct',lastCorrect);
            appendNewEvent_local(theseEvents,'RT',lastRT);
            appendNewEvent_local(theseEvents,'resp_word',lastRecWord);
            

            % Fill in posthoc but only when REC_END
            if strcmp(thisTYPE,'REC_START')
                % also fill in the study_pair with the probepos
                thisPair = [events.list]==trialNum & [events.serialpos]==serialPos+1;
                thisStudy = thisPair & strcmp({events.type},'STUDY_PAIR_START');
                appendNewEvent_local(thisStudy,'probepos',probePos);
                appendNewEvent_local(thisStudy,'RT',thisRT);
                appendNewEvent_local(thisStudy,'probe_word',probe);
                appendNewEvent_local(thisStudy,'expecting_word',expecting);
                appendNewEvent_local(thisStudy,'cue_direction',cue_direction);
                appendNewEvent_local(thisStudy,'isStim',isStim);
                appendNewEvent_local(thisStudy,'resp_word',thisRecWord);
                
                % if there was a stim, fill that in too
                thisStim = thisPair & strcmp({events.type},'STIM_ON');
                if any(thisStim)
                    appendNewEvent_local(thisStim,'probepos',probePos);
                    appendNewEvent_local(thisStim,'serialpos',serialPos+1);
                    appendNewEvent_local(thisStim,'probe_word',probe);
                    appendNewEvent_local(thisStim,'expecting_word',expecting);
                    appendNewEvent_local(thisStim,'cue_direction',cue_direction);
                    appendNewEvent_local(thisStim,'RT',thisRT);
                end
            end
            
            
    end
    

end



function mkNewEvent_local(evCounter,mstime,msoffset)
global SUBJECT SESSION STIMTYPE STIMTRIAL events versionNum stimParams currParamSet

events(evCounter).subject       = SUBJECT;
events(evCounter).session       = SESSION;
events(evCounter).stimType      = STIMTYPE;
events(evCounter).stimTrial     = STIMTRIAL;
events(evCounter).type          = -999;
events(evCounter).list          = -999;
events(evCounter).serialpos     = -999;
events(evCounter).probepos      = -999;
events(evCounter).study_1       = -999;
events(evCounter).study_2       = -999;
events(evCounter).cue_direction = -999;
events(evCounter).probe_word    = -999;
events(evCounter).expecting_word = -999;
events(evCounter).resp_word     = -999;
% events(evCounter).lag           = -999;
events(evCounter).correct       = -999;
events(evCounter).intrusion     = -999;
events(evCounter).pass          = -999;
events(evCounter).vocalization  = -999;
events(evCounter).RT            = -999;
events(evCounter).mstime        = mstime;
events(evCounter).msoffset      = msoffset;
events(evCounter).isStim        = -999;
events(evCounter).expVersion    = versionNum;
if currParamSet==0
    events(evCounter).stimAnode=nan;
    events(evCounter).stimAnodeTag = '';
    events(evCounter).stimCathode=nan;
    events(evCounter).stimCathodeTag = '';
    events(evCounter).stimAmp=nan;
else
    events(evCounter).stimAnode = stimParams(currParamSet).Anode;
    events(evCounter).stimAnodeTag = stimParams(currParamSet).AnodeTag;
    events(evCounter).stimCathode = stimParams(currParamSet).Cathode;
    events(evCounter).stimCathodeTag = stimParams(currParamSet).CathodeTag;
    events(evCounter).stimAmp = stimParams(currParamSet).Amp;
end

function appendNewEvent_local(evCounter,varargin)
global events
nVar = length(varargin)/2;
for v=1:nVar
    thisVarField = varargin{2*(v-1)+1};
    thisVarData  = varargin{2*(v-1)+2};
    [events(evCounter).(thisVarField)] = deal(thisVarData);
end

function [out] = getExpInfo_local(expDir,str2get)
fid_foo1 = fopen(fullfile(expDir,'../../config.py'),'r');
while true
    thisLine = fgetl(fid_foo1);
    if ~ischar(thisLine);break;end
    if numel(thisLine)==0;continue;end
    if strcmp(thisLine(1),'#');continue;end
    possible_str=textscan(thisLine,'%s%f','Delimiter','=');
    if strcmp(possible_str{1}{1},str2get)
        out=possible_str{2};
        break
    end
end
fclose (fid_foo1);

function [words] = getNounPool_local(expDir)
fid_foo1 = fopen(fullfile(expDir,'nounpool.txt'),'r');
X=textscan(fid_foo1,'%s');
words=X{1};
fclose (fid_foo1);
if ismember(upper('pass'),words)
    error('pass is in the nounpool!')
end

function thisStr = lastPart(thisStr)
splitStr = regexp(thisStr,'_','split');
thisStr = splitStr{end};


function [AnodeNum, AnodeTag, CathodeNum, CathodeTag] = getStimParams(expDir, subject, AnodeStr, CathodeStr, startingElec)
% Get the Tag/Num that was recorded:
if all(isstrprop(AnodeStr,'digit')) && all(isstrprop(CathodeStr,'digit'))
    AnodeNum = str2double(AnodeStr) + (1-startingElec);
    CathodeNum = str2double(CathodeStr) + (1-startingElec);
    AnodeTag = '';
    CathodeTag = '';
    isNum = true;
    isTag = false;
else
    AnodeTag = AnodeStr;
    CathodeTag = CathodeStr;
    AnodeNum = nan;
    CathodeNum = nan;
    isTag = true;
    isNum = false;
end

% Now, to assign the Tag/Num that was not recorded:

% First, try getting the talairach structure and getting the number from there
try
    monoTal = getBipolarSubjElecs(subject, false);
    foundTal = true;
catch e
    foundTal = false;
    warning('could not retrieve tal struct');
end

if foundTal
    if isTag
        anodeMask = strcmp({monoTal.tagName},AnodeStr);
        cathodeMask = strcmp({monoTal.tagName},CathodeStr);
    else
        anodeMask = [monoTal.channel]==AnodeNum;
        cathodeMask = [monoTal.channel]==CathodeNum;
    end
    
    if any(anodeMask)
        AnodeNum = monoTal(anodeMask).channel;
        AnodeTag = monoTal(anodeMask).tagName;
    else
        warning('COULD NOT ASSIGN ANODE NUMBER. CHECK THAT STIM PARAMS ARE CORRECT');
    end
    
    if any(cathodeMask)
        CathodeNum = monoTal(cathodeMask).channel;
        CathodeTag = monoTal(cathodeMask).tagName;
    else
        warning('COULD NOT ASSIGN CATHODE NUMBER. CHECK THAT STIM PARAMS ARE CORRECT');
    end
else % Otherwise, try the jacksheet
    
    jackFile    = fullfile('/data/eeg',subject,'docs','jacksheet.txt');
    subjDir = expDir(1:strfind(expDir,subject)-1);
    jackFileNIH = fullfile(subjDir,subject,'docs','jacksheetMaster.txt');  % NIH modification
    if ~exist(jackFile,'file') & exist(jackFileNIH,'file'),
        
        fid = fopen(jackFileNIH);
        jackOut = textscan(fid, '%d%s%s');
        fclose(fid);
        
        jackFile = fullfile(subjDir,subject,'docs','jacksheet.txt');
        fid = fopen(jackFile,'w');
        for iJack=1:length(jackOut{1}), fprintf(fid,'%d\t%s\n',jackOut{1}(iJack),jackOut{2}{iJack}); end
        fclose(fid);
    end
    if exist(jackFile,'file')
        fid = fopen(jackFile);
        jackOut = textscan(fid, '%d%s');
        nums = jackOut{1};
        tags = jackOut{2};
        if isTag
            anodeMask = strcmp(tags, AnodeStr);
            cathodeMask = strcmp(tags, CathodeStr);
        else
            anodeMask = nums==AnodeNum;
            cathodeMask = nums==CathodeNum;
        end
        
        if any(anodeMask)
            AnodeNum = nums(anodeMask);
            AnodeTag = tags{anodeMask};
        else
            warning('COULD NOT ASSIGN ANODE NUMBER. CHECK THAT STIM PARAMS ARE CORRECT');
        end
        
        if any(cathodeMask)
            CathodeNum = nums(cathodeMask);
            CathodeTag = tags{cathodeMask};
        else
            warning('COULD NOT ASSIGN CATHODE NUMBER. CHECK THAT STIM PARAMS ARE CORRECT');
        end
    else
        warning('COULD NOT GET TALAIRACH STRUCTURE OR JACKSHEET!!!!');
    end
end




