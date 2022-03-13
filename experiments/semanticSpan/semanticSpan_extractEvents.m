function [events] = semanticSpan_extractEvents(sessLogFile, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from semanticSpan %%%%

% % % %% Temporary: just load the one session
% % % clear;
% % % sessLogFile = '/Volumes/Shares/FRNU/dataWorking/eeg/NIH047/behavioral/semanticSpan/session_1/session.log';
% % % subject = 'NIH047';
% % % sessionName = 'session_1';

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1
    iKeep=1:find(diff(strNumeric)>1,1,'first');
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep);
end
sessionNum = str2num( sessionName(strNumeric) );
if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...


%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 1;
moreToRead = true;

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

sessNum = strsplit(sessionDir,'session_');
sessNum = str2double(sessNum{2});

% Initialize empty event structure
eventVars = {'sessNum','trialTypeStr','isWordOnset','blockNumChrono','blockNum','listNum','serialPos','word','categoryStr',...
    'categoryNum','isTarget','isPseudo','resp','rt','isPractice'};
emptyEvent = struct;
for var = 1:length(eventVars)
    emptyEvent.(eventVars{var}) = nan;
end
index = 1;

while moreToRead
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); moreToRead = false; end
    
    if moreToRead
        thisEvent = emptyEvent;
        
        % Read this line separated by spaces
        xTOT = textscan(thisLine,'%s',20);
        thisLineStr = {};
        for k = 1:length(xTOT)
            if ~isempty(xTOT{k})
                thisLineStr = cat(2,thisLineStr,xTOT{k});
            end
        end
        
        % Third column identifies trial type
        trialTypeStr = thisLineStr{3};
        includeEvent = false;
        switch trialTypeStr
            case 'SESSION_TYPE'   % Grab session type (practice or test?)
                if strcmp(thisLineStr{4},'TRAIN'), isPractice = 1;
                elseif strcmp(thisLineStr{4},'TEST'), isPractice = 0;
                else, error('Can''t find session type');
                end
            case 'BLANK_SCREEN'
                includeEvent = true;
            case {'TAR_ON_SCREEN','TXT_ON_SCREEN'}
                % Block number in chronological order
                blockNumChrono = strsplit(thisLineStr{4},'_');
                thisEvent.blockNumChrono = str2double(blockNumChrono{1});
                % Category string
                categoryStr = strsplit(thisLineStr{5},'Category_');
                thisEvent.categoryStr = categoryStr{2};
                % Category number
                categoryNum = strsplit(thisLineStr{6},'CategoryNum_');
                thisEvent.categoryNum = str2double(categoryNum{2});
                % Is target?
                isTarget = strsplit(thisLineStr{7},'IsTarget_');
                thisEvent.isTarget = str2double(isTarget{2});
                % Is pseudo word?
                isPseudo = strsplit(thisLineStr{8},'IsPseudo_');
                thisEvent.isPseudo = str2double(isPseudo{2});
                % Block number (fixed for every session)
                blockNum = strsplit(thisLineStr{9},'BlockNum_');
                thisEvent.blockNum = str2double(blockNum{2});
                % List number (for each block)
                listNum = strsplit(thisLineStr{10},'ListNum_');
                thisEvent.listNum = str2double(listNum{2});
                % Serial position
                serialPos = strsplit(thisLineStr{11},'SerialPos_');
                thisEvent.serialPos = str2double(serialPos{2});
                % Word
                thisEvent.word = thisLineStr{13};
                includeEvent = true;
            case 'RESPONSE_KEY'
                resp = thisLineStr{5};
                if strcmp(resp,'NO_RESP')
                    thisEvent.resp = 0;
                else
                    resp = strsplit(resp,'=');
                    thisEvent.resp = str2double(resp{1});
                    rtStr = thisLineStr{~cellfun('isempty',strfind(thisLineStr,'RxnTime'))};
                    rtStr = strsplit(rtStr,'=');
                    thisEvent.rt = str2double(rtStr{2});
                end
                includeEvent = true;
            case 'SHOW_FEEDBACK'
                includeEvent = true;
        end
        
        if includeEvent
            % Add variables that are common to all events
            thisEvent.sessNum = sessNum;
            thisEvent.mstime = str2double(thisLineStr{1});
            thisEvent.isPractice = isPractice;
            thisEvent.trialTypeStr = trialTypeStr;
            if (index==1)
                events        = thisEvent; %- before events defined must convert to structure
            else
                events(index) = thisEvent;
            end
            index = index+1;
        end
    end
end
fclose(fid);  % close session.log



% Clean up events (fill in NaNs)
ev = 1;
events_clean = [];
while true
    % Every trial starts with BLANK_SCREEN
    if strcmp(events(ev).trialTypeStr,'BLANK_SCREEN')
        eventChunk = events(ev);
        chunkDone = false;
        ev = ev+1;
        while ~chunkDone
            if strcmp(events(ev).trialTypeStr,'BLANK_SCREEN')
                chunkDone = true;
            else
                eventChunk = cat(2,eventChunk,events(ev));
                ev = ev+1;
            end
        end
        % For each chunk (trial-related events), make all variables
        % consistent
        for f = fieldnames(eventChunk)'
            thisFieldAllVals = {eventChunk.(char(f))};
            nanIndex = false(length(thisFieldAllVals),1);
            for ff = 1:length(thisFieldAllVals), if isnan(thisFieldAllVals{ff}), nanIndex(ff) = true; end, end
            if sum(~nanIndex) == 1  % If there's only one unique value
                for v = 1:length(eventChunk)
                    eventChunk(v).(char(f)) = thisFieldAllVals{~nanIndex};
                end
            end
        end
        % Every block, there's a blank screen at the end of the block that get's separated out as an event chunk. Exclude these
        if length(eventChunk) > 1
            if isempty(events_clean), events_clean = eventChunk;
            else, events_clean = cat(2,events_clean,eventChunk);
            end
        end
        
        if ev == length(events), break; end
    end

end

% Fill in isWordOnset
for t = 1:length(events_clean)
    if ismember(events_clean(t).trialTypeStr,{'TAR_ON_SCREEN','TXT_ON_SCREEN'})
        events_clean(t).isWordOnset = 1;
    else
        events_clean(t).isWordOnset = 0;
    end
end

events = events_clean;



