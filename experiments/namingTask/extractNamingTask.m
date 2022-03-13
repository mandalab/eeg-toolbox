function events = extractNamingTask(sessLogFile, subject, sessionName)
% Function for extracting behavioral data from namingTask
% Inputs:
% sessLogFile: directory of session.log file
% subject: e.g. 'NIH036'
% sessionName: e.g. 'session_0'


trialMod = 'IMAGE'; %keyword that let's us know this line has a stimulus

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: ' sessionDir]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1
    iKeep=1:find(diff(strNumeric)>1,1,'first');
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep);
end
sessionNum = str2double( sessionName(strNumeric) );
if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...


%- Read session.log line-by-line and convert to events structure
events(120)      = struct('experiment',[],'subject',[],'sessionName',[],'sessionNum',[],'mstime',[],'stimID',[],'reactionTime',[]);

index       = 1;
moreToRead = true;
while moreToRead
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); moreToRead = false; end
    
    if moreToRead
        % Check whether this line is for a trial
        xTOT = textscan(thisLine,'%s',20);
        thisLineStr = {};
        for k = 1:length(xTOT)
            if ~isempty(xTOT{k})
                thisLineStr = cat(2,thisLineStr,xTOT{k});
            end
        end
        
        if size(thisLineStr,1) >= 7 && strcmp(thisLineStr{7}, trialMod)
            mstime = str2double(thisLineStr{1});
            reactionTime = nan;
            
            stimID = thisLineStr{8};
            fileType = strfind(stimID,'_1.bmp');
            stimID = stimID(1:fileType-1);

            clear thisEvent
            thisEvent.experiment        = 'namingTask';
            thisEvent.subject           = subject;
            thisEvent.sessionName       = sessionName;
            thisEvent.sessionNum        = sessionNum;
            thisEvent.mstime            = mstime;
            thisEvent.stimID            = stimID;
            thisEvent.reactionTime     = reactionTime;

            events(index) = thisEvent;
                
            index = index+1;
        end
    end
end
fclose(fid);  % close session.log

%clean up structure in case it has empty elements.
%this will be the case when only 1 or 2 trials was presented (i.e. training task)
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)), events);
events(empty_elems) = [];


