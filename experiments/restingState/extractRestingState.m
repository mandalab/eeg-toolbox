function [events] = extractRestingState(sessLogFile, subject, sessionName)
% Function for extracting behavioral data from paRemap
% Inputs:
% sessLogFile: directory of session.log file
% subject: e.g. 'NIH036'
% sessionName: e.g. 'session_0'

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1, iKeep=[1:find(diff(strNumeric)>1,1,'first')]; fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep); end;
sessionNum = str2num( sessionName(strNumeric) );               if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...



%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 1;
while true
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); break; end
    
    
    %- Generic text scan to get time, offset, and type
    xTOT                = textscan(thisLine,'%f%d%s');
    mstime              = xTOT{1}(1);   %- must add (1) because numbers after the string in the above line cause overflow to first %f
    msoffset            = xTOT{2}(1);
    type                = xTOT{3}{1};

    %- create dummy event structure that is upddated below based on type
    clear thisEvent
    thisEvent.experiment        = 'restingState';
    thisEvent.subject           = subject     ;
    thisEvent.sessionName       = sessionName ;
    thisEvent.sessionNum        = sessionNum  ;   % store in state var so all events are assigned a sessionNum  %% JW updated 2/2015
    thisEvent.type              = type        ;
    thisEvent.msoffset          = msoffset    ;
    thisEvent.mstime            = mstime      ;
    
    if (index==1)
        events        = thisEvent; %- before events defined must convert to structure
    else
        events(index) = thisEvent;
    end
    
    index = index+1;
end
fclose(fid);  % close session.log