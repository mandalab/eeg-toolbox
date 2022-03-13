function events=RAM_PAL_CallCreateTaskFunctions(subject,expDir,session,sessionDir,forceSESSION, startingElec)
%
% FUNCTION:
%  Calls RAM_PAL_CreateTASKEvents.m for version 1 and
%  RAM_PAL_CreateTASKEVENTS_v2.m for later versions.
%
% DESCRIPTION:
%  extracts the events associated with palRAM
%
% INPUTS:
%  subject.......... 'UP020'
%  expDir........... '/data/eeg/UP020/behavioral/pa3/'
%  session.......... 0
%  forceSESSIONS.... [optional] 1 = forces session to this number
%  startingElec..... [optional] for systems with recording systems that
%                    start at 0, enter 0. Defaults to 1.
% OUTPUTS:
%  events=the events structure
%
clear global
global SUBJECT SESSION events versionNum
SUBJECT = subject;
SESSION = session;
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
else
    forceSESSION = [];
end
if ~exist('startingElec','var') || isempty(startingElec)
    startingElec = 1;
end

% should only be sys3  (subj 48 to 58)
while true
    thisLine = fgetl(fid);
    if ~ischar(thisLine); fclose(fid); return;end
    
    % get the third string before the underscore
    xTOT=textscan(thisLine,'%f\t%f\t%s');
    if isempty(xTOT{1})
        xTOT= textscan(thisLine, '(%fL, %f)\t%*s\t%s','delimiter','\t');
    end
    thisTYPE   = xTOT{3}{1};
    thisMSTIME = xTOT{1}(1);
    thisMSOFF  = xTOT{2};
    
    % based on the type write different fields for this event
    switch upper(thisTYPE)
        case 'INSTRUCT_START' % this should be for post-NIH048 (sys 3 RAM)
            versionNo = 3;
            versionNum = '3.0';
            events=RAM_PAL_CreateTASKEvents_sys3(subject,expDir,session,sessionDir,forceSESSION, startingElec);
    end
end


end %function
