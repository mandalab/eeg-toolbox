function [events] = RAM_FR_CreateAllEvents(subject, expDir, session, forceSession, startingElec)
%
% FUNCTION:
%   RAM_FR_CreateAllEvents(subject,expDir,session,forceSession)
% 
% DESCRIPTION:
%   Wrapper function that makes and saves free recall and math events.
%   This 'wrapper' function calls RAM_FR_CreateTASKEvents and
%   RAM_FR_CreateMATHEvents to create events for the individual components
%   of the experiment.
%
% INPUTS:
%   SUBJECT.........'TJ038_1'
%   EXPDIR..........path to 'session_['SESSION']' directory.  Examples:
%                   '/data/eeg/TJ083/behavioral/RAM_FR/nonstim_data/TJ083/'
%                   '/data/eeg/TJ082/behavioral/RAM_FR_nostim_Olivia/'
%   SESSION.........0 = looks in 'session_0' in EXPDIR
%   FORCESESSION....1 = [optional] sets session to 1 (despite the
%                       fact that behavioral data are in session_0)
%                       Leave blank or empty if session number is
%                       same as SESSION
%   STARTINGELEC....1 = [optional] for systems with recording systems that
%                       start at 0, enter 0. Defaults to 1.
% OUTPUTS:
%   Return and save events:
%     events.mat......... contains 'events'
%
% LAST UPDATED:
%    09/03/14 YE    created from extractPYFR_allEVENTS originally by JFB
%    07/08/16 MST   Updated to work with Zaghloul lab file structure
%
events = [];

if ~exist('startingElec','var') || isempty(startingElec)
    startingElec = 1;
end

% check tp see if subject names match
fprintf('\n')
if isempty(regexp(expDir,subject, 'once'))
  fprintf('  WARNNG: %s not found in %s\n',upper(subject),upper(expDir))
  fprintf('          you might be making an error.\n')
  fprintf('          please check this before making events. EXITING\n\n')
  fprintf('               !!! NO EVENTS SAVED !!! \n\n')
  return
end

% set defaults
if ~exist('forceSession','var')
  forceSession = [];
end

% get the directories
thisSessDirNAME = sprintf('session_%d',session);
thisSessDir     = fullfile(expDir,thisSessDirNAME);
evFile          = fullfile(thisSessDir,'events.mat');
%mevFile         = fullfile(thisSessDir,'MATH_events.mat');

% print opening line
fprintf('  Making FREE RECALL events for ')
fprintf('%s, session %d: \n',subject,session)

%--------------------------------------------
% fprintf('    %-15.15s','EEGLOG:')
% sourceName = 'eeg.eeglog';
% targetName = 'eeg.eeglog.up';
% eeglogFile_source = fullfile(thisSessDir,sourceName);
% eeglogFile_target = fullfile(thisSessDir,targetName);
% if exist(eeglogFile_target,'file')
%   fprintf('%s exists.\n',targetName)  
% else
%   if exist(eeglogFile_source,'file')
%     fixEEGLog(eeglogFile_source,eeglogFile_target);
%     fprintf('DONE.\n')  
%   else
%     error(sprintf('%s does not exist.\n',sourceName))
%   end
% end

%--------------------------------------------
fprintf('    %-15.15s','FREE RECALL: ')
events = RAM_FR_CreateTASKEvents(subject,expDir,session,forceSession, startingElec);
if isempty(events)
return
end
save(evFile,'events');

fprintf('DONE.\n')


%--------------------------------------------
% fprintf('    %-15.15s','MATH:')
% if ~exist(mevFile,'file')
%   [events MATHcfg]=RAM_FR_CreateMATHEvents(subject,expDir,session,forceSession);
%   save(mevFile,'events','MATHcfg');
%   cd(thisSessDir); unix('chgrp RAM_eeg MATH_events.mat');
%   clear events
%   fprintf('DONE.\n')  
% else
%   fprintf('SKIPPING (events exist).\n')  
% end

fprintf('\n\n')
