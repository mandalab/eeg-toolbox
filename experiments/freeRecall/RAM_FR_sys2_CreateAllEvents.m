function RAM_FR_sys2_CreateAllEvents(subject,exp,session,forceSession, diagnostic_plots, eegRoot, force_redo)
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
%
% OUTPUTS:
%   Makes and save three events:
%     (1) events.mat......... contains 'events'
%     (2) MATH_events.mat.... contains 'events' and 'MATHcfg
%
% LAST UPDATED:
%    09/03/14 YE    created from extractPYFR_allEVENTS originally by JFB
%
%

if ~exist('force_redo','var') || isempty(force_redo)
    force_redo = false;
end

if ~exist('eegRoot','var') || isempty(eegRoot)
    eegRoot = '/data/eeg';
end

if ~exist('diagnostic_plots', 'var') || isempty(diagnostic_plots)
    diagnostic_plots = true;
end

expDir = fullfile(eegRoot, subject, 'behavioral', exp);

% set defaults
if ~exist('forceSession','var') 
  forceSession = [];
end

% get the directories
thisSessDirNAME = sprintf('session_%d',session);
thisSessDir     = fullfile(expDir,thisSessDirNAME);
 evFile         = fullfile(thisSessDir,'events.mat');
mevFile         = fullfile(thisSessDir,'MATH_events.mat');

% print opening line
fprintf('  Making FREE RECALL and MATH events for ')
fprintf('%s, session %d: \n',subject,session)

%--------------------------------------------
fprintf('    %-15.15s','FREE RECALL: ')
if ~exist(evFile,'file') || force_redo
  events=RAM_FR_sys2_CreateTASKEvents(subject,expDir,session,forceSession);
  if isempty(events)
    return
  end
  
  events=RAM_sys2_alignEvents(events, subject, exp, session, diagnostic_plots, eegRoot); 
  
  save(evFile,'events');
  clear events
  fprintf('DONE.\n')
else
  fprintf('SKIPPING (events exist).\n')
end


%--------------------------------------------
fprintf('    %-15.15s','MATH:')
if ~exist(mevFile,'file') || force_redo
  [events MATHcfg]=RAM_FR_sys2_CreateMATHEvents(subject,expDir,session,forceSession);
  
  events=RAM_sys2_alignEvents(events, subject, exp, session, diagnostic_plots, fullfile(eegRoot)); 
  save(mevFile,'events','MATHcfg');
  
  
  clear events
  fprintf('DONE.\n')  
else
  fprintf('SKIPPING (events exist).\n')  
end

fprintf(' Adding to database...\n')

RAM_sys2_AddEventsToDatabase(subject, ['RAM_' exp], session);
fprintf('\n\n')