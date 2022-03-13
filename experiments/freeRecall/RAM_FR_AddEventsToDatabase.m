function [] = RAM_FR_AddEventsToDatabase(subject,expDir,session,RAMFRver, eventsDir)
% function [] = RAM_FR_AddEventsToDatabase(subject,expDir,session,RAMFRver)
%
% FUNCTION:
%   RAM_FR_AddEventsToDatabase
%
% DESCRIPTION:
%   Add events in a RAM_FR folder to the events database. performs
%   all the checks to make sure that you do not screw up the events
%   database. 
%
% INPUT:
%   subject........ 'TJ057'
%   expDir......... path to the events that you want to add.  Examples:
%                   '/data/eeg/TJ057/behavioral/pyFR/TJ057_Olivia/session_0'
%                   '/data/eeg/TJ057/behavioral/pyFR/TJ057a_Olivia/session_0' 
%   session........ the session number that should be in the events dir
%   RAMFRver....... the version of RAM_FR that the events belong to
%                   e.g. 'RAM_FR1' or 'RAM_FR2',...
%   eventsDir...... (Optional) The directory in which the events are stored. 
%                   Defaults to: /data/events/
%
% OUTPUT:
%   saves the pyFR and the recognition and math events
%
% LAST UPDATED:
%   09/03/14 YE     created function from addPyFREventsToDatase by JFB
%

% get the location of the events database
if ~exist('eventsDir','var') || isempty(eventsDir)
    if isdir('/Volumes/rhino/data/events/')
        baseDir = '/Volumes/rhino';
    elseif isdir('/data/events/')
        baseDir = '';
    else
        error('can''t identify connection to rhino');
    end
    % get directories. locate the given version of RAM_FR
    eventsDir = fullfile(baseDir,'/data/events',RAMFRver);
else
    eventsDir = fullfile(eventsDir, RAMFRver);
end
if ~exist(eventsDir,'dir');
  fprintf('EXITING....Events directory does not exist\n\n')
  return
end

% check tp see if subject names match
fprintf('\n')
if isempty(regexp(expDir,subject))
  fprintf('  WARNNG: %s not found in %s\n',upper(subject),upper(expDir))
  fprintf('          you might be making an error.\n')
  fprintf('          please check this before making events. EXITING\n\n')
  fprintf('               !!! NO EVENTS SAVED !!! \n\n')
  return
end

% get the directories
 evFile         = fullfile(expDir,'events.mat');
mevFile         = fullfile(expDir,'MATH_events.mat');

% load the events
ev_RAM_FR = loadEvents_local(evFile,session); 
ev_math = loadEvents_local(mevFile,session);

% load the events in the events database
ev_RAM_FR_db = loadEventsDB_local(eventsDir,subject,session,'events'); 
ev_math_db = loadEventsDB_local(eventsDir,subject,session,'math');

% save the concatenated events
fprintf('\n\n')
fprintf([RAMFRver '  events: '])
if ~isempty(ev_RAM_FR)&&~isempty(ev_RAM_FR_db)  
  fprintf('adding new session\n')
  ev_RAM_FR_new = [ev_RAM_FR_db ev_RAM_FR];
  saveEventsDB_local(eventsDir,subject,'events',ev_RAM_FR_new);
elseif ~isempty(ev_RAM_FR)&&isempty(ev_RAM_FR_db)
    fprintf('no database events found for this subject...\ncreating new events struct in database from the current session\n');
    ev_RAM_FR_new = ev_RAM_FR;
    saveEventsDB_local(eventsDir,subject,'events',ev_RAM_FR_new);
else
  fprintf('not adding session\n')
end

fprintf('MATH  events: ')
if ~isempty(ev_math)&&~isempty(ev_math_db)
  fprintf('adding new session\n')
  ev_math_new = [ev_math_db ev_math];
  saveEventsDB_local(eventsDir,subject,'math',ev_math_new);
elseif ~isempty(ev_RAM_FR)&&isempty(ev_math_db)
  fprintf('no database events found for this subject...\ncreating new events struct in database from the current session\n');
  ev_math_new = ev_math;
  saveEventsDB_local(eventsDir,subject,'math',ev_math_new);
end

fprintf('\n\n')

%----------------------------------------------------
function saveEventsDB_local(eDir,subj,evStr,events);
  thisFileName = sprintf('%s_%s.mat',subj,evStr);
  thisFile     = fullfile(eDir,thisFileName);
  save(thisFile,'events')
  
%----------------------------------------------------
function ev = loadEventsDB_local(eDir,subj,sess,evStr); 
  thisFileName = sprintf('%s_%s.mat',subj,evStr);
  thisFile     = fullfile(eDir,thisFileName);
  if ~exist(thisFile,'file')
    ev=[];
    fprintf('%s does not exist in %s\n',thisFileName,eDir)
    return
  end
  
  % load the events
  ev_tmp = load(thisFile);
  ev     = ev_tmp.events;
  
  % check to make sure that the session doed not already exist in
  % the events on the database
  unSess  = unique([ev.session]);
  if sum(ismember(sess,unSess))>0    
    fprintf('%s session is already loaded in the events database\n',thisFile);
    ev = [];
    return
  end
    
%-----------------------------------------------
function ev = loadEvents_local(evFil,sessNum)
  [~,thisFile] = fileparts(evFil);
  if ~exist(evFil,'file')
    fprintf('%s does not exist\n',thisFile)
    ev = [];
    return
  end  
  
  % load the events
  ev_tmp = load(evFil);
  ev     = ev_tmp.events;
  
  % check to make sure that the events have the session that you
  % think that thet have
  unSess  = unique([ev.session]);
  numSess = length(unSess);
  if numSess~=1
    error(sprintf('I found more than one session in %s'),thisFile)
  end
  if unSess~=sessNum
    error(sprintf('%s: expecting session %d, found session %d'),...
	  sessNum,unSess)
  end
 
