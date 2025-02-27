function subj = DBSalign_subj(subj, session, varargin)
%ALIGN_SUBJ   Align a subject's events to EEG data.
%
%  subj = align_subj(subj, ...)
%
%  INPUTS:
%    subj:  a subject structure, where each sess subfield has the
%           following fields:
%            dir     - directory where behavioral data is stored
%            eegfile - path to the EEG files, including the filename,
%                      without the .XXX suffix that indicates the
%                      channel
%
%  OUTPUTS:
%     subj:  the subject structure, unmodified. This is only output for
%            purposes of compatibility with apply_to_subj.
%
%  PARAMS:
%  These options may be specified using parameter, value pairs or by
%  passing a structure. Defaults are shown in parentheses.
%   eventfile - path (relative to each sess.dir) to the events
%               structure. ('events.mat')
%   pulse_ext - file extension for the pulse files, to be appended to
%               the EEG fileroot for each session. May contain wildcards
%               (*). ('*.sync.txt')
%   pulse_dir - directory where pulse files are stored, relative to
%               subj.dir. ('eeg.noreref')
%
%  EXAMPLE:
%   % path to directory with behavioral data
%   clear subj
%   subj.sess.dir = '/data/eeg/TJ003/behavioral/iCatFR/session_0';
%
%   % path to an EEG channel file, minus the .XXX extension
%   subj.sess.eegfile = '/data/eeg/TJ003/eeg.reref/TJ003_18Feb09_1335';
%   pulse_dir = '/data/eeg/TJ003/eeg.noreref';
%
%   % load [subj.sess.dir]/events.mat, align, add eegoffset and eegfile
%   % fields, and resave
%   align_subj(subj, 'pulse_dir', pulse_dir);
%
%  NOTES:
%   This script makes various assumptions about directory structure.
%   Each sess.dir must contain:
%    eeg.eeglog.up OR eeg.eeglog
%
%   The directory containing each EEG file must contain a file called
%   params.txt which has information about samplerate.
%
%   There must be a EEG channel file [sess.eegfile].001 for each
%   session.

% subj.sess.dir = '/Users/dongj3/Jian/data/DBS/TTL/behavioral';
% subj.sess.eegfile ='/Users/dongj3/Jian/data/DBS/TTL/eeg.noreref/TTL_090911_1200'; 
% subj.sess.pulse_dir = '/Users/dongj3/Jian/data/DBS/TTL/eeg.noreref/';

% input checks
if ~exist('subj', 'var')
  error('You must pass a subject structure.')
end

% options
defaults.eventfile = 'events.mat';
defaults.pulse_ext = '*.sync.txt';
%defaults.pulse_dir = 'eeg.noreref';
params = propval(varargin, defaults);

for sess=subj.sess
  % read the eegfile(s) from the sess structure
  if ~isfield(sess, 'eegfile')
    error('Each session must have an "eegfile" field.')
  end
  eegfiles = sess.eegfile;
  if ~iscell(eegfiles)
    eegfiles = {eegfiles};
  end

  % get the events structure: events structure has to be in behavioral
  % folder
  eventfile = fullfile(sess.dir,['session_' session], params.eventfile);
  if ~exist(eventfile, 'file')
    error('Events file not found: %s\n', eventfile)
  end
  
  % get sync files
  eegsyncfiles = cell(1, length(eegfiles));
  for i=1:length(eegfiles)
    [pathstr, basename, ext] = fileparts(eegfiles{i});
    basename=[basename ext];

    % get the EEG sync pulse file
    pulse_path = fullfile(sess.pulse_dir, ...
                          [basename params.pulse_ext]);
    temp = dir(pulse_path);
    if length(temp) == 0
      warning('No EEG sync pulse files found that match: %s', ...
              pulse_path);
      return;
    elseif length(temp) > 1
      fprintf('Multiple EEG sync pulse files found that match: %s\n', pulse_path)
      % use the first one that matches this pattern: '###.sync.txt'
      for j=1:length(temp)
        if length(temp) > 1
          if ~isempty(regexp(temp(j).name,'\w*\d\d\d.sync.txt'))
            temp = temp(j);
            fprintf('Using: %s\n', temp.name);
          end
        end
      end
      if length(temp) > 1
        error('No valid sync files found.');
      end
    end

    eegsyncfiles{i} = fullfile(sess.pulse_dir, temp.name);
    % for runAlign, make eegfile point to a specific channel
    eegfiles{i} = [eegfiles{i} '.001'];
  end

  % there should be only one behavioral sync pulse file
  behsyncfile = fullfile(sess.dir,['session_' session],'eeg.eeglog.up');
  if ~exist(behsyncfile, 'file')
    % if we haven't already, extract the UP pulses
    raw_behsyncfile = fullfile(sess.dir,sprintf('session_%s/eeg.eeglog',session));
    if ~exist(raw_behsyncfile,'file')
      error('Behavioral pulse file not found: %s\n', raw_behsyncfile)
    end
    fixEEGLog(raw_behsyncfile, behsyncfile);
  end
  
  % get the samplerate
  samplerate = GetRateAndFormat(fileparts(eegfiles{1}))
  % run the alignment
  
  try
    runAlign(samplerate, {behsyncfile}, eegsyncfiles, eegfiles, {eventfile}, ...
             'mstime', 0, 0);
  catch err
    fprintf('runAlign threw an error.\n');
    getReport(err)
  end
end
