function [eeg_avg, eeg_all] = reref(fileroot, refChans, masterLeadNames, eeg_all)
%REREF - Rereference EEG channels with respect to average of given channels
%
% FUNCTION:
%   [eeg_avg, eeg_all] = reref( fileroot, refChans [, masterLeadNames, eeg_all] )
%
% INPUTS:
%   fileroot:  timestamped directory that contains raw channels [rootEEG]/[subj]/eeg.processed/[timestamp]
%   refChans   :  channels to reference to
%
% OPTIONAL_INPUT
%   masterLeadNames: call getLeads(subj, rootEEGdir, chanType','PHYS', 'leadsType','all') ahead of time to save file I/O time
%   eeg_all :  if you have the matrix of all lead's timeseries already, you can pass it to avoid reading files
%
% OUTPUTS:
%   eeg_avg : timeseries of average of refChans
%   eeg_all : chan x nsamples matrix. Here length(chans)=length(masterLeadNames)
%
% Revision History:
%   05/17 MST - Created
%
% SEE ALSO: reref_legacy, rerefWrapper

REREF_DIR = 'eeg.reref';

subjDir = fileparts(fileparts(fileroot));
[rootEEGdir,subj] = fileparts(subjDir);
[~, eegFileStem] = fileparts(fileroot);

% Get channel names
if nargin < 3, masterLeadNames = getLeads(subj, rootEEGdir, 'leadsType','all','chanType','PHYS'); end
if nargin < 4, eeg_all = []; end

% create dir
if ~exist(REREF_DIR,'dir')
    mkdir(REREF_DIR);
end

fprintf('Processing %s...\n',fileroot);

% Instead of using all electrodes, use only those that were recorded for a
% specific session (from that session's jacksheet.txt)
% *note: these numbers are based on jacksheetMaster.csv
[sessionChans] = nkElectrodeFilt(rootEEGdir, subj, eegFileStem);
%if isempty(sessionChans)
%    sessionChans = masterLeadNames;  %- always returns something as of 12/2018
%end

if ~isempty(eeg_all)
    assert(size(eeg_all,1) == length(sessionChans), 'rows in eeg_all must match # session channels (%d)', length(sessionChans));
end

% filter based on given ref chans
chan_mask       = ismember(upper(sessionChans), upper(refChans));
if ~all(ismember(refChans, sessionChans))
    error('reref: channels were passed that were not found by getLeads (masterLeadNames)');
end

% Load all channel eeg's from split files in fileroot
if isempty(eeg_all)
    
    % get data info
    [~,~,dataformat] = GetRateAndFormat(fileroot);

    % load channels
    malloc = @(m,n,format) repmat(eval(sprintf('%s(0)', format)), [m n]); % memory allocation
    nsamples   = [];
    for i = 1:length(sessionChans)
        chan = sessionChans{i};
        filename = fullfile(fileroot, chan);
        assert(exist(filename, 'file') > 0, 'File not found: %s', filename);
        
        fprintf('%s ',chan);
        fid = fopen(filename, 'r');
        eeg = fread(fid, inf, dataformat); 
        try fclose(fid); catch, end;
        if mod(i,20)==0, fprintf('\n'); end

        % allocate on first loop when we know nsamples
        if isempty(eeg_all)
            nsamples = length(eeg);
            eeg_all = malloc(length(sessionChans), nsamples, dataformat);
        else
            assert(length(eeg) == nsamples, 'Not all channels have same # samples');
        end

        eeg_all(i,:) = eeg;
    end
end
fprintf('\n');
eeg_avg = mean(eeg_all(chan_mask,:), 1); % mean along rows (chans)

    
end %reref
