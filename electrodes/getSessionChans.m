function chans = getSessionChans(subj, rootEEGdir, timestamp, task, session)
% chans = getSessionChans(subj, rootEEGdir, timestamp, task, session)
% gets the channel names that were recorded for the given session based on
% raw_info.csv
%
% Ex. 1: chans = getSessionChans(subj, rootEEGdir, '160614_1008')
% Ex. 2: chans = getSessionChans(subj, rootEEGdir, [], 'paRemap',0)
%
% revision history
%   06/17 MST - Created
%
% see also: getSessionDate

    filename = fullfile(rootEEGdir, subj, 'docs', 'raw_info.csv');
    assert(exist(filename, 'file') > 0, 'File not found: %s\n', filename);
    rawinfo = readtableSafe(filename);
    
    if isempty(timestamp)
        assert(nargin >= 5, 'If timestamp is empty, must pass task and session params');
        timestamp = getSessionDate(subj, rootEEGdir, task, session);
    end
    col_name = ['raw_' timestamp];
    recorded = rawinfo.(col_name);
    chans    = rawinfo.chanName(find(recorded));
    

end