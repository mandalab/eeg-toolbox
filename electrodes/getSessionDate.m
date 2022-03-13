function timestamp = getSessionDate(subj, rootEEGdir, task, session)
% timestamp = getSessionDate(subj, rootEEGdir, task, session)
% Gets the timestamp YYmmDD_hhMM of a task session
%
% data source: events.mat, alignment
%
% Inputs
%   subj    -
%   rootEEGdir -
%   task    - e.g. 'moveTask'
%   session - e.g. 0
%
% % revision history
%   06/17 MST - Created
%
% see also: getSessionChans

    % turn on verbose for print statements
    VERBOSE = 0;
    print = @(x,varargin) fprintf(turnary(VERBOSE,x,''),varargin{:});

    if isnumeric(session), session = sprintf('session_%d', session); end
    
    d = fullfile(rootEEGdir, subj, 'behavioral', task, session);
    assert(exist(d, 'dir') > 0, 'Not found: %s\n', d);
    
    align_timestamp = [];
    event_timestamp = [];
    
    fname_events = fullfile(d, 'events.mat');
    if exist(fname_events, 'file')
        data = load(fname_events);
        event = data.events(1);
        if isfield(event, 'eegfile')
            eegfile = event.eegfile;
            [~,event_timestamp] = fileparts(eegfile);
        else
            print('eegfile is not a field of events\n');
        end
        
        
    else
        print('No events file: %s\n', fname_events);
    end
        
    fname_align = fullfile(rootEEGdir, subj, 'behavioral', 'alignmentSummary.txt');
    if exist(fname_align, 'file')
        line = get_line_w_text(fname_align, sprintf('%s/%s',task,session));
        words = strsplit(line);
        words = words(~cellfun('isempty',words));
        % 1st word = task/session
        % 2nd word = date
        % 3rd word = time
        % 4th word = AM/PM
        try
            ndate = datenum([words{2} ' ' words{3} ' ' words{4}]);
            align_timestamp = datestr(ndate, 'YYmmDD_hhMM');
        catch e
            disp(e.message);
            print('Bad date for line: %s\n', line);
        end        
        
    else
        print('file not found: %s\n', fname_align);
    end
    
    if ~isempty(event_timestamp)
        timestamp = event_timestamp;
    end
    
    print('timestamp (events.mat):    %s\n', event_timestamp);
    print('timestamp (align_summary): %s\n', align_timestamp);
    print('timestamp returned       : %s\n', timestamp);
    
    if ~isempty(event_timestamp) && ~isempty(align_timestamp) && ~strcmpi(event_timestamp, align_timestamp)
        warning('Event and align timestamps do not match (%s vs %s)', event_timestamp, align_timestamp);
    end
    
end

function line = get_line_w_text(filename, srch)
    % returns first line with given srch strings on it
    if ischar(srch) && ~iscellstr(srch), srch={srch}; end

    fid = fopen(filename, 'r');

    while 1
        line = fgetl(fid);
        if ~ischar(line), break; end;
        
        % search this line for all words
        found = 1;
        for i = 1:length(srch)
            if isempty(strfind(line, srch{i}))
                found = 0;
            end
        end
        
        if found, break; end % return current line if all found
    end
end

% function line = get_line_w_text(filename, srch)
%     % returns first line with given srch strings on it
%     if ischar(srch) && ~iscellstr(srch), srch={srch}; end
% 
%     fid = fopen(filename, 'r');
% 
%     while 1
%         line = fgetl(fid);
%         if ~ischar(line), break; end;
%         %search_str = sprintf('%s/%s', task,session);a
%         fnd_ndx = strfind(line, srch);
%         if fnd_ndx > 0
%             break;
%             words = strsplit(line);
%             words = words(~cellfun('isempty',words));
%             % 1st word = task/session
%             % 2nd word = date
%             % 3rd word = time
%             % 4th word = AM/PM
%             try
%                 ndate = datenum([words{2} ' ' words{3} ' ' words{4}]);
%                 align_timestamp = datestr(ndate, 'YYmmDD_hhMM');
%                 break;
%             catch e
%                 disp(e.message);
%                 print('Bad date for line: %s\n', line);
%             end
% 
%         end
%     end
% end
