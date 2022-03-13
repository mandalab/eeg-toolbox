function leads = getLeads(subj, rootEEGdir, varargin)
% GETLEADS Returns the leads (channels) for a subject.
%
% Default behavior: PHYS channels, all hardwareTypes EXCEPT utah
%
% Lead definition: A "lead," also called a "channel," is a row from
%   jacksheetMaster. Leads are ordered in the order of elements in
%   element_invo.csv and are numbered 1-n, where n is the total number of
%   non-synch. The lead number is given in the recordIdx column of jacksheetMaster.csv.
%
% Required Inputs
%   subj - e.g. NIH042
%   rootEEGdir - e.g. /Volumes/Shares/FRNU/data/eeg
% 
% Name-Value Pair Inputs  
%
%   jackTable - For speed, pass in an existing reference to jacksheetMaster table
%               to avoid re-loading it
%
%   leadsType - Pass one of [good,bad,all] to get a specific lead type.
%               Default is all. If passing good/bad, also must pass timestamp
%               todo: intersection across all sessions
%
%   isBP    -   Pass true to get bipolar leads. False (default) for monopolar.
%
%   chanType -  Pass PHYS, SYNC, or CLIN to filter by chanType. Default is
%               {PHYS}. You can also pass 'ALL'
%
%   hardwareType - filter by hardwareType column. Default is 
%               {'subdural','micro-subdural','depth'} (all execpt 'utah','microwire').
%               You can also pass 'ALL'
%
%   spacing - filter by spacing column. Default is {'5','10'} (mm)
%
%   whichHemi - Default is {'rh','lh'}
%
%   timestamp - session timestamp (e.g. 120907_1342)
%
%   sessionCommon - if 1, return those channels which are common across ALL sessions
%               if you pass good, good across all
%               if you pass bad, bad across all
%
%   implant -   can be a dateIn or an integer (1 for first, 2 for 2nd) to filter on implant (e.g. to exclude shift)
%               default no filter. Implant filters on dateIn, so if you pass in 2, you will get *only* those channels
%               *added* for the 2nd implant (not those that may have persisted from the first)
%
%   markedAs -  Pass one of "isIctal, isInterictal, isResected" to get those channels matching one of these columns, as
%               marked in element_info.csv
%
% Outputs
%   leads   -   1 x n cell array of *strings* of channel names (or numbers if useNumeric).
%               For BP output, expect arrays to contain entries like '12 19' or
%               'TG2 TG3'
%
% Revision History
%   07/16 MST - Created
%   06/17 MST - Remove locType filter
%   08/17 JHW - remove "useNumeric" option... only using chanName indexing now
%   04/18 MST - Add "markedAs" for getting isIctal, isInterictal, isResected 
%   02/20 SNJ - replace 'subduralHD' with 'subdural' in jacktable when using it to filter by hardware type
%   04/21 SNJ - added extra to implant: check dateOut in addition to dateIn so you don't miss cases where
%               the dateIn is the same for everyone but dateOut is different (ie if 1 electrode gets
%               removed and all others remain in the same location). This may only be the case for the
%               1st run of tal before the shift is calculated (ex. NIH087)
%
% See Also: getGrids, getElementInfo, getJackTable
% 

% CONSTANTS
chanTypeAll = {'PHYS','CLIN','SYNC','','NaN'};
hardwareTypeAll = {'subdural','micro-subdural','depth','NaN','utah','microwire',''};
DATE_FORMAT = 'M/d/yy'; % e.g. read 5/4/14 as May 4 from element_info

% input parsing
p = inputParser;
p.addParameter('jackTable', []);
p.addParameter('leadsType',     'all');
p.addParameter('isBP',          false);
p.addParameter('stopOnWarning', true);
p.addParameter('chanType', {'PHYS'});
p.addParameter('hardwareType', {'subdural','micro-subdural','depth','NaN',''}); % 'utah','microwire'
p.addParameter('spacing', [], @isnumeric);
p.addParameter('whichHemi', {'rh','lh',''});
p.addParameter('timestamp', [], @ischar);
p.addParameter('sessionCommon', 0);
p.addParameter('implant',[]);
p.addParameter('markedAs', [], @ischar);

parse(p, varargin{:});

jacktable   = p.Results.jackTable; 
leadsType   = p.Results.leadsType;
isBP        = p.Results.isBP;
chanType    = cellstr(p.Results.chanType);
hardwareType= cellstr(p.Results.hardwareType);
spacing     = p.Results.spacing;
whichHemi   = p.Results.whichHemi;
timestamp   = p.Results.timestamp;
sessionCommon = p.Results.sessionCommon;
implant     = p.Results.implant;
markedAs    = p.Results.markedAs;

if strcmpi(chanType, 'all'), chanType = chanTypeAll; end
if strcmpi(hardwareType, 'all'), hardwareType = hardwareTypeAll; end

assert(all(ismember(hardwareType, hardwareTypeAll)), 'hardwareType not found in list of all types');
assert(all(ismember(chanType, chanTypeAll)), 'chanType not found in list of all types');

% variable declaration / setup
leads   = {};

% load from jasheetMaster.csv

if isBP
    %% biploar section
    filename = fullfile(rootEEGdir, subj, 'docs', 'leads_bp.txt'); %- 4/2018, leads_bp moved to docs
    if exist(filename, 'file')
        fd = fopen(filename);
        assert(fd > 0, 'leads_bp failed to open: %s', filename);
        leads = textscan(fd,'%d%d', 'delimiter','-');
        if isempty(leads{1}) || isempty(leads{2})
            fseek(fd, 0, 'bof');
            leads = textscan(fd,'%s%s', 'delimiter','-');
            leads = [leads{1} leads{2}];
        end
        fclose(fd);
    else
        % create leads_bp.txt
        leads = extractBipolarity(subj, rootEEGdir);
    end
    
    if strcmpi(leadsType, 'good') || strcmpi(leadsType, 'bad')
        assert(~isempty(timestamp) || sessionCommon, 'GetLeads: good/bad leads decided on a session-by-session basis; must pass timestamp parameter');
        fname_variance = fullfile(rootEEGdir, subj, 'eeg.processed', timestamp, 'variance.csv');
 
        assert(exist(fname_variance, 'file') > 0 || sessionCommon, 'Variance file not found: %s\n', fname_variance);
        t = readtable(fname_variance);
        bad_chans = t.chanName(~logical(t.is_good));
        if strcmpi(leadsType, 'bad')
            leads = leads(any(ismember(leads, bad_chans), 2), :); % either pair can be bad
        elseif strcmpi(leadsType, 'good')
            leads = leads(all(~ismember(leads, bad_chans), 2), :); % both pairs must be good
        end
        fprintf('** Important Note: nonsubdurals are not included in good/bad **\n');
    end
    
    if ~isempty(timestamp)
        fraw = fullfile(rootEEGdir, subj, 'docs/raw_info.csv');
        assert(exist(fraw, 'file') > 0, 'Cannot find raw_info: %s\n', fraw);
        traw = readtableSafe(fraw);
        cols = traw.Properties.VariableNames;
        raw_col_mask = strncmpi(cols, 'raw_', length('raw_')) & ~strcmp(cols, 'raw_info');
        raw_cols = cols(raw_col_mask);
        assert(ismember(['raw_' timestamp], traw.Properties.VariableNames), 'Timestamp not in raw info: %s', timestamp);
        mask = traw{:,['raw_' timestamp]}; % matrix
        sessionLeads = traw.chanName(logical(mask));
        leads = leads(ismember(leads(:,1),sessionLeads) & ismember(leads(:,2),sessionLeads), :);        
    end
    
    if ~isempty(whichHemi)
        if isempty(jacktable)
            jackfile = fullfile(rootEEGdir, subj, 'docs', 'jacksheetMaster.csv');
        if ~exist(jackfile, 'file')
            createMasterJack(subj, rootEEGdir);
        end

        jacktable = getJackTable(subj, rootEEGdir);
        end
    
        try
            jacktable = jacktable(ismember(jacktable.whichHemi,whichHemi), :);
        catch er
            fprintf('getLeads error filtering: %s\nLikely cause: column has bad values for whichHemi\nNo whichHemi filter applied!\n\n', er.message);
        end
        leads = leads(all(ismember(leads,jacktable.chanName),2),:);
    end
else
    %% Monopolar section
    % get the jacktable if it wasn't passed in
    if isempty(jacktable)
        jackfile = fullfile(rootEEGdir, subj, 'docs', 'jacksheetMaster.csv');
        if ~exist(jackfile, 'file')
            createMasterJack(subj, rootEEGdir);
        end

        jacktable = getJackTable(subj, rootEEGdir);
    end
    
    % Change 'subduralHD' to 'subdural' in jacktable  -SJ
    jacktable.hardwareType(ismember(upper(jacktable.hardwareType),upper('subduralHD'))) = {'subdural'};
    
    %== Filters ==%
    jacktable = jacktable(ismember(upper(jacktable.chanType), upper(chanType)), :);
    
    try
        jacktable = jacktable(ismember(upper(jacktable.hardwareType),upper(hardwareType)), :);
    catch er
        fprintf('getLeads error filtering: %s\nLikely cause: column has bad values for hardwareType\nNo hardwareType filter applied!\n\n', er.message);
    end

    if ~isempty(spacing)
        try
            jacktable = jacktable(ismember(jacktable.spacing,spacing), :);
        catch er
            fprintf('getLeads error filtering: %s\nLikely cause: column has bad values for spacing\nNo spacing filter applied!\n\n', er.message);
        end
    end
    
    if ~isempty(whichHemi)
        try
            jacktable = jacktable(ismember(jacktable.whichHemi,whichHemi), :);
        catch er
            fprintf('getLeads error filtering: %s\nLikely cause: column has bad values for whichHemi\nNo whichHemi filter applied!\n\n', er.message);
        end
    end
    
    if ~isempty(implant) && ~isempty(jacktable) && ~all(cellfun('isempty',jacktable.dateIn))
        dates_in = [];
        date_str_in = jacktable.dateIn;
        for i = 1:numel(jacktable.dateIn)
            try 
                dates_in = [dates_in, datetime(date_str_in{i}, 'format',DATE_FORMAT)]; 
            catch
            end
        end
        
        if numel(dates_in) > 0
            dates_in = dates_in(~isnat(dates_in));
            dates_in = sort(unique(dates_in));
            assert(numel(dates_in) <= 2, 'getLeads needs to be rewritten in the case that >= 2 dateIn values');
            if numel(dates_in) > 0
                if ischar(implant)
                    d = datetime(implant, 'format',DATE_FORMAT);
                    implant = find(d == dates_in);
                    mask = datetime(jacktable.dateIn, 'format',DATE_FORMAT) == d;
                elseif numel(dates_in) < implant
                    %SJ adding: check dateOut
                    % If you only have 1 dateIn, then you need to check whether there is a second set
                    % from the dateOuts
                    dates_out = [];
                    date_str_out = jacktable.dateOut;
                    for i = 1:numel(jacktable.dateOut)
                        try
                            dates_out = [dates_out, datetime(date_str_out{i}, 'format',DATE_FORMAT)];
                        catch
                        end
                    end
                    dates_out = dates_out(~isnat(dates_out));
                    dates_out = sort(unique(dates_out));
                    if numel(dates_out) < implant
                        warning('getLeads: You wanted implant #%d but there is only %d unique date in the jacksheet--update your element_info date fields', implant, numel(dates_in));
                        fprintf('Proceeding using the only date in element_info: %s\n', datestr(dates_in));
                        d = dates_in;
                        mask = datetime(jacktable.dateIn, 'format',DATE_FORMAT) == d;
                    elseif numel(dates_out) == implant % Only grab the electrodes with the second dateOut 
                        d = dates_out(implant);
                        mask = datetime(jacktable.dateOut, 'format',DATE_FORMAT) == d;
                    else
                        keyboard % Do we get here?
                    end
                else
                    d = dates_in(implant);
                    mask = datetime(jacktable.dateIn, 'format',DATE_FORMAT) == d;
                end

                if implant == 2
                    dout = jacktable.dateOut;
                    mask = mask | (isempty(dout) | dout > d);
                end

                jacktable = jacktable(mask, :);
            end
        end
        
    end
    
    % MarkedAs
    if ~isempty(markedAs)
        if ~contains(lower(markedAs),'is')
            markedAs = strcat('is', markedAs);
        end

        assert(ismember(lower(markedAs), {'isinterictal', 'isictal', 'isresected'}), 'markedAs must be one of isInterictal, isIctal, isResected');
        
        if iscellstr(jacktable.isIctal)
            toLogical = @(c) logical(str2num(cell2mat(c)));
        else
            toLogical = @(c) logical(cell2mat(c));
        end
        
        switch lower(markedAs)
            case 'isictal'
                mask = toLogical(jacktable.isIctal);
            case 'isinterictal'
                mask = toLogical(jacktable.isInterictal);
            case 'isresected'
                mask = toLogical(jacktable.isResected);
            otherwise
                error('Not handled %s', markedAs);
        end
        assert(numel(mask) == height(jacktable));
        jacktable = jacktable(mask, :);
    end
    
    
    %== End Filters ==%
    
    
    leads = jacktable.chanName;
    
    
    if sessionCommon
        %get only those leads which are common across all sessions' noreref folder
        fraw = fullfile(rootEEGdir, subj, 'docs/raw_info.csv');
        assert(exist(fraw, 'file') > 0, 'Cannot find raw_info: %s\n', fraw);
        traw = readtableSafe(fraw);
        cols = traw.Properties.VariableNames;
        raw_col_mask = strncmpi(cols, 'raw_', length('raw_')) & ~strcmp(cols, 'raw_info');
        raw_cols = cols(raw_col_mask);
        nraws = numel(raw_cols);
        mask = traw{:,raw_cols}; % matrix
        common_mask = sum(mask,2) == nraws;
        common = traw.chanName(common_mask);
        leads = intersect(leads, common, 'stable');
        
        if strcmpi(leadsType, 'good') || strcmpi(leadsType, 'bad')
            procDir = fullfile(rootEEGdir, subj, 'eeg.processed');
            assert(exist(procDir, 'dir') > 0, 'getLeads: eeg.processed directory must exist to use good/bad. not found: %s\n', procDir);
            contents =  lsCell(procDir);
            timestamps = contents(cellfun(@isdir, fullfile(procDir,contents)));
            
            for i = 1 : length(timestamps)
                fname_variance = fullfile(rootEEGdir, subj, 'eeg.processed', timestamps{i}, 'variance.csv');
                assert(exist(fname_variance,'file') > 0 , 'variance file does not exist: %s\n', fname_variance);
                t = readtable(fname_variance);
                bad_chans = t.chanName(~logical(t.is_good));
                if i == 1
                    if strcmpi(leadsType, 'good')
                        leads = setdiff(leads, bad_chans, 'stable');
                    else
                        % throw away all good channels
                        leads = intersect(leads, bad_chans, 'stable');
                        all_bad = bad_chans;
                    end
                else
                    if strcmpi(leadsType, 'good')
                        % Gives us channels that are good in all sessions by removing any bad channels at each session
                        leads = setdiff(leads, bad_chans, 'stable');
                    else
                        % Gives us channesl that are bad in any session
                        all_bad = union(all_bad, bad_chans, 'stable');
                        if i == length(timestamps)
                            leads = intersect(leads, all_bad, 'stable');
                        end
                    end
                end
                fprintf('** getLeads note: non-subduralshave no good/bad distinction **\n');
            end
            
        end
        
        
    else
        % timestamp-specific: get only those leads which are in a session's noreref folder
        if ~isempty(timestamp)
            splitDir = fullfile(rootEEGdir, subj, 'eeg.noreref');
            assert(exist(splitDir, 'dir') > 0, 'getLeads: Noreref directory must exist to use timestamp. not found: %s\n', splitDir);
            sessionDir = fullfile(splitDir, timestamp);
            assert(exist(sessionDir, 'dir') > 0, 'getLeads: timestamp given does not exist: %s\n', sessionDir);
            rawFilenames = lsCell(sessionDir);
            leads = intersect(leads, rawFilenames, 'stable');
        end

        if strcmpi(leadsType, 'good') || strcmpi(leadsType, 'bad')
            assert(~isempty(timestamp), 'GetLeads: good/bad leads decided on a session-by-session basis; must pass timestamp parameter');
            fname_variance = fullfile(rootEEGdir, subj, 'eeg.processed', timestamp, 'variance.csv');

            assert(exist(fname_variance, 'file') > 0, 'Variance file not found: %s\n', fname_variance);
            t = readtable(fname_variance);
            bad_chans = t.chanName(~logical(t.is_good));
            if strcmpi(leadsType, 'good')
                leads = setdiff(leads, bad_chans, 'stable');
            elseif strcmpi(leadsType, 'bad')
                leads = intersect(leads, bad_chans, 'stable');
            end
            fprintf('** getLeads note: non-subduralshave no good/bad distinction **\n');
        end
    end
end
    


end % end function
