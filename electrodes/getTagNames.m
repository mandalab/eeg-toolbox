function tagNames = getTagNames(subj, rootEEGdir, varargin)
% GETTAGNAMES Returns the tagNames for a subject.
%
% Inputs
%   subj - e.g. NIH042
%   rootEEGdir - e.g. /Volumes/Shares/FRNU/data/eeg
% 
% Name-Value Pair Inputs  
%   chanType      - Filter by given chanTypes (default: {'PHYS'})
%                   can pass {PHYS,CLIN,SYNC} or 'ALL'
%
%   hardwareType - filter by hardwareType column. Default is 
%           {'subdural','micro-subdural','depth','microwire} (all execpt 'utah').
%           You can also pass 'ALL'
%           **NOTE: 'subdural' includes subduralHD
%
%   elementInfo - For speed, pass in an existing reference to elementInfo table
%           to avoid re-loading it
%
% Outputs
%   tagNames - 1 x n cell array of tag names. This is identical to what you
%           would find in the old tagNames.txt file.
%
% See also: getLeads
%
%
% 2/11/2020 - SJ: added functionality for recognizing subduralHD as subdural
%


chanTypeAll = {'PHYS','CLIN','SYNC','','NaN'};
hardwareTypeAll = {'subdural','micro-subdural','depth','microwire','NaN','utah',''};

% input parsing
p = inputParser;
p.addParameter('chanType', {'PHYS'}); %{'CLIN','PHYS','SYNC'}
p.addParameter('hardwareType', {'subdural','micro-subdural','depth','NaN',''}); % 'utah', 'microwire'
p.addParameter('elementInfo',[]);
parse(p, varargin{:});

chanType = cellstr(p.Results.chanType);
hardwareType = cellstr(p.Results.hardwareType);
info = p.Results.elementInfo;

if strcmpi(chanType, 'all'),     chanType     = chanTypeAll;     end
if strcmpi(hardwareType, 'all'), hardwareType = hardwareTypeAll; end


% variable declaration / setup
tagNames = {};
docsDir = fullfile(rootEEGdir, subj, 'docs');

if isempty(info)
    % load from element_info.csv
    filename = fullfile(docsDir, 'element_info.csv');
    if ~exist(filename, 'file')
       error('element_info.csv was not found here: %s\n', filename);
    end
    %verifyElementInfo(subj, rootEEGdir);
    info = readtableSafe(filename);
end

% Change 'subduralHD' to 'subdural' in element info  -SJ
if any(ismember(upper(info.hardwareType),upper('subduralHD')))
    info.hardwareType(ismember(upper(info.hardwareType),upper('subduralHD'))) = {'subdural'}; 
end

tagNames = info.tagName;
mask_chanType = ismember(upper(info.chanType), upper(chanType));
mask_hw = ismember(upper(info.hardwareType), upper(hardwareType)); %Cant use contains because that will return '' for everything
tagNames = tagNames(mask_chanType & mask_hw, :);

end % end function
