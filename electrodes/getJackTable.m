function table = getJackTable(subj, rootEEGdir, varargin)
% GETJACKTABLE gets the jacksheetMaster.csv table
%
% Required Inputs
%   subj - e.g. NIH042
%   rootEEGdir - e.g. /Volumes/Shares/FRNU/data/eeg
%
% Optional key-value inputs:
%   createOnTheFly - If 1, create and return the jacktable instead of reading/writing from/to file
%
% Outputs
%   table
%
% Notes
%   Uses element_info.csv-type files
%   Side Effect: creates jacksheetMaster.csv if not found
%   
% Revision History
%   08/16 MST - Created
%   05/18 MST - Add createOnTheFly parameter

jackFile = fullfile(rootEEGdir, subj, 'docs', 'jacksheetMaster.csv');
useFormat = false;

ip = inputParser;
ip.addParameter('createOnTheFly', 0);
ip.parse(varargin{:});
createOnTheFly = ip.Results.createOnTheFly;

if ~exist(jackFile, 'file') || createOnTheFly
    table = createMasterJack(subj, rootEEGdir, 'skipFileCreationAndRawCheck',1);
    return
end

try
    % which columns in key are numeric format?
    key = readtableSafe(which('element_info_key.csv'),'readrownames',true);
    keyCols = key.Properties.VariableNames;
    keyFormat = key{'Format',:};
    numericMask = strcmpi(keyFormat, 'numeric');
    numericCols = keyCols(numericMask);
    
    % custom add "chanNum" as numeric
    numericCols = [numericCols {'chanNum'}];
    
    % which of of those numeric columns are in jacksheet?
    rawJack = readtableSafe(jackFile);
    jackCols = rawJack.Properties.VariableNames;
    numericCols = intersect(numericCols, jackCols);
    numericColMask = ismember(jackCols, numericCols);
    
    nCols = length(jackCols);
    charFormat = repmat({'%s'}, 1, nCols);
    numFormat = repmat({'%d'}, 1, nCols);
    
    % start format with all char %s. Then change numeric to %d
    parseFormat = charFormat;
    parseFormat(numericColMask) = numFormat(numericColMask);
    parseFormat = strjoin(parseFormat, '');
    
    useFormat = true;
catch e
end

if useFormat
    table = readtableSafe(jackFile, 'format',parseFormat);
    
else    
    table = readtableSafe(jackFile);
end

if numel(unique(table.subjId)) ~= 1
    disp(table);
    warning('subjId column is not unique. This is usually an indication that jacksheetMaster.csv''s row''s weres parsed incorrectly');
    if inputYN('Do you want to re-createMasterJack now?')
        createMasterJack(subj, rootEEGdir);
        table = getJackTable(subj, rootEEGdir);
        disp(table);
        return;
    end
end

end % getJackTable
