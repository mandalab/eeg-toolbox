function [stemsChecked, stemsRerefed] = rerefWrapper(subj,eegDir, varargin)
% REREFWRAPPER is DEPRECATED
%
% Function rerefWrapper(subj,eegDir)
%
%   Description: This function implements both bipolar and average
%   re-referencing schemes. It gets gridLayout and
%   missingElecs from the electrodes.m file to be created in the patient's docs directory.
%
%   Input:
%         --subj=subject name (ex. 'NIH001')
%         --eegDir=directory where data is held
%                  (ex. '/Users/damerasr/Sri/data/eeg/')
%
%   Key-Value Pair Input
%         --useOldSource, [true/false] = whether to use old leads.txt files
%         (old) or new csv files (true, default)
%
%   Output:
%         -- number of EEG stems checked, and number of those rereferenced
%
% Revision History
%   05/17 MST - Deprecated, see PROCESSANDREREF
%
% See Also: PROCESSANDREREF

% read input
p = inputParser;
p.addParameter('useOldSource', false, @islogical);
parse(p, varargin{:});
useOldSource = p.Results.useOldSource;

warning('rerefWrapper is DEPRECATED. processAndReref is its replacement');

% sets up input and output directories
subjDir    = fullfile(eegDir,subj); % where subject data is located
taldir     = fullfile(subjDir,'tal'); % where tailrach info is located
noRerefDir = fullfile(subjDir,'eeg.noreref'); % where no-reref data is located
rerefDir   = fullfile(subjDir,'eeg.reref'); % where reref data will be saved to

%if the reref direcetory does not exist it is now made
if ~exist(rerefDir,'dir')
    mkdir(rerefDir)
end


% gets unique fileStems from the files in the no-reref direcetory
noreref_dates = lsCell(noRerefDir);
noreref_dates = noreref_dates(cellfun(@isdir, fullfile(noRerefDir, noreref_dates)));
reref_dates = lsCell(rerefDir);
reref_dates = reref_dates(cellfun(@isdir, fullfile(noRerefDir, reref_dates)));


rerefCount = 0;

% for each fileStem run bipolar and 'laplacian' re-referencing schemes and save out the data
for i = 1:length(noreref_dates)
    noreref_date = noreref_dates{i};
    
    if ~ismember(noreref_date, reref_dates) % if filestem hasn't already been rereferenced then do it now
        mkdir(fullfile(rerefDir, noreref_date));
        
        % bipolar rereferencing
        grids = reref_Bipolarity(subjDir, noreref_date, 'useOldSource', useOldSource);
            
        % global weighted average rerefrencing
        fileroots = {fullfile(noRerefDir,noreref_date)};
        reref(fileroots,grids,rerefDir,taldir, 'useOldSource', useOldSource);
        
        rerefCount = rerefCount + 1;
    end
end


stemsChecked = length(noreref_dates);
stemsRerefed = rerefCount;

end % rerefWrapper

function yes = isExpectedFilename(filename)
% check if file has subj_date_time format (look for 2 '_' chars)
    yes = (2 == length(strfind(filename, '_')));
end