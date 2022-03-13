function extract(subj, rootEEGdir)
% EXTRACT begins the extraction process from element_info.cvs.
%
% REQUIRED INPUT
%   subj - e.g. NIH042
%   rootEEGdir - e.g. /Volumes/Shares/FRNU/data/eeg
% 
% NAME-VALUE PAIR OPTIONAL INPUT
%         
%
% OUTPUT
%   
%   Files created:
%    - docs/jacksheetMaster.csv
%    - eeg.noreref/(extracted channel data)
%    - eeg.reref/(rereferenced channel data)
%
% NOTES
%   
%   
% REVISIONS
%   07/16 MST - Created
% 
% See also

% CONSTANTS

% input parsing


% variable declaration / setup
patientDir = fullfile(rootEEGdir, subj);

% Check file exists
filename = fullfile(patientDir, 'docs', 'element_info.csv');
if ~exist(filename, 'file')
    fprintf(' ERROR: element_info.csv not found here: %s\n', filename);
    keyboard;
    return;
end

verifyElementInfo(subj, rootEEGdir);
eegPrepAndAlign(subj, rootEEGdir);

end % end function