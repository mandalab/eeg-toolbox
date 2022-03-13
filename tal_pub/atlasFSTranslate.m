function expanded_name = atlasFSTranslate(atlas, fs_name, varargin)
% ATLASFSTRANSLATE translates the Freesurfer atlas to its full name
%
% expanded_name = atlasFSTranslate('destrieux', 'calcarine', 'translation','pretty')
%
% Freesurfer outputs an atlas which is abbreviated. This function
% translates the abbreviated form to its full name
%
% INPUTS:
%   atlas - One of: {'destrieux','desikan','aseg'}
%   fs_name - atlas label name to look up
%
% OPTIONAL KEY-VALUE INPUTS:
%   silent - if 0 (default), show error for labels not found
%   translation - One of: {'full', 'pretty'} (default 'full')
%
% Currently implemented translations:
%   Destrieux Atlas (surface parcellation, "aparc.a2009s")
%   Desikan-killiany atlas (surface parcellation, "aparc")
%   Aseg Atlas (subcortical segmentation, "aseg")
%
% HOW IT WORKS:
%   *Translations are looked up from CSV files in: eeg_toolbox/tal_pub/
%       desikan_lookup.csv
%       destrieux_lookup.csv
%       aseg_lookup.csv
%
%   Lookup key is in a column named "short", value column has the name that user gives in the translation parameter
%       
%
% Resources:
%   Freesurfer atlas outputs: 
%       https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
%   Destrieux parc table:
%       "Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature"
%       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2937159/
%       $FREESURFER_HOME/Simple_surface_labels2009.txt
%   Destrieux table:
%   Aseg table:
%   
%       
%
% REVISION HISTORY
%   ??/?? MST - Created
%   07/18 MST - Added pretty, changed the parameters
%
ip = inputParser();
ip.addParameter('silent', 0);
ip.addParameter('translation', 'full');
ip.parse(varargin{:});
silent = ip.Results.silent;
translation = ip.Results.translation;

expanded_name = [];

if ~ismember(atlas, {'destrieux','desikan','aseg'})
    error('Atlas %s not available. Try destrieux, desikan, or aseg', atlas);
end

% Find look-up table
fname = sprintf('%s_lookup.csv', atlas);
fpath = which(fname);
if isempty(fpath)
    tbd = getToolboxDir();
    fpath = fullfile(tbd, 'tal_pub', fname);
end
assert(exist(fpath, 'file') > 0, 'Cannot locate atlas: %s', fpath);

% Load table
t_lookup = readtableSafe(fpath);

assert(ismember(translation, t_lookup.Properties.VariableNames), ...
    'Translation %s is not a valid column in %s', ...
    translation, fpath);

% Find given atlas label row
match_ndx = find(strcmpi(fs_name, t_lookup.short), 1);
hem = 'lr';
ihem = 1;
while isempty(match_ndx) && ihem <= 2
    ndx = max([strfind(fs_name, [hem(ihem) 'h_']), strfind(fs_name, [hem(ihem) 'h-'])]);
    if ~isempty(ndx)
        fs_name = fs_name(ndx(1) + length([hem(ihem) 'h_']) : end);
        match_ndx = find(strcmpi(fs_name, t_lookup.short), 1);
    end
    ihem = ihem + 1;
end
if isempty(match_ndx)
    
end

if isempty(match_ndx)
    if ~silent
        fprintf('atlasFSTranslate Error: atlas %s not found in lookup table (%s)\n', fs_name, fpath);
    end
else
	expanded_name = strtrim(t_lookup.(translation){match_ndx});
end


end % end function 