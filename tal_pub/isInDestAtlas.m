function [found, dest_labels] = isInDestAtlas(label, dest_labels)
    % isInDestAtlas(label) returns True if the given label is in the Destrieux atlas
    %
    % Description:
    %   We store a file in the toolbox called eeg_toolbox/tal_pub/atlas_labels.mat
    %   This contains a struct atlas_labels, with fields for each atlas. For example:
    %       .destrieux.labels 
    %   
    %   This labels fields contain all labels plus an "unknown"
    %
    %
    % More Info:
    %   https://projects.ninds.nih.gov:8090/display/Home/For+Users
    %   https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
    %
    % Revision History
    %   07/18 MST   - Created
    %
    % See Also: localizer_atlas, atlas_simplify, atlasFSTranslate, aparc2lobe, isInDKAtlas
   
    
    if nargin < 1, label = ''; end
    if nargin < 2, dest_labels = []; end
    
    if ~isempty(dest_labels)
        % You already have the labels!
        % This parameter, which is returned and accepted, is for the calling
        % function to cache the value so as to avoid repeated file I/O
        found = ismember(lower(label), lower(dest_labels));
        return
    end
    
    toolbox_dir = [];
    
    % Get definition file
    fname_atlas_labels = which('atlas_labels.mat');
    if isempty(fname_atlas_labels)
        toolbox_dir = getToolboxDir();
        fname_atlas_labels = fullfile(toolbox_dir, 'tal_pub/atlas_labels.mat');
    end
    
    % Check the file for the dest atlas
    good_file = 0;
    if exist(fname_atlas_labels, 'file')
        data = load(fname_atlas_labels);
        if isfield(data, 'atlas_labels')
            atlas_labels = data.atlas_labels;
            if isfield(atlas_labels, 'destrieux') && ~isempty(atlas_labels.destrieux)
                good_file = 1;
            end
        end
    end
    
    if ~good_file
        atlas_labels = create_dest_label(toolbox_dir);
    end
    
    % do the lookup
    found = ismember(label, lower(atlas_labels.destrieux));
    dest_labels = atlas_labels.destrieux;

end


function atlas_labels = create_dest_label(toolbox_dir)
    % manually created csv file to our dest atlas file
    if nargin < 1, toolbox_dir = []; end
    
    fname_dest_colorLUT = which('destrieux_lookup_surf.csv');
    fname_atlas_labels  = which('atlas_labels.mat');
    if isempty(fname_dest_colorLUT) || isempty(fname_atlas_labels)
        if isempty(toolbox_dir)
            toolbox_dir = getToolboxDir();
        end
        fname_dest_colorLUT = fullfile(toolbox_dir, 'tal_pub/destrieux_lookup_surf.csv');
        fname_atlas_labels  = fullfile(toolbox_dir, 'tal_pub/atlas_labels.mat');
    end
    
    assert(exist(fname_dest_colorLUT, 'file') > 0, 'Cannot find colortable: %s', fname_dest_colorLUT);
        
    % Get 2nd column of text file
    t = readtableSafe(fname_dest_colorLUT);

    dest_labels = t.short;
    dest_labels = dest_labels(~cellfun(@isempty, dest_labels));
    dest_labels = lower(dest_labels(:));
    
    % Add labels to our definitions file
    atlas_labels = struct();
    if exist(fname_atlas_labels, 'file')
        data = load(fname_atlas_labels);
        if isfield(data, 'atlas_labels')
            atlas_labels = data.atlas_labels;
        end
    end
    atlas_labels.destrieux = dest_labels;
    save(fname_atlas_labels,'atlas_labels');  %#ok<*STRNU>
end






