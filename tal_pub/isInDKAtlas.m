function [found, dk_labels] = isInDKAtlas(label, dk_labels)
    % isInDKAtlas(label) returns True if the given label is in the Desikan-Killiany atlas
    %
    % Description:
    %   We store a file in the toolbox called eeg_toolbox/tal_pub/atlas_labels.mat
    %   This contains a struct atlas_labels, with fields for each atlas. For example:
    %       .desikan_killiany.labels 
    %   
    %   This labels fields contain all labels (e.g. 35 for desikan-killiany, plus an "unknown")
    %
    %   Creation note:
    %       This file also contains the code that generates atlas_labels.mat from the freesurfer-distributed atlas.
    %       The DK atlas raw data from Freesufer is eeg_toolbox/visualize/data/colortable_desikan_killiany.txt
    %       atlas_labels.csv's DK column was genereated from this
    %
    %
    % More Info:
    %   https://projects.ninds.nih.gov:8090/display/Home/For+Users
    %   https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
    %
    % Revision History
    %   03/18 MST   - Created
    %
    % See Also: localizer_atlas, atlas_simplify, atlasFSTranslate, aparc2lobe, isInDestAtlas
   
    FORCE_RECREATE_DK_LABELS = 0;
    
    if FORCE_RECREATE_DK_LABELS
        create_dk_label();
    end
    
    if nargin < 1, label = ''; end
    if nargin < 2, dk_labels = []; end
    
    if ~isempty(dk_labels)
        % You already have the labels!
        % This parameter, which is returned and accepted, is for the calling
        % function to cache the value so as to avoid repeated file I/O
        found = ismember(lower(label), lower(dk_labels));
        return
    end
    
    
    toolbox_dir = [];
    
    % Get definition file
    fname_atlas_labels = which('atlas_labels.mat');
    if isempty(fname_atlas_labels)
        toolbox_dir = getToolboxDir();
        fname_atlas_labels = fullfile(toolbox_dir, 'tal_pub/atlas_labels.mat');
    end
    
    % Check the file for the DK atlas
    good_file = 0;
    if exist(fname_atlas_labels, 'file')
        data = load(fname_atlas_labels);
        if isfield(data, 'atlas_labels')
            atlas_labels = data.atlas_labels;
            if isfield(atlas_labels, 'desikan_killiany') && ~isempty(atlas_labels.desikan_killiany)
                good_file = 1;
            end
        end
    end
    
    if ~good_file
        atlas_labels = create_dk_label(toolbox_dir);
    end
    
    % do the lookup
    found = ismember(label, lower(atlas_labels.desikan_killiany));
    dk_labels = atlas_labels.desikan_killiany;

end


function atlas_labels = create_dk_label(toolbox_dir)
    % freesurfer colortable text file to our DK atlas file
    if nargin < 1, toolbox_dir = []; end
    
    fname_dk_colorLUT = which('colortable_desikan_killiany.txt');
    fname_atlas_labels  = which('atlas_labels.mat');
    if isempty(fname_dk_colorLUT) || isempty(fname_atlas_labels)
        if isempty(toolbox_dir)
            toolbox_dir = getToolboxDir();
        end
        fname_dk_colorLUT = fullfile(toolbox_dir, 'visualize/data/colortable_desikan_killiany.txt');
        fname_atlas_labels  = fullfile(toolbox_dir, 'tal_pub/atlas_labels.mat');
    end
    
    assert(exist(fname_dk_colorLUT, 'file') > 0, 'Cannot find colortable: %s', fname_dk_colorLUT);
        
    % Get 2nd column of text file
    cmd = sprintf('awk ''{print $2}'' %s', fname_dk_colorLUT);
    [status,stdout] = unix(cmd);
    assert(status==0, 'Unix read command failed: %s', cmd);

    dk_labels = strsplit(stdout);
    dk_labels = dk_labels(~cellfun(@isempty, dk_labels));
    dk_labels = lower(dk_labels(:));
    
    % Add labels to our definitions file
    atlas_labels = struct();
    if exist(fname_atlas_labels, 'file')
        data = load(fname_atlas_labels);
        if isfield(data, 'atlas_labels')
            atlas_labels = data.atlas_labels;
        end
    end
    atlas_labels.desikan_killiany = dk_labels;
    save(fname_atlas_labels,'atlas_labels');  %#ok<*STRNU>
end






