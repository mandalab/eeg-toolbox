function lobe = aparc2lobe(aparc)
% lobe = aparc2lobe(aparc) Maps a Desikan-Killiany atlas location to a lobe
% Overview: https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation 
% Labels  : $FREESURFER_HOME/average/colortable_deskan_killiany.txt

    if contains(aparc, '-')
        % strip the prefix: "ctx-lh-" for example
        parts = strsplit(aparc, '-');
        if numel(parts) >= 3
            aparc = strjoin(parts(3:end),'-');
        end
    end

    FRONTAL = {     'superiorfrontal'
                    'rostralmiddlefrontal'
                    'parsopercularis'
                    'parstriangularis'
                    'parsorbitalis'
                    'caudalmiddlefrontal'
                    'lateralorbitofrontal'
                    'medialorbitofrontal'
                    'precentral'
                    'paracentral'
                    'frontalpole'
                    'rostralanteriorcingulate'
                    'caudalanteriorcingulate'
                };     
    PARIETAL = {    'superiorparietal'
                    'inferiorparietal'
                    'supramarginal'
                    'postcentral'
                    'precuneus'
                    'isthmuscingulate'
                    'posteriorcingulate'
                };   
    TEMPORAL = {    'superiortemporal'
                    'bankssts'
                    'middletemporal'
                    'fusiform'
                    'transversetemporal'
                    'temporalpole'
                    'inferiortemporal'
                    'entorhinal'
                    'parahippocampal'
                };
    OCCIPITAL = {   'lateraloccipital'
                    'cuneus'
                    'lingual'
                    'pericalcarine'
                };        
    INSULA    = { 'insula' };

    aparc = lower(aparc);

    switch aparc
        case FRONTAL,   lobe = 'frontal';
        case PARIETAL,  lobe = 'parietal';
        case TEMPORAL,  lobe = 'temporal';
        case OCCIPITAL, lobe = 'occipital';
        case INSULA,    lobe = 'insula';
        otherwise,      lobe = 'unknown';
    end

end