function [r, gridLayout, missingElecs] = electrodesDotM(subj, rootEEGdir)
% ELECTRODESDOTM returns [r, gridLayout, missingElecs] from electrodes.m
%
% This function can be used to read electrodes.m for now, but electrodes.m
% should be retired since all info is in element_info.csv now

useElementInfo = false;

filename = fullfile(rootEEGdir, subj, 'docs', 'electrodes.m');
clear r gridLayout missingElecs;

if ~useElementInfo && ~exist(filename, 'file')
    %fprintf('Could not find filename %s. Using element_info.csv\n', filename);
    useElementInfo = true;
end

if useElementInfo
    % read from element_info.csv
    info = getElementInfo(subj, rootEEGdir);
    phys = strcmpi(info.chanType, 'PHYS');
    noUtah = ~strcmpi(info.hardwareType, 'UTAH');
    info = info(phys & noUtah,:);
    gridLayout = [info.bigDim, info.smallDim];
    missingElecs = cellfun(@eval, info.isCut, 'uniformOutput',false);

    % make r
    numTags = length(gridLayout);
    r = zeros(numTags, 2);
    chan = 0;
    for i = 1 : numTags
        fullSize = gridLayout(i, 1) * gridLayout(i, 2);
        realSize = fullSize - length(missingElecs{i});
        r(i, 1) = chan + 1;
        r(i, 2) = chan + realSize;
        chan = chan + realSize;
    end
    
else % use electrodes.m
    run(filename);

    for var = {'r','gridLayout','missingElecs'};
        if ~exist(var{1}, 'var')
            error('Could not read %s from %s', var{1}, filename);
        end
    end
end




end % electrodesDotM function