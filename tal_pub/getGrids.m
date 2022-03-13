function grids = getGrids(subj, rootEEGdir)
% getGrids returns channel names of grid elements
%
% Outputs a cell array of channel name strings 
%
% Author: MST
    info = getElementInfo(subj, rootEEGdir);
    row_mask = info.smallDim > 1;
    grids = info.tagName(row_mask);
%     if isscalar(grids) && iscellstr(grids)
%         grids = grids{1};
%     end
end