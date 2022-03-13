function pairs = findNeighbors(subj, rootEEGdir, tagName, distance, varargin)
    % findNeighbors finds neighboring electrodes within a hardware piece
    %
    % Inputs:
    %   subj
    %   rootEEGdir
    %   tagName - Name of hardware to look at (ex: 'G')
    %   distance - number of units away (ex: 2 finds 2-neighbors *only*)
    %
    %
    % Optional Key-Value Input:
    %   useDiagonal, true/[false] - if true, distance is on diagonal
    %   element_info, info - can pass in element_info to avoid loading it (you
    %                       can skip subj/rootEEGdir params in this case)
    %
    % Output:
    %   pairs - n by 2 array of neighbor pairs. E.g. If a row is [1, 9], then G1 and G9 
    %           elements are neighbors
    %
    % Revision History
    %   7/17 MST - Updated
    %  10/17 MST - Use recDontUse column data
    %
    % See Also: createBipolarPairs

    %  Copyright (C) 2017 Mike Trotta

    % parse input
    ip = inputParser;
    ip.addParameter('useDiagonal', false);
    ip.addParameter('elementInfo', [], @istable);
    ip.addParameter('element_info', [], @istable);
    ip.parse(varargin{:});
    info = ip.Results.element_info;
    if ~isempty(ip.Results.elementInfo), info=ip.Results.elementInfo; end
    useDiagonal = ip.Results.useDiagonal;

    % var setup
    if isempty(info), info = getelement_info(subj, rootEEGdir); end
    row = info(strcmpi(info.tagName, tagName), :);
    assert(~isempty(row), 'tagName not found in element_info: %s\n', tagName);

    pairs = [];
    x_max = row.smallDim;
    y_max = row.bigDim;

    isCut = [];
    recDontUse = [];
    try
        temp = char(row.isCut);         if ~isempty(temp), isCut = eval(temp); end
        if ismember('recDontUse', row.Properties.VariableNames)
            temp = char(row.recDontUse);    if ~isempty(temp), recDontUse = eval(temp); end
        end
        
        missing = union(isCut, recDontUse);
        
    catch e
        fprintf('Error: element_info isCut/recDontUse column entries must be valid matlab code: %s\n', char(row.isCut));
        error(e.message);
    end

    % set dx/dy to neighbors
    if useDiagonal
        step1 = distance * [1,1];
        step2 = distance * [1,-1];
    else
        step1 = distance * [0,1];
        step2 = distance * [1,0];
    end

    % row/col to index
    getNum = @(x,y) (x-1) * row.bigDim + y;

    % for each electrode, find 2 neighbors and add their edge
    for x = 1 : row.smallDim
        for y = 1 : row.bigDim

            n = getNum(x,y);
            if ismember(n, missing), continue; end

            x1 = x + step1(1);
            y1 = y + step1(2);

            x2 = x + step2(1);
            y2 = y + step2(2);

            n1 = getNum(x1, y1);
            n2 = getNum(x2, y2);

            if ~ismember(n1, missing) &&  (1 <= x1 && x1 <= x_max) && (1 <= y1 && y1 <= y_max)
                pairs = [pairs; n, n1];
            end

            if ~ismember(n2, missing) &&  (1 <= x2 && x2 <= x_max) && (1 <= y2 && y2 <= y_max)
                pairs = [pairs; [n, n2]];
            end
        end
    end
end
