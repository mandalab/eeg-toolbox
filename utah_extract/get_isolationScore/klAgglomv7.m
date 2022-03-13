%% klAgglom takes an mxn matrix obs where rows are observations and
%  columns are data points. Note that in comments, the recursion is such
%  that I don't want to distinguish between individuals and clusters, so I
%  consider individuals as one-member-clusters.

function [clusters, clusMembers, distRaw, linkMat] = klAgglomv7(obs,varargin)

global catVars

% Set defaults
type        = 'corr';
combType    = 'all'; % use "mean" or "all"
k           = 6;
nMembers    = 10;
refine      = 1;
refType     = 'same';
kType       = 'correlation';
goodType    = 'sumsquares';
nPrint      = 100;
print       = 0;
catVars     = zeros(1,size(obs,2));

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'},
            type = varargin{varStrInd(iv)+1};
        case {'-k','k'},
            k = varargin{varStrInd(iv)+1};
        case {'-r','refine'},
            refine = varargin{varStrInd(iv)+1};
        case {'-p','print'},
            print = 1;
        case {'-c','comb'},
            combType = varargin{varStrInd(iv)+1};
        case {'-n'}
            nMembers = varargin{varStrInd(iv)+1};
        case {'catvars','catVars','-cv'}
            catVars = varargin{varStrInd(iv)+1};
        case {'rtype','reftype','refType'}
            refType = varargin{varStrInd(iv)+1};
        case {'ktype','kType'}
            kType   = varargin{varStrInd(iv)+1};
        case {'nprint'}
            nPrint = varargin{varStrInd(iv)+1};
    end
end

% Make sure catVars is logical...
catVars = logical(catVars);
% obsRaw  = obs;

% If computing prod, normalize obs
if ismember(type,{'prod'}),
    maxObs = nanmax(obs,[],1);
    minObs = nanmin(obs,[],1);
    normObs = (obs-repmat(minObs,size(obs,1),1))./repmat((maxObs-minObs),size(obs,1),1);
    normRaw = normObs;
%     obs = normObs;
end

% Figure out range of observations: Works for corr, prod
switch type
    case 'corr'
        simRange = [0, 2];
    case 'prod'
        simRange = [0, 1];
    otherwise
        simRange = [0, 1];
end


% Start by computing pair-wise distances
fprintf('Initializing pairwise distances...')
if ismember(type,{'corr'}),
    % Clean up obs, if "corr"
    goodCols = zeros(1,size(obs,2));
    for ic = 1:size(obs,2),
        goodCols(ic) = sum(~isnan(obs(:,ic))) == size(obs,1);
    end
    obs = obs(:,logical(goodCols));
    distMat = 1-corr(obs',obs');
    for ir = 1:size(obs,1),
        for ic = ir:size(obs,1),
            distMat(ic,ir) = nan;
        end
    end
else
    distMat = nan(size(obs,1),size(obs,1));
    for im = 1:size(obs,1),
        for in = (im+1):size(obs,1),
            if ismember(type,{'prod','exprod'}),
                distMat(im,in) = klGetDist(normObs(im,:),normObs(in,:),'-t',type,'cat',catVars);
            else
                distMat(im,in) = klGetDist(obs(im,:),obs(in,:),'-t',type,'cat',catVars);
            end
        end
    end
end
fprintf('Done!\n');

% Save distance and observations in raw form in case I manipulate them in
% the next step...
distRaw = distMat;
obsRaw  = obs;

% Initialize clusMembers
clusMembers = num2cell((1:size(distRaw,1))');
linkMat = [];
linkIDs(1:size(distRaw,1),1) = 1:size(distRaw,1);
nLoops = 0;
while size(distMat,1) > 1,
    if mod(size(distMat,1),nPrint) == 0 && print,
        fprintf('%d members remaining...\n',size(distMat,1));
    end
    nLoops = nLoops+1;
    % Get the row and column corresponding to the minimum distance
    % remaining. If there are multiple, take just the first index
    % because other indices should be captured in the next loop
    % through.

    [minRow,minCol] = find(distMat == min(distMat(:)),1);
    if isnan(min(distMat(:))), break; end
    nRow = length(clusMembers{minRow,end});
    nCol = length(clusMembers{minCol,end});

    % Combine observations (weighted mean)
    combObs = sum([obs(minRow,:).*nRow; obs(minCol,:).*nCol],1)./(nRow+nCol);

    % Get the indices of clusMembers to keep
    nonClust = 1:size(distMat,1); nonClust(ismember(nonClust,[minRow,minCol])) = [];
    newMems = clusMembers(nonClust,end);
    newMems{end+1} = [clusMembers{minRow,end},clusMembers{minCol,end}];
    theseMems = newMems{end};

    % Let's make a linkage matrix. Remember:
    % For example, suppose there are 30 initial nodes and at 
    % step 12 cluster 5 and cluster 7 are combined. Suppose their 
    % distance at that time is 1.5. Then Z(12,:) will be [5, 7, 1.5]. 
    % The newly formed cluster will have index 12 + 30 = 42. If 
    % cluster 42 appears in a later row, it means the cluster created 
    % at step 12 is being combined into some larger cluster.
    linkMat = cat(1,linkMat,[linkIDs(minRow),linkIDs(minCol),min(distMat(:))]);
    newLink = nLoops + size(distRaw,1);
    linkIDs = linkIDs(nonClust);
    linkIDs(end+1) = newLink;

    % Debugging to be sure we're not overcounting, should be worked
    % out...
    if any(ismember(clusMembers{minRow,end},clusMembers{minCol,end})),
        keyboard
    end

    % Save combined members in clusMembers
    clusMembers(1:length(newMems),size(clusMembers,2)+1) = newMems;

    % Put new cluster distances at the end of "obs" after removing the
    % members of that cluster
    obs = obs(nonClust,:);
    obs = cat(1,obs,combObs);

    % Cut the cluster members from distMat
    distMat = distMat(nonClust,nonClust);

    % Calculate distance between other clusters and the most recent
    % combination
    distMat = cat(1,distMat,nan(1,size(distMat,2)));
    distMat = cat(2,distMat,nan(size(distMat,1),1));
    for im = 1:(size(distMat,1)-1),
        switch combType
            case 'mean'
                thisDist = klGetDist(obs(im,:),obs(end,:),'-t',type,'cat',catVars);
            case 'all'
                tstMems = clusMembers{im,end};
                allCorrs = distRaw(tstMems,theseMems);
                allCorrsTrans = distRaw(theseMems,tstMems)';
                noNanCorrs = nan(size(allCorrs));
                noNanCorrs(~isnan(allCorrs)) = allCorrs(~isnan(allCorrs));
                noNanCorrs(~isnan(allCorrsTrans)) = allCorrsTrans(~isnan(allCorrsTrans));

                thisDist = nanmean(noNanCorrs(:));
        end
        distMat(end,im) = thisDist;
    end    
end

% Let's put the next part in a loop so we can ask for multiple "k" values
% without having to rerun the time consuming correlation part above.
clusters = nan(size(obsRaw,1),length(k));
icdInd = nan(1,length(k));
for ik = 1:length(k),
    % Go backwards and find the last instance with "k" clusters of "nMembers"
    % or more
    nClusts = 0;
    nLoops = 0;
    premExit = 0;

    while nClusts < k(ik)
        nLoops = nLoops + 1;
        clustLens = cellfun(@length,clusMembers(:,end-nLoops));
        nClusts = sum(clustLens >= nMembers);
        if nClusts == 0,
            fprintf('Criteria not reached before waveforms are too scattered. Exiting loop\n');
            premExit = 1;
            break;
        end
        
    end

    if premExit,
        clusters(:,ik) = nan(size(obsRaw,1),1);
        return;
    end

    % Get mean data points by cluster
    clustInds   = find(clustLens >= nMembers);
    clustersPre = nan(size(obsRaw,1),1);
    for ic = 1:length(clustInds),
        memInds = clusMembers{clustInds(ic),end-nLoops};
        clustMeans(ic,:) = nanmean(obsRaw(memInds,:),1);
        if exist('normRaw','var'), normMeans(ic,:)  = nanmean(normRaw(memInds,:),1); else normRaw = nan(size(obsRaw)); end
        clustersPre(memInds) = ic;
    end

    % If asked to, loop through all observations and compare to the mean data
    % points and reassign each observation based on the nearest distance
    if refine
        switch refType
            case 'kmeans'
                if ismember(type,{'prod','exprod'}),
                    [clusters(:,ik),~,~] = kmeans(normRaw,k(ik),'Start',clustMeans,'Distance',kType);
                else
                    [clusters(:,ik),~,~] = kmeans(obsRaw,k(ik),'Start',clustMeans,'Distance',kType);
                end
                keyboard
            otherwise
                refRestart = 1;
                nTries = 0;
                while refRestart == 1 && nTries < 10,
                    nTries = nTries + 1;
                    if nTries == 8,
                        keyboard
                    end
                    for io = 1:size(obsRaw,1),
                        testDist = nan(length(clustInds),1);
                        for ic = 1:length(clustInds),
                            if ismember(type,{'prod','exprod'}),
                                testDist(ic) = klGetDist(normRaw(io,:),normMeans(ic,:),'-t',type,'cat',catVars);
                            else
                                testDist(ic) = klGetDist(obsRaw(io,:),clustMeans(ic,:),'-t',type,'cat',catVars);
                            end
                        end
                        if sum(isnan(testDist)) == numel(testDist),
                            tempClust = nan;
                        else
                            % I haven't figured out how to deal with identical distances yet...
                            % put in a catch for debugging this when necessary
        %                     if sum(testDist == min(testDist)) > 1,
        %                         keyboard
        %                     end

                            % Here, put in a catch to make sure that the distance
                            % isn't identical.
                            tempClust = find(testDist == min(testDist),1);
                            breakIO = 0; reRefine = 0;
                            if ismember(type,{'prod','exprod'}),
                                if any(abs(normRaw(io,catVars) - normMeans(tempClust,catVars)) > .0001),
                                    diffVar = find(abs(normRaw(io,catVars) - normMeans(tempClust,catVars)) > .0001,1);
                                    % Combine closest two mean clusters, then start
                                    % over.
                                    % Get the distances
                                    for ic = 1:length(clustInds),
                                        clustDists(ic,ic) = nan;
                                        for iic = (ic+1):length(clustInds)
                                            clustDists(ic,iic) = klGetDist(normMeans(iic,:),normMeans(ic,:),'-t',type,'cat',catVars);
                                            clustDists(iic,ic) = clustDists(ic,iic);
                                        end
                                    end
                                    % Find the closest
                                    [closeRow,closeCol] = find(clustDists == min(clustDists(:)),1);
                                    normMeans(closeRow,~catVars) = (normMeans(closeRow,~catVars).*sum(clustersPre == closeRow) + normMeans(closeCol,~catVars).*sum(clustersPre == closeCol))./sum(ismember(clustersPre,[closeCol,closeRow]));
                                    normMeans(closeCol,~catVars) = nanmean(normRaw(normRaw(:,catVars) == normRaw(io,catVars),~catVars),1);
                                    catVarInd = find(catVars);
                                    uVar = unique(normRaw(abs(normRaw(:,catVarInd(diffVar)) - normRaw(io,catVarInd(diffVar))) < .0001,1));
                                    normMeans(closeCol,catVarInd(diffVar)) = uVar(1);
                                    reRefine = 1;

                                end
                            end
                        end
                        if reRefine, 
                            break; 
                        end
                        clusters(io,ik) = tempClust;
                        clear testDist;
                        
                    end
                    if reRefine
                        refRestart = 1;
                    else
                        refRestart = 0;
                    end
                    
                end
        end
    else
        clusters(:,ik) = clustersPre;
    end

end

end