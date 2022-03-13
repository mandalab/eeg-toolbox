function [SNR,isoScore,fnScore,fpScore,outClustInfo] = klUnitIsolation_AIJ(spkWaves,nzWaves,quickQuality)
%- Inputs: spkWaves, 
%          nzWaves, 
%          quickQuality [optional]:  0=10,000 snippets max per calculation;  1=1,000 snippets.  Use quickQuality for time series
%
%-
if nargin==2, quickQuality = 0; end

% keyboard
if ~isempty(spkWaves) & ~isempty(nzWaves)
    noiseC = 5;
    kPerc = .05;
    faCut = .5;
    distType = 'euc';
    % distType = 'corr';
    lamb = 10;
    
    maxWaves = 10000;    %- paper uses 10000
    if quickQuality,
        maxWaves = 2000; %- 2000, faster but less accurate.  OK for time series
    end
    
    silent = 1;
    
    % Check if anyi nzWaves cross bottomThresh
    % nzCross = abs(nzWaves(:,alTimes == 0)) > abs(bottomThresh);
    nzTemp = nzWaves;%(nzCross,:);
    
    % Now we have Sclust and noise clust (spkClust, nzClust)
    Savg = nanmean(spkWaves,1);
    p2p  = max(Savg)-min(Savg);
    
    %% Start by getting SNR
    % Get noise type 1
    residK = spkWaves - repmat(Savg,size(spkWaves,1),1);
    noiseSpk = nanstd(residK(:));
    
    % Noise type 2 seems unavailable from our current data...
    
    % Get SNR:
    SNR = p2p./(noiseSpk*noiseC);
    
    %% Now let's get the isolation score
    % Subsample
%     maxSpk = ceil(maxWaves*(size(spkWaves,1)./(size(spkWaves,1)+size(nzTemp,1))));
%     maxNz = maxWaves-maxSpk;
%     spkRand = randperm(size(spkWaves,1));
%     spkInds = spkRand(1:min([maxSpk,size(spkWaves,1)]));
%     spkClust = spkWaves(spkInds,:);
%     nzRand = randperm(size(nzTemp,1));
%     nzInds = nzRand(1:min([maxNz,size(nzTemp,1)]));
%     nzClust  = nzTemp(nzInds,:);
    
    spkClust = spkWaves;
    nzClust = nzWaves;
    
%     outClustInfo.spkInds = spkInds;
%     outClustInfo.nzInds = nzInds;
    outClustInfo.spkClust = spkClust;
    outClustInfo.nzClust = nzClust;
    
    % Concatenate all events
    allEvents = [spkClust;nzClust];
    
    % Get d0 = average Euclidean distance in spike cluster
    if ~silent,
        fprintf('Initializing similarity matrix\n'); distTic = tic;
    end
    if strcmp(distType,'corr'),
        for i = 1:size(allEvents,2),
            goodCols(i) = sum(~isnan(allEvents(:,i))) == size(allEvents,1);
        end
        
        distMat = 1-corr(allEvents(:,goodCols)',allEvents(:,goodCols)');
        for i = 1:size(distMat,1),
            distMat(i:end,i) = deal(nan);
        end
    else
        %     keyboard
        distMat = squareform(pdist(allEvents));
        for i = 1:size(distMat,1),
            distMat(i:end,i) = deal(nan);
        end
    end
    if ~silent, fprintf('Distmat created in %s\n',printTiming(distTic)); end
    
    spkDist = distMat(1:size(spkClust,1),1:size(spkClust,1));
    d0 = nanmean(spkDist(:));
    simMat = getSim(distMat,lamb,d0);
    
    if ~silent, fprintf('Getting all P(X)...'); pxTic = tic; end
    % Now get all P(X) where X is an event in the spike cluster
    PX = zeros(size(spkClust,1),1);
    for ix = 1:size(spkClust,1),
        % Get sum distance to other clusters
        % 	sumSim = nansum(simMat(ix,:));
        sumSim = nansum([simMat(1:ix,ix)',simMat(ix,(ix+1):end)]);
        
        % Now for all Y's, get PxY
        %     PX(ix) = nansum(simMat(ix,1:size(spkClust,1))./sumSim);
        PX(ix) = nansum([simMat(1:ix,ix)',simMat(ix,(ix+1):size(spkClust,1))]./sumSim);
    end
    if ~silent, fprintf('Done in %s\n',printTiming(pxTic)); end
    
    isoScore = nanmean(PX);
    
    %% Let's work on false negative scores
    
    % Get k nearest neighbors for each noise event
    k = ceil(kPerc.*size(spkClust,1));
    isNFA = nan(size(nzClust,1),1);
    for inz = 1:size(nzClust,1),
        simNz = [simMat(1:(inz+size(spkClust,1)),(inz+size(spkClust,1)))',simMat((inz+size(spkClust,1)),(inz+size(spkClust,1)+1):end)];
        [sortNz,sortInd] = sort(simNz,'descend');
        knnInds = sortInd(1:k);
        if any(isnan(simNz(knnInds))),
            tmpInds = nan(1,k);
            tmpInds(1:sum(~isnan(simNz(knnInds))))     = knnInds((sum(isnan(simNz(knnInds)))+1):end);
            tmpInds((sum(~isnan(simNz(knnInds)))+1):k) = sortInd((k+1):(k+sum(isnan(simNz(knnInds)))));
            knnInds = tmpInds;
        end
        nSimSpks(inz) = sum(ismember(knnInds,[1:size(spkClust,1)]));
        isNFA(inz) = nSimSpks(inz)/k >= faCut;
    end
    nfa = sum(isNFA);
    fnScore = nfa/(nfa+size(spkClust,1));
    
    %% Now false positive scores
    for ispk = 1:size(spkClust,1),
        simSpk = [simMat(1:ispk,ispk)',simMat(ispk,(ispk+1):end)];
        [sortSpk,sortInd] = sort(simSpk,'descend');
        knnInds = sortInd(1:k);
        if any(isnan(simSpk(knnInds))),
            tmpInds = nan(1,k);
            tmpInds(1:sum(~isnan(simSpk(knnInds)))) = knnInds((sum(isnan(simSpk(knnInds)))+1):end);
            tmpInds((sum(~isnan(simSpk(knnInds)))+1):k) = sortInd((k+1):(k+sum(isnan(simSpk(knnInds)))));
            knnInds = tmpInds;
        end
        fpK(ispk,:) = knnInds;
        nSimNz = sum(~ismember(knnInds,[1:size(spkClust,1)]));
        isNFP(ispk) = nSimNz/k >= faCut;
    end
    nfp = sum(isNFP);
    fpScore = nfp/size(spkClust,1);
else
    SNR = nan;
    isoScore = nan;
    fnScore = nan;
    fpScore = nan;
    outClustInfo = nan;
end



%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
thisSim = exp(-(dist.*lamb./dNorm));



