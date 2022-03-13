function [SNR,isoScore,fnScore,fpScore,outClustInfo] = klUnitIsolationv2(waves,varargin)

upSamp = 0;
times = 1:size(waves,2);
noiseC = 5;
kPerc = .05;
faCut = .5;
sortType = 'classify';
distType = 'euc';
lamb = 10;
maxWaves = 10000;
silent = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','times'},
            times = varargin{varStrInd(iv)+1};
        case {'-c','c'},
            noiseC = varargin{varStrInd(iv)+1};
        case {'-u','upsamp'},
            upSamp = varargin{varStrInd(iv)+1};
        case {'-k','kperc'}
            kPerc = varargin{varStrInd(iv)+1};
        case {'type'},
            sortType = varargin{varStrInd(iv)+1};
        case {'-d'},
            distType = varargin{varStrInd(iv)+1};
        case {'-l'},
            lamb = varargin{varStrInd(iv)+1};
                
    end
end

normWaves = waves-repmat(nanmean(waves,2),1,size(waves,2));
if upSamp,
    warning off
    if size(waves,1) <= maxWaves,
        diffWaves = [nan(size(normWaves,1),5),diff(normWaves(:,5:20),[],2),nan(size(normWaves,1),length(21:32))];
        cutWaves=normWaves;
        cutWaves(abs(diffWaves) < 2 & abs(normWaves) > nanstd(nanmean(normWaves,1))*2) = nan;
        for ii = 1:size(waves,1),
            smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
        end
        smoothTimes = spline(1:32,times,1:.1:32);
        alignWind = 10;
    else
        if size(waves,1) <= (3*maxWaves),
            [isAP,isNoise, apMean, apTimes] = klSortSpikesv2(normWaves,'-u',0,'-a',1,'-t',times,'type',sortType);
            spkWaves = normWaves(isAP,:); 
            nzTemp   = normWaves(~isAP,:);
        else
            [sortThresh,sortInd] = sort(abs(nanmean(normWaves(:,times > -50 & times < 50),2)));
            nzTemp = normWaves(sortInd(1:floor(size(normWaves,1)/2)),:);
            spkWaves = normWaves(sortInd((floor(size(normWaves,1)/2)+1):end),:);
        end
        maxSpk = ceil(maxWaves*(size(spkWaves,1)./(size(spkWaves,1)+size(nzTemp,1))));
        maxNz = maxWaves-maxSpk;
        
        spkRand = randperm(size(spkWaves,1));
        spkInds = spkRand(1:min([maxSpk,size(spkWaves,1)]));
        spkClust = spkWaves(spkInds,:);
        nzRand = randperm(size(nzTemp,1));
        nzInds = nzRand(1:min([maxNz,size(nzTemp,1)]));
        nzClust  = nzTemp(nzInds,:);

        subWaves = [spkClust;nzClust];
        diffWaves = [nan(size(subWaves,1),5),diff(subWaves(:,5:20),[],2),nan(size(subWaves,1),length(21:32))];
        cutWaves=subWaves;
        cutWaves(abs(diffWaves) < 2 & abs(subWaves) > nanstd(nanmean(subWaves,1))*2.5) = nan;
        smoothWaves = nan(size(subWaves,1),length(1:.1:32));
        for ii = 1:size(subWaves,1),
            smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
        end
        smoothTimes = spline(1:32,times,1:.1:32);
        alignWind = 10;
        
        spkClust = smoothWaves(1:size(spkClust,1),:);
        nzClust = smoothWaves((size(spkClust,1)+1):end,:);
%         clear spkWaves nzTemp maxSpk maxNz
        
    end
    warning on
else
    smoothWaves = normWaves;
    smoothTimes = times;
    alignWind = 2;
end

% Align these waves on their troughs
[alWaves, alTimes] = klTroughAlignv3(smoothWaves,smoothTimes,0,'-w',alignWind);

inputWaves = alWaves;
if strcmp(distType,'corr'),
    for ic = 1:size(inputWaves,2),
        if sum(~isnan(inputWaves(:,ic))) ~= size(inputWaves,1);
            inputWaves(:,ic) = deal(nan);
        end
    end
end
 

% Sort units
[isAP,isNoise, apMean, apTimes] = klSortSpikesv2(inputWaves,'-u',0,'-a',0,'-t',alTimes,'type',sortType);


spkWaves = alWaves(isAP,:);
nzWaves  = alWaves(~isAP,:);

% Get bottom 2% of spike events
sortSpkZero = sort(abs(spkWaves(:,alTimes == 0)));
isNeg = 1;
if spkWaves(abs(spkWaves(:,alTimes == 0)) == max(sortSpkZero),alTimes == 0) >= 0,
    isNeg = 0;
end
bottom2 = nanmean(sortSpkZero(1:ceil((.02*size(sortSpkZero,1)))));

if isNeg,
    bottomThresh = -bottom2/2;
else
    bottomThresh = bottom2/2;
end

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

if ~upSamp,
    maxSpk = ceil(maxWaves*(size(spkWaves,1)./(size(spkWaves,1)+size(nzTemp,1))));
    maxNz = maxWaves-maxSpk;
    spkRand = randperm(size(spkWaves,1));
    spkInds = spkRand(1:min([maxSpk,size(spkWaves,1)]));
    spkClust = spkWaves(spkInds,:);
    nzRand = randperm(size(nzTemp,1));
    nzInds = nzRand(1:min([maxNz,size(nzTemp,1)]));
    nzClust  = nzTemp(nzInds,:);
else
    if size(waves,1) <= maxWaves,
        allInds = 1:size(waves,1);
    else
        allInds = [spkInds,nzInds];
    end
    spkInds = allInds(isAP);
    nzInds = allInds(~isAP);
    spkClust = spkWaves;
    nzClust = nzWaves;
end
  
outClustInfo.spkInds = spkInds;
outClustInfo.nzInds = nzInds;
outClustInfo.spkClust = spkClust;
outClustInfo.nzClust = nzClust;
outClustInfo.times = alTimes;

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
    distMat = nan(size(allEvents,1));
    for is = 1:size(allEvents,1),
%         for iis = (is):size(allEvents,1),
%             if iis ~= is,
%                 distMat(is,iis) = klGetDist(allEvents(is,:),allEvents(iis,:),'-t',distType);
%             end
%         end
        thisMult = allEvents-repmat(allEvents(is,:),size(allEvents,1),1);
        thisMult(isnan(thisMult)) = 0;
        distMat(:,is) = diag(thisMult*thisMult');
        
        
    end
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
        tmpInds(1:sum(~isnan(simNz(knnInds)))) = knnInds((sum(isnan(simNz(knnInds)))+1):end);
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

%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
    thisSim = exp(-(dist.*lamb./dNorm));
end


end
        