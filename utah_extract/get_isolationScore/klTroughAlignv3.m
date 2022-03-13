function [outMat, alignedTimes] = klTroughAlignv3(inMat,times,tCrit,varargin)

wind = 10;
tInds = find(times == tCrit);
tInc  = nanmean(unique(diff(times)));
minVal = .15;
visualize = 0;
tMax = 200;
tMin = -50;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'wind','-w'},
            wind = varargin{varStrInd(iv)+1};
        case {'-v','vis'},
            visualize = varargin{varStrInd(iv)+1};
    end
end

% Loop through rows of inMat
for ir = 1:size(inMat,1),
    myWind = 1;
    nTries = 0;
    [xmax,imax,xmin,imin] = extrema(inMat(ir,:));
    imax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(ir,:)))) = [];
    xmax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(ir,:)))) = [];
    imin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(ir,:)))) = [];
    xmin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(ir,:)))) = [];
    
    diffInds = [imax,imin] - tInds;
    nearExtreme = diffInds(abs(diffInds) == min(abs(diffInds)))+tInds;
    if length(nearExtreme) > 1,
        myT(ir) = nearExtreme((nearExtreme - tInds) > 0);
    else
        myT(ir) = nearExtreme;
    end
    
    % Redo without the amplitude criterion if myT isn't within the correct range 
%     if inMat(ir,myT(ir)) > 0 && min(inMat(ir,:)) < -2,
%         keyboard
%     end
    if times(myT(ir)) < tMin || times(myT(ir)) > tMax,
        myT(ir) = tInds;
    end
    tempCell{ir,1} = inMat(ir,:);
    
    if visualize,
        figure()
        plot(times,inMat(ir,:));
        vline(times(myT(ir)));
        keyboard
    end
end

[outMat, tZero] = klAlignv2(tempCell,myT);
tMin = (tZero-1).*(-tInc); tMax = (size(outMat,2)-tZero).*tInc;
alignedTimes = tMin:tInc:tMax;    
    