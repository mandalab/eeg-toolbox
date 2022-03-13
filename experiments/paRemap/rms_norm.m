function [r,x] = rms_norm(x)
    %takes the root mean square of a vector or a matrix along the first
    %dimension
    meanx = mean(x,1);
    x = x - repmat(meanx,size(x,1),1);
    r = sqrt(mean(x.^2,1));
end