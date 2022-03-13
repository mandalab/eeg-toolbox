function [mat1_subsampled,mat2_subsampled] = subsampleWaves(mat1,mat2,maxNum)
%Get subsamples of rows of two matrices so their count is maxNum
%   
% keyboard

for k = 1:2
    if k == 1, thisRawMat = mat1;
    else thisRawMat = mat2;
    end

    if size(thisRawMat,1) > maxNum
        i_rand = randperm(size(thisRawMat,1));
        subsampled = thisRawMat(i_rand(1:maxNum),:);
    else
        % Keep the number of spike waves
        subsampled = thisRawMat;
    end
    
    if k == 1, mat1_subsampled = subsampled;
    else mat2_subsampled = subsampled;
    end
end



end

