function [chanGroups fileSizeGB]= getChanGroupsToLoad(NSxFileDir,chans,compMemory)
%getChanGroupsToLoad
%Written by: Anthony I. Jang
%
%When grabbing ephys data from NSx files (utah arrays), the computer's
%memory may not be large enough to load all channels at once. This function
%gives the optimal number of channels to load at one time given the
%computer's memory capacity
%
% Input:
% NSxFileDir: directory of the NSx file (.ns5 or .ns6).
% chans: channel numbers to load
% compMemory: memory of the computer running the function (Gb)
%
% Output:
% chanGroups: cell array of channels numbers to be grouped


%giga         = 1000000000;
giga  = 2^30;
if ~exist(NSxFileDir,'file'),
    fprintf('\n ERROR: %s is not found and thus cant be split', NSxFileDir);
    keyboard;
end
fileProperty = dir(NSxFileDir);
fileSizeGB     = fileProperty.bytes/giga;


% Use only ~25% of total memory capacity to give some leeway
%memoryThresh = compMemory * 0.25; %- anthonys setting
memoryThresh = compMemory * 0.7;  %- JWs update... now we only ever convert a single channel to double, and LFPs are saved as int16 from the beginning.  So filesize*1.20 should be ample

%- turns out that the peak memory usage is identicle when pulling out all chanenls or a subset.  
%   This is because blackrock code actually pulls all of them, then trunkates after alloting all the memory, which is about 1.8x the size of the file
%   now JW is very efficient with not allocating more memory after that... only allocating about 25% more than file itself
%   question is whether that will cause a crash...

howManyBins = 1;
% If file size is big, load it in segments (divide total number of channels to appropriately sized bins)
if fileSizeGB > memoryThresh
    foundOptimalBins = false;
    while ~foundOptimalBins
        howManyBins = howManyBins+1;
        if fileSizeGB/howManyBins < memoryThresh, foundOptimalBins = true; end
    end
end


chanNumBoundary = 1:floor(length(chans)/howManyBins):length(chans);
chanGroups = {};
for b = 1:length(chanNumBoundary)
    if b < length(chanNumBoundary)
        chanGroups{b} = chans(chanNumBoundary(b):chanNumBoundary(b+1)-1);
    else
        chanGroups{b} = chans(chanNumBoundary(b):end);
    end
end



end

