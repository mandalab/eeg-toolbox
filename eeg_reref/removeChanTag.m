function fileStem = removeChanTag(tag)
% function fileStems = removeChanTag(tag)
%   Description: removes the channel information and returns only filestems
%   used in conjunction with rerefWrapper
%
%   Input:
%          --tag: NIH011_170113_1602.061
%
%   Output:
%          -- fileStem: NIH011_170113_1602
if iscellstr(tag), tag = char(tag); end
idx = strfind(tag,'.');
fileStem = tag(1:idx-1);
