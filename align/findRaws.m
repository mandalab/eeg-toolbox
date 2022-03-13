function eegRootList = findRaws(rawPath, eegRootListIn)
% eegRootList = findRaws(rawPath, eegRootListIn)
%
%   Search all subdirectories and create list of raws
%     does not search "_*" directories like _grabbed
%
%   Inputs
%    rawPath - root under which to search
%    eegRootListIn - existing list of raws from previous searches (cell
%    array)
%
%   Outputs
%    eegRootList - list of filepaths of all raws (without extensions)
%
%
% Revision History
%   12/17 melkalliny - handles either NK or BR in raw folders


%- first check to see if raw file (instead of raw path)
if nargin < 2, eegRootListIn = {}; end


%if nihon-kohden raw
% gets filepath for .EEG file
listEEGNK = dir( fullfile(rawPath, '*.EEG'));
listEEGBR = dir( fullfile(rawPath, '*.ns2') );
listEEGCV = dir( fullfile(rawPath, '*.TRC') );
listEEG = [listEEGBR;listEEGNK;listEEGCV];

% check if multiple BR files are present

if (size(listEEG,1) > 2)
    fprintf('\nWARNING: raws from multiple filesystems are present\n');
    keyboard;
end

for iList=1:length(listEEG)
    rootEEG = listEEG(iList).name(1:end-4);
    listHeaderNK = dir( fullfile(rawPath, [rootEEG '.21E']) );
    listHeaderBR_one = dir( fullfile(rawPath, [rootEEG '.ns3']) );
    listHeaderBR_two = dir( fullfile(rawPath, [rootEEG '.nev']) );
    listHeaderCV = dir( fullfile(rawPath, [rootEEG '.TRC']) );
    
    list21E = [size(listHeaderNK,1) size(listHeaderBR_one,1) size(listHeaderBR_two,1) size(listHeaderCV,1)];
    if sum(list21E) > 2;
        fprintf('\nWARNING: headers from multiple filesystems in %s\n',rootEEG);
        keyboard;
    end
    
    if sum(list21E)==0;
        fprintf('>>>>>>>\n> in %s MISSING MATCHED %s .header .signal files!!\n>>>>>>>\n', rawPath, listEEG(iList).name);
    else
        eegRootListIn{end+1} = fullfile(rawPath,rootEEG);
        %fprintf('> in %s found %s .EEG and .21E\n', rawPath, rootEEG)
    end
end

dir_list = lsCell(rawPath);
% exclude _ directories
dir_list = dir_list(~cellfun(@(s) s(1)=='_', dir_list));

for iList=1:length(dir_list),
    eegRootListIn = findRaws( fullfile(rawPath,dir_list{iList}), eegRootListIn );
end

eegRootList = unique(eegRootListIn);