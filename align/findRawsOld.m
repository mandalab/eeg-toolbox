function eegRootList = findRaws(rawPath, eegRootListIn)
    % eegRootList = findRaws(rawPath, eegRootListIn)
    %
    %   Search all subdirectories and create list of raws
    %   *does not search _* directories like _grabbed
    %
    %   Inputs
    %    rawPath - root under which to search
    %    eegRootListIn - existing list of raws from previous searches (cell
    %    array)
    %
    %   Outputs
    %    eegRootList - list of filepaths of all raws (without extensions)

%- first check to see if raw file (instead of raw path)
if nargin < 2, eegRootListIn = {}; end

if ~isempty(strfind(rawPath,'.EEG')) || ~isempty(strfind(rawPath,'.21E')),
    rootEEG = rawPath(1:end-4);
    listEEG = dir([rootEEG '.EEG']);
    list21E = dir([rootEEG '.21E']);
    
    if length(listEEG)+length(list21E)==0, 
        % do nothing
    elseif length(list21E)~=1 || length(listEEG)~=1,
        fprintf('>>> %s missing matched .EEG and .21E files\n', rawPath);
    else
        eegRootListIn{end+1} = fullfile(rootEEG);
        %fprintf('> found %s .EEG and .21E\n', rootEEG)
    end
    %- its not an .EEG or .21E file (should be a directory).  recursively search subdirectories
else
    % gets filepath for .EEG file
    listEEG = dir( fullfile(rawPath, '*.EEG') );
    
    for iList=1:length(listEEG)
        rootEEG = listEEG(iList).name(1:end-4);
        list21E = dir( fullfile(rawPath, [rootEEG '.21E']) );
        
        if length(list21E)~=1
            fprintf('>>>>>>>\n> in %s MISSING MATCHED %s .21E .EEG files!!\n>>>>>>>\n', rawPath, listEEG(iList).name);
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
end
eegRootList = unique(eegRootListIn);