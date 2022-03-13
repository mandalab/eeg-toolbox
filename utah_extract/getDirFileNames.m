function folders = getDirFileNames(path,includeFolders,includeFiles)
% function getDirFileNames
%   input:  path
%       optional: includeFolders (1; default True)
%       optional: includeFiles   (1; default True)

if nargin<2,
    includeFolders = 1;
    includeFiles   = 1;
end

folders = dir(path);
if includeFolders==0,  folders = folders(~[folders.isdir]); end
if includeFiles==0,    folders = folders([folders.isdir]); end
    
    

for k = length(folders):-1:1

    % remove folders starting with .
    fname = folders(k).name;
    if fname(1) == '.'
        folders(k) = [ ];
    end
    
    
end
% Just get the names
folders = {folders.name};