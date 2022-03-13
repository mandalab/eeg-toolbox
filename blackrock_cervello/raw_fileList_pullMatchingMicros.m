function [outputArg1] = raw_fileList_pullMatchingMicros(subj, ecogDir, microDir)

if ~exist(sprintf('%s/%s/behavioral/',ecogDir,subj))
    fprintf('\ninput subject directory doesnt exist.. aborting\n')
    return
end

temp_dir = dir(sprintf('%s/%s/raw/',ecogDir,subj));
rawFileFound = (find(contains({temp_dir.name},'matchWithMicro')));
if isempty(rawFileFound)
    fprintf('\nraw file list doesnt exist.. aborting\n')
    
    [resp] = (input(sprintf('\nShould we go ahead and copy over all the data in the micro folder? Y or N?\n '),'s'));
    if (resp=='Y') || (resp=='y')
        fprintf('\ncopying all of %s data from %s into %s\n',subj,microDir,ecogDir)
        [out] = copyAllMicroIntoSubjFolder(subj,ecogDir,microDir);
        fprintf('\nall done\n')
        return
    elseif (resp=='N') || (resp=='n')
        return
    end
end

if ~exist(sprintf('%s/%s/manual_sorts/',microDir,subj))
    fprintf('\nmicro subject directory doesnt exist.. aborting\n')
    return
end

rawFileList = readtable(fullfile(temp_dir(rawFileFound).folder,temp_dir(rawFileFound).name));

% for each of these files, lets go through and pull it into the subject
% folder
for file=1:size(rawFileList,1)
    grabThis = cell2mat(table2cell(rawFileList(file,{'containsBehavior'})));
    if grabThis==1
        fileToGrab = table2cell(rawFileList(file,{'folderName'}));
        grabFolder = sprintf('%s/%s/manual_sorts/data_processed/%s',microDir,subj,fileToGrab{1});
        if ~exist(grabFolder) %#ok<*EXIST>
            fprintf('\nneed a micro file that is not in the micro dir!\n-- %s --\n',fileToGrab{1})
            continue % this should be a keyboard
        end
        
        fprintf('\nbeginning to copy %s\n',fileToGrab{1})
        
        % now lets copy things in to the subject folder
        % but lets be sure not to overwrite anything that may have already
        % been accomplished in the subject folder
        grabFolderDir = dir(grabFolder);
        for i=1:length(grabFolderDir)
            
            if contains(grabFolderDir(i).name,'sortNotes')
                copyfile(fullfile(grabFolderDir(i).folder,grabFolderDir(i).name),sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,fileToGrab{1}))
            elseif contains(grabFolderDir(i).name,'lfp') || contains(grabFolderDir(i).name,'other') || contains(grabFolderDir(i).name,'sorted')
            %elseif contains(grabFolderDir(i).name,'sorted')

                if ~exist(sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,fileToGrab{1},grabFolderDir(i).name))
                    mkdir(sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,fileToGrab{1},grabFolderDir(i).name))
                end
                copyfile(fullfile(grabFolderDir(i).folder,grabFolderDir(i).name),sprintf('%s/%s/micro/manualsort/%s/%s/',ecogDir,subj,fileToGrab{1},grabFolderDir(i).name))
            end
            
        end
        
        fprintf('\ndone copying %s\n\n',fileToGrab{1})
        

    end
    
end

outputArg1 = [];

end









function [out] = copyAllMicroIntoSubjFolder(subj,ecogDir,microDir)

micro_dirs = dir(sprintf('%s/%s/manual_sorts/data_processed/',microDir,subj));
micro_dirs(find(strcmpi({micro_dirs.name},'.'))) = [];
micro_dirs(find(strcmpi({micro_dirs.name},'..'))) = [];
for j=1:length(micro_dirs);
    fprintf('\nbeginning to copy %s\n',micro_dirs(j).name)
    
    grabFolder = sprintf('%s/%s/manual_sorts/data_processed/%s/',microDir,subj,micro_dirs(j).name);
    % now lets copy things in to the subject folder
    % but lets be sure not to overwrite anything that may have already
    % been accomplished in the subject folder
    grabFolderDir = dir(grabFolder);
    for i=1:length(grabFolderDir)
        
        if contains(grabFolderDir(i).name,'sortNotes')
            copyfile(fullfile(grabFolderDir(i).folder,grabFolderDir(i).name),sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,micro_dirs(j).name))
        elseif contains(grabFolderDir(i).name,'lfp') || contains(grabFolderDir(i).name,'other') || contains(grabFolderDir(i).name,'sorted')
            if ~exist(sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,micro_dirs(j).name,grabFolderDir(i).name))
                mkdir(sprintf('%s/%s/micro/manualsort/%s/',ecogDir,subj,micro_dirs(j).name,grabFolderDir(i).name))
            end
            copyfile(fullfile(grabFolderDir(i).folder,grabFolderDir(i).name),sprintf('%s/%s/micro/manualsort/%s/%s/',ecogDir,subj,micro_dirs(j).name,grabFolderDir(i).name))
        end
    end
    fprintf('done copying %s\n\n',micro_dirs(j).name)
end
out = [];



end

