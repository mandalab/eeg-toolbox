function table = getElementInfo(subj, rootEEGdir, pathToElmtInfo)
% GETELEMENTINFO returns a patient's element_info.csv table
%
% Looks at element_info_key.csv to know the column types of element_info.
% Also addresses a common readtable error
%
% Required Inputs
%   subj - e.g. NIH042
%   rootEEGdir - e.g. /Volumes/Shares/FRNU/data/eeg
%
% Optional Inputs
%   pathToElmtInfo - If this is passed, fn will use this full path instead of
%                   the default (which is [rootEEGdir]/[subj]/docs/)
%
% Revision History
%   03/2018 MST - Fix parsing error (e.g. commas in note) by calling readtableSafe with format string

KEY_FILEPATH_OVERRIDE = '';

% get key path
if ~isempty(KEY_FILEPATH_OVERRIDE)
    filename = KEY_FILEPATH_OVERRIDE;
else
    filename = which('element_info_key.csv');
end

% try to open key
if exist(filename, 'file')
   key = readtableSafe(filename, 'readRowNames',true);
   % get format from key
   isNumericCol = strcmpi(key{'Format', :}, 'numeric');
else
    key = [];
    warning('\nERROR: element_info_key.csv not found here %s.\n', filename);
    keyboard;
end


% open element_info
if exist('pathToElmtInfo', 'var')
    filename = pathToElmtInfo;
    clear('pathToElmtInfo');
else
    filename = fullfile(rootEEGdir, subj, 'docs', 'element_info.csv');
end

% Check that the file can be opened
% Developer Note: If it CAN be opened but readtableSafe throws an error, that is a table-parsing error and its 
%       handling should be done in readtableSafe.
%       If readtableSafe doesn't throw an error, but you want to check the table's format, do that check
%       after the table is read.
%
%       * Do NOT try to read the table without a format string *
%       If there is an error with the format string, update element_info_key.csv to address the issue

[fd, msg] = fopen(filename);
if fd <= 0
    fprintf('\n ERROR: element_info.csv could not be opened; you probably have it opened in an application... save and close and try again\n');
    fprintf('Message: %s\n', msg);
    reply = input('\n >>>> CLOSE ELEMENT_INFO.CSV and then press [RETURN] to proceed, or ''q'' to break. <<<< \n','s');
    if ~isempty(reply) && (reply=='Q' || reply=='q')
        fprintf('\n keyboard here so you can see what up... ');
        keyboard;
    end
    table = getElementInfo(subj, rootEEGdir); %- recursive call... force users to fix the OFFENDING channels to move forward to extractBipolarity!!!
    return
else
    fclose(fd);
end

% File can be opened, next try to open it as a table
if isempty(key)
    table = readtableSafe(filename);
else
    %- key found, so read element_info again with formats (assuming dimensions match up)
    
    % use formats from key to properly read this table
    formats = cell(width(key), 1);
    formats(isNumericCol) = {'%d'};
    formats(~isNumericCol) = {'%s'};
        
    try
        table = readtableSafe(filename, 'format',[formats{:}]);
    catch e_readtableSafe
        
        % Any error not handled by readtableSafe is serious and unexpected
        fprintf('Error:getElementInfo::readtableSafe: %s\n', e_readtableSafe.message);
        
        try
            % it's possible that the table can be read without the key's format
            % if so, the table will not parse correctly, but this code will help the user
            % update the element_info_key.csv file
            
            table = readtable(filename);
            
            % list of columns from key and this table.. do they match?
            keyCol  = key.Properties.VariableNames;
            thisCol = table.Properties.VariableNames;
            notMatching = 0;

            %- while we are looking for matches, make a nice output to help user find any inconsistencies
            keyColOut = {'**KEY.CSV**    --> **ELEMENT_INFO.CSV**'};
            for ii=1:max([length(keyCol) length(thisCol)])
                if ii<=length(keyCol)
                    keyStr = keyCol{ii};
                else
                    keyStr = '<empty>';
                    notMatching = 1;
                end
                if ii<=length(thisCol),
                    if ii<=length(keyCol) && strcmp(keyCol(ii),thisCol(ii)),
                        eleStr = '<match>';
                    else
                        eleStr = thisCol{ii};
                        eleStr(end+1:12) = ' ';
                        eleStr = sprintf('%s  << must fix!!!', eleStr);
                        notMatching = 1;
                    end
                else
                    eleStr = '<empty>';
                    notMatching = 1;
                end
                keyStr(end+1:12) = ' ';
                keyColOut{end+1} = sprintf('%s   -->   %s', keyStr, eleStr);
            end



            %- this condition should be true unless element_info_key has been updated with new columns
            if notMatching==0
                table = getElementInfo(subj, rootEEGdir);
                return;

            else
                fprintf('\n\n ERROR: element_info.csv not matching format of element_info_key.csv:\n ');
                disp(keyColOut')
                reply = input('\n >>>> UPDATE element_info.csv NOW, **CLOSE IT**, and then press [RETURN] to proceed, or ''q'' to break. <<<< \n','s');
                if ~isempty(reply) && (reply=='Q' | reply=='q'),  
                    fprintf('\n keyboard here so you can see what up... ');
                    keyboard;    
                end
                table = getElementInfo(subj, rootEEGdir); %- recursive call... force users to fix the OFFENDING channels to move forward to extractBipolarity!!!
                return
            end
            
        catch e_readtable
 
            fprintf('\n ERROR::getElementInfo::readtable: Probably element_info.csv is open now or not saved... save and close and try again\n');
            fprintf('       readtableSafe error: %s\n', e_readtableSafe.message);
            fprintf('subsequent readtable error: %s\n', e_readtable.message);

            reply = input('\n >>>> CLOSE ELEMENT_INFO.CSV and then press [RETURN] to proceed, or ''q'' to break. <<<< \n','s');
            if ~isempty(reply) && (reply=='Q' || reply=='q')
                fprintf('\n keyboard here so you can see what up... ');
                keyboard;
            end
            table = getElementInfo(subj, rootEEGdir); %- recursive call... force users to fix the OFFENDING channels to move forward to extractBipolarity!!!
            return
        end
    end
end % isempty(key)
    
end % getElementInfo function
