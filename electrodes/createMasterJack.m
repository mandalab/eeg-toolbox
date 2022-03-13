function jacktable = createMasterJack(subj, rootEEGdir, varargin)
% Creates jacksheetmaster.csv based on element_info and raws
%
% createMasterJack.m
%
%
% FUNCTION:
%    createMasterJack(subj,rootEEGdir)
%
% INPUT ARGs:
%   subj = 'TJ022'
%   rootEEGdir
%
% OPTIONAL KEY-VALUE PAIRS:
%   skipFileCreationAndRawCheck - If 1, don't write to file or cross-check raws. Default 0. Should only be used
%                                 by getJackTable.m
%
% Output
%   jacktable
%   writes to file (if not skipFileCreation)
%
% Notes
%   If raw_info.csv format is changed, also change createRawInfo
%
% REVISION HISTORY
%   08/16 MST Created
%   10/16 MST Tweak to give more user flexibility (dont auto add cut chans)
%   5/17 TCS auto-add DC channels
%   07/17 MST - Creates tal/leads_bp.txt
%   07/17 MST - Don't add UTAH/Microwires to jacksheet
%   04/18 MST - Fix the parsing of arrays in element_info to use eval
%   05/18 MST - Add skipFileCreation param
%
% See also: createRawInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JACKSHEET_FILEPATH  = fullfile(rootEEGdir, subj, 'docs', 'jacksheetMaster.csv');
RAWINFO_FILEPATH    = fullfile(rootEEGdir, subj, 'docs', 'raw_info.csv');
EL_INFO_FILEPATH    = fullfile(rootEEGdir, subj, 'docs', 'element_info.csv');
HW_EXCLUDE          = {'utah'};

ip = inputParser;
ip.addParameter('skipFileCreationAndRawCheck', 0);
ip.parse(varargin{:});
skipFileCreationAndRawCheck = ip.Results.skipFileCreationAndRawCheck;

%%-----------------------------------------------------------------------------------------------------------%%
%%- Create "jacktable" expected jackSheetMaster based on element_info
temp = dbstack();
callers = {temp.file};
if ~any(contains(callers, 'verifyElementInfo.m'))
    verifyElementInfo(subj, rootEEGdir);
end


elFilename    = fullfile(rootEEGdir, subj, 'docs', 'element_info.csv');
elInfo        = getElementInfo(subj, rootEEGdir);
key           = readtableSafe(which('element_info_key.csv'), 'readRowNames', true);
jacktable     = table();
numChans      = 0;
allRecNotUsed = {};
hw_exl_tags   = {};

% get channels from element_info, push to jacktable
for i = 1 : height(elInfo)
    hardwareSize = elInfo.bigDim(i) * elInfo.smallDim(i);
    cutStr     = elInfo.isCut{i};
    try
        dontUseStr = elInfo.recDontUse{i};
    catch
        dontUseStr = [];
    end
    hwType = elInfo.hardwareType{i};
    if ~ismember(hwType, HW_EXCLUDE)
        if ~isempty(cutStr)
            try
                cut = eval(cutStr);
            catch e
                fprintf('%s:%s\n', e.message, cutStr);
                system(sprintf('open %s', EL_INFO_FILEPATH));
                keyboard;
                rethrow(e);
            end
        end
        dontUse = [];
        %- dontUse electrodes were recorded but should not be extractd or split, e.g., because physio issues, or DC channels without info we use
        if ~isempty(dontUseStr)
            try
                dontUse = eval(dontUseStr);
            catch e
                fprintf('%s:%s\n', e.message, dontUseStr);
                system(sprintf('open %s', EL_INFO_FILEPATH));
                keyboard;
                rethrow(e);
            end
        end
        for hardwareNdx = 1 : hardwareSize
            
            if ismember(hwType, {'microwire'}) && isempty(dontUse)
                continue;
            end
            if ismember(hardwareNdx, [dontUse])
                if (strcmp('DC',elInfo.tagName{i}) == 1) %If EKG or R, don't use 'O#', use '#'
                    allRecNotUsed{end+1}=sprintf('%s%02d',elInfo.tagName{i},hardwareNdx); continue; else
                    allRecNotUsed{end+1}=sprintf('%s%01d',elInfo.tagName{i},hardwareNdx); continue;
                end
            end % dontUse should get processed before cut so they can be excluded below
            if ismember(hardwareNdx, [cut]), continue; end % skip cut
            %if strcmpi(elInfo.chanType, 'CLIN'), continue; end % skip clin *MST 07/17 keep CLIN channels (so EKG/REF will be extracted!)
            numChans  = numChans + 1;
            tagName   = strtrim(elInfo.tagName{i});
            row       = makeJacksheetRow(subj, elInfo, tagName, jacktable, hardwareNdx, key);
            jacktable = [jacktable; struct2table(row)]; %#ok<AGROW>
        end
    else
        hw_exl_tags = [hw_exl_tags elInfo.tagName(i)];
    end
end

if isempty(jacktable)
    fprintf('elementInfo: %s\n', EL_INFO_FILEPATH);
    error('Jacktable is empty; are you sure elementInfo.csv is filled out correctly?');
end

assert(numel(unique(jacktable.chanName)) == numel(jacktable.chanName), 'Channel names are not all unique');




%%-----------------------------------------------------------------------------------------------------------%%
%%- Now loop over Raw files and confirm that the saved channels match jackSheetMaster
% get channels from raws
skip_raws = 0;
rawDir = fullfile(rootEEGdir, subj, 'raw');
if exist('findRaws','file')
    raws = findRaws(rawDir, {});
    if isempty(raws)
        warning('No raws for %s. Exiting createMasterJack without checking raws', subj);
        skip_raws = 1;
    end
else
    warning('MATLAB function findRaws not found on path (check eeg_toolbox/align). Exiting createMasterJack without checking raws');
    skip_raws = 1;
end

if ~skipFileCreationAndRawCheck
    writetable(jacktable, JACKSHEET_FILEPATH);  % write the jackSheetMaster before calling createRawInfo
end

if ~skip_raws && ~skipFileCreationAndRawCheck
    fprintf('\nChecking element_info-derived jacksheet against raw files...');
    
    % generate raw_info file... lists every possible channel from each raw
    [rawInfo, rawChanNotInJacksheet, jackChanNotInRaws] = createRawInfo(subj, rootEEGdir);
    rawChanNotSpecified =  setdiff(rawChanNotInJacksheet, allRecNotUsed);     % remove rawChanNotInJack that were exclued using "recNotUsed" in element_info.csv
    excludedRawNotRec   =  setdiff(allRecNotUsed, rawChanNotInJacksheet);
    
    if ~isempty(rawChanNotSpecified) || ~isempty(jackChanNotInRaws) || ~isempty(excludedRawNotRec)
        fprintf(' \n');
        fprintf('\n ************** ERROR: CHANNELS FOUND IN RAWs NOT SPECIFIED IN ELEMENT_INFO.CSV **************'); %following channels were found in RAWs but not specified in element_info.csv:\n');
        fprintf('\n  ALL recorded channels (including REFs, EKG, DC, etc) must be accounted for in element_info.csv, or removed from .21Es');
        fprintf('\n   -- add all CLIN electrodes with correct range of bigDim and smallDim (e.g., EKG, bigDim=2, smallDim=1;  REF, bigDim=2, smallDim=1');
        fprintf('\n   -- add all SYNC or DC electrodes, exluding any below DC09 (e.g., DC, bigDim=12, smallDim=1, recNotUsed=[1 2 3 4])');
        fprintf('\n   -- if channel existed but should NOT be used or split, include that channel in ''recNotUsed'' column (e.g., DC01,DC02)');
        fprintf('\n   -- if channel was mistakenly recorded but doesnt exist, delete that string from all .21E files (e.g., ''G33'' --> '''')');
        fprintf('\n   -- if a tagName was mispelled, do NOT include, but QUIT and FIX the ChanNames manually NOW');
        fprintf('\n');
        fprintf('\t *Hint: If a tagName was mispelled, check/correct ALL raws for the mispelling.\n');
        fprintf('\t\t To open affected raw files, you can use this terminal search command:\n');
        fprintf('\t\t open `grep -irl "[mispelled_search_string]" %s/*/ChanNames.txt`\n', rawDir);
 
        %- tell the user which channels are not matching up
        if ~isempty(rawChanNotSpecified)
            fprintf('\n >>>> LIST OF OFFENDING CHANNELS (in RAWs, not element_Info):   ');         % found in raw, but not element.info (check docs/raw_info.csv to see which raws)
            fprintf('%s  ',rawChanNotSpecified{:});
        end
        if ~isempty(jackChanNotInRaws)
            fprintf('\n >>>> LIST OF OFFENDING CHANNELS (in element_Info, not RAWs):   ');         % found in raw, but not element.info (check docs/raw_info.csv to see which raws)
            fprintf('%s  ',jackChanNotInRaws{:});
        end
        if ~isempty(excludedRawNotRec)
            fprintf('\n >>>> LIST OF OFFENDING CHANNELS (excluded in element_info, but not found in RAWs):   ');   % listed as "recNotUsed", but not found in anyraw, but not element.info (check docs/raw_info.csv to see which raws)
            fprintf('%s  ',excludedRawNotRec{:});
        end
        
        system(sprintf('open %s', elFilename));
        reply = input('\n >>>> Manually UPDATE element_info.csv and/or ChanNames NOW and then PRESS [RETURN] TO RECREATE JACKSHEET AND RAW_INFO, \n    or ENTER ''M'' to interactively Modify all 21Es,  or ''q'' to break. <<<< \n','s');
        if ~isempty(reply) && (reply=='Q' || reply=='q')
            fprintf('\n keyboard here so you can see what up... ');
            keyboard;
        elseif reply=='M' | reply=='m' %#ok<OR2>
            eegModifyAll21Es(subj,rootEEGdir);
        end
        createMasterJack(subj, rootEEGdir); %- recursive call... force users to fix the OFFENDING channels to move forward to extractBipolarity!!!
        return; 
    else
        fprintf('\n SUCCESS!  Element_info accounts for all channels recorded in Raws!');
        fprintf('\n  --> docs/jacksheetMaster.csv and docs/rawInfo.csv are valid and have been created\n\n');
    end
    
    
    
    % recreate leads_bp.txt ONLY if raw's check out
    extractBipolarity(subj, rootEEGdir, 1);
    
end % end if ~skip_raws



end % createMasterJack2





% private functions
function [row, foundInElInfo] = makeJacksheetRow(subj, elInfo, tagName, jacktable, chanNumPart, key)
    % Makes a blank jacksheet table row by reading rows to include from key
    % and pulling values from element_info.csv
    if isnumeric(chanNumPart)
        chanNumPart = num2str(chanNumPart);
    end

    if ~istable(key)
        try
            keyFile = which('element_info_key.csv');
        catch ME
            fprintf('ERROR: no key file found. Should be on path and in: ');
            fprintf('eeg_toolbox/eeg_extract/element_info_key.csv\n');
            keyboard;
            error('Cannot make jacksheet row: No element_info_key');
        end
        key = readtableSafe(keyFile, 'readRowNames',true);
    end

    elRowNdx = find(strcmpi(elInfo.tagName, tagName), 1); % matching row
    foundInElInfo = ~isempty(elRowNdx);

    % Blank new row - standard columns
    row = [];
    row.subjId = {subj};

    format = '%s%s';    % tag# for everything else
    if strcmp(tagName,'DC'), format = '%s%02d'; chanNumPart = str2num(chanNumPart); end  %- JW added 8/17 now that DC's are getting split here
    row.chanName = {sprintf(format, tagName, chanNumPart)};

    % For each column marked for jacksheet inclusion in element_info_key.csv
    for j = 1 : width(key)
        incl = key{'inJacksheet',j}{:};

        if ~isempty(incl) && str2double(incl) > 0 
            % Marked for inclusion

            fieldName = key.Properties.VariableNames{j};

            % Initialize new fields as 1-element cells strings
            if ~isfield(row, fieldName)
                row.(fieldName) = {''};
            end

            if ~isempty(elRowNdx)

                % get row from element_info
                elRow = elInfo(elRowNdx, :);

                % Make non-strings into cells for consistence
                if ~iscell(elRow.(fieldName))
                    elCell = { elRow.(fieldName) };
                    %keyboard;
                else
                    elCell = elRow.(fieldName);
                end

                el = elCell{1};

                if strcmpi(key{'Format', j}{:}, 'array')
                     % Handle arrays (e.g. "[1 2:4]")
                    if isempty(el)
                        is_chan_included = 0;
                    else
                        % parse '[1 4 5]' --> 0 or 1 (membership boolean)
                        include_chans = eval(el);
                        is_chan_included = ismember(str2double(chanNumPart), include_chans);
                    end
                    row.(fieldName) = is_chan_included;
                else
                    % Handle normal values
                    row.(fieldName) = el;
                end
            end
        end
    end

    % strip numbers off utah's
    if strcmpi(row.hardwareType, 'utah')
        row.chanName = util_split_stringnum(row.chanName);
    end
    
    % struct2table requires that every field be a cell
    fns = fieldnames(row);
    for i=1:numel(fns)
        fn = fns{i};
        if ~iscell(row.(fn))
            row.(fn) = { row.(fn) };
        end
    end
    
    if ~isempty(jacktable) && istable(jacktable) && height(jacktable) > 0
        columns = jacktable.Properties.VariableNames;
        assert(numel(intersect(columns, fieldnames(row))) == numel(columns), 'Row does not match table');
        
    end

end