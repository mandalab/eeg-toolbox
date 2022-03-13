function verifyElementInfo(subj, rootEEGdir, varargin)
% VERIFYELEMENTINFO Checks element_info.csv to make sure it's setup correctly
%
%   DESCRIPTION
%    Examines the element_info_key.csv to verify correct setup
%
%    Current checks are:
%       column headers - warnings for unexpected/missing
%       Category check - checks that entry is member of accepted category
%                       given in key's Format row in [choiceA/B/C] type format
%       Synch check -   Check that at least one Synch channel is present
%       Empty columns - check that no columns are empty
%       unmarked _sh resected - if there are _sh channels, but only pre-_sh are marked resected, error.
%
%    Current corrections:
%       Left/right --> lh/rh.
%       Auto-correct column header capitalization.
%       For array columns (e.g. isCut), map '' --> '[]'.
%       If only _sh channels (implant-2) marked resected, mark nearby implant-1 channels resected.
%       Enforce isIcatal and isInterictal attributes match between _sh and pre-_sh channels
%
%
%   INTPUT
%    subj       - subject
%    rootEEGdir - eeg root folder
%
%   NAME-VALUE PAIR OPTIONAL INPUT
%    pathToElmtInfo, filepath - full filepath to element_info.csv instead of
%                   inferring from eeg/subj/docs/ path
%
%   OUTPUT
%    Returns error or alters element_info.csv to be correct
%
%   REVISION HISTORY:
%    07/2016 MST - Created
%    10/2016 MST - Add array column and auto add Sync checks
%    09/2017 JHW - modified a bit... removed nkRecOrder column and check, added check for non-zero bigDim/smallDim
%              - added some checks/safty for saving element_info.csv
%    04/2018 MST - Added an error to enforce no "," in array columns
%    05/2018 MST - Adding attribute (e.g. isIctal) checks/corrections for shifted (_sh) channels
%    03/2020 SNJ - Changed search between shift electrodes in ictal and interictal to isequal instead of comparing all elements, which didn't work for differing lengths
%                - Switched input reply == letter to strcmpi
%                - General performance improvement
%
%   TODO: verify that cut/icta/interictal fall within big*small range


% columns that require '[]' like array entries
arrayCols = {'isCut', 'recDontUse', 'isIctal', 'isInterictal', 'isResected'};


% INPUT PARSING
p = inputParser;
p.addParameter('pathToElmtInfo', '', @ischar);
parse(p, varargin{:});  %- JW removed 'silent' input... YOU CANT SILENCE ME, BITCH!


% VARIABLE SETUP
keyFilename = which('element_info_key.csv');
if ~isempty(p.Results.pathToElmtInfo)
    elFilename = fullfile(p.Results.pathToElmtInfo);
else
    elFilename  = fullfile(rootEEGdir, subj, 'docs', 'element_info.csv');
end

fileAltered = false;


% file checks
if ~exist(elFilename, 'file')
    error('ERROR: element_info.csv not found: %s', elFilename);
elseif ~exist(keyFilename, 'file')
    error('ERROR: key not found here: %s\n', keyFilename);
end

info = getElementInfo(subj, rootEEGdir);
key  = readtableSafe(keyFilename, 'ReadRowNames', true);
info_orig = info;

fname_jacksheet = fullfileEEG(rootEEGdir, subj, 'docs', 'jacksheetMaster.csv');

% empty columns
EMPTY_COLUMN_WARNING = 0;
if EMPTY_COLUMN_WARNING
    for i = 1 : length(info.Properties.VariableNames)
        colData = info{:, info.Properties.VariableNames{i}};
        fun = turnary(iscell(colData), @cellfun, @arrayfun);
        emptyMask = fun(@isempty, colData);
        if sum(emptyMask) == height(info)
            if inputYN(sprintf('Column %s is completely empty. Populate now? (see verifyElementInfo to silence this warning) ', info.Properties.VariableNames{i}))
                system(sprintf('open %s', elFilename));
                keyboard; % allow user to save
                info = readtableSafe(elFilename);
            end
        end
    end
end %- if EMPTY COLUMN WARNING


% Column Name checks
infoCols    = info.Properties.VariableNames;
keyCols     = key.Properties.VariableNames;
extraCols   = setdiff(infoCols, keyCols);
missingCols = setdiff(keyCols, infoCols);


% check for capitalization errors and correct them
capErrorColNames = intersect(upper(missingCols), upper(extraCols));
if ~isempty(capErrorColNames)
    % make corrections
    [~,infoNdx] = ismember(capErrorColNames, upper(infoCols));
    [~, keyNdx] = ismember(capErrorColNames, upper(keyCols));
    infoCols(infoNdx) = keyCols(keyNdx);
    info.Properties.VariableNames = infoCols;
    fileAltered = true;
    % redo missing columns checks
    infoCols = info.Properties.VariableNames;
    keyCols = key.Properties.VariableNames;
    extraCols = setdiff(infoCols, keyCols);
    missingCols = setdiff(keyCols, infoCols);
end


% Show column name mismatch errors
if ~isempty(missingCols) && ~contains(subj, 'DBS')
    fprintf(' WARNING: The following columns were missing in element_info:\n');
    for i = 1:length(missingCols)
        fprintf('\t%s: %s\n', missingCols{i}, key{'Format', missingCols{i}}{:});
    end
    fprintf('You should fix this now.\n');
    system(sprintf('open %s', elFilename));
    keyboard;
    verifyElementInfo(subj, rootEEGdir, varargin{:});
    return;
end
if ~isempty(extraCols) && ~contains(subj, 'DBS')
    fprintf(' WARNING: The following columns were unexpected in element_info:\n');
    for i = 1:length(extraCols)
        fprintf('\t%s\n', extraCols{i});
    end
    keyboard;
    return;
end


% map left --> lh and right --> rh
newCol = lower(info.whichHemi);
newCol = strrep(newCol, 'left', 'lh');
newCol = strrep(newCol, 'right','rh');
if ~isempty(find(~strcmp(newCol, info.whichHemi), 1))
    fileAltered = true;
    info.whichHemi = newCol;
end


%%%  Enforce non-numeric columns %%%%
% charNdxs = find(~strcmpi(formats, 'numeric'));
% for i = charNdxs
%    vals = info{:,i};
%    if isnumeric(vals) && sum(isnan(vals)) == length(vals)
%        colName = info.Properties.VariableNames{i};
%        newColVal = repmat({''}, length(vals), 1);
%        eval(sprintf('info.%s = newColVal;', colName));
%        fileAltered = true;
%    end
% end


% For columsn with choices in [A/B/c] format, make sure entries are one of the choices
formats = strtrim(key{'Format',:});
formatCols = find(cellfun(@(x) ~isempty(x) && x(1)=='[' && x(end)==']', formats));
for i = formatCols
    format = formats{i};
    format = strrep(format, '[', '');
    format = strrep(format, ']', '');
    choices = strsplit(format, '/');
    vals = info{:,i}';
    if ~isnumeric(vals) && ~contains(subj, 'DBS')
        badValNdx = find(~ismember(upper(vals), upper(choices)) & ~cellfun('isempty', vals));
        if ~isempty(badValNdx)
            fprintf('\n\n ********* ERROR: found values in column %s which have unexpected values. ***** \n', infoCols{i});
            fprintf('\tExpected values are:'); disp(choices);
            fprintf('\tFound values:'); disp(vals(badValNdx));
            fprintf('   Correct the error or edit format row in key file to add an expected value (%s).\n', keyFilename);
            
            system(sprintf('open %s', elFilename));
            reply = input('\n >>>> UPDATE element_info.csv NOW and then PRESS [RETURN] TO RE-CHECK IT, or ''q'' to break. <<<< \n','s');
            if ~isempty(reply) && strcmpi(reply,'Q') %SJ changed from (reply=='Q' || reply=='q')
                fprintf('\n keyboard here so you can see what up... ');
                keyboard
            end
            verifyElementInfo(subj, rootEEGdir);
            return
        end 
    end
end


% Synch check -- should be exactly 1 SYNC row
if sum(strcmpi(info.chanType,'SYNC'))~=1
    fprintf('\n\n ********* ERROR: There must be exactly 1 row labeled with chanType==SYNC in elemnent_info.csv');
    fprintf('  \n *********       but %d were found in %s *********\n', sum(strcmpi(info.chanType,'SYNC')), elFilename);
    if sum(strcmpi(info.chanType,'SYNC'))==0
        fprintf('\tTo fix this, identify the actual channel used for SYNC (probably DC, maybe EKG or R) \n');
        fprintf('\tand add it to element_info.csv (if not there) with chanType==SYNC\n');
    else   
        fprintf('\tTo fix this, make sure only the actual channel used for SYNC is labeled SYNC.\n');
        fprintf('\t(Probably that means making DC into CLIN, or EKG into CLIN.)\n');
    end
    fprintf('  Note: SYNC should always be the last row of element_info.csv');
    
    system(sprintf('open %s', elFilename));
    reply = input('\n >>>> UPDATE element_info.csv NOW and then PRESS [RETURN] TO RE-CHECK IT, or ''q'' to break. <<<< \n','s');
    if ~isempty(reply) && strcmpi(reply,'Q') %SJ changed from (reply=='Q' || reply=='q')
        fprintf('\n keyboard here so you can see what up... ');
        keyboard
    end
    verifyElementInfo(subj, rootEEGdir);
    return
end


% Force order of chanType:  PHYS-->CLIN->SYNC
iPHYS = find(strcmpi(info.chanType,'PHYS'));
iCLIN = find(strcmpi(info.chanType,'CLIN'));
iSYNC = find(strcmpi(info.chanType,'SYNC'));
if length([iPHYS; iCLIN; iSYNC])<size(info,1),
    fprintf('\n ERROR: chanType column must be PHYS, CLIN, or SYNC for each row (no blanks!)');
    system(sprintf('open %s', elFilename));
    reply = input('\n >>>> UPDATE element_info.csv NOW and then PRESS [RETURN] TO RE-CHECK IT, or ''q'' to break. <<<< \n','s');
    if ~isempty(reply) && strcmpi(reply,'Q') %SJ changed from (reply=='Q' || reply=='q')
        fprintf('\n keyboard here so you can see what up... ');
        keyboard
    end
    verifyElementInfo(subj, rootEEGdir);
    return
end
if max(iPHYS)>min(iCLIN) || max(iPHYS)>iSYNC || max(iCLIN)>iSYNC,
    fprintf('\n  ********* ERROR: element_info.csv rows should be ordered with all PHYS, then CLIN, then SYNC last *****');
    
    system(sprintf('open %s', elFilename));
    reply = input('\n >>>> UPDATE element_info.csv NOW and then PRESS [RETURN] TO RE-CHECK IT, or ''q'' to break. <<<< \n','s');
    if ~isempty(reply) && strcmpi(reply,'Q') %SJ changed from (reply=='Q' || reply=='q')
        fprintf('\n keyboard here so you can see what up... ');
        keyboard
    end
    verifyElementInfo(subj, rootEEGdir);
    return
end
    

    
% map '' --> '[]'
arrayCols = intersect(arrayCols, info.Properties.VariableNames);
for i = 1 : length(arrayCols)
    colName = arrayCols{i};
    rowMask = cellfun('isempty', info{:, colName});
    numEmpty = sum(rowMask);
    if numEmpty > 0
        info{rowMask, colName} = repmat({'[]'}, sum(rowMask), 1);
        fileAltered = true;
    end
end

% Alert user that there are commas that need to be removed from array columns
array_subt = info(:, arrayCols);
cells = array_subt{:,:};
good_cells = contains(cells, '[') & contains(cells, ']');
bad_columns = ~all(good_cells, 1);
if any(bad_columns)
    system(sprintf('open %s', elFilename));
    error('** Due to potential CSV-parsing errors, array-type columns may not contain commas.\nPlease edit element_info column(s) {%s} to use spaces instead of commas. **', strjoin(arrayCols(bad_columns), ' and '));
end



%%%%  check nkRecOrder --- JW removed nkRecOrder Column 8/2017  %%%%
% This column is NOT meant to accurately reflect the eeg montage order, because that order
% may change over time or for each session. It does incidentally usually roughly correspond with that
% order. The index is meant to be just a convenient numerical index for PHYS channels only.
% UTAH channels get a 0 for nkRecOrder
% ndx_utah_mask = strncmpi(info.tagName, 'utah', 4);
% phys_mask = strcmpi(info.chanType, 'PHYS') & ~ndx_utah_mask;
% num_phy = sum(phys_mask);
% ndx_pos_mask = info.nkRecOrder > 0;
% if ~all(ndx_pos_mask == phys_mask)
%     % set non-phys channel ndx to 0
%     info{~phys_mask, 'nkRecOrder'} = zeros(sum(~phys_mask),1);
%     info(phys_mask,:).nkRecOrder = [1:sum(phys_mask)]';
%     fileAltered = true;
% end
% if max(info.nkRecOrder) > num_phy
%     warning('element_info.csv has an nkRecOrder ndx higher than the number of PHYS channels (Okay for multi-implant)');
% end

%% _sh-related attribute checks
% Check _sh (shift, "multi-implant") subjects the following:
%   1) if only pre-_sh channels marked resected, error
%   2) if implant-2 channels are marked resected, show user additional implant-1
%      channels that should be marked as resected based on their location
%   3) Check that ictal and inter-ictal attributes of _sh channels and their originals exactly match
jack  = getJackTable(subj, rootEEGdir, 'createOnTheFly', 1); % don't read from file. We just do this so that we can use getLeads functions
ictal = getLeads(subj, rootEEGdir, 'markedAs', 'ictal',      'jackTable',jack);
inter = getLeads(subj, rootEEGdir, 'markedAs', 'interictal', 'jackTable',jack);
resec = getLeads(subj, rootEEGdir, 'markedAs', 'resected',   'jackTable',jack);
mask_sh = contains(jack.chanName, '_sh');
isShiftSubj = any(mask_sh);
if isShiftSubj
    preshift = @(chan) strrep(chan, '_sh', ''); % 'G_sh' -> 'G'
    addSh = @(chan) strcat(chan, '_sh'); % 'G' -> 'G_sh'
    addShC = @(chans) cellfun(addSh, chans, 'uniformOutput',0); % cells -> cells_sh
    preshiftC = @(chans) cellfun(preshift, chans, 'uniformOutput',0); % cells -> pre_sh
    toTagC = @(chans) unique(util_split_stringnum(chans)); % chans -> tags
    AUTO_RESECT_STR = 'resected channels identified programmatically based on resected shift channels.';
    
    %- use the existance of leads.csv to indicate that tal has been processed
    leadsCSV = fullfileEEG(rootEEGdir,subj,'tal','leads.csv');
    
    % 1) For each implant-1 element that is marked resected and that ends up shifting, verify that 
    %   its implant-2 _sh element has some resected leads marked. If not, this doesn't make sense:
    %   At the time of the resection, the _sh version of an element was most recent, so those xyz-
    %   locations are most accurate, but Kareem marked the older, pre-sh channels. Have Kareem resolve
    %   the discrepency by updating elementInfo.
    if ~isempty(resec) & exist(leadsCSV,'file')
        sh_tags = toTagC(jack.chanName(mask_sh));
        resec_tags = toTagC(resec);
        auto_resect_tags = toTagC(info.tagName(contains(info.notes, AUTO_RESECT_STR)));
        resec_presh_tags = preshiftC(toTagC(resec));
        resec_presh_tags = setdiff(resec_presh_tags, auto_resect_tags);
        resec_presh_tags = unique(resec_presh_tags(ismember(addShC(resec_presh_tags), sh_tags)));
        for itag = 1:numel(resec_presh_tags)
            tag = resec_presh_tags{itag};
            if ismember(tag, resec_tags) && ~ismember(addSh(tag), resec_tags)
                fprintf('\n\nverifyElementInfo error: Notify Kareem that channels have shifted; _sh channels (most recent) should likely be marked resected in elementInfo instead of original chans\n');
                system(sprintf('open %s', elFilename));
                error('%s shifted between implants but has pre-shifted channels marked resected instead of post-shift (more recent) channels', tag);
            end
        end
        
        % 2) if implant-2 channels are marked resected, identify additional implant-1  channels that should be marked as resected based on their location (only if not marked yet)
        DIST_WITHIN = 5; % mm
        SRC_IMPLANT = 2;
        
        % get channels ("cp2Chan") in first implant that should be resected based on distance to those marked resected in 2nd implant:
        [~, cp2Chan, ~] = elecAttribsCrossImplant(subj, rootEEGdir, SRC_IMPLANT, DIST_WITHIN, 'resected', 1);
        [tags, chanNums] = util_split_stringnum(unique(cp2Chan));
        tags_unique = unique(tags);
        % For each tag with exactly 0 channels marked resected, update its resected list
        for i = 1:numel(tags_unique)
            irow = strcmpi(info.tagName, tags_unique{i});
            resectArray = eval(info.isResected{irow}); 
            if isempty(resectArray)
                resectStr = '';
                for j = 1:numel(chanNums)
                    if strcmpi(tags_unique{i}, tags{j})
                        resectStr = sprintf('%s %d', resectStr, chanNums(j));
                    end
                end
                resectStr = ['[', strtrim(resectStr), ']'];
                info.isResected{irow} = resectStr;
                fileAltered = 1;
                info.notes{irow} = strtrim([info.notes{irow}, ' ', AUTO_RESECT_STR]);
            end 
        end
    elseif ~isempty(resec) & ~exist(leadsCSV,'file')
        fprintf('\n WARNING: element info has resection info with multiple implants surgeries, \n but automatic correction of element info cannot be made before tal is complete.\n Rerun createMasterJack once tal is done. skipping step for now.');
        pause(2);
    end
    
    % 3) Check that ictal and inter-ictal attributes of _sh channels and their pre-_sh originals exactly match
    ictal_tags = toTagC(ictal);
    inter_tags = toTagC(inter);
    sh_tags = toTagC(jack.chanName(mask_sh));
    attrib_tags = union(ictal_tags, inter_tags);    
    pre_sh_attrib_tags = unique(preshiftC(attrib_tags)); % root (sh-stripped) of all tags with ictal/interictal marked
    OVERRIDE_STR_ICTAL ='Ictal mismatch acknowledged.';
    OVERRIDE_STR_INTER ='Interictal mismatch acknowledged.';
    
    
    for itag = 1:numel(pre_sh_attrib_tags)
        tag = pre_sh_attrib_tags{itag};
        mask_row_pre_sh = strcmpi(info.tagName, tag);
        mask_row_sh = strcmpi(info.tagName, addSh(tag));
        if sum(mask_row_sh) == 0
            continue;
        end
        
        % Ictal
        isIctalStrings(1) = info.isIctal(mask_row_pre_sh);
        isIctalStrings(2) = info.isIctal(mask_row_sh);
        isIctalArrays{1} = eval(isIctalStrings{1});
        isIctalArrays{2} = eval(isIctalStrings{2});
        if isempty(isIctalArrays{1}) && ~isempty(isIctalArrays{2})
            info.isIctal(mask_row_pre_sh) = isIctalStrings(2);
            fileAltered = 1;
            
        elseif ~isempty(isIctalArrays{1}) && isempty(isIctalArrays{2})
            info.isIctal(mask_row_sh) = isIctalStrings(1);
            fileAltered = 1;
            
        elseif ~isempty(isIctalArrays{1}) && ~isempty(isIctalArrays{2}) && ~isequal(isIctalArrays{1},isIctalArrays{2}) %SJ change from ~all(isIctalArrays{1} == isIctalArrays{2})
            fprintf('\nError: %s isIctal mismatch between shifted (_sh) and original. These are expected to match!\n', tag);
            if any(contains(info.notes(mask_row_sh | mask_row_sh), OVERRIDE_STR_ICTAL))
                fprintf('\tAcknowledged.\n');
            else
                if inputYN('Do you want to open/update element_info to correct discrepency now?\nSelect "n" to *allow* the mismatch (unusual)');
                    system(sprintf('open %s', elFilename));
                    input('Press Enter to recheck');
                    verifyElementInfo(subj, rootEEGdir, varargin{:})
                    return;
                else
                    info.notes{mask_row_pre_sh} = addNote(info.notes{mask_row_pre_sh}, OVERRIDE_STR_ICTAL);
                    info.notes{mask_row_sh} = addNote(info.notes{mask_row_sh}, OVERRIDE_STR_ICTAL);
                    fileAltered = 1;
                end
            end
        end
        
        % Same for interictal
        isInterStrings(1) = info.isInterictal(mask_row_pre_sh);
        isInterStrings(2) = info.isInterictal(mask_row_sh);
        isInterArrays{1} = eval(isInterStrings{1});
        isInterArrays{2} = eval(isInterStrings{2});
        if isempty(isInterArrays{1}) && ~isempty(isInterArrays{2})
            info.isInterictal(mask_row_pre_sh) = isInterStrings(2);
            fileAltered = 1;
            
        elseif ~isempty(isInterArrays{1}) && isempty(isInterArrays{2})
            info.isInterictal(mask_row_sh) = isInterStrings(1);
            fileAltered = 1;
            
        elseif ~isempty(isInterArrays{1}) && ~isempty(isInterArrays{2}) && ~isequal(isInterArrays{1},isInterArrays{2}) %SJ Change from ~all(isInterArrays{1} == isInterArrays{2})
            fprintf('\nError: %s isInterictal mismatch between shifted (_sh) and original. These are expected to match!\n', tag);
            if any(contains(info.notes(mask_row_pre_sh | mask_row_sh), OVERRIDE_STR_INTER))
                fprintf('\tAcknowledged.\n');
            else
                if inputYN('Do you want to open/update element_info to correct discrepency now?\nSelect "n" to *allow* the mismatch (unusual)');
                    system(sprintf('open %s', elFilename));
                    input('Press Enter to recheck');
                    verifyElementInfo(subj, rootEEGdir, varargin{:})
                    return
                else
                    info.notes{mask_row_pre_sh} = addNote(info.notes{mask_row_pre_sh}, OVERRIDE_STR_INTER);
                    info.notes{mask_row_sh} = addNote(info.notes{mask_row_sh}, OVERRIDE_STR_INTER);
                    fileAltered = 1;
                end
            end
        end
    end
end


%% display any changes made and save if user chooses to do so
if fileAltered
    disp('ORIGINAL:');
    disp(info_orig);
    disp('CORRECTED:');
    disp(info);
    fprintf('verifyElementInfo has made minor corrections to element_info.csv.\n');
    fprintf('Usual correction: capitalization of column titles or added [].\n');
    fprintf('Before overwriting, please verify that things look okay and nothing crazy happened.\n\n');
    choice = input('Save corrected element_info.csv (displayed above)? [y/n] ', 's');
    if strcmpi(choice, 'y')
        
        %writetable(info, elFilename);
        elFilenameTemp = [elFilename '.temp.csv'];
        writetable(info, elFilenameTemp);
        
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(elFilenameTemp,elFilename);
        if SUCCESS,
            fprintf('\n element_info.csv successfully updated');
        else
            fprintf('\n Uh Oh, element_info.csv NOT successfully updated... maybe its open? \n Manually delete and rename element_info.csv.temp.csv --> element_info.csv\n');
            keyboard
        end
        
        
        if exist(fname_jacksheet,'file')
            delete(fname_jacksheet);
            if ~exist(fname_jacksheet,'file')
                fprintf('\n since changes were made, jacksheetMaster was just deleted and needs to be recreated\n');
            else
                fprintf('\n counldnt delete jacksheetmaster.csv, but it should be recreated now');
                keyboard
            end
        end
    end
end

end % verifyElementInfo

function note = addNote(note, txt)
    if ~ischar(note) || strcmpi(note, 'NaN')
        note = txt;
    else
        note = [note, ' ', txt];
    end        
    note = strtrim(note);
end

