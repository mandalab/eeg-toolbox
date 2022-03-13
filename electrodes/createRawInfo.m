function [rawInfo, rawChanNotInJacksheet, jackChanNotInRaws] = createRawInfo(subj, rootEEGdir)
% createRawInfo writes raw_info.csv based on a subject's raws
%   check for, and facilitate correction of, duplicate channel names and shift chanenls "_sh"
%
% Usage: rawInfo = createRawInfo(subj, rootEEGdir)
%
% Input:
%   subj
%   rootEEGdir
%
% Output:
%   rawInfo - table
%   rawChanNotInJacksheet -- list of channels found in raws but not specified in element_info
%
% Output Files:
%   docs/raw_info.csv
%
% Description:
%   This is designed to replicate what is done originally in createMasterJack.
%   Its purpose was to attempt to provide an updated sheet in case more raws
%   were added after createMasterJack was done. It was first called in
%   dungeon_toolbox/GetFilesAndChannels
%
%   It is important that this gets updated when createMasterJack does
%
% Dependencies:
%   align/findRaws
%   getChanNames
%   getJackTable
%
%
% Revision History
%   07/16 MST - Created
%   08/17 JHW - output all channels, even those not specified in element_info.csv
%             - add rawChanNotInJacksheet output to be used by createMasterJack
%   09/17 JHW - Move MT's checks from nk_split to here:  duplicate entries detected and fixed; _sh channel pushed to .21E's
%   12/17 melkalliny - code now handles either NK or BR raws
%   7/2020 SNJ- functionality/clarity improvements
%
% See also: createMasterJack

RAWINFO_FILEPATH = fullfile(rootEEGdir, subj, 'docs', 'raw_info.csv');


%-  grab the jackSheetMaster list of split channels
jacktable  = getJackTable(subj, rootEEGdir);
cJackChans = jacktable.chanName;


% get channels from raws
rawDir = fullfile(rootEEGdir, subj, 'raw');
raws   = findRaws(rawDir, {});
if isempty(raws)
    return;
end

%- location of file that will store raw change history
rawChangeFile = fullfile(rootEEGdir, subj, 'raw/raw21E_changeHistory.txt');


%- Loop over all Raw's and grab all saved channels from each raw
c21Es           = strcat(raws, '.21E');
chansPerRaw     = cell(1, length(c21Es));
codesPerRaw     = cell(1, length(c21Es)); %- string of codes in .21E indicating each row
cAllRawNewChans = {};
foundDuplicate  = zeros(length(c21Es),1);
for i = 1:length(c21Es)
    d21E = c21Es{i};
    [parent,filestem] = fileparts(d21E);
    if     contains(filestem,'ieeg')     | contains(filestem,'INST') | contains(filestem,'ns2')
        dEEG = fullfile(parent, [filestem '.ns2']);
        c21Es{i} = strrep(c21Es{i},'21E','ns2');
    elseif contains(filestem,'cervello') | contains(filestem,'EEG')
        dEEG = fullfile(parent, [filestem '.TRC']);
        c21Es{i} = strrep(c21Es{i},'21E','TRC');
    else
        dEEG = fullfile(parent, [filestem '.EEG']);
    end
    % the "21E" could actually be
    assert(exist(dEEG, 'file') > 0, 'No EEG file for %s', d21E);
    
    % Get channel names and codes from EEG and 21E file
    [cRawChans  cRawCodes] = getChanNames(d21E, dEEG); % channels as specified in this .21E (can return duplicates if exist, doesn't return empty strings ''
    cRawChans = deblank(cRawChans); 
    chansPerRaw(i) = {cRawChans};
    codesPerRaw(i) = {cRawCodes};
    
    %- look for duplicates
    if length(unique(cRawChans))<length(cRawChans)
        foundDuplicate(i)=1;
    end
   
    %- find any channels NOT in jacksheetMaster... these will be added to bottom of raw_info
    cThisRawNewChans = setdiff(cRawChans, cJackChans);               % in raw but not in jacksheet
    cAllRawNewChans  = union(cAllRawNewChans, cThisRawNewChans);     % all new chans from raws
end
rawChanNotInJacksheet = cAllRawNewChans;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- Check for a shift channel --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- look for _sh shifted channel specification in element_info (and thus jacksheetMaster)
% actualName_ALL is updated to _sh along with .21E if user specifies to do so
shift_mark = '_sh';
if any(strfound(cJackChans, shift_mark))
    info = getElementInfo(subj, rootEEGdir);
    shiftChangeTotal = 0; fprintf('\n');
    for i = 1:length(c21Es)
        if  ~(contains(c21Es{i},'ns2') || contains(c21Es{i},'TRC'))
        shiftChangeTotal = updateShiftChannels(c21Es{i}, info, cJackChans, rawChangeFile) + shiftChangeTotal;
        end
    end
    if shiftChangeTotal>0
        reply = input(sprintf('\n HEADS UP: %d changes made to .21Es to account for electrode_info.csv shift channels. Press [RETURN] to continue or ''q'' to break',shiftChangeTotal),'s');
        if ~isempty(reply) && (reply=='Q' | reply=='q')
            fprintf('\n keyboard here so you can see what up... ');
            keyboard
        end
        %- recursive call here so updated .21E is incorporarted into createMasterJack when corrected
        [rawInfo, rawChanNotInJacksheet, jackChanNotInRaws] = createRawInfo(subj, rootEEGdir);
        return;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-- Fix duplicate entries --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- if duplicate channels detected there is an error in the .21E... help to fix it
%- found one or more duplicates in this .21E... fix this on and all others then rerun
if sum(foundDuplicate)>0
    fprintf('\n  ERROR: %d of %d .21E files have at least 1 dupicate entry. \n  Focus on first erroneous .21E and try to fix it: \n', sum(foundDuplicate), length(foundDuplicate));
    
    %- Idenify the duplicates
    iFirstDupFile = find(foundDuplicate,1);
    cRawChans = chansPerRaw{iFirstDupFile};
    cRawCodes = codesPerRaw{iFirstDupFile};
    [~,uniqueIdx,~] = unique(cRawChans,'stable');
    dupIdx = setdiff([1:length(cRawChans)],uniqueIdx);
    
    
    %- list all duplicates so user knows whats going on
    fprintf('\n Copy of .21E contents (%s)', c21Es{iFirstDupFile});
    for ii=1:length(cRawChans)
        fprintf('\n   %s = %s',cRawCodes{ii},cRawChans{ii});
        if sum(strcmp(cRawChans,cRawChans{ii}))>1, fprintf('  << duplicate entry'); end
    end
    fprintf('\nDuplicate channel(s) in %s: ',  c21Es{iFirstDupFile});
    for iDup=1:length(dupIdx),  fprintf(' %s (%d instances), ', cRawChans{dupIdx(iDup)}, sum(strcmp(cRawChans,cRawChans{dupIdx(iDup)}))); end
    fprintf('\n');
    
    
    %- now attempt to correct the first one (hopefully only 1 present!!!)
    dupChan   = cRawChans{dupIdx(1)};
    dupInd    = find(strcmp(cRawChans,dupChan));
    ELEC_file = c21Es{iFirstDupFile};
    correctRepeatedChannel(dupInd, codesPerRaw{iFirstDupFile}, ELEC_file, dupChan, c21Es,rawChangeFile);
    
    
    %- should be fixed: output the revised .21E to show
    fprintf('\n Updated .21E contents (%s)', c21Es{iFirstDupFile});
    d21E = c21Es{iFirstDupFile};
    [parent,filestem] = fileparts(d21E);
    dEEG = fullfile(parent, [filestem '.EEG']);
    [cRawChans  cRawCodes] = getChanNames(d21E, dEEG); % channels as specified in this .21E (can return duplicates if exist, doesn't return empty strings ''
    
    %- if all duplicates removed, show the resultant list so user knows its all good.  else the list will be generated with next duplicate in recursive call below
    if length(unique(cRawChans))==length(cRawChans)
        for ii=1:length(cRawChans)
            fprintf('\n   %s = %s',cRawCodes{ii},cRawChans{ii});
            if sum(strcmp(cRawChans,cRawChans{ii}))>1, fprintf('  << duplicate entry'); foundDup=foundDup+1; end
        end
        fprintf('\n GOOD WORK, all duplicates removed.  Take a look above');
    end
    fprintf('\n');
    
    %- recursive call here so updated .21E is incorporarted into createMasterJack when corrected
    [rawInfo rawChanNotInJacksheet jackChanNotInRaws] = createRawInfo(subj, rootEEGdir); %#ok<*NCOMMA>
    return;
end


% generate raw_info file -- ALL channels included (those that are in jacksheetMaster and those that arent)
rawInfo = table();
numRaws = length(chansPerRaw);
cAllChans = [cJackChans; cAllRawNewChans(:)];
jackChanNotInRaws = {};
for i = 1 : length(cAllChans)
    chan = cAllChans(i);
    rawInfoRow.chanName = chan;
    rawInfoRow.in_jackSheetMaster = double(ismember(chan, cJackChans));
    rawInfoRow.raw_info = '';
    rawHasChanCnt = 0;
    for j = 1 : numRaws
        [~,timestamp] = fileparts(fileparts(c21Es{j})); % should be timestamp
        % if doubleNSP-blackrock, name the raws differently
        % edit- actually, ^^ caused errors during step 5. commenting it out
        %if     contains(c21Es{j},'ieeg1')|contains(c21Es{j},'INST0'),
        %    rawName = ['raw_' timestamp];
        %elseif contains(c21Es{j},'ieeg2')|contains(c21Es{j},'INST1'),
        %    rawName = ['raw_' timestamp];
        %else
        %    rawName = ['raw_' timestamp];
        %end
        rawName = ['raw_' timestamp];
        cChansInRaw = chansPerRaw{j};
        rawHasChan = double(length(find(ismember(cChansInRaw, chan), 1)));
        rawInfoRow.(rawName) = rawHasChan;
        rawHasChanCnt = rawHasChanCnt + rawHasChan;
    end
    if rawHasChanCnt==0
        rawInfoRow.raw_info = {sprintf('MISSING_FROM_ALL_RAWS')};
        jackChanNotInRaws{end+1} = chan{1}; %- convert from cell to string
    else
        rawInfoRow.raw_info = {sprintf('FOUND_IN_%d_OF_%d_RAWS', rawHasChanCnt, numRaws)};
        
    end
    [rawInfo; struct2table(rawInfoRow)];
    rawInfo = [rawInfo; struct2table(rawInfoRow)]; %#ok<AGROW>
end

writetable(rawInfo, RAWINFO_FILEPATH);
end




%---------------------------------------------------------------------------------------------------------------------------%
%------------------------------------------------- Duplicate-Fix Functions -------------------------------------------------------%
%---------------------------------------------------------------------------------------------------------------------------%


%---------------------------------------------------------------------------------------------------------------------------%
function correctRepeatedChannel(actualChan, actualCode_ALL, ELEC_file, masterName, c21Es,rawChangeFile)
newNames = cell(length(actualChan), 1);
chanNums = cell(length(actualChan), 1);


% Give user a pulse visualization of the repeated channel to allow them to examine
fprintf('\n\nGenerating visualization for each repeated channel... (this takes up to a minute, so hold on)\n');
for i = 1:length(actualChan)
    fprintf('\n==============================--> Instance %d: %s=%s <--==============================', i, actualCode_ALL{actualChan(i)}, masterName);
    hFig = eegPulseVisualize(fileparts(ELEC_file), {masterName}, i);
    set(hFig(1),'name',sprintf('CHAN %s, INSTANCE #%d',masterName,i));
end


% User Identifies which is the true channel to keep
success = false;
while success==false
    [iGood, success] = str2num(input(sprintf('\nWhich repeated instance is the TRUE %s channel? Look at the figure titles. (Enter instance number 1-%d): ', masterName, length(actualChan)), 's'));
    if ~(success && 1 <= iGood && iGood <= length(actualChan))
        success = false;
        fprintf(' -- input no good... try again idiot');
    end
end
fprintf('\nInstance %d: %s=%s marked correct\n', iGood, actualCode_ALL{actualChan(iGood)}, masterName);
fprintf('\nThe other instances will be renamed (makes .21E file changes but backs up original file).\n');


% Rename the other channels
for i = 1:length(actualChan)
    if i == iGood, continue; end % good name will have empty newName entry
    newNames{i} = input(sprintf('Rename instance %d (channel %s) to: ', i, actualCode_ALL{actualChan(i)}), 's');
    chanNums{i} = actualCode_ALL{actualChan(i)};
end

% Create new file(s)


% rename all raws (this used to be an option, now its the default
for i = 1:length(c21Es)
    %if rewriteDupChannel([eegRawList{i} '.21E'], newNames, chanNums, backupAffix, masterName)
    rewriteDupChannel(c21Es{i}, newNames, chanNums, masterName,rawChangeFile);
end

end % correctRepeatedChannel function



%---------------------------------------------------------------------------------------------------------------------------%
%- helper function that loops through .21E files and updates/corrects channel names
function rewriteDupChannel(fileName, newNames, chanNums, masterName,rawChangeFile)
iRepeat = 0;
temp_name = fullfile(tempdir, 'matlab.nk_split.rewriteChannel');


fid_raw = fopen(fileName, 'r');
fid_tmp = fopen(temp_name, 'w');
outStr  = '';

line = 'enterLoop';
while ischar(line)
    line = fgets(fid_raw);
    % check line has repeated name with the right channel number and associated renaming
    if strfind(line, masterName)
        iRepeat = iRepeat + 1;
        if ~isempty(newNames{iRepeat}) && ~isempty(strfind(line, chanNums{iRepeat}))
            % copy a modified line
            fprintf(fid_tmp, '%s', strrep(line, masterName, newNames{iRepeat}));
            outStr = sprintf('%s  -->  %s', strtrim(line), strtrim(strrep(line, masterName, newNames{iRepeat})));
            continue;
        end
    end
    fprintf(fid_tmp, '%s', line);
end

fclose(fid_tmp);
fclose(fid_raw);

if     iRepeat==0
    %fprintf('\n  NOTE: didnt find ANY instances of %s in %s. \nMaybe duplicate happened after rearrangement. No change made to .21E', masterName, fileName);
elseif iRepeat==1
    %fprintf('\n  NOTE: found only 1 instance of %s in %s. \nPossibly duplicates already fixed. No change made to .21E', masterName, fileName);
else
    %fprintf('\n  NOTE: found %s instances of %s in %s. Correcting .21E file now (saving backup of original if not backuped up already)', iRepeat, masterName, fileName);
    [path, name, ext] = fileparts(fileName);
    backupAffix = '.original';
    if ~exist( fullfile(path, [name ext backupAffix]),'file' )
        copyfile(fileName, fullfile(path, [name ext backupAffix])); % create backup only ONCE (so original version can still be accessed after multiple fixes
    end
    copyfile(temp_name, fileName); % create new
    delete(temp_name); % remove temp
    
    %- keep a record of changes in a single file
    fid_rawRecord = fopen(rawChangeFile,'a+');
    fprintf(fid_rawRecord,'\n%s  (%s)',outStr,fileName);
    fclose(fid_rawRecord);
    fprintf('\n NOTE: duplicate entries in .21E corrected:  %s  (%s)',outStr,fileName);
    fprintf('\n  updated %s',fileName);
end


end %function success = rewriteChannel




%---------------------------------------------------------------------------------------------------------------------------%
%------------------------------------------------- Shift_fix functions -------------------------------------------------------%
%---------------------------------------------------------------------------------------------------------------------------%



%---------------------------------------------------------------------------------------------------------------------------%
function numChangesMade = updateShiftChannels(ELEC_file, info, jackmaster_names, rawChangeFile)

shift_mark  = '_sh';
shift_names = jackmaster_names(strfound(jackmaster_names, shift_mark));
[shift_tags,~] = util_split_stringnum(shift_names);
shift_tags  = unique(shift_tags);

% defaults to 0 in case early session
numChangesMade   = 0;

for i = 1:length(shift_tags) % tag loop
    shift_tag = shift_tags{i};
    old_tag   = strsplit(shift_tag, shift_mark);
    old_tag   = old_tag{1}; % 1st peice (not the _sh)
    
    %- how many times do we expect to see this shift tag (if past the date of shift)
    numExpected = sum(strfound(jackmaster_names,shift_tag));
    
    %- this is a shift channel, but has the shift change already been made to .21E? Check with grep
    [grepOut,grepStr] = unix(sprintf('grep %s %s', shift_tag, ELEC_file));
    %isInFile    = ~grepOut; % grep return 0 when found
    numFound    = length(strfind(grepStr,'=')); %- each line has an equal sign in it
    
    %if ~isInFile    %- grep version was missing case in NIH049 where single channel label was mis-spelled and then shift wasn't being applied after correction
    if numFound<numExpected
        
        % get date
        info_row       = info(strcmp(info.tagName, shift_tag), :);
        shift_date     = info_row.dateIn;
        try shift_date = datetime(shift_date,'format','MM/dd/yy');
        catch e
            fprintf('%s\n', e.message);
            error('nk_split cannot continue: Please fill in the dateIn column in element_info.csv for all shifted (_sh) channels');
        end
        
        sameDay  = @(d1,d2) (d1.Year == d2.Year) && (d1.Month == d2.Month) && (d1.Day == d2.Day);
        aBeforeB = @(d1,d2) datenum(d1) < datenum(d2);
        
        % check to see if this session is after the shift
        [~,sessionDate] = fileparts(fileparts(ELEC_file));
        sessionDate = datetime(sessionDate(1:11),'format','yyMMdd_HHmm');
        if sameDay(shift_date, sessionDate) || aBeforeB(shift_date, sessionDate)
            
            
            % just this session, check if same day
            if sameDay(shift_date, sessionDate)
                fprintf('The session is on the same date as the shift (%s)\n', datestr(sessionDate));
                shift_date = addUserTimeToDate(shift_date);
            end
            
            if ~aBeforeB(shift_date, sessionDate)
                return; % nothing need change!
            end
            
            % Create new file(s)
            numChangesMade = rewriteTagShift(ELEC_file, old_tag, shift_tag,rawChangeFile)+numChangesMade;
            
        end % after shift
    end % if ~inFile
end % tag loop
end % function



%---------------------------------------------------------------------------------------------------------------------------%
function success = rewriteTagShift(fileName, oldTag, newTag,rawChangeFile)
success = false;
temp_name = fullfile(tempdir, 'matlab.nk_split.rewriteChannel');

fid_raw = fopen(fileName, 'r');
fid_tmp = fopen(temp_name, 'w');
outStr  = '';


line = 'enterLoop';
anyFound = 0;
while ischar(line)
    line = fgets(fid_raw);
    lineFound = 0;
    % check line has repeated name with the right channel number and associated renaming
    
    % in general, look for old tag and not new tag
    isFound = @(s) strfound(s, oldTag) && ~strfound(s, newTag);
    if isFound(line)
        if strfound(line, '=')
            % double check it's not a superstring (check length)
            temp = strsplit(line, '=');
            rhs = temp{2};
            rhs = util_split_stringnum(rhs); % <-- look at tag only
            if strcmpi(rhs, oldTag)
                lineFound = 1;
            end
        else
            lineFound = 1;
        end
    end
    
    if lineFound
        anyFound = 1;
        fprintf(fid_tmp, '%s', strrep(line, oldTag, newTag));
        outStr = sprintf('%s  -->  %s', strtrim(line),strtrim(strrep(line, oldTag, newTag)));
    else
        fprintf(fid_tmp, '%s', line);
    end
end
success = anyFound;

fclose(fid_tmp);
fclose(fid_raw);
if ~success, return; end % didn't find name

[path, name, ext] = fileparts(fileName);
backupAffix = '.original';
backupFile = fullfile(path, [name ext backupAffix]);
if ~exist(backupFile,'file')
    copyfile(fileName, backupFile); % create backup
end
copyfile(temp_name, fileName); % create new
delete(temp_name); % remove temp


%- keep a record of changes in a single file
fid_rawRecord = fopen(rawChangeFile,'a+');
fprintf(fid_rawRecord,'\n%s  (%s)',outStr,fileName);
fclose(fid_rawRecord);
fprintf('\n NOTE: element_info.csv specified shift ''_sh'' channel corrected:  %s  (%s)',outStr,fileName);
fprintf('\n  updated %s',fileName);

end %function success = rewriteChannel



%---------------------------------------------------------------------------------------------------------------------------%
function dt = addUserTimeToDate(day_date)
format = 'MM-dd-yy HHmm';
text = 'Enter estimated time of shift (HH:mm) eg 14:00 for 2pm: ';
valid = 0;
while ~valid
    date_str = [datestr(day_date) ' ' input(text,'s')];
    try
        dt = datetime(date_str,'format',format);
        valid = 1;
    catch
    end
end

end %function dt = addUserTimeToDate(day_date)

