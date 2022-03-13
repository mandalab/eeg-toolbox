function name = phys2name(subj, rootEEGdir, timestamp, physNum, montageLabels, montageJack)
% PHYS2NAME translates the physical channel number to a channel name
%
% name = phys2name(subj, rootEEGdir, timestamp, physNum [,montageLabels, montageJack])
%
% Input:
%   subj       - subject name
%   rootEEGdir - contains subject folder
%   timestamp  - name of folder which contains session file
%   physNum    - numeric physical channel you want (e.g. 128)
%
% Optional Input:
%   montageLabels - if you have read the montage, pass in the labels
%   montageJack   - if you have read the montage, the jacks (e.g. A13)
% Output:
%   name
%
%
% This function is based on the same hard-coded scheme that nk_split uses.
% Presumably the scheme that nk_split uses is based on the idiosyncracies of
% the Nihon Khoden's .21E channel numbering scheme. In this scheme, each line of
% the .21E file receives a number (a "code") and a label (the form is: 004=G1). Presumably
% these lines represent electrode channels. The catch is that not all lines represent
% physical channels. By "Physical channels," we mean those lines/channels which
% could receive signals from real channels, and these should be in the same order as
% the .21E, but non-physical channels are skipped. Note in this code, the term "good
% code" refers to lines which represent physical channels.
%
%
% Dependencies:
%   basic/lsCell
%
% Called by extractStimEventsFromAll and extractStimEventsFromStimLog_v10 (in behavioralProcessing)
%
% Revision History:
%   06/01/17 MST - Created
%   03/30/18 MST - Try fixing blackrock non-found errors
%   01/02/19 JW  - JW making this work correctly for blackrock and cervello files
%   09/29/20 SJ  - Getting rid of call to chanNames.txt, formatting changes, switch to readtableSafe
%   10/12/20 SK  - Adding case for isCervello to use jacksheetCV
%
% See also: nk_split, readMontage

% -------------------------------------------------------------------------

%- initialized the return variable
name = [];




%% Figure out if this is Blackrock or Cervello.. if so, use "jacksheetBR.csv" to figure out the phys2chan mapping
rawLocationRoot = fullfile(rootEEGdir,subj,'raw',timestamp);
rawLocationStim = fullfile(rootEEGdir,subj,'raw','STIM_MAP',timestamp);
if     exist(rawLocationRoot,'dir')
    rawLocation = rawLocationRoot;
elseif exist(rawLocationStim,'dir')
    rawLocation = rawLocationStim;
else
    fprintf('\n ERROR: cant find raw directory %s', rawLocationRoot); 
    keyboard
    return 
end

rawDirList = lsCell(rawLocation);
isBlackrock = any(contains(rawDirList, '.nev')) || any(contains(rawDirList, '.ns')) || any(contains(rawDirList, '.TRC'));
isCervello  = any(contains(rawDirList, '.TRC')); %- subset of "isBlackrock"... hardware is blackrock but recorded file is different type
if isBlackrock
    %- open the BR jacksheet to get true mapping between physical and name
    %brJacksheet = fullfileEEG(rawLocation,'jacksheetBR_complete.csv');
    brJacksheet = fullfileEEG(rawLocation,'jacksheetBR_local.csv');
    if ~exist(brJacksheet,'file')
        if any(contains(rawDirList,'jacksheetBR'))
            %- alternate exists, so that that
            brJacksheet = fullfileEEG(rawLocation,rawDirList{find(contains(rawDirList,'jacksheetBR'),1,'first')});
        else
            %- create a local br jacksheet
            if isCervello
                fprintf('\n ERROR: cant make a phys2chan mapping with cervello files.  \n  Copy jacksheetBR_complete.csv made from a blackrock session with same recording settings \n   to %s and rerun. \n  Talk to JW if stuck',rawLocation);
                keyboard
                error('\n fix it!');
            else
                fprintf('\n HEADS UP: about to make a quick local jacksheetBR for phys2chan for %s, \n  but it would be better to grab the "jacksheetBR_complete.csv" from FRNU56.  \n  Ask John what this means', rawLocation);
                keyboard
                jt = makeJacksheetBR(rawLocation,'',0,1);
                brJacksheet = fullfileEEG(rawLocation,'jacksheetBR_local.csv');
            end
        end
    end %- ~exist brJacksheet)
    
    %- table found or just created... read it now
    jackTable = readtableSafe(brJacksheet); %SJ: changed from readtable
    
    iRow = find(jackTable.PhysicalChan == physNum);
    if length(iRow)~=1
        fprintf('\n\n Uh oh... this shouldnt happen. physical channel (%d) found %d times in jacksheetBR (should find exactly 1)', physNum, length(iRow));
        disp(iRow);
        keyboard
    else
        % SJ: Do we need to change this to jackTable.ChanNameNew if applicable?
        name = deblank(jackTable.ChanName{iRow}); %- for some reason jackTable.ChanName is returning trailing zeros... strtrim wont cut them but deblank will
        
        %SJ replacing the following commented code with this (Centralized pulling of chanNames from jacksheetBR. chanNames.txt is now obsolete):
        if isCervello
            % If Cervello, use jackSheetCV_local
            if contains(brJacksheet,'jacksheetBR_complete')
                cvJacksheet = strrep(brJacksheet,'jacksheetBR_complete','jacksheetCV_local');
                if exist(cvJacksheet,'file')
                    jackTable = readtableSafe(cvJacksheet);
                else
                    fprintf('%s\n','jacksheetCV_local does not exist. Go create that now');
                    keyboard
                end
            else
                keyboard % Usually if isCervello, it should have grabbed the jacksheetBR_complete...
            end
        end
            
        chanNames = getChansFromJack(jackTable);
        
%         %- with blackrock and cervello, there is a CHANCE that the name stored in the ns2/nev/trc is NOT what was used
%         %   in prepAndAlign... confirm the channel name is in "chanNames.txt", which is how we rename in prepAndAlign
%         fid = fopen(fullfile(rawLocation,'ChanNames.txt'));
%         chanNames = textscan(fid,'%s'); chanNames = chanNames{1};
%         fclose(fid);
        
        %- did we find it in ChanNames?
        if ~any(strcmp(chanNames,name))
            
            %- didn't find it, but one possibility is that this became a shift channel...
            tagEnd = find(name >= '0' & name <= '9',1,'first')-1;
            strNum = name(tagEnd+1:end); 
            strNumNoLeadZero = strNum(find(strNum > '0',1,'first'):end);  %- shift code cuts out leading zeros... so do that here
            nameShift = [name(1:tagEnd) '_sh' strNumNoLeadZero];
            
            if any(strcmp(chanNames,nameShift))
                fprintf('\n phys2name dynamically converted to a shift name to make a match (%s -> %s)', name,nameShift);
                name = nameShift;
            else
                %- wasnt a shift channel... not sure how to deal with this
                fprintf('\n\n uh oh... chan name found in jacksheetBR does NOT exist in "ChanNames.txt", so it was probably renamed.\n get john');
                keyboard;
                error('\n i said get john');
                %- somehow map back to ChanNames....
            end
        end
    end
    
    return;
    
end  %- if isBlackrock





%%%%%%%% FOLLOWING CODE IS FOR Nihon Khoden RAW FILES

% Constants and declarations
DC_CODES   = [42:73]; %#ok<*NBRAK>
GOOD_CODES = [0:36, 74:75, 100:253];    % classic good codes... does not include DC
BAD_NAMES  = {'E'};
BANKS      = 'ABC????';

%GOOD_CODES = union(GOOD_CODES, DC_CODES);




DBG = 0;

phys_codes = {};
phys_names = [];

if nargin < 5 || isempty(montageLabels)
    montageLabels = [];
    montage_col_name = '';
else
    montage_col_name = ' | MONTAGE';
end
if nargin < 6
    montageJack = [];
end

% -------------------------------------------------------------------------
% Get 21E file
subj_dir = fullfile(rootEEGdir, subj);
assert(exist(subj_dir,'dir') > 0, 'File not found: %s\n', subj_dir);

%- first check in STIM_MAP folder... if not there check in raw folder and spit and error if not there
session_dir = fullfile(subj_dir, 'raw/STIM_MAP', timestamp);
if ~exist(session_dir,'dir')
    session_dir     = fullfile(subj_dir, 'raw', timestamp);
    assert(exist(session_dir,'dir') > 0, 'File not found: %s\n', session_dir);
end

name_21e = char(lsCell([session_dir '/*.21E']));
filename_21e = fullfile(session_dir, name_21e);
assert(exist(filename_21e, 'file') > 0, 'File not found: %s\n', filename_21e);
fd = fopen(filename_21e, 'r');
temp_cells = textscan(fd, '%s%s', 'delimiter','=');
fclose(fd);

% -------------------------------------------------------------------------
% Read channel codes from 21E
[allCodes, allNames] = deal(temp_cells{:});
endRange = find(strcmpi(allCodes, '[SD_DEF]')); % all electrodes are before this index
allCodes  = allCodes(1:endRange-1);  %stores the codes
allNames  = allNames(1:endRange-1);  %stores the names
phys_chan_ndx = 0;

if DBG
    row_header = 'CODE | CHANNEL | PHYS | LABEL';
    %if isempty(montageJack), row_header = [row_header '_GUESS']; end
    row_header = [row_header montage_col_name];
    disp(row_header);
end


for i = 1 : length(allCodes)
    code_str = allCodes{i};
    code_num = str2double(code_str);
    is_numeric_code = (~isempty(code_num) && ~isnan(code_num));
    i_name = allNames{i};
    
    if is_numeric_code && ismember(code_num, GOOD_CODES) && ...
            ~ismember(i_name, BAD_NAMES) % && ~isempty(name)
        
        % then this channel is a physical channel
        phys_chan_ndx = phys_chan_ndx + 1;
        phys_names = [phys_names i_name];
        phys_codes = [phys_codes code_num];
        
        if phys_chan_ndx == physNum
            % this is the one we're looking for
            name = i_name;
        end
        
        if DBG
            fprintf('%03d  = %4s\t%3d', code_num, i_name, phys_chan_ndx);
            marker = '';
            if phys_chan_ndx == physNum
                marker = ' <-- ';
            end
            fprintf('%-6s', marker);
            
            if ~isempty(montageJack) && phys_chan_ndx < length(montageJack)
                fprintf('%-8s', montageJack{phys_chan_ndx})
            else
                within_bank_ndx = mod(phys_chan_ndx - 1, 64) + 1; % +/- 1 ensures 0-->64
                bank_ndx        = floor((phys_chan_ndx-1) / 64) + 1;
                bank_name       = BANKS(bank_ndx);
                bank_guess = sprintf('%s%d(?)', bank_name, within_bank_ndx);
                fprintf('%-8s', bank_guess);
            end
            
            if ~isempty(montageLabels)
                if phys_chan_ndx < length(montageLabels)
                    fprintf('%-8s', montageLabels{phys_chan_ndx})
                end
            end
            
            fprintf('\n');
        end
    end
    
end

if DBG
    fprintf('(?) - denotes guess\n');
    fprintf('<-- - denotes phys2name returned channel\n\n');
    fprintf('%d physical channels found\n', phys_chan_ndx);
end











