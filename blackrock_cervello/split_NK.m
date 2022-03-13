function [allTags] = split_NK(subj, rootEEGdir, EEG_file, varargin)
% nk_split - Splits an nk .EEG datafile into separate channels into
% the specified directory.
%
% FUNCTION:
%    nk_split(subj,nk_dir,EEG_file,...)
%
% INPUT ARGs:
% subj = 'TJ022'
% rootEEGdir = '/Users/damerasr/damerasr/data/eeg/'
% EEG_file = '/Users/damerasr/damerasr/data/eeg/160614_1008/x.EEG'
%
% Output
%	output is saved out to: output_dir = '/Users/damerasr/damerasr/data/eeg/[SUBJ]/eeg.noreref/' set in line 49
%
%
% Edited from previous version so that manual input of the ouput_dir would not be necessary
%
% 12/2013 ... now uses jacksheetMaster.txt to guide channel number outputs
% 10/2015 ... now can handel "new" EEG-1200 extended file formats
% 08/2016 MST outputs using channel names instead of numbers
% 07/2017 MST Check for shift channels, modify .21E's


VERBOSE = 0 ; %[0,1] = output info about each channel's remapping


tagNameOrder = getTagNames(subj, rootEEGdir);


%%- Load the jacksheetMaster, or create it, or give a warning that it can't be created...
subjDir             = fullfile(rootEEGdir, subj);
rawDir              = fullfile(subjDir, 'raw');
[rawSubDir,~,~]     = fileparts(EEG_file); %- path to the raw file being processed, so FRNU/data/eeg/NIHXXX/raw/121401_1415/
jackMaster_file_new = fullfile(subjDir, 'docs/jacksheetMaster.csv');


% create jacksheetMaster if it doesn't exist or hasn't looked at raws
if ~exist(jackMaster_file_new, 'file')
    createMasterJack(subj, rootEEGdir);
end

jacktable = getJackTable(subj, rootEEGdir);
jackMaster_names = jacktable.chanName;
jackMaster_chans = [1:length(jackMaster_names)];


%- read the raw_info file to identify channels intentionally EXCLUDED from jacksheetmaster... no reason to give a warning if those are found below
raw_info_file = fullfile(subjDir, 'docs/raw_info.csv');
raw_info      = readtableSafe(raw_info_file);
chanDontSplit = raw_info.chanName(find(raw_info.in_jackSheetMaster==0));
chanDontSplit{end+1} = ''; %- dont report error or try to split channels with no name


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the EEG file based on the timestamp given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(exist(EEG_file, 'file') > 0);
[parentDir, filestem] = fileparts(EEG_file);

% Same as above for .21E file
ELEC_file = fullfile(parentDir, [filestem, '.21E']);
assert(exist(ELEC_file, 'file') > 0, '.21E file not found');

%- remove the patients name from the raw data (only happens when extracted)
nknih_anon(EEG_file);  % make sure the following works still... used to be in prepAndAlign
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies filepath for output directory and then checks to see if the
% directory already exists. If it doesn't exist it then creates it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_dir = fullfile(subjDir,'eeg.noreref');
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

%open and obtain *EEG* file information
fid = fopen(EEG_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skipping EEG device block
%%%%%%%%%%%%%%%%%%%%%%%%%%%
deviceBlockLen=128; %skips the first 128 bytes
fseek(fid,deviceBlockLen,'bof');  %fseek(fileID, offset, origin) moves to specified position in file. bof=beginning of file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading EEG1 control Block (contains names and addresses for EEG2 blocks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fread(fileID, sizeA, precision) reads data from a binary file.
%sizeA=output array size. uint8=unsigned integer with 8 bits.
%precision=string that specifies the form and size of the values to read
x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('block ID: %d\n',x); end; %#ok<*NASGU,*UNRCH>
x=fread(fid,16,'*char');                if VERBOSE, fprintf('device type: %s\n',x); end;   if strcmp(x(1:9)','EEG-1200A'), NEW_FORMAT=1; else NEW_FORMAT=0; end;
x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('number of EEG2 control blocks: %d\n',x); end

numberOfBlocks=x;
if numberOfBlocks > 1
    % we think we will never have this
    % throw an error for now and re-write code if necessary
    fprintf('ERROR: %d EEG2 control blocks detected (only expecting 1).\n', numberOfBlocks);
    return
end
% if numberOfBlocks is ever > 1, the following should be a for loop
blockAddress=fread(fid,1,'*int32');     if VERBOSE, fprintf('address of block %d: %d\n',i,blockAddress); end
x=fread(fid,16,'*char');                if VERBOSE, fprintf('name of EEG2 block: %s\n',x); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading EEG2m control block (contains names and addresses for waveform blocks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fseek(fid,blockAddress,'bof');          if VERBOSE, fprintf('\nin EEG21 block!\n'); end
x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('block ID: %d\n',x); end
x=fread(fid,16,'*char');                if VERBOSE, fprintf('data format: %s\n',x); end
numberOfBlocks=fread(fid,1,'*uint8');   if VERBOSE, fprintf('number of waveform blocks: %d\n',numberOfBlocks); end
if numberOfBlocks > 1,
    % we think we will never have this
    % throw an error for now and re-write code if necessary
    fprintf('ERROR: %d waveform blocks detected (only expecting 1).\n', numberOfBlocks);
    return
end
% if numberOfBlocks is ever > 1, the following should be a for loop
blockAddress=fread(fid,1,'*int32');     if VERBOSE, fprintf('address of block %d: %d\n',i,blockAddress); end
x=fread(fid,16,'*char');                if VERBOSE, fprintf('name of waveform block: %s\n',x); end


%%%%%%%%%%%%%%%%%%%%%%%
%Reading waveform block
%%%%%%%%%%%%%%%%%%%%%%%
fseek(fid,blockAddress,'bof'); %fprintf('\nin EEG waveform block!\n')
x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('block ID: %d\n',x); end
x=fread(fid,16,'*char');                if VERBOSE, fprintf('data format: %s\n',x); end
x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('data type: %d\n',x); end
L=fread(fid,1,'*uint8');                if VERBOSE, fprintf('byte length of one data: %d\n',L); end
M=fread(fid,1,'*uint8');                if VERBOSE, fprintf('mark/event flag: %d\n',M); end


%%- annonomous function to convert binary to decimal.  input is binary string created with dec2bin
bcdConverter2 = @(strDec2bin)  10*bin2dec(strDec2bin(1:4)) + bin2dec(strDec2bin(5:8));

% get the start time
T_year   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
T_month  = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
T_day    = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
T_hour   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
T_minute = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
T_second = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
strTime  = sprintf('%d/%d/%d %02d:%02d:%02d',T_month,T_day,T_year,T_hour,T_minute,T_second); %
fprintf(' Date of session: %d/%d/%d\n',T_month,T_day,T_year)
fprintf(' Time at start: %02d:%02d:%02d\n',T_hour,T_minute,T_second)
fileStemDate = sprintf('%02d%02d%02d_%02d%02d',T_year,T_month,T_day,T_hour,T_minute);     % new version: file stem of extracted channels (new version.... YYMMDD_HHMM -- JHW 11/2013


% get the sampling rate
x=fread(fid,1,'*uint16');  %fprintf('sample rate (coded): %d\n',x);
switch(x)
    case hex2dec('C064'),
        actSamplerate=100;
    case hex2dec('C0C8'),
        actSamplerate=200;
    case hex2dec('C1F4'),
        actSamplerate=500;
    case hex2dec('C3E8'),
        actSamplerate=1000;
    case hex2dec('C7D0'),
        actSamplerate=2000;
    case hex2dec('D388'),
        actSamplerate=5000;
    case hex2dec('E710'),
        actSamplerate=10000;
    otherwise
        fprintf('UNKNOWN SAMPLING RATE\n');
end

fprintf(' Sampling rate: %d Hz\n',actSamplerate);

% get the number of 100 msec block
num100msBlocks=fread(fid,1,'*uint32');      if NEW_FORMAT==0, fprintf(' Length of Session: %2.2f hours\n',double(num100msBlocks)/10/3600); end %- num100msBlocks=10 if new format
numSamples=actSamplerate*num100msBlocks/10; if VERBOSE, fprintf('number of samples: %d\n',numSamples); end
AD_off=fread(fid,1,'*int16');               if VERBOSE, fprintf('AD offset at 0 volt: %d\n',AD_off); end
AD_val=fread(fid,1,'*uint16');              if VERBOSE, fprintf('AD val for 1 division: %d\n',AD_val); end
bitLen=fread(fid,1,'*uint8');               if VERBOSE, fprintf('bit length of one sample: %d\n',x); end
comFlag=fread(fid,1,'*uint8');              if VERBOSE, fprintf('data compression: %d\n',x); end
numChannels=fread(fid,1,'*uint8');          if VERBOSE, fprintf('number of RAW recordings: %d\n',numChannels); end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ****** EXTENDED FORMAT (NEW) .EEG FILE                          ******
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (numChannels==1 && numSamples-actSamplerate==0) && NEW_FORMAT==0, fprintf('\n expecting old format, but 1 channel for 1 second'); keyboard; end
if (numChannels>1)                                && NEW_FORMAT==1, fprintf('\n expecting new format, but >1 channel ');           keyboard; end

if NEW_FORMAT || (numChannels==1 && numSamples-actSamplerate==0),
    
    if VERBOSE, fprintf('** New File Format **'); end
    
    %- seek the file location of the new wave data... need to make a pit stop in the new EEG2 header, which will provide the direct address
    waveformBlockOldFormat = 39 + 10 + 2*actSamplerate + double(M)*actSamplerate; %- with new format the initial waveform block contains 1 channel (10 bytes of info) for 1 second (2bytes x 1000)
    controlBlockEEG1new    = 1072;
    blockAddressEEG2       = blockAddress + waveformBlockOldFormat + controlBlockEEG1new;% + controlBlockEEG1new + controlBlockEEG2new;
    
    
    
    %- EEG2' format
    addTry = blockAddressEEG2;
    fseek(fid,addTry,'bof');                if VERBOSE, fprintf('--EEG2-prime format--\n'); end
    x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('block ID: %d\n',x); end
    x=fread(fid,16,'*char');                if VERBOSE==2, fprintf('data format: %s\n',x); end
    x=fread(fid,1,'*uint16');               if VERBOSE==2, fprintf('number of waveform blocks: %d\n',x); end;
    x=fread(fid,1,'*char');                 if VERBOSE==2, fprintf('reserved: %s\n',x); end
    x=fread(fid,1,'*int64');ii=1;           if VERBOSE==2, fprintf('address of block %d: %d\n',ii,x);    end; waveBlockNew = x;
    
    
    %- EEG2' waveform format
    fseek(fid,waveBlockNew,'bof');          if VERBOSE, fprintf('--EEG2-prime WAVE format--\n'); end
    x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('block ID: %d\n',x); end
    x=fread(fid,16,'*char');                if VERBOSE==2, fprintf('data format: %s\n',x); end
    x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('data type: %d\n',x); end
    L=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('byte length of one data: %d\n',L); end
    M=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('mark/event flag: %d\n',M); end
    
    %- now things get a little different with the new header
    x=fread(fid,20,'*char');                if VERBOSE==2, fprintf('start time string: %s\n',x); end
    x=fread(fid,1,'*uint32');               if VERBOSE==2, fprintf('data interval (sample rate): %d\n',x); end; actSamplerate  = double(x);
    x=fread(fid,1,'*uint64');               fprintf('Length of Session: %2.2f hours\n',double(x)/10/3600); num100msBlocks = double(x);
    
    numSamples  = actSamplerate*num100msBlocks/10; if VERBOSE==2, fprintf('number of samples: %d\n',numSamples); end
    AD_off      = fread(fid,1,'*int16');           if VERBOSE==2, fprintf('AD offset at 0 volt: %d\n',AD_off); end
    AD_val      = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('AD val for 1 division: %d\n',AD_val); end
    bitLen      = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('bit length of one sample: %d\n',bitLen); end
    comFlag     = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('data compression: %d\n',comFlag); end
    reserveL    = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('reserve length: %d\n',reserveL); end
    x           = fread(fid,reserveL,'*char');     if VERBOSE==2, fprintf('reserve data: %s\n',x); end
    
    numChannels = fread(fid,1,'*uint32');          if VERBOSE==2, fprintf('number of RAW recordings: %d\n',numChannels); end
    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ******.21E FILE******set the look-up tables to get the electrode names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fd_elecFile = fopen(ELEC_file, 'r');
elecCells = textscan(fd_elecFile,'%s%s','delimiter','='); % *MST textread -> textscan
[allCodes, allNames] = deal(elecCells{:});
% ndx_empty = cellfun('isempty', allNames);
% allCodes = allCodes(~ndx_empty);
% allNames = allNames(~ndx_empty);
fclose(fd_elecFile);

endRange  = find(strcmp(allCodes,'[SD_DEF]')); %finds the range of the electrodes
allCodes  = allCodes(1:endRange-1);  %stores the codes
allNames  = allNames(1:endRange-1);  %stores the names
info      = getElementInfo(subj, rootEEGdir);
syncTags  = info{strcmpi(info.chanType, 'SYNC'), 'tagName'};


goodCodes = [0:36 74 75 100:253];    % classic good codes... does not include DC
if ~isempty(tagNameOrder)
    if ismember('DC', syncTags)  % if DC is sync, then extract the DC channels
        goodCodes = [0:36 42:73 74:75 100:253];  %include DC channels, 42-73, but not mark channels 76-77  (JHW 10/2013)
    end
end

badNames  = {'E'};
actualName_ALL = {};
actualCode_ALL = {};


%the following for loop iterates through all the channels, which each
%stores a channel code. For each channel, it reads in the channel code,
%changes it to string format, and matches it to "allCodes" to find a match
%and stores that in matchingRow. It stores the matchingRow-th element in
%allNames into actualName, and checks to see if that is a good or bad
%electrode code. If good, it appends it to the list of good electrodes
%(actualName_ALL).

for k=1:numChannels
    x=fread(fid,1,'*int16');  %reads in 1 byte every time you iterate the loop
    chanCode(k)=x; %and stores it in chanCode(k)
    if (VERBOSE), fprintf(' Index %d ''name'': Channel %d\n',k,x); end
    chanCodeString=sprintf('%04d',x); %format data into string. Same as chanCode except the format is string.
    matchingRow=find(strcmp(chanCodeString,allCodes)); %looks for this particular string in allCodes and stores its locations in matchingRow
    actualName = allNames{matchingRow};
    
    if ~ismember(chanCode(k),goodCodes) %if not a member of goodCodes
        if VERBOSE, fprintf(' chan %d (%s) is a bad channel code and excluded\n',chanCode(k),actualName); end
        goodElec(k)=false;
        badElec(k) = true;
    elseif any(strcmp(actualName,badNames)) %or if it's part of badNames
        %fprintf(' chan %d (%s) is a bad address\n',chanCode(k),actualName);
        goodElec(k)=false;
        badElec(k) = true;
    else
        if (VERBOSE), fprintf(' chan %d (%s) is good!\n',chanCode(k),actualName); end
        goodElec(k)=true;
        badElec(k) = false;
    end
    
    % save out the names for the jacksheet
    if goodElec(k)
        %if it is a good electrode, append it to the jacksheet
        actualName_ALL(end+1)=allNames(matchingRow);
        actualCode_ALL(end+1)=allCodes(matchingRow);
        
        if isempty(actualName_ALL{end}) & (VERBOSE | ~any(strcmp(chanDontSplit,''))), %
            fprintf(' WARNING: Channel code %s maps to empty name (index %d)\n', chanCodeString, matchingRow);
        end
    end
    
    
    fseek(fid,6,'cof'); %skipping the six most sig. bits of 'name'
    
    %finds the channel sensitivity
    chan_sensitivity=fread(fid,1,'*uint8');
    if (VERBOSE) fprintf('channel sensitivity: %d\n',chan_sensitivity); end
    switch fread(fid,1,'*uint8');  %fprintf('         unit: %d\n',chan_unit);
        case 0; CAL=1000;%microvolt
        case 1; CAL=2;%microvolt
        case 2; CAL=5;%microvolt
        case 3; CAL=10;%microvolt
        case 4; CAL=20;%microvolt
        case 5; CAL=50;%microvolt
        case 6; CAL=100;%microvolt
        case 7; CAL=200;%microvolt
        case 8; CAL=500;%microvolt
        case 9; CAL=1000;%microvolt
    end
    GAIN(k)=CAL/double(AD_val); %OK TO ASSUME THIS IS CONSTANT FOR ALL ELECS?
end  % for for k=1:numChannels

% ADD THIS CODE TO CLEAN UP EEG ERROR WITH WHITESPACE AFTER SOME TAGS (NIH018, channel "LF2 ") -- JW 10/2013
actualName_ALL = strtrim(actualName_ALL);


% CHECKS TO SEE IF USER HAS SPECIFIED TAG NAME ORDER OTHERWISE PROGRAM ENDS AND TAG NAME ORDER IS RETURNED
if isempty(tagNameOrder)
    allTags = actualName_ALL;
else
    allTags = {};
    %the function 'assert' generates error when condition is violated. 'unique'
    %returns an array with the same value but no repetition
    assert(length(unique(GAIN))==1,'All channels do not have the same gain!');
    
    
    %%%%%%%%%%%%%%
    % get the data
    %%%%%%%%%%%%%%
    fprintf('Reading Data...');
    
    d=fread(fid,[double(numChannels+1) double(numSamples)],'*uint16'); %reads the content into an array
    if VERBOSE, fprintf(' done\n'); end
    
    %- pull the trigger bits out of the event/mark entry (JHW 10/2013)
    dEvntMrk = d((numChannels+1),:);  %additional element in time series (chan+1) is 16-bit event/mark data, where bits 7-14 encode DC09-DC13 triggers
    trigDC09 = bitget(dEvntMrk,7);    %trigDC09 encodes the aligment pulses
    trigDC10 = bitget(dEvntMrk,8);    %trigDC10 encodes stimulation pulses
    trigDC11 = bitget(dEvntMrk,9);
    trigDC12 = bitget(dEvntMrk,10);
    trigDC13 = bitget(dEvntMrk,11);
    trigDC14 = bitget(dEvntMrk,12);
    trigDC15 = bitget(dEvntMrk,13);
    trigDC16 = bitget(dEvntMrk,14);
    
    d=d([goodElec false],:);  %removing bad electrodes (including 16-bit event/mark data)
    
    fprintf('Removing offset...');
    d_int16=int16(int32(d)+int32(AD_off)); %convert to int16
    
    %the line below proves the above is lossless... it eats up ram, so only do it on a beefy mac
    if ismac,
        fprintf('Validating conversion...');
        assert(isequal(d,uint16(double(d_int16)-double(AD_off))))
        fprintf('.');
    end
    
    % scale the DC input lines (different scaling than EEG signals)  (JHW 10/2013)
    iDC = find(~cellfun('isempty',strfind(actualName_ALL,'DC')));
    if isempty(iDC) & VERBOSE, fprintf('\nWARNING: DC chan not found in nk_split, probably need to do manual pulse extraction from EKG channels'); end
    
    %GAIN_DC = 500 / 2730 ; % E11FFmt.pdf says "Ox0555 corresponds to 500 mV"; Ox555 = 1365;  looks like actually OxAAA-->500mV (2730)... should confirm on multiple machines though
    GAIN_DC = 500 / 1365 ; % E11FFmt.pdf says "Ox0555 corresponds to 500 mV"; Ox555 = 1365;  looks like actually OxAAA-->500mV (2730)... should confirm on multiple machines though
    for thisDC = iDC,
        d_int16(thisDC,:) = d_int16(thisDC,:)*GAIN_DC;  %possibly this should not be scaled here... could loose data with the int16 transform
    end
    
    % clear the temp variable and close the file
    clear d
    fprintf(' done\n');
    fprintf(' Number of channels to split from raw file: %d\n',size(d_int16,1));
    fclose(fid);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-- Check for a shift channel --%  %% this is now down in createRawInfo
    % actualName_ALL is updated to _sh along with .21E if user specifies to do so
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %shift_mark = '_sh';
    %if any(strfound(jackMaster_names, shift_mark))
    %    actualName_ALL = updateShiftChannels(ELEC_file, info, jackMaster_names, actualName_ALL);
    %end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reorder the jacksheet and the electrodes... two options: sort using tagNames, or sort using jacksheetMaster.  Always use master if possible
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    chanOut = nan(1, size(d_int16,1));  % chanOut is same length as d_int16... maps raw channel order to jackMaster_chans
    newElectrodeNames = {};             % newElectrodeNames is same size as jackMaster_names... modified to indicate whether some channels are missing
    newElectrodeChans = [];
    
    %- loop through all jackMaster channels... if jackMaster name is found, map the channel,
    %  if not found, modify the channel name for this session's jacksheet.
    %  if channel is shifted based on localization pull/rearrangement, modify channel name
    for iChan = 1:length(jackMaster_chans),
        
        masterChan = jackMaster_chans(iChan);
        masterName = jackMaster_names{iChan};
        
        actualChan = find(strcmp(actualName_ALL,masterName));
        
        if length(actualChan) == 1, % Good case
            fail = false;
            
        elseif isempty(actualChan), % missing
            fail = true;
            
            if VERBOSE,
                fprintf('\n  WARNING: raw channel list missing %s, which is found in jacksheetMaster.txt.  Channel will be skipped for this eeg file stem.\n', masterName, masterChan);
            end
            
            newElectrodeNames{end+1} = sprintf('<%s_missing_from_raw>', masterName);
            newElectrodeChans(end+1) = masterChan;
            
        else % duplicates
            fail = true; % (but we might correct it!)
            fprintf('\n  You should not be able to get here anymore.  Call Wittig if you are here');
            keyboard;
            
            %fprintf('\n  ERROR: more than 1 instance of %s found in raw channel list...\n  Press "Y" to examine/correct the discrepancy: ', masterName);
            %if strcmpi(input('', 's'), 'Y')
            %    % user wants to correct
            %    [corrected, actualChan] = correctRepeatedChannel(actualChan, actualCode_ALL, ELEC_file, masterName);
            %    if corrected
            %        fail = false;
            %    end
            %end
            
            if fail
                error('More than one instance of %s found in raw channel list', masterName);
            end
        end
        
        if ~fail
            chanOut(actualChan)      = masterChan;
            newElectrodeNames{end+1} = masterName;
            newElectrodeChans(end+1) = masterChan;
        end
        
    end % channel loop
    
    %- any raw channels not accounted for?
    notInMaster = find(isnan(chanOut));
    if ~isempty(notInMaster),
        
        for iChan = 1:length(notInMaster),
            rawName = actualName_ALL{notInMaster(iChan)};
            if any(strcmp(chanDontSplit, rawName))==0,
                fprintf(['\n  WARNING: raw channel list contains %s, but jacksheetMaster.csv does not.'...
                    '\n Channel will not be extracted.  Correct 21E, or add it to element_info.csv to split, '...
                    '\n or supress in element_info with recNotUsed column and delete jacksheetMaster' ...
                    '\n if just added a new raw delete jackSHeetMaster and re-run prepAndAlign so the channel can be fixed'], rawName) ;
                fprintf('\n NOTE: probably shouldnt be possible to get here. so call wittig if you get here');
                keyboard;
            end
            newElectrodeNames{end+1} = sprintf('<%s_missing_from_jacksheetMaster.txt>', rawName);
            newElectrodeChans(end+1) = nan;
            
        end
        
        %- cut the raw channels that were not specified in the master
        inMaster = find(~isnan(chanOut));
        chanOut  = chanOut(inMaster);
        d_int16  = d_int16(inMaster,:);
    end
    
    %- sort so output in order (makes directory sorted by date look same as sorted by name)
    [chanOut_sorted sortIdx] = sort(chanOut);
    chanOut = chanOut(sortIdx);
    d_int16 = d_int16(sortIdx,:);
    
    
    %%%%%%%%%%%%%%%%%%%%
    % make the jacksheet
    %%%%%%%%%%%%%%%%%%%%
    %keyboard
    fileRoot = sprintf('%s',fileStemDate);     % new version uses "fileStemDate".... YYMMDD_HHMM... defined above where date is extracted
    fileDir = fullfile(output_dir,fileRoot);
    if ~exist(fileDir, 'dir'), mkdir(fileDir); end
    
    %- make the raw file specific jacksheet (with file stem)
    jackFile = fullfile(fileDir,'jacksheet.txt');
    fout3    = fopen(jackFile,'w','l');
    for iChan = 1:length(newElectrodeChans)
        thisChan = newElectrodeChans(iChan);
        thisName = newElectrodeNames{iChan};
        fprintf(fout3,'%d %s\n',thisChan,thisName);  %JHW... mod to all for jacksheetMaster implementation
    end
    fclose(fout3);
    
    %- check for global jacksheet (if no jacksheetMaster)
    if isempty(jackMaster_chans),
        fprintf('\n you should never get here.  call wittig.  seriously');
        keyboard;
    end
    pause(.5);%JFB: pause to smooth output. I am slow and so I like slow output!
    % the real JFB: above comment editorialized by JJ
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the params.txt file
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    paramsFile = fullfile(fileDir, 'params.txt');
    fout4=fopen(paramsFile,'w','l');
    fprintf(fout4,'samplerate %d\n',actSamplerate);
    fprintf(fout4,'dataformat ''int16''\n');
    fprintf(fout4,'gain %d\n',GAIN(1));
    fclose(fout4);
    % copy over the most params.txt file
    if 0 %- dont make params.txt at root level anymore because can be different for different sessions
        if ismac,
            system(sprintf('cp "%s" "%s"',paramsFile,fullfile(output_dir,'params.txt')));
            fprintf(' params.txt is made\n')
        else
            [success, ~, ~] = copyfile(paramsFile,fullfile(output_dir,'params.txt'));
            if success==1
                fprintf(' params.txt is made\n')
            else
                fprintf('\n WARNING: running nk_split on PC instead of MAC... params.txt not copied sucessfully');
            end
        end
        filename = fullfile(output_dir, 'params.txt');
        if ~exist(filename, 'file')
            copyfile(paramsFile, filename);
        end
    end
    
    
    % make sourcetype.txt (which system recorded raw)
    sourceFile=fullfile(fileDir,'sourcetype.txt');
    fileOut = fopen(sourceFile,'w','l');
    fprintf(fileOut,'nihon khoden raw\n');
    fclose(fileOut);


    pause(.5);%JFB: pause to smooth output
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write the electrodes to file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pause(.5);%JFB: pause to smooth output
    fprintf('\nwriting files:')
    ticker=0;
    tick_inc=10;
    for thisChan=chanOut
        if thisChan/length(chanOut)*100>=ticker
            fprintf(' %2.0f%%',ticker)
            ticker=ticker+tick_inc;
        end
        thisName = newElectrodeNames{find(newElectrodeChans == thisChan, 1)};
        chanfile = fullfile(fileDir, thisName);
        [fchan,msg] = fopen(chanfile,'w','l');
        assert(fchan > 0, 'Could not open file %s for writing (error: %s)', chanfile, msg);
        fwrite(fchan,d_int16(find(chanOut==thisChan),:),'int16');
        fclose(fchan);
        if ismac
            try
                % JHW - change files to "executable"... helps for sorting in mac finder
                % MST - avoid permission error (only owner can change permission)
                fileattrib(chanfile, '+x', 'a');
            catch
            end
        end
    end
    
    
    %%- DIGIAL SYNC CHANNELS... OUTPUT THE TRIGGERED VERSIONS OF DC09-DC12
    %- output a sync for DC10 as well, this will be useful for newer stimMapping sessions
    numTrigDC09 = length(find([trigDC09(1) diff(trigDC09)]==1));
    numTrigDC10 = length(find([trigDC10(1) diff(trigDC10)]==1));
    numTrigDC11 = length(find([trigDC11(1) diff(trigDC11)]==1));
    numTrigDC12 = length(find([trigDC12(1) diff(trigDC12)]==1));
    
    syncMade = 0;
    if numTrigDC09>10 & find(trigDC09>0) %JHW - output sync file if trigger bit is non-zero 11/2013
        updowns = [trigDC09(1) diff(trigDC09)];
        uptimes = find(updowns==1);
        chanfile = fullfileEEG(fileDir, 'trigDC09.sync.txt');
        fchan = fopen(chanfile,'w','l');
        fprintf(fchan,'%d\n',uptimes);
        fclose(fchan);
        fprintf('\nDigital pulse trigger found... extracted %d pulse up times to:\n   %s', length(uptimes), chanfile); syncMade=1;
    end
    
    if length(find(trigDC10>0))>0 & numTrigDC10>10, %JHW - only output DC10 sync if >20 pulses (really should have >200, so lower bound)
        updowns = [trigDC10(1) diff(trigDC10)];
        uptimes = find(updowns==1);
        chanfile = fullfile(fileDir, 'trigDC10.syncStim.txt');
        fchan = fopen(chanfile,'w','l');
        fprintf(fchan,'%d\n',uptimes);
        fclose(fchan);
        fprintf('\n Digital pulse trigger found for DC10 (stim)... extracted %d pulse up times to:\n   %s', length(uptimes), chanfile); syncMade=1;
    end
    
    if length(find(trigDC11>0))>0 & numTrigDC11>10, %JHW - only output if >20 pulses (really should have >200, so lower bound)
        updowns = [trigDC11(1) diff(trigDC11)];
        uptimes = find(updowns==1);
        chanfile = fullfile(fileDir, 'trigDC11.syncStim.txt');
        fchan = fopen(chanfile,'w','l');
        fprintf(fchan,'%d\n',uptimes);
        fclose(fchan);
        fprintf('\n Digital pulse trigger found for DC11 (stim)... extracted %d pulse up times to:\n   %s', length(uptimes), chanfile); syncMade=1;
    end
    
    if length(find(trigDC12>0))>0 & numTrigDC12>10, %JHW - only output if >20 pulses (really should have >200, so lower bound)
        updowns = [trigDC12(1) diff(trigDC12)];
        uptimes = find(updowns==1);
        chanfile = fullfile(fileDir, 'trigDC12.syncBR.txt');
        fchan = fopen(chanfile,'w','l');
        fprintf(fchan,'%d\n',uptimes);
        fclose(fchan);
        fprintf('\n Digital pulse trigger found for DC12 (Blackrock)... extracted %d pulse up times to:\n   %s', length(uptimes), chanfile); %syncMade=1; %this shouldn't count as a sync
    end
    
    %- OPTIONAL: create a more complete output for the DC channels, showing both up and down times.
    %            for now only do this for DC12 because the duration of the up period from the blackrock might be important for aligning during non-task periods
    %
    %- loop through trig arrays and see whether they have any non-zero entries... if so, export a up-down file
    for trigOut = [12],  % used to do this for 10, 11, 12... now just make sync files for those
        trigStr  = sprintf('DC%02d',trigOut); %-create variable representing trigDC09, 10, 11, etc
        thisTrig = eval(sprintf('trig%s',trigStr)); %-create variable representing trigDC09, 10, 11, etc
        
        %- additional outputs if DC10, 11, or 12 had trigger events: this is used to determine stimulation timining 2/2014
        if sum(thisTrig)>0,
            trigDCout   = double(thisTrig);
            
            %-create list of pulse start (up) and stop (down) times (in units of sample), and a string for each event
            strUpDown   = {sprintf('PULSE_HI \t %s',trigStr),sprintf('PULSE_LO \t %s',trigStr)};
            updowns     = [trigDCout(1) diff(trigDCout)];  %diff requires double input for proper functionality
            updownTimes = find(updowns==1  |  updowns==-1);
            updownStr   = {strUpDown{ ((updowns(updownTimes)-1)*-.5)+1 }} ; %-convert -1-->2 and 1->1
            
            %-output updown file (annotation not incorporated here, only in versuion used to make behavioral/stimMapAnn folder
            chanfile    = fullfile(fileDir,sprintf('trig%s.updown.txt', trigStr));;
            fchan       = fopen(chanfile,'w','l');
            fprintf(fchan,'1 \t FILENAME \t %s \n1 \t EEGSTEM \t %s \n', chanfile(min(strfind(chanfile, fileRoot)):end),fileRoot); %- make first entry "FILENAME"; trim off the path info
            for iOut = 1:length(updownTimes),
                fprintf(fchan,'%d \t %s\n',updownTimes(iOut), updownStr{iOut});
            end
            fprintf(fchan,'%d \t SESS_END \n', size(d_int16,2)); %- make first entry "FILENAME"; trim off the path info
            fclose(fchan);
            if ismac, fileattrib(chanfile, '+x', 'a'); end %JW - change files to "executable"... helps for sorting in mac finder
            fprintf('\nDigital trigger events found on %s... extracted %d pulse up and down times to:\n   %s', trigStr, length(updownTimes), chanfile);
        end
    end
    
    
    
    %%- MANUAL SYNCS (old subjects, <NIH020) were previously processed and are now saved in the raw folder so they can be automatically copied to eeg.noreref
    %- if no tripDC09, then check to see if sync file exists in the raw directory... if so, make a copy in the noreref split folder
    oldSync  = dir(fullfileEEG(rawSubDir,'*.sync.txt'));
    for iSync=1:length(oldSync),
        if ~strcmp(oldSync(iSync).name(1:3),'NIH') & ~strcmp(oldSync(iSync).name(1:2),'TJ') & ~strcmp(oldSync(iSync).name(1:2),'UP'),
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfileEEG(rawSubDir,oldSync(iSync).name), fullfileEEG(fileDir, oldSync(iSync).name), 'f');
            if SUCCESS, fprintf('\n Manually generated SYNC file found in raw folder and copied to: %s',fullfileEEG(fileDir, oldSync(iSync).name));  syncMade=1;
            else        fprintf('\n Uh Oh.'); keyboard; end
        elseif length(oldSync)==1,
            fprintf('\n ERROR: manually generated sync file present, but not named correctly. Use jacksheet to rename.  \n Ex) NIHXXX_XXXX_XXX.083.084.sync.txt --> EKG1.EKG2.sync.txt');
            keyboard;
        end
    end
    
    
    %%- ANOTATION FILES
    %-  If annotation file exists and contains any notes from clinicians save it to the split folder
    annStruct = nk_parseAnnotation(rawDir);
    annTimes = [];
    annStr   = {};
    for iAnn=1:length(annStruct)
        annTimes(iAnn) = (annStruct(iAnn).timeSec+1)*actSamplerate;  %-add 1 sample to avoid indexing 0 (could cause problem for end of file?)
        annStr{iAnn}   = sprintf('ANNOTATE \t %s',annStruct(iAnn).str);
    end
    %-output annotation file to split directory
    if ~isempty(annTimes)
        chanfile    = fullfile(fileDir, 'annotation.txt');
        fchan       = fopen(chanfile,'w','l');
        fprintf(fchan,'1 \t FILENAME \t %s \n1 \t EEGSTEM \t %s \n', chanfile(min(strfind(chanfile, fileRoot)):end),fileRoot); %- make first entry "FILENAME"; trim off the path info
        for iOut = 1:length(annTimes),
            fprintf(fchan,'%d \t %s\n',annTimes(iOut), annStr{iOut});
        end
        fclose(fchan);
        fprintf('\n Annotation entries found (%d total) and extracted to: %s', length(annTimes), chanfile);
    end
    
    % JC ignore the bit about manualy finding pulses if this function is being called from 'eeg_BL_SZ_extract', since there are no pulses for these sessions
    parent_funcs=dbstack;
    parent_funcs={parent_funcs.name};
    batch_funcs={'eeg_BL_SZ_extract'};
    bl_sz_flag = any(ismember(batch_funcs,parent_funcs));
    
    if syncMade==0 && ~bl_sz_flag
        fprintf('\n WARNING:  NO SYNC FILE MADE... might have to manualy find pulses for this raw');
        keyboard
    end
    
    
    fprintf('\nExtraction complete\n\n\n')
    
end % end tagNameOrder empty conditiona


end

