function timestamps = eegRawServerGrabFiles( subj, rootEEGdir, rawFilePrefix, rawFileList, typeListOptional, nktdir, stimFlag, szFlag )
% eegRawServerGrabFiles( subj, rootEEGdir, rawFilePrefix, rawFileList, typeListOptional, nktdir, stimFlag, szFlag )
%
%   NOTE: ONLY WORKS IF YOU HAVE ACCESS TO and have mounted '/Volumes/Shares/EEG/...' -- confirm this first!
%
%   INPUTS:
%      -subject         -- 'NIH024'
%      -rootEEGdir      -- e.g. '/Volumes/Shares/FRNU/dataWorking/eeg'
%      -rawFileNaPrefix -- 'DA866'
%      -rawFileNameList -- cell array of strings that complete the prefix:  e.g.  {?JD?,?JS?,?JW?} will grab DA8662JD.21E DA8662JD.EEG DA8662JD.LOG,  and DA8662JS.21E...
%      -typeListOptional-- optional cell array of strings defining the file types to grab. If not passed , grab {'21E','EEG','LOG'}
%      -nktdir [optional]- directory of Dr. Inati's EEG server. If unpassed, default directories are searched.
%      -stimFlag [optional] -- if raw files should be placed in a STIM/ subdirectory
%      -szFlag [optional]   -- if raw files should be placed in a sz/ seizure subdirectory
%
%   OUTPUTS:
%       timestamps: list of timestamps grabbed
%
%   FILE OUTUT:
%       copies .21E, .EEG, and .LOG (or specified typeListOptional files) to [rootEEGdir]/[subj]/raw/[timestamp]/
%
%   example call:   eegRawServerGrabFiles('NIH025', '/Volumes/Shares/FRNU/dataWorking/eeg', 'DA8662', {'GH','GI','GY','H2','H6'}, {'21E','EEG','LOG'});
%

% sometimes this code hangs/crashes, probably due to slow interaction with the server.  Forcing all files closed seems to help
fclose all;
timestamps = {};

nkt_root = '/Volumes/Shares/EEG/LTVEEG_DATA'; % <-- Update this path if EEG is mounted under a different name in your machine
all_nkt_dirs = fullfile(nkt_root,...
    {'';
    'nkt/EEG2100';
    'Archive';
    'Archive2';
    'Archive3';
    'Archive4';
    'Archive5'}, 'NKT/EEG2100');

% -------------- Parameter handling ------------------- %
if exist('nktdir','var') && ~isempty(nktdir) && ~exist(nktdir, 'dir')
    % user specified
    error('Given nktdir not found: %s', nktdir);
    
elseif ~exist(nkt_root, 'dir')
    fprintf(['\n\nDefault nkt root directory not found....\n\tDid you mount EEG?\n\tDo you have permission to access?\n\t'...
        'You may need to update hard-coded nkt_root in eegRawServerGrabFiles\n\n']);
    error('Directory not found: %s', nktdir);
    
elseif ~exist('nktdir', 'var')
    nktdir = '';
end

if exist('typeListOptional','var') && ~isempty(typeListOptional)
    typeList = typeListOptional;
else
    typeList = {'21E','EEG','LOG'};  %- the file types that will be copied.  Log is only used for stimulation, but easier to get it for all of them
end

if  ~exist('stimFlag','var') || isempty(stimFlag)
    stimFlag = false; %- copies files to STIM folder, creates if it does not exist
end

if  ~exist('szFlag','var') || isempty(szFlag)
    szFlag = false; %- copies files to sz folder, creates if it does not exist
end

if ischar(rawFileList)
    rawFileList = {rawFileList};
end

if ischar(rawFileList)
    rawFileList = {rawFileList};
end


% -------------- end Parameter handling ------------------- %


tryCopy = 0;
numCopy = 0;
for iRaw = 1:length(rawFileList)
    
    % no matter the type, have to look at EEG file for timestamp
    eegSourceName = sprintf('%s%s.EEG', rawFilePrefix, rawFileList{iRaw});
    eegSource = fullfile(nktdir, eegSourceName);
    ndx_nkt_dir = 0;
    
    % loop through list of nkt dirs    
    while ~exist(eegSource, 'file')
        ndx_nkt_dir = ndx_nkt_dir + 1;
        assert(ndx_nkt_dir <= length(all_nkt_dirs), 'File %s not found in any nkt_dir', eegSourceName);
        nktdir = all_nkt_dirs{ndx_nkt_dir};
        eegSource = fullfile(nktdir, eegSourceName);
    end
    fprintf('%s found in %s\n', eegSourceName, nktdir);
    
    timestamp = readEEG(eegSource);
    timestamps = [timestamps {timestamp}];
    
    %%- destination directory: in subjects/raw/[timestamp]/
    grabbedSubDir = '_grabbed/';;
    %grabbedSubDir = '';
    if ~stimFlag && ~szFlag
        dest = sprintf('%s/%s/raw/%s%s/',      rootEEGdir, subj, grabbedSubDir, timestamp);
    elseif stimFlag
        dest = sprintf('%s/%s/raw/%sSTIM_MAP/%s/', rootEEGdir, subj, grabbedSubDir, timestamp);
    elseif szFlag
        dest = sprintf('%s/%s/raw/%ssz/%s/',   rootEEGdir, subj, grabbedSubDir, timestamp);
    end

    if ~exist(dest,'dir'), mkdir(dest); end
    if ~exist(dest,'dir'), error(' ERROR: destination directory cant be created ');  end
    
    for iType=1:length(typeList)
        
        sourceName = sprintf('%s%s.%s', rawFilePrefix, rawFileList{iRaw},typeList{iType});
        source = sprintf('%s/%s', nktdir, sourceName);

        fprintf('\n copying: %s ', sourceName); tic;
        [status,message,messageid] = copyfile(source,dest,'f');
        fprintf('%s [%.1f s]',message, toc);
        
        tryCopy = tryCopy+1;
        numCopy = numCopy+status;
        
    end
    fprintf('\n');
    
end

fprintf('\n  %d of %d files successfully copied to %s \n', numCopy, tryCopy, fileparts(dest)); 

end % function

function timestamp = readEEG(EEG_file)
    % The following code was taken from nk_split and reads EEG binary until
    % the timestamp is gotten
    VERBOSE = 0;
    
    fid = fopen(EEG_file);
    if fid < 0 % note: move this outside loop for every archive
        error('Cannot find file: %s. \nYou might try another directory.\n', EEG_file);
    end

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
    x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('block ID: %d\n',x); end;
    x=fread(fid,16,'*char');                if VERBOSE, fprintf('device type: %s\n',x); end;   if strcmp(x(1:9)','EEG-1200A'),NEW_FORMAT=1;else NEW_FORMAT=0; end; 
    x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('number of EEG2 control blocks: %d\n',x); end

    numberOfBlocks=x;
    if numberOfBlocks > 1
        % we think we will never have this
        % throw an error for now and re-write code if necessary
        fprintf('ERROR: %d EEG2 control blocks detected (only expecting 1).\n', numberOfBlocks);
        keyboard
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
    %strTime  = sprintf('%d/%d/%d %02d:%02d:%02d',T_month,T_day,T_year,T_hour,T_minute,T_second); % 
    %fprintf(' Date of session: %d/%d/%d\n',T_month,T_day,T_year)
    %fprintf(' Time at start: %02d:%02d:%02d\n',T_hour,T_minute,T_second)
    timestamp = sprintf('%02d%02d%02d_%02d%02d', T_year, T_month, T_day, T_hour, T_minute);
    fclose(fid);
end