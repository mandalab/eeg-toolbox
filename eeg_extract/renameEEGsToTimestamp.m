function renameEEGsToTimestamp(subj, rootEEGdir)
% .EEG reading to get time/date
% taken from nk_split

raw = fullfile(rootEEGdir, subj, 'raw');
d = dir([raw '/*.EEG']); 
cNames = {d.name};

for i = 1 : length(cNames)
    eegFile = fullfile(raw, cNames{i});
    timestamp = readEEG(eegFile);
    makeNewFile(subj, timestamp, raw, cNames{i});
end

end % main function

function makeNewFile(subj, timestamp, rawPath, eegFile)
    % designed to rename a folder in raw/ with timestamp
    newFilename = eegFile;
    d = dir(rawPath);
    cFiles = {d.name};
    
    prefix = eegFile(1 : strfind(eegFile,'.') - 1);
    cFilesToCpy = cFiles(~cellfun('isempty', strfind(cFiles, prefix)));
    
    for i = 1 : length(cFilesToCpy)
        filename = sprintf('%s_%s', subj, cFilesToCpy{i});
        newDir = fullfile(rawPath, timestamp);
        if ~exist(newDir, 'dir')
            mkdir(newDir);
        end
        from = fullfile(rawPath, filename);
        x = strfind(filename,'_');
        if length(x) == 3
            filename = filename(x(end)+1:end);
        end
            
        to = fullfile(newDir, filename);
        %fprintf('mv %s %s\n', from, to);
        movefile(from, to);
    end
end


function timestamp = readEEG(EEG_file)
    % taken from nk_split
    VERBOSE = 0;
    
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
    x=fread(fid,1,'*uint8');                if VERBOSE, fprintf('block ID: %d\n',x); end;
    x=fread(fid,16,'*char');                if VERBOSE, fprintf('device type: %s\n',x); end;   if strcmp(x(1:9)','EEG-1200A'),NEW_FORMAT=1;else NEW_FORMAT=0; end; 
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
    %strTime  = sprintf('%d/%d/%d %02d:%02d:%02d',T_month,T_day,T_year,T_hour,T_minute,T_second); % 
    %fprintf(' Date of session: %d/%d/%d\n',T_month,T_day,T_year)
    %fprintf(' Time at start: %02d:%02d:%02d\n',T_hour,T_minute,T_second)
    timestamp = sprintf('%02d%02d%02d_%02d%02d', T_year, T_month, T_day, T_hour, T_minute);
    fclose(fid);
end