function [EEG resampleFreq] = gete_ms(channel,events,DurationMS,OffsetMS,BufferMS,filtfreq,filttype,filtorder,resampleFreq,RelativeMS)
%GETE_MS - Get EEG event data based on MSec ranges instead of samples.
%
% Returns data from an eeg file.  User specifies the channel,
% duration, and offset along with an event.  The event struct MUST
% contain both 'eegfile' and 'eegoffset' members.
%
% You can optionally resample the data with the Signal Toolbox's
% resample function.  The resampling occurs following the filter.
%
% The distinctive feature of this function is that it can handle
% files that have different sampling rates. If this occurs, then it
% will resample all EEG files to have the same sampling rate.  If
% you specify a sampling rate, then it uses that one. Otherwise, it
% takes the sampling rate from the first eegfile that it looks at.
%
% FUNCTION:
%   [EEG resampleFreq] = gete_ms(channel,events,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampleFreq,RelativeMS)
%
% INPUT ARGS:

%   channel = 'G28';        % channel name (string) or cell string (for bipolar)
%   events = events(8:12);  % event struct to extract [eegfile eegoffset]
%   durationMS = 2000;      % signal time length in milliseconds
%   offsetMS = 0;           % offset at which to start in milliseconds 
%   bufferMS = 1000;        % buffer (needed for filtering or resampling)
%                           %   default is 0
%   filtfreq = [58 62];     % Filter freq (depends on type, see buttfilt)
%                           %   default is []
%   filttype = 'stop';      % Filter type (see buttfilt)
%   filtorder = 1;          % Filter order (see buttfilt)
%   resampleFreq = 200;     % Sample rate in Hz of the returned data
%   RelativeMS = [-200 0];  % Range for use with the relative subtraction
%
% OUTPUT ARGS:
%   EEG(Trials,Time) - The data from the file
%   resampleFreq     - The sampling rate used

%CHANGE LOG:
% 10/21/19 -  SJ   -  convert int16max to NaN
% 9/01/11  -  JFB  -  made the sample rate a second output
% 7/27/11  -  EH   -  completely rewrote function for performance
% 9/27/10  -  JRM  -  re-wrote gete_ms as a wrapper for gete.m.
% 7/2/10   -  JRM  -  return nans when eegoffset for an event occurs after the
%                     end of the file
% 12/18/07 -  MvV  -  changed the indices into readbytes when saving
%                     to EEG, such that it was always fit and not be affected by
%                     rounding differences.
% 11/29/04 -  PBS  -  Changed round to fix to fix range problem
% 4/20/04  -  PBS  -  Added Relative Range subtraction

%events(1).eegoffset = 1; %SJ- to go to beginning of file (test)

% check the arg
if ~exist('OffsetMS','var') || isempty(OffsetMS)
    OffsetMS = 0;
end
if ~exist('BufferMS','var') || isempty(BufferMS)
    BufferMS = 0;
end
if ~exist('filtfreq','var')
    filtfreq = [];
end
if ~exist('filttype','var') || isempty(filttype)
    filttype = 'stop';
end
if ~exist('filtorder','var') || isempty(filtorder)
    filtorder = 1;
end
if ~exist('resampleFreq','var')
    resampleFreq = [];
end
if ~exist('RelativeMS','var')
    RelativeMS = [];
end

% dereference char cells
if iscellstr(channel) && isscalar(channel)
    channel = char(channel);
end

% global IS_MEF
% if isempty(IS_MEF)
%     IS_MEF = false;
% end
IS_MEF = false;

fileList = {events.eegfile};
fileList(~cellfun( @ischar, fileList )) = {''}; % set any [] filenames to '' for unique to work
fileSet = unique(fileList);

if length(fileSet) <= 1
    if isempty(fileSet) || isempty(fileSet{1})
        EEG = nan(length(events), 0);
        return;
    end
elseif isempty(fileSet{1})
    fileSet = fileSet(2:end);
end

% store the butterowrth filter to save computation
bFilterData = ~isempty(filtfreq);
butterFilt = {};
butterFiltSamplerate = -1;

% store the resample filter to save computation
resampleFilt = [];
resampleFiltSamplerate = -1;
resampleExcess = -1;    % extra samples created during resampling
  
if isempty(resampleFreq)
    % if resampleFreq is not specified, always re-sample to the first
    % session's sampling rate to ensure the same number of samples for all
    % events, across sessions
    resampleFreq = round(GetRateAndFormat( fileSet{1} ));
else
    resampleFreq = round(resampleFreq);
end
        
% base final datasize on resampled data
buffDur  = ceil( double((DurationMS+2*BufferMS))*resampleFreq/1000 );
buffer   = fix( (BufferMS)*resampleFreq/1000 );
duration = buffDur - 2*buffer;

% allocate space for EEG data
EEG = nan(length(events),duration);

% loop through the files
for s = 1:length(fileSet)
    
    fileMask = strcmp(fileSet{s}, fileList);
    
    [sampleFreq,nBytes,dataformat,gain] = GetRateAndFormat( fileSet{s} );
    sampleFreq = round(sampleFreq);
    
    if isempty(sampleFreq)
        % cannot find params.txt file, cannot determine sampling rate, so error
        error('EEGTOOLBOX:GETE_MS:NOSAMPLERATE','The sample rate for file %s could not be found',fileSet{s});
    end
    
    %compute input params for gete based on next session's sampling rate
    readDuration = ceil( (DurationMS+2*BufferMS)*sampleFreq/1000 );
    readOffset   = round( (OffsetMS-BufferMS)*sampleFreq/1000 );
    
    %%% open file and read in EEG
    if iscellstr(channel) && numel(channel) == 2
        % Bipolar
        eegfname = fullfile(fileSet{s}, sprintf('%s-%s', channel{:}));
    else
        % Monopolar
       eegfname = fullfile(fileSet{s}, channel);
       
       %- JW needed support for old format
       if ~exist(eegfname,'file'),
           oldStyleName = sprintf('%s.%03d',fileSet{s}, channel);
           %- only do the rename if rename actually exists
           if exist(oldStyleName,'file'), 
               eegfname = oldStyleName;
           end
       end
    end
    
    startIndices = [events(fileMask).eegoffset] + readOffset;
    
    fileEEG = nan(length(startIndices),max(readDuration,buffDur));
    %fileEEG = int16(zeros(length(startIndices),max(readDuration,buffDur))); %SJ test

    if IS_MEF
        [readEEG,fileEEGlen] = decomp_mef_events(eegfname,startIndices,readDuration,'');
        fileEEG(:,1:readDuration) = readEEG';
        % EH - need to add code (either here or to decomp_mef_events) to
        % check the 'unpadded lead' filename (see a few lines below)
    else
        fileEEGlen = zeros(length(startIndices),1,'uint32');
        
        eegfile = fopen(eegfname,'r','l'); % NOTE: the 'l' means that it came from a PC!

        % if the file still can't be opened, throw an error
        if eegfile == -1                
            error('Missing EEG file: %s\n',eegfname);
        end
        
        for i = 1:length(startIndices)
            status = fseek(eegfile,nBytes*startIndices(i),-1);
            if status == 0
                readbytes = fread(eegfile,readDuration,dataformat)';
                %readbytes = fread(eegfile,readDuration,'int16=>int16')';%dataformat)'; %SJ to keep int16
                fileEEG(i,1:length(readbytes)) = double(readbytes); %  fread converts int16 to double already, but ok
                %fileEEG(i,1:length(readbytes)) = readbytes; %SJ to keep int16
                fileEEGlen(i) = length(readbytes);
            end
        end
        
        fclose(eegfile);
    end
    
    max16 = max(fileEEG(:)); %SJ
    
    % Check for max int16 values (32767) which were assigned before writing
    % processed data as int16 -SJ
    if max16 > intmax('int16') %SJ
        fprintf('\t%s\n',['ERROR!!! Somehow maximum EEG value (' num2str(max16) ') is greater than max int16 (32767). Error when converting to double?? Ask SJ!']);
        keyboard;
    elseif max16 == intmax('int16') % Change max int16 values back to double NaN
        fprintf('\t%s\n',['Stored EEG data contains values at int16 max. Max EEG value = ' num2str(max16) '; Converting int16 max to NaN ... '])
        fileEEG_orig = fileEEG;
        fileEEG(fileEEG==intmax('int16'))= NaN;
        if ~(sum(fileEEG_orig == intmax('int16'),'all') == sum(isnan(fileEEG),'all'))
            fprintf('\t%s\n',['ERROR!!! ' num2str(sum(fileEEG_orig == intmax('int16'),'all')) ' instances of int16 max found, but ' num2str(sum(isnan(fileEEG),'all')) ' replaced with NaN (Should be equal!) Ask SJ!']);
        end
    else
        %fprintf('\t%s\n',['Max stored EEG value = ' num2str(max16) '. No int16 max found, no conversion necessary.']);
    end
    % Create figure
%     perm_fileEEG = permute(fileEEG_orig, [2,1]);
%     yvals = perm_fileEEG(:);
%     figure; plot(yvals);
%     line(xlim,[intmax('int16'),intmax('int16')],'Color','r','LineStyle','--');
%     title(['Session ' num2str(events(1).sessionNum) ', Channel ' channel ', all events in order']);
    
    noEEGmask   = (fileEEGlen == 0);
    fullEEGmask = (fileEEGlen == readDuration);
    shortEEGind = find( ~(fullEEGmask | noEEGmask) );
    
    if any(noEEGmask)
        warning('EEGTOOLBOX:GETE_MS:NODATA','%s: EEG data were not found for %d event(s) -- setting NaN',...
                eegfname,sum(noEEGmask));
    end
    
    if ~isempty(shortEEGind)
        warning('EEGTOOLBOX:GETE_MS:INCOMPLETEDATA','%s: not all samples read for %d event(s) -- appending NaN',...
                eegfname,length(shortEEGind));
    end
    
    if IS_MEF
        % set unread EEG data to NaN
        for i = 1:length(shortEEGind)
            ind = shortEEGind(i);
            fileEEG(ind,fileEEGlen(ind)+1:end) = NaN;
        end
        fileEEG(noEEGmask,:) = NaN;
    end

    if bFilterData
        if( butterFiltSamplerate ~= sampleFreq )
            [fileEEG(fullEEGmask,1:readDuration),butterFilt] = buttfilt(fileEEG(fullEEGmask,1:readDuration),filtfreq,sampleFreq,filttype,filtorder);
            butterFiltSamplerate = sampleFreq;
        else
            fileEEG(fullEEGmask,1:readDuration) = buttfilt(fileEEG(fullEEGmask,1:readDuration),butterFilt);
        end
        
        for i = 1:length(shortEEGind)
            ind = shortEEGind(i);
            % no need to check that butterFilt is the correct sampling rate, since that was just done above
            fileEEG(ind,1:fileEEGlen(ind)) = buttfilt(fileEEG(ind,1:fileEEGlen(ind)),butterFilt);
        end
    end
    
    if( resampleFreq ~= sampleFreq )
        if( resampleFiltSamplerate ~= sampleFreq )
            [resamplebytes,resampleFilt] = resample(fileEEG(fullEEGmask,1:readDuration)',resampleFreq,sampleFreq);
            resampleFiltSamplerate = sampleFreq;
            resampleExcess = fix( (size(resamplebytes,1) - buffDur)/2 );
        else
            resamplebytes = resample(fileEEG(fullEEGmask,1:readDuration)',resampleFreq,sampleFreq,resampleFilt);
        end
        
        % often, resamplebytes will be longer than buffDur
        % b/c resample uses: length(output) = ceil(length(input)*resampleFreq/sampleFreq)
        fileEEG(fullEEGmask,1:buffDur) = resamplebytes(resampleExcess+1:resampleExcess+buffDur,:)';

        for i = 1:length(shortEEGind)
            ind = shortEEGind(i);
            % no need to check that resampleFilt is the correct sampling rate, since that was just done above
            resamplebytes = resample(fileEEG(ind,1:fileEEGlen(ind))',resampleFreq,sampleFreq,resampleFilt);
            fileEEG(ind,1:length(resamplebytes)) = resamplebytes';
        end
    end
    
    EEG(fileMask,:) = gain .* fileEEG(:,buffer+1:duration+buffer);
end

%JRM NOTE: this is old code -- I left it in from the old version of gete_ms
%take relative baseline correction
if length(RelativeMS) == 2
    % get the average for the range
    offset = fix((OffsetMS)*resampleFreq/1000);
    relative = fix((RelativeMS)*resampleFreq/1000);
    relative = relative - offset + 1;
    relative(2) = relative(2) - 1;

    % calculate the relative
    releeg = mean(EEG(:,relative(1):relative(2)),2);

    % subtract the baseline
    EEG = bsxfun(@minus, EEG, releeg);
end
