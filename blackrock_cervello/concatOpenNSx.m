function nsx_data = concatOpenNSx(nsx_filepath, applyGain, skipfactor, chanToOpen, nonClockResetTimeStampThresh)
%this function will read an NSx file using openNSx
%it will remove the first data segment if it's a junk segment
%and concatenate the remaining data segments and other segmented struct fields
%
%  skipfactor - (1) dont skip any sample, >1 skip samples
%  chanToOpen - [] for all channels, else specify channels


VERBOSE = 0; %- 0:  set to 1 to get more outputs about cells existing, etc

if ~exist(nsx_filepath, 'file')
    error('%s is not a valid filepath', nsx_filepath);
end



%- options for how openNSx is called
if nargin < 2,
    applyGain = 0;
end
if nargin < 3,
    skipfactor = 1;
end
if nargin < 4,
    chanToOpen = [];
end
if nargin<5;
    nonClockResetTimeStampThresh = 200;
end

if isempty(chanToOpen)
    nsx_data = openNSx(nsx_filepath,'skipfactor',skipfactor);
else
    nsx_data = openNSx(nsx_filepath,'skipfactor',skipfactor,'channels',chanToOpen); %- JW tweaked so chanToOpen can actually be passed 12/2018
end

%- add catch for CZ's new escape from openNSx (trying to avoid those infinite loops)
if ~isstruct(nsx_data) && nsx_data==-99,
    fprintf('\n caught infinite loop error from openNSx... pushing error up');
    return
end

nsx_data.postProc = struct;
nsx_data.postProc.samplingFreq = nsx_data.MetaTags.SamplingFreq;
nsx_data.postProc.removedJunk = 0;
nsx_data.postProc.appliedGain = 0;

% disp(nsx_data);

if iscell(nsx_data.Data)
    
    if VERBOSE, fprintf('\n in concatOpenNSx, data contains %d segments ', length(nsx_data.Data)); end 
    nsx_data.postProc.numRawSegments     = length(nsx_data.Data);
    nsx_data.postProc.numNonJunkSegments = length(nsx_data.Data);
    
    %multi-segment data
    if length(nsx_data.Data) > 1
        
        %if first segment is less than 5 sec, treat it as a pre-sync junk segment
        if nsx_data.MetaTags.DataDurationSec(1) < 5
            
            removedTimeSec = nsx_data.MetaTags.DataDurationSec(1);
            
            if VERBOSE, fprintf('first data segment only %.1f seconds long, removing it\n', removedTimeSec); end
            
            nsx_data.postProc.numJunkSamplesRemoved = nsx_data.MetaTags.DataPoints(1);
            
            nsx_data.MetaTags.Timestamp(1)  = [];
            nsx_data.MetaTags.DataPoints(1) = [];
            nsx_data.MetaTags.DataDurationSec(1) = [];
            nsx_data.Data(1) = [];
            
            nsx_data.postProc.removedJunk = 1;
            nsx_data.postProc.numNonJunkSegments = nsx_data.postProc.numRawSegments - 1;
            
        end
        
        %if there are still multiple segments after the junk segment is
        %possible removed, we have run into the blackrock clock reset issue
        if length(nsx_data.Data) > 1
            
            if VERBOSE, fprintf('concatenating remaining %d data segments\n\n', length(nsx_data.Data)); end
            
            
            original_timestamps = nsx_data.MetaTags.Timestamp;
            original_cell_lengths = cellfun(@(x) size(x,2), nsx_data.Data);
            
            nsx_data.postProc.nonJunkTimeStamps = original_timestamps;
            nsx_data.postProc.nonJunkCellLengths = original_cell_lengths;
            nsx_data.postProc.samplesAdded = [];
            nsx_data.postProc.nonJunkSegConcatIdx = [1 size(nsx_data.Data{1,1}, 2)];
            nsx_data.postProc.fillerSegConcatIdx = [];
            
            postProcIter = 1;
            
            while length(nsx_data.Data)~=1
                
                
                timeStampDifference = original_timestamps(2) - original_timestamps(1);
                
                samplingFreq = nsx_data.MetaTags.SamplingFreq;
                samplingRateOfTimestamps = 30000;
                segmentLength = original_cell_lengths(1);
                samplesInGap = round(timeStampDifference * (samplingFreq/samplingRateOfTimestamps));
                
                numSamplesToInsert = floor((samplesInGap - segmentLength) / skipfactor);  
                    
                % first check for a packet loss, in which both timestamps
                % are large (first segment was not after clock reset)
                check1 = (original_timestamps(1) > nonClockResetTimeStampThresh) && (original_timestamps(2) > nonClockResetTimeStampThresh) && (numSamplesToInsert>=1);
                % check for a packet loss, in which first segment was
                % result of clock reset
                check2 = (original_timestamps(1) < nonClockResetTimeStampThresh) && (original_timestamps(2) > nonClockResetTimeStampThresh) && (numSamplesToInsert>=1);
                if check1 || check2 % checks for non-clock reset split
                    
                    if numSamplesToInsert==1
                        samplesToInsert = nsx_data.Data{1,1}(:,end);
                        
                    elseif mod(numSamplesToInsert,2)==1
                        
                        takeFromBefore = grabFillerSamples(nsx_data.Data{1,1}, round(numSamplesToInsert/2), 'end');
                        takeFromAfter = grabFillerSamples(nsx_data.Data{1,2}, round(numSamplesToInsert/2)-1, 'start');
                       
                        samplesToInsert = cat(2,fliplr(takeFromBefore),fliplr(takeFromAfter));
                        
                    elseif mod(numSamplesToInsert,2)==0
                        
                        takeFromBefore = grabFillerSamples(nsx_data.Data{1,1}, round(numSamplesToInsert/2), 'end');
                        takeFromAfter = grabFillerSamples(nsx_data.Data{1,2}, round(numSamplesToInsert/2), 'start');          
                        
                        samplesToInsert = cat(2,fliplr(takeFromBefore),fliplr(takeFromAfter));
                        
                    end
                    
                    current_fillerSegConcatIdx_start = size(nsx_data.Data{1}, 2) + 1; %next index after current data
                    current_fillerSegConcatIdx_end = current_fillerSegConcatIdx_start + size(samplesToInsert, 2)-1; % last index of filler samples added
                    nsx_data.postProc.fillerSegConcatIdx = [nsx_data.postProc.fillerSegConcatIdx; current_fillerSegConcatIdx_start current_fillerSegConcatIdx_end];

                    
                    current_nonJunkSegConcatIdx_start = current_fillerSegConcatIdx_end + 1; % next index after filler samples
                    current_nonJunkSegConcatIdx_end = current_nonJunkSegConcatIdx_start + size(nsx_data.Data{2}, 2)-1; % last index of data segment added
                    nsx_data.postProc.nonJunkSegConcatIdx = [nsx_data.postProc.nonJunkSegConcatIdx ; current_nonJunkSegConcatIdx_start current_nonJunkSegConcatIdx_end];
                    
                    nsx_data.Data{1} = horzcat(nsx_data.Data{1},samplesToInsert,nsx_data.Data{2});
                    
                else
                    
                    numSamplesToInsert = 0;
                    
                    current_nonJunkSegConcatIdx_start = size(nsx_data.Data{1}, 2) + 1;
                    current_nonJunkSegConcatIdx_end = current_nonJunkSegConcatIdx_start + size(nsx_data.Data{2}, 2)-1;
                    nsx_data.postProc.nonJunkSegConcatIdx = [nsx_data.postProc.nonJunkSegConcatIdx ; current_nonJunkSegConcatIdx_start current_nonJunkSegConcatIdx_end];
                    
                    nsx_data.Data{1} = horzcat(nsx_data.Data{1},nsx_data.Data{2});
                    
                end
                
                nsx_data.postProc.samplesAdded(postProcIter) = numSamplesToInsert;
                postProcIter = postProcIter + 1;
                
                nsx_data.Data(2) = [];
                original_timestamps(1) = [];
                original_cell_lengths(1) = [];
            end

            
        end   % real data segment end

    end   % junk segment end
    
    nsx_data.Data = nsx_data.Data{1};
    
    nsx_data.MetaTags.DataPoints = size(nsx_data.Data,2);
    nsx_data.MetaTags.DataDurationSec = size(nsx_data.Data,2) / nsx_data.MetaTags.SamplingFreq;
    
end


if applyGain == 1
    
    fprintf('apply 1/4 gain to data\n');
    nsx_data.postProc.appliedGain = 1;
    
    nsx_data.Data = nsx_data.Data/4;
end

end


function samps = grabFillerSamples(data_mat, numSamples, locationStr)

    if isequal(locationStr, 'start')
        
        if size(data_mat, 2) >= numSamples
            
            samps = data_mat(:,1:numSamples);

        else % need to multiply samples
            
            samps = replicate_data(data_mat, numSamples); 

        end
        
    elseif isequal(locationStr, 'end')
        
        if size(data_mat, 2) >= numSamples
            
            samps = data_mat(:,(end-numSamples+1):end);
            
        else % need to multiply samples
            
            samps = replicate_data(data_mat, numSamples);
           
        end
        
    else
        error('locationStr: %s is not valid. Must be either "start" or "end"\n', locationStr);
    end

end

function samps = replicate_data(data_mat, numSamples)
            
    % grab the samples we have
    temp_samps = data_mat;

    % replicate the data until the needed sample count is reached, mirrored repeats

    increment_idx = size(data_mat, 2);
    incrementer = 1; % 1 or -1

    while size(temp_samps, 2) < numSamples

        temp_samps = [ temp_samps data_mat(:,increment_idx) ];

        % if the data contains more than a single point, move an incrementer
        if size(data_mat, 2) > 1

            if increment_idx == size(data_mat, 2) && incrementer == 1 %reached the end of the data while moving forward

                % switch the incrementer direction
                incrementer = -1;

            elseif increment_idx == 1 && incrementer == -1 %reached the beginning of the data while moving backward

               % switch the incrementer direction
               incrementer = 1;

            end

            increment_idx = increment_idx + incrementer;

        end
    end

    samps = temp_samps;
    
end