function [timeSeries, mask_stats] = createNoiseMaskSPK(noiseUnits,lengthMicro,varargin)
%
%  Creates a combined timeseries mask based on the timestamps of the provided units. For each separate
%  list of units provided, createNoiseMask will create separate noise masks for each list. To create a
%  noise mask, the function creates a timeseries mask for each separate unit, applies a convolution with 
%  a specified kernel length to each individual unit mask, and then keeps any values that are present in 
%  a specified number of units. Each timeseries mask will be the length of the micro file and contain 1s 
%  where the mask exists (0 otherwise).
% 
%  [timeSeries, mask_stats] = createNoiseMaskSPK(noiseUnits,lengthMicro,'windowDur',windowDur,'numChanReq',numChanReq)
% 
%  Called by extractSpikeInfo in which noiseUnits = {globalUnits1, globalUnits2, stimChannelUnits};
%
%   Inputs (Required):
%       noiseUnits      = (cell) array of units you want to apply mask to, in which each element
%                         contains 1+ units that will used to generate a mask for all units in that
%                         element. For example, noiseUnits = {list1, list 2, list3, ..., listn}, where
%                         listn = {unit1, unit2, unit3, ..., unitn}, and unitn is a cell of timestamps
%                         ex) {{120,9387,90003,...},{2,239,11223,...},{876,3456,9276,...}},
%                              {{287,1678,19876,...},{98,8269,12876,...}}}
%       lengthMicro     = (double) the length of the micro file in ms
%   Inputs (Optional):
%       windowDur       = (double) length of convolution kernel (Default is 5000 ms)
%       numChanReq      = (double) the number of units required to have noise in order for the mask to
%                         count it as noise (Default is 2)
%  
%    Outputs:
%       timeSeries      = (cell) 1xnumel(noiseUnits) array containing the timeseries (1xlengthMicro) mask 
%                         for each element in noiseUnits: {mask1_ts,mask2_ts,mask3_ts,...,maskn_ts}
%                         1 means part of the mask
%       mask_stats      = (double) 1xnumel(noiseUnits) array containing the percentage of the micro file
%                         duration that is covered by each noise mask
%
% Example:
%       noiseUnits = {{{120,9387,90003,...},{2,239,11223,...},{876,3456,9276,...}},
%                     {{287,1678,19876,...},{98,8269,12876,...}}}
%           will create a mask for each set of units (2 sets here, the 1st with 3 units, the 2nd with 2
%           units). timeSeries will be a 1x2 cell, each cell containing a 1xlengthMicro mask of 0s and
%           1s. mask_stats will also be a 1x2 cell.
%
%
%   Created by Ai Phuong Tong and Srijan Bhasin on 03/10/2020, edited by Samantha Jackson
%   
% 3/10/2020 - Created by APT and SB
% 3/11/2020 - Edited by SJ: changed input variables to 1 cell array, added numChanReq, changed
%             definition, added varargin parse inputs
% 4/23/2020 - SJ: added 'round' to account for decimal spike timestamps
% 8/24/2020 - SJ: if empty array given, skip and count timeseries and mask_stats as []
% 1/14/2021 - SJ: Get rid of any spiketimes that occur after the end of the micro file (due to trimming)
% 2/26/2021 - SJ: Fix numChanReq- beforehand it was ignoring this and keeping all

% parse inputs
inp_pars = inputParser;
defaultwindowDur = 5000; %Default size of convolution kernel set to 5000ms
defaultnumChanReq = 2; %Default is to require 2 channels

addParameter(inp_pars,'windowDur',defaultwindowDur,@isnumeric);
addParameter(inp_pars,'numChanReq',defaultnumChanReq,@isnumeric);

parse(inp_pars,varargin{:})
windowDur = inp_pars.Results.windowDur;
numChanReq = inp_pars.Results.numChanReq;

arguments = noiseUnits;

% initialize variables
inputTemp = cell(1,numel(arguments));
outputCell = cell(1,numel(arguments));
mask_stats = cell(1,numel(arguments));

for ii = 1:numel(arguments)
    if isempty(arguments{ii})
        outputCell{ii} = [];
        mask_stats{ii} = [];
        continue
    end
    for units = 1:numel(arguments{ii})
        inputTemp{ii}{units} = zeros(1,lengthMicro); % Initialize with zeros
        % SJ: if any of the timestamps in arguments are > lengthMicro, this is because there was trimming from the start or end that made so that some timestamps are not included.
        % This happens because when timestamps are aligned, they are added or subtracted to after drift correction, but only to account for the beginning add/cut, not for the end
        if any(round(arguments{ii}{units}) > lengthMicro)
            arguments{ii}{units}(round(arguments{ii}{units}) > lengthMicro) = [];
        end
        inputTemp{ii}{units}(round(arguments{ii}{units})) = 1; % Assign all of the timestamps for this unit to 1; SJ: added 'round' for decimal cases
        temp = conv(inputTemp{ii}{units}, ones(1,windowDur)); % Convolution
        inputTemp{ii}{units} = temp(round(windowDur/2):end-round(windowDur/2)); % Trimming ends
    end
    outputCell{ii} = cat(1, inputTemp{ii}{:}); % Concatenate each unit 
    outputCell{ii}(outputCell{ii}>0) = 1; % If any decimals, assign to 1
    outputCell{ii} = sum(outputCell{ii},1); % Sum to see how many units have noise at each timepoint
    outputCell{ii}(outputCell{ii}<numChanReq) = 0; % Get rid of
    outputCell{ii}(outputCell{ii}>=numChanReq) = 1; % Only keep noise when the sum > numChanReq
    mask_stats{ii} = sum(outputCell{ii})/length(outputCell{ii}); % Calculate the percentage of the micro file that the mask covers
end

%Final Output Cells
timeSeries = outputCell;
%mask_stats = {sum(timeSeries{1})/length(timeSeries{1}), sum(timeSeries{2})/length(timeSeries{2}), sum(timeSeries{3})/length(timeSeries{3})};

end