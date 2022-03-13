function [output_br_timeStamp_up, output_ts, IPIViolationFlag] = getBlackRockPulsesDC(NEVdata,whichDC,postProc)
%Given the NEV data, outputs the pulse information from a certain DC
%channel
%
% whichDC = 9, 11, 12, etc.

% keyboard

output_br_timeStamp_up = [];
output_ts = [];
IPIViolationFlag = 0;

if nargin < 3,
    postProc = [];
end

% 2. Get raw pulse information
rawUnparsedData = NEVdata.Data.SerialDigitalIO.UnparsedData;
raw_br_timeStamp = double(NEVdata.Data.SerialDigitalIO.TimeStamp);

binaryChar = dec2bin(rawUnparsedData);
binaryMat = double(binaryChar);
binaryMat(binaryMat==49) = 1; binaryMat(binaryMat==48) = 0;

diffMat = nan(size(binaryMat));
for instance = 2:size(binaryMat,1)
    diffMat(instance,:) = binaryMat(instance,:) - binaryMat(instance-1,:);
end

timeStamps = cell(1,size(binaryMat,2));
for loc = 1:size(binaryMat,2)
    timeStamps{loc} = cat(2,diffMat(diffMat(:,loc)~=0,loc),raw_br_timeStamp(diffMat(:,loc)~=0)');
end

br_timeStamp_up = [];
br_timeStamp_down = [];

if ~isempty(timeStamps)
    
    % DC09 is D3, which is first from the back (D0,D1,D2,D3)
    % DC11 is D1, which is third from the back (D0,D1,D2,D3)
    % DC12 is D0, which is fourth from the back (D0,D1,D2,D3)
    
    DCposition = length(timeStamps) - (whichDC - 9);
    
    br_timeStamp_up = timeStamps{DCposition}(timeStamps{DCposition}(:,1)==1,2);
    br_timeStamp_down = timeStamps{DCposition}(timeStamps{DCposition}(:,1)==-1,2);
        
end


%- cut initial pulses; important for downsampling and removing crap
numSampCutFromStart = 30000; %- 30k sample rate, so this is 1 second
br_timeStamp_up = br_timeStamp_up(br_timeStamp_up>numSampCutFromStart); 


output_br_timeStamp_up = br_timeStamp_up;   


%- if any uptimes were recorded after the first second, and none of them were nan
if ~isempty(br_timeStamp_up) && ~any(isnan(br_timeStamp_up))
    
    if ~isempty(postProc) && isfield(postProc,'samplesAdded')
        
        [br_timeStamp_up, IPIViolationFlag] = correctSplitNEV(br_timeStamp_up, whichDC, postProc);
        [br_timeStamp_down, ~] = correctSplitNEV(br_timeStamp_down, whichDC, postProc);
    end
    
    %set output value
    output_br_timeStamp_up = br_timeStamp_up;

    %combine to produce ts
    [br_timeStamp_combo, sorted_idx] = sort([br_timeStamp_up ; br_timeStamp_down]);
    
    up_down_indicator = [ repmat(1, length(br_timeStamp_up), 1) ; repmat(-1, length( br_timeStamp_down), 1) ];
    sorted_up_down_indicator = up_down_indicator(sorted_idx);
    
    %remove pulses that happened before 30000 samples
    
    br_timeStamp_combo = br_timeStamp_combo(br_timeStamp_combo > numSampCutFromStart);
    sorted_up_down_indicator = sorted_up_down_indicator(br_timeStamp_combo > numSampCutFromStart);
    
    br_timeStamp_combo_ms = floor(br_timeStamp_combo ./ 30);

    i = 1;

    %skip any initial pulse-down indicators
    while sorted_up_down_indicator(i) == -1
       i = i + 1; 
    end

    last_ts = 1;
    while i <= length(br_timeStamp_combo_ms)

        current_ts = br_timeStamp_combo_ms(i);

        if sorted_up_down_indicator(i) == 1 % reach a pulse-up

            %preceeding interval was pulse-down
            output_ts(last_ts:current_ts-1) = 0;

            last_ts = current_ts;

        elseif sorted_up_down_indicator(i) == -1 % reach a pulse-down

            %preceeding interval was pulse-up
            output_ts(last_ts:current_ts-1) = 1;

            last_ts = current_ts;

        end

        i = i + 1;
    end

end


output_br_timeStamp_up = output_br_timeStamp_up(output_br_timeStamp_up > numSampCutFromStart);

end

