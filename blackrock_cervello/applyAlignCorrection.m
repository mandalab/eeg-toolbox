function [signal1Fixed, signal2Fixed, SPKtimestampsFixed, transformsFixed, unstretchedOffset] = applyAlignCorrection(signal1, signal2, iMod, driftRate, driftFull, removeDrift, keepFirstSignal, varargin)
%applyCorrectDriftFull applies the drift correction specified in the inputs
%   This is an adaptation of the function correctDriftFull that was originally in align_nsps, but it now
%   replaces most of its original functionality
%   
%   Called by align_nsps and align_SpikesWithEcog_ManualORAuto
%
% Inputs:
%   signal1
%   signal2
%   seglength        = 1x2 vector with the lengths of the segment1 and segment2 (in that order). The
%                       segment lengths are the lengths between the riding and falling edges
%   signal2type      = The type of data signal2 is. Options: 'ECOG', 'MICRO', 'SPK', 'LFP'
%   keepFirstSignal  = 1 if you do not want signal1 to change
%                      0 if you do not care
%
% 9/18/2020 SJ (committed with JW's account)- added floor to idx initialization for CV cases because SRmatchECOG = 1024/1000 = 1.024
% 4/02/2021 SJ: get rid of SRmatchECOG

inp_pars = inputParser;
defaultSPKtimestamps = {};
defaultalignType = 'ECOG';
defaultSRmatchECOG = 1;
defaultSPKtrimStart = 0;
defaultidxRE = [];
defaultapply_theseChans = [];

addParameter(inp_pars,'SPKtimestamps',defaultSPKtimestamps,@iscell);
addParameter(inp_pars,'alignType',defaultalignType,@ischar);
addParameter(inp_pars,'SRmatchECOG',defaultSRmatchECOG,@isnumeric);
addParameter(inp_pars,'SPKtrimStart',defaultSPKtrimStart,@isnumeric);
addParameter(inp_pars,'idxRE',defaultidxRE,@isnumeric);
addParameter(inp_pars,'apply_theseChans',defaultapply_theseChans,@isnumeric);

parse(inp_pars,varargin{:})
SPKtimestamps = inp_pars.Results.SPKtimestamps;
alignType = inp_pars.Results.alignType;
SRmatchECOG = inp_pars.Results.SRmatchECOG;
SPKtrimStart = inp_pars.Results.SPKtrimStart;
idxRE = inp_pars.Results.idxRE;
apply_theseChans = inp_pars.Results.apply_theseChans;

%- shouldn't be too much... 30 is most I've seen for 1hour session
if abs(driftFull)>100
    fprintf('\n CORRECTING DRIFT: heads up, extrapolated drift is >50 samples (%d)',driftFull);
    keyboard;
end

if iMod == 1
    sigToMod = signal1;
    whichNSPnotmod = 2;
elseif iMod == 2
    sigToMod = signal2;
    whichNSPnotmod = 1;
else
    fprintf('%s\n','ERROR!!! sigToMod must be equal to 1 or 2');
    keyboard
end

if SRmatchECOG ~=1
    fprintf('%s\n','SRmatchECOG is not equal to 1. Are you sure you want to continue?');
    keyboard
end

sigLength = [size(signal1,2), size(signal2,2)];

timeStamp = SPKtimestamps;


%- prepare to stretch the slower signal
numChans      = size(sigToMod,1);

%For this, driftFull will not be negative yet!!!
if strcmpi(alignType,'ECOG') || strcmpi(alignType,'MICRO')    
    if ~removeDrift % If you are already changing signal2 or you don't care which gets changed, proceed normally
        modSignal = zeros(numChans,sigLength(iMod) + driftFull);
        %- where to add in the extra samples
        correctionInterval = floor(sigLength(iMod)/driftFull);
        if correctionInterval<1000
            fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
            keyboard;
        end
        %- indices that will get a second copy (steps of correction interval, starting 1/2 way in correction interval
        idx = zeros(1,driftFull);
        for ii=1:driftFull
            %-  put it in the middle of the correction interval because that makes sense
            %       and so an extra sample isn't added to the end on even drifts
            idx(ii) = floor(correctionInterval/2) + correctionInterval*(ii-1);
        end
        %- setup template for how extra samples are filled.. just a copy of the immediately preceeding sample
        repN = ones(1,sigLength(iMod));
        repN(idx) = 2;
        %- make the drift correction
        for jj=1:numChans
            modSignal(jj,:) = repelem(sigToMod(jj,1:sigLength(iMod)),repN); %- fill gap with a copy of the preceeding sample
        end
    else %Need to delete from signal2 instead
        modSignal = sigToMod;
        % Where to remove samples
        correctionInterval = floor(sigLength(iMod)/driftFull);
        if correctionInterval<1000
            fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
            keyboard;
        end    
        %- indices that will be removed (steps of correction interval, starting 1/2 way in correction interval
        idx = zeros(length(driftFull));    
        for ii=1:driftFull
            %-  remove from the middle of the correction interval because that makes sense
            %       and so an extra sample isn't added to the end on even drifts
            idx(ii) = floor(correctionInterval/2) + correctionInterval*(ii-1); % Don't have the +1 anymore
        end
        modSignal(:,idx) = [];
    end

    %- tweak the rising edge index of the streched signal to account for the stretch
    %if iShort==1, %iMod?
    if ~removeDrift  % Go ahead and add Drift to signal1 or 2
        if iMod == 1
            idxRE1shift = idxRE(1) + sum(idx<idxRE(1)); %- any correction intervals before idxRE1 will offset it
            deltaRE = idxRE1shift-idxRE(2);
            signal1Fixed = modSignal;
            signal2Fixed = signal2;
        else
            idxRE2shift = idxRE(2) + sum(idx<idxRE(2));
            deltaRE = idxRE(1)-idxRE2shift;
            signal1Fixed = signal1;
            signal2Fixed = modSignal;
        end
    else %removing Drift from signal2
        if iMod ~= 2
            fprintf('%s\n','ERROR!!!!!!! removeDrift flagged, but trying to delete from Signal1 (Should never happen!). Ask SJ!');
            keyboard
        else
            idxRE2shift = idxRE(2) - sum(idx<idxRE(2)); %subtract the number of idxs that were removed before the rising edge, from the rising edge
            deltaRE = idxRE(1)-idxRE2shift; %new delta Rising edge is probably a little larger now
            signal1Fixed = signal1;
            signal2Fixed = modSignal;   
        end
    end

    %- trim the time series so there is no delta rising edge
    trimStart(1) = 0; trimStart(2) = 0;
    if deltaRE>0
        %- idx1 greater than idx2, so trim 1 (cant add imaginary stuff to 2 to match) -SJ: well now we can
        if keepFirstSignal ~= 1 %okay, go ahead and trim from start of signal1
            signal1Fixed = signal1Fixed(:,deltaRE+1:end);
            trimStart(1) = deltaRE;
        else %Not okay, need to pad signal2
            signal2Fixed = cat(2,NaN(size(signal2Fixed,1),abs(deltaRE)),signal2Fixed);
            trimStart(2) = -deltaRE;
        end
    elseif deltaRE<0
        signal2Fixed = signal2Fixed(:,abs(deltaRE)+1:end);
        trimStart(2) = abs(deltaRE);
    end

    %- and last trim the end so they are the same length
    minDur = min([size(signal1Fixed,2) size(signal2Fixed,2)]);
    trimEnd1 = size(signal1Fixed,2)-minDur;  
    trimEnd2 = size(signal2Fixed,2)-minDur;

    if trimEnd1 < 0
       fprintf('%s\n','ERROR!!! Why is the trim end already negative??? Ask SJ!');
       keyboard
    elseif trimEnd1 > 0 %Going to trim signal1 - is that okay?
        if trimEnd2 ~= 0
            fprintf('%s\n','ERROR!!! Why are we trying to trim both??? Ask SJ!');
            keyboard
        end
        if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
            signal1Fixed = signal1Fixed(:,1:minDur); % Trimming end of signal1
        else %Not okay, need to pad end of signal2
            trimEnd2 = -1*trimEnd1;
            trimEnd1 = 0;
            signal2Fixed = cat(2,signal2Fixed,NaN(size(signal2Fixed,1),abs(trimEnd2))); %adding NaN values to end
        end
    elseif trimEnd2 > 0 % Trim the end of signal2
        signal2Fixed = signal2Fixed(:,1:minDur);
    elseif trimEnd1 == 0 && trimEnd2 == 0
        % Not trimming enything
    else
        fprintf('%s\n','ERROR!!! How do we get here??? Ask SJ!');
        keyboard
    end

    %- update NSP sync selection (non-stretched NSP should be used for getting sync data, with the offset here)
    % If you have keepFirstSignal, then you will always be modifying only signal2 and whichNSPnotmod will always
    % be 1. Because Sig1 was not modified ever, trimStart(1) will always be 0
    % If you don't have keepFirstSignal, then sig1 or sig2 will be modified, but you will never end up with a
    % negative trimStart(whichNSPnotmod). Do a check just to make sure:
    % (I don't think this matters for micro but make sure it always works for the other cases!)

    if keepFirstSignal == 1 && whichNSPnotmod==1 && trimStart(whichNSPnotmod)==0
        % Good
    elseif keepFirstSignal ~= 1 && trimStart(whichNSPnotmod) >= 0
        % Good
    else
        fprintf('%s\n','ERROR!!!!!! You should never get here!! Ask SJ!');
        keyboard;
    end

    %Assign output values
    if removeDrift
        driftAdded = -driftFull;
    else
        driftAdded = driftFull;
    end
    if iMod == 1
        driftNSP1 = driftAdded;
        driftNSP2 = 0;
    else
        driftNSP1 = 0;
        driftNSP2 = driftAdded;
    end

    %- update the transform outputs
    transformsFixed{1,1} = [trimStart(1) trimEnd1  driftNSP1]; %- NSP1: trim start, trim end, add drift
    transformsFixed{1,2} = [trimStart(2) trimEnd2  driftNSP2]; %- NSP2: trim start, trim end, add drift
    transformsFixed{1,3} = [ driftRate driftAdded iMod];                %- DRIFT: drift rate, drift added, whichNSP
    SPKtimestampsFixed = {};
    unstretchedOffset = 1+trimStart(whichNSPnotmod); %- this will be used to offset the DC channels from the unstreched NSP in split_BR

elseif strcmpi(alignType,'SPK')
    
    if iMod ~= 2 || keepFirstSignal ~= 1 || isempty(timeStamp)
        fprtinf('%s\n','ERROR!!! Attempting to align Spikes but incorrect input settings!!! Ask SJ!');
        keyboard
    end
    
    sizeMicro = size(signal2,2);
    
    correctionInterval = floor(sizeMicro/abs(driftFull));
    if correctionInterval<1000
        fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
        keyboard;
    end

    timeStamp_aligned = cell(size(timeStamp,1),1);
    for unit_num = 1:length(apply_theseChans)
        unit = apply_theseChans(unit_num);
        if driftFull ~= 0
            %- indicies that will get a second copy (steps of correction interval, starting 1/2 way in correction interval
            idx = zeros(1,floor(abs(driftFull)));%SJ: Adding floor to SRmatchECOG for when it is a fraction for CV
            for ii=1:numel(idx)
                %-  put it in the middle of the correction interval because that makes sense
                %       and so an extra sample isn't added to the end on even drifts
                idx(ii) = floor(correctionInterval/2) + correctionInterval*(ii-1);
            end

            % given those indices, what changes are we going to apply to each spike
            applyChange = zeros(1,length(timeStamp{unit,1}));
            for sample = 1:length(idx)
                if sample==length(idx) % then only search for spikes beyond this time
                    temp_match_spikes = find(timeStamp{unit,1} > idx(sample)); 
                else %look for spikes that fall in the range of indices
                    temp_match_spikes = find((timeStamp{unit,1} > idx(sample)) & (timeStamp{unit,1} <= idx(sample+1)));
                end
                %if we find a match, assign this change to it
                if ~isempty(temp_match_spikes)
                    applyChange(temp_match_spikes) = sample;
                end
            end
            % apply the changes to spikes, whether subtraction or removal
            %keyboard %is applyChange correct?
            if driftFull < 0
                timeStamp_aligned{unit,1} = timeStamp{unit,1}  - applyChange;
            elseif driftFull > 0
                timeStamp_aligned{unit,1} = timeStamp{unit,1}  + applyChange;
            end
        else
            timeStamp_aligned{unit,1} = timeStamp{unit,1}; % Leave as-is
        end
        % now, apply offset correction  - (the transforms
        % outputted by align_nsps describe a drift that
        % should be applied to the full signal, and
        % the start/end cuts describe operations that
        % should be completed after drift correction)
        % and since transforms is in 1kHz indices but
        % timestamps are in 1kHz indices
        % Note: if the 'trim start' from transforms is negative, then this will add that amount instead
        % of subtract.
        timeStamp_aligned{unit,1} = timeStamp_aligned{unit,1} - (SPKtrimStart);
        % Now, remove any negative values (ie occurring before the start of the eCOG file)
        if any(timeStamp_aligned{unit,1}<0)
            % First verify that they are all at the start and in a row
            mask_neg = timeStamp_aligned{unit,1}<0;
            if sum(diff(mask_neg(1:find(mask_neg,1,'last')))) ~= 0
                fprintf('%s\n','ERROR!!!! Detecting negative spike times but they are either not in a row or not starting at the beginning!! What is up?? Ask SJ!');
                keyboard
            else
                timeStamp_aligned{unit,1}(mask_neg) = [];
            end
        end
    end

    %transforms{1,2}(1,1)
    signal1Fixed = signal1;
    signal2Fixed = signal2;
    transformsFixed = {};
    SPKtimestampsFixed = timeStamp_aligned;
    unstretchedOffset = [];
    
else
    fprintf('%s\n',['Not equipped to align this type (' alignType ')']);
    keyboard
end
    
    
    

    
    %- prepare to stretch the slower signal
%numChans      = size(sigToMod,1);

    
end

