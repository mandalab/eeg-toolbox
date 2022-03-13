function [signal, syncPulse, syncStartIndex, errorStr, transforms] = align_nsps(sync, signal1, signal2,sampRate,systemType,pulseType,writeFigsHere,sessTitle,differentFileSystems,varargin)

%% Drift correction/alignment script for non-synchronized datafiles
% Purpose:          Determine offset/alignment off of sync pulses into two
%                   NSPs, then align and correct drift in data files
%
% adapted by M.Elkalliny from Blackrock code, to interface with eegPrepAndAlign
% called by split_BR (blackrock) and split_CV (cervello). December 2018
%
% Inputs:
%     Required:
%         sync                  = 1x4 cell array, containing NSP1-DC12-analog, NSP2-DC12-analog,
%                                 NSP1-DC09-analog, and NSP2-DC09-analog    
%         signal1               = signal from EEG file for system 1
%         signal2               = signal from EEG file for system 2
%         sampRate              = sample rate of data (ecog) which should be the same as the second data type
%         systemType            = 'CV' or 'BR'
%         pulseType             = 9 or 12
%         writeFigsHere         = directory of raw file, to write summary figs into
%         sessTitle             = session name
%                                 ex) '180331_0415'
%         differentFileSystems  = 1 or 0
%                                 1: different file systems, such as BR NSP + micro
%                                 0: same file systems, such as BR NSP to BR NSP
%
%     Optional: (name-value pairs)
%         keepFirstSignal       = 1 or 0
%                                 1: if you only want to apply changes to the second type of system (ie
%                                    if you are doing micro alignment, in which case signal2 should be
%                                    micro data, and you only want the micro data to change)
%                                 0: if you don't care which changes (ie for NSP-NSP alignment) -
%                                    eventually, all should probably run with 1
%         isMicro               = 1 or 0
%                                 1: micro alignment to eCOG
%                                 0: any other type of alignment
%         overwriteLabels       = 1x3 cell array used for micro alignment containing the signal 1 type, 
%                                 signal 2 type, and signal 2 type channel names
%                                 ex) {'ECOG','MICRO','utah'};
%         zeroNeeded            = 1x4 numerical array defining which portions of the signals were zeroed out if the start times were not within 5 min
%                                 zeroNeeded = [EcogBeginning, EcogEnd, MicroBeginning, MicroEnd]
%           
% outputs: signal (1x2 cell if blackrock, 1x1 cell if cervello)
% syncPulse - which NSP to use for behavioral alignment
% syncStartIndex - index of analog timeseries to use for pulses
% errorStr. these are predefined: 0 = NSPs already synced. 1 = no pulses
%   present on either DC09 or DC12. 2 = NSPs were not synced, but have been
%   synced using DC12. 3 = NSPs were not synced, but have now been synced
%   using DC09. 4 = pulses are there, but they dont seem to be matching
% transforms - 1x5 cell
%    {1} = NSP1 add/remove start, end, drift,  {2} = NSP2 add/remove start, end, drift
%           ^^ positive value = trim
%    {3} = DRIFT: drift rate, drift added, whichNSP
%    {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
%    {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
%    {6} = XCORR correction on indices: rising diff, falling diff, used correction
%
% 2/2020 SJ: Overhaul 
%           - clean up function definition, input, output
%           - Add 'isMicro' flag
%           - Add options to keep first signal 
%           - Use applyAlignCorrection to centralize all alignment
%           - Change sample rates to be consistent -> just one sample rate, which should be same for ecog
%             and micro
%           - Change labels to reflect data types
%           - fix a lot of alignment steps
%           - report statistics here instead of transformSync
%           - Include many checks for data type and inputs
% 6/2020 SJ: fixed fprintf statement
% 7/2020 SJ: (committed with JW account, combined with next) - added in case where one micro din channel does not have any pulses
% 8/2020 SJ: (committed with JW account) - Added in info about zeroing; moved first call to pulsealign/alignment.txt creation to inside if statement for BR or micro files only (not for CV)
% 9/2020 SJ: (committed with JW acocunt) - added ~CV requirement in a second call to pulsealign, also changed length of segLength_xcorr to end of micro file - rising edge if the rising edge + segLength_xcorr was longer than the micro file 

fprintf('Aligning pulses ...\n')

% parse inputs
inp_pars = inputParser;
defaultkeepFirstSignal = 0;
defaultisMicro = 0;
defaultoverwriteLabels = {};
defaultzeroNeeded = [0 0 0 0];

addParameter(inp_pars,'keepFirstSignal',defaultkeepFirstSignal,@isnumeric);
addParameter(inp_pars,'isMicro',defaultisMicro,@isnumeric);
addParameter(inp_pars,'overwriteLabels',defaultoverwriteLabels,@iscell);
addParameter(inp_pars,'zeroNeeded',defaultzeroNeeded,@isnumeric);

parse(inp_pars,varargin{:})
keepFirstSignal = inp_pars.Results.keepFirstSignal;
isMicro = inp_pars.Results.isMicro;
overwriteLabels = inp_pars.Results.overwriteLabels;
zeroNeeded = inp_pars.Results.zeroNeeded;


if ~isempty(overwriteLabels) %#ok<*EXIST>
    topType = overwriteLabels{1,1};
    bottomType = overwriteLabels{1,2};
    suffixChans = overwriteLabels{1,3};   
    if ~contains(bottomType,'micro','IgnoreCase',true) || contains(topType,'micro','IgnoreCase',true)
        fprintf('%s\n\t%s\n\t%s\n','ERROR!!! Please check your optionalOverwriteLabels! As of now, signal 1 should be from eCOG and signal 2 should be from micro!', ...
            ['topType = ' topType], ['bottomType = ' bottomType]);
        keyboard
    end
    if keepFirstSignal == 0
        fprintf('%s\n','ERROR!!!!! You have specified overwrite labels but keepFirstSignel is not flagged. Ask SJ!')
        keyboard
    end
    if isMicro == 0
        fprintf('%s\n','ERROR!!!!! You have specified overwrite labels but isMicro is not flagged. Ask SJ!')
        keyboard
    end
else
    topType = 'NSP1';
    bottomType = 'NSP2';
    suffixChans = '';
end

if keepFirstSignal ~= isMicro
    fprintf('%s\n','ERROR!!!!! keepFirstSignal and isMicro are not equal! Did you make a mistake? Ask SJ!')
    keyboard
end



%% CONSTANTS FOR THIS FUNCTION (could tweak but better if not)

minGap = 5; % minimum number of ms offset b/w two files required to say they are aligned already


% 2V will capture all pulses and no noise for DC12 (usually), blackrock or cervello
% but DC09s go higher and sometimes there are aberrant pulses, so for those
% we'll use 5V default threshold and later increase if needed
MIN_RANGE_PULSE_mV = 1000; %- JW: all signals are transformedto mV before getting here.  1V peak-to-peak is a good signal

MAX_THRESH_PULSES  = 6000; % if default threshold leads to picking up pulses that
% seem to be noise, then we'll increase the threshold to this, and see if
% we can pick up real pulses and not noise ones. no need to change unless
% pulses someday go beyond 6V



%% initialize some variables
sessTitle(find(sessTitle=='_'))='-'; %- remove underscore for figures

FigAndAbort = 0; % dont change this value from 0. it will be overwritten if needed

originalSig1Length = size(signal1,2);
originalSig2Length = size(signal2,2);

%- Transforms holds information about how aligment changed the files
transforms = cell(1,6);
% {1} = NSP1 add/remove start, end, drift,  {2} = NSP2 add/remove start, end, drift
% {3} = DRIFT: drift rate, drift added, whichNSP
% {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
% {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
% {6} = XCORR correction on indices: rising diff, falling diff, used correction
transforms{1,1} = [nan nan nan]; %- NSP1: trim start, trim end, add drift
transforms{1,2} = [nan nan nan]; %- NSP2: trim start, trim end, add drift
transforms{1,3} = [nan nan nan]; %- DRIFT: drift rate, drift added, whichNSP
transforms{1,4} = [originalSig1Length nan nan nan nan]; % {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
transforms{1,5} = [originalSig2Length nan nan nan nan]; % {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
transforms{1,6} = [nan nan nan]; % {6} = XCORR correction on indices: rising diff, falling diff, used correction



%%------------------------------------------------------------------------------------------------%
%%      START ALIGNMENT PROCEEDURE: CHECK PULSE QUALITY, TWEAK, ALIGN, RETURN
%%------------------------------------------------------------------------------------------------%

%%% SJ added %%%%%%%

if max(sync{2}) == 1 || max(sync{4}) == 1 %DC12 and DC09
    % Check to see if the max sync value is 1 (din) or 0 (no pulses?)
    if max(sync{2}) <= 1 || max(sync{4}) <= 1
        if isMicro == 1
            if (max(sync{2}) == 0 || max(sync{4}) == 0)
                fprintf('%s\n','Note: One of the micro sync channels (DC09 or DC12) has a max value of 0 (no pulses?) Ask SJ!');
            end
            % This is a digital pulse. Need to multiply by the max of sync1 in order to get it to work with
            % the thresholds!
            sync_orig = sync;
            sync{2} = sync{2}*max(sync{1});
            sync{4} = sync{4}*max(sync{3});
        else
            fprintf('%s\n','ERROR!!!! your max for sync2 is 1 but this is not micro data? Ask SJ!');
            keyboard
        end
    else
        fprintf('%s\n','ERROR!!!! Both channels do not have a max of 1 or less?? Ask SJ!');
        keyboard
    end
elseif max(sync{2}) == 0 && max(sync{4}) == 0
    fprintf('%s\n','ERROR!!! max(sync2) = 0 for DC09 AND DC12. Why is this?? Ask SJ!')
    keyboard
end

%%%%%%%%%%%%%%%%%%%%

DC09or12 = pulseType;  %- what is being passed in... might change below based on pulse quality
if DC09or12==12
    sync1 = sync{1}; sync2 = sync{2};
elseif DC09or12==9
    sync1 = sync{3}; sync2 = sync{4};
end
    
% first pass... just see if you can find a rising edge
[edgeThresh]   = startingThresh(sync1, sync2); %- starts 500mV above 5th percentile of sync1... should protect against outliers and give a range for increasing if necesssary
idxRisingEdge1 = checkEdge(sync1,edgeThresh,'rising');
idxRisingEdge2 = checkEdge(sync2,edgeThresh,'rising');


if (isnan(idxRisingEdge1) || isnan(idxRisingEdge2)) && (pulseType==9)
    % we started off with DC09, which means there were no DC12 chans recorded.
    % at this point we are dead. no alignment possible
    
    % does this happen?  We should probably not output 2nd NSP if so (unless DRIFT is zero because of hardware SYNC)
    %keyboard;
    
    
    syncPulse      = 1;
    syncStartIndex = 1;
    errorStr       = 1;
    FigAndAbort    = 1;
    
    
    
elseif (isnan(idxRisingEdge1) || isnan(idxRisingEdge2)) && (pulseType==12)
    % lets see if DC09 is intact
    sync1 = sync{3}; sync2 = sync{4};
    
    % find first and last pulses present in the sync timeseries
    [edgeThresh]   = startingThresh(sync1, sync2);
    idxRisingEdge1 = checkEdge(sync1,edgeThresh,'rising');
    idxRisingEdge2 = checkEdge(sync2,edgeThresh,'rising');
    
    if (isnan(idxRisingEdge1) || isnan(idxRisingEdge2))
        
        % does this happen?  We should probably not output 2nd NSP if so (unless DRIFT is zero because of hardware SYNC)
        fprintf('\npulses might totally be missing on one machine. maybe try using DC12? to visualize, continue\n')
        %keyboard;
        
        
        % ok, DC09 doesnt work either. no alignment possible
        syncPulse      = 1;
        syncStartIndex = 1;
        errorStr       = 1;
        FigAndAbort    = 1;
        
    else
        % DC09 does work, so lets switch them
        DC09or12 = 9;
    end
    
end %- check edges for DC12 and DC09


% series of more checks. above we were checking for presence of pulses at
% all. now lets check that the selected sync channels show the range we expect
% if not, we'll switch things around and see if we can use DC09
if DC09or12==9
    if range(sync{1,3}) < MIN_RANGE_PULSE_mV || range(sync{1,4}) < MIN_RANGE_PULSE_mV
        % there were no DC12 chans at all. and these are no good. so, we're canned
        FigAndAbort = 1;
        errorStr = 1;
        
    end
    
elseif DC09or12==12
    if range(sync{1,1}) < MIN_RANGE_PULSE_mV || range(sync{1,2}) < MIN_RANGE_PULSE_mV
        
        if range(sync{1,3}) < MIN_RANGE_PULSE_mV || range(sync{1,4}) < MIN_RANGE_PULSE_mV
            % both DC09 and DC12 are bad
            errorStr = 1;
            FigAndAbort = 1;
            
        else
            % DC12 is bad but DC09 seems usable, so switch them around
            sync1 = sync{1,3};
            sync2 = sync{1,4};
            DC09or12 = 9;
        end
    else
        % DC12 is good, dont do anything
    end
    
end



%%-----------------------------------------------------------------------------------------------------%
%%--------ERROR CONDITION: where NSPs cannot be synced (no valid pulses)   ----------------------------%
%%------- nsps cant be synced, so lets just write some figs out into the raw dir ----------------------%
%%-----------------------------------------------------------------------------------------------------%
if FigAndAbort==1
    % write some figs out into the raw dir to help understand why this failed
    
    %% generate figure summarizing alignment
    fS = 12;
    xTimeMin_nsp1 = [1:length(sync1)]/(sampRate*60);  %- time in minutes %SJ: used to be sampRates(2)
    xTimeMin_nsp2 = [1:length(sync2)]/(sampRate*60);  %- time in minutes %SJ: used to be sampRates(2)
    yLimits = [-500 MAX_THRESH_PULSES+1000];
    %yLimits = [min([0 cell2mat(sync)])-500 max([maxThresh cell2mat(sync)])+500];
    
    %- visualize sync signals
    figure(999); clf; set(gcf,'color','w','name','DEBUG SYNC ISSUE','units','normalized','outerposition',[0 0 1 1]);
    subplot(421)
    plot(xTimeMin_nsp1(1:length(sync{1})),sync{1});
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    ylim([yLimits(1) yLimits(2)]);
    if DC09or12==12; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
        set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    end
    title(sprintf('<%s DC 12>',topType));
    ylabel('DC Signal (mV)')
    
    subplot(423)
    plot(xTimeMin_nsp1(1:length(sync{3})),sync{3})
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    ylim([yLimits(1) yLimits(2)]);
    if DC09or12==9; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
        set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    end
    xlabel('Time (min)');
    title(sprintf('<%s DC 09>',topType));
    ylabel('DC Signal (mV)')
    
    subplot(422)
    plot(xTimeMin_nsp2(1:length(sync{2})),sync{2})
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    ylim([yLimits(1) yLimits(2)]);
    if DC09or12==12; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
        set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    end
    title(sprintf('<%s DC 12>',bottomType));
    
    subplot(424)
    plot(xTimeMin_nsp2(1:length(sync{4})),sync{4})
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    ylim([yLimits(1) yLimits(2)]);
    if DC09or12==9; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
        set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    end
    xlabel('Time (min)');
    title(sprintf('<%s DC 09>',bottomType));
    
    [ax, h] = suplabel(sprintf('%s\nNo NSP Alignment Possible',sessTitle),'t');
    set(h,'FontSize',fS+4)
    
    subplot(425)
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    h = title(sprintf('no xcorr done on first pulse, upstream errors'));
    
    subplot(426)
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    h = title(sprintf('no xcorr done, upstream errors'));
    
    subplot(427)
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    h = title(sprintf('no xcorr done on last pulse, upstream errors'));
    
    subplot(428)
    set(gca,'fontsize',fS,'box','off','tickdir','out')
    h = title(sprintf('no xcorr done, upstream errors'));
    
    fig2pngSimple(gcf,sprintf('%s%s',writeFigsHere,'/NSP_alignment_plot.png'));
    
    
    fprintf('\n Do we actually get here?  If so, we should consider tossing the 2nd NSP data for this session instead of trimming ');
    % as of Apr 2019, we only get here if these different file systems and
    % pulses were totally missing on one file system. in that case, we wont
    % apply any transforms and errorStr will make it clear to the user that
    % this data shouldnt be used. keeping the keyboard here though, in case
    % we hit this condition while doing double-NSP alignment
    % keyboard
    
    
    if differentFileSystems==1
        trimEnd1 = 0; trimEnd2 = 0; transforms{1,1} = [0 trimEnd1 0];  %- NSP1: add/remove start, end, drift
        transforms{1,2} = [0 trimEnd2 0];
        signal = signal2; return;
    end
    
    %- Do a minimum cleanup of the signals -- make them the same length
    if strcmp(systemType,'CV')
        signal = cat(1,signal1,signal2);
        trimEnd1 = 0; trimEnd2 = 0;
    else % -with SJ edits
        signal = cell(1,2); % force them to be equal sizes anyway (necessary for step 5)
        minSize = min(size(signal2,2),size(signal1,2));
        trimEnd1 = size(signal1,2) - minSize; 
        trimEnd2 = size(signal2,2) - minSize;
        if trimEnd1 < 0
           fprintf('%s\n','ERROR!!! Why is the trim end already negative??? Ask SJ!');
           keyboard
        elseif trimEnd1 > 0 %Going to trim signal1 - is that okay?
            if trimEnd2 ~= 0
                fprintf('%s\n','ERROR!!! Why are we trying to trim both??? Ask SJ!');
                keyboard
            end
            if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
                signal1 = signal1(:,1:minSize); % Trimming end of signal1
            else %Not okay, need to pad signal2
                trimEnd2 = -1*trimEnd1;
                trimEnd1 = 0;
                signal2 = cat(2,signal2,NaN(size(signal2,1),abs(trimEnd2))); %adding NaN values to end
            end
        elseif trimEnd2 > 0 % Trim the end of signal2
            signal2 = signal2(:,1:minSize);
        elseif trimEnd1 == 0 && trimEnd2 == 0
            % Not trimming enything
        else
            fprintf('%s\n','ERROR!!! How do we get here??? Ask SJ!');
            keyboard
        end
        signal{1,1} = signal1;
        signal{1,2} = signal2;
    end
    transforms{1,1} = [0 trimEnd1 0];  %- NSP1: add/remove start, end, drift
    transforms{1,2} = [0 trimEnd2 0];  %- NSP2: add/remove start, end, drift
    
    [done] = writeAlignmentSummary(transforms, errorStr, writeFigsHere, sessTitle, keepFirstSignal, isMicro, topType, bottomType, suffixChans, zeroNeeded);
    
    close gcf
    return % problem, so break out early
    
end %
%%--------END ERROR CONDITION where NSPs cannot be synced (no valid pulses) ---------------------------%
%%-----------------------------------------------------------------------------------------------------%


%%-----------------------------------------------------------------------------------------------------%
%% TO GET HERE PULSES EXIST AND APPEAR TO BE ALIGNABLE... RUN MORE DETAILED CHECKS THEN DO IT %%
%%-----------------------------------------------------------------------------------------------------%

% find first and last pulses present in the sync timeseries again, given potential corrections above
[edgeThresh]    = startingThresh(sync1, sync2);
idxRisingEdge1  = checkEdge(sync1,edgeThresh,'rising');
idxRisingEdge2  = checkEdge(sync2,edgeThresh,'rising');
idxFallingEdge1 = checkEdge(sync1,edgeThresh,'falling');
idxFallingEdge2 = checkEdge(sync2,edgeThresh,'falling');


%- one more round of checks for DC12, if not within a minute of the start and end, try tweaking the threshold
if differentFileSystems ~= 1
    if DC09or12==12,
        tempEdgeThresh = edgeThresh;
        edgeTurnedNan  = 0;
        while tempEdgeThresh<MAX_THRESH_PULSES & edgeTurnedNan==0 & (idxRisingEdge1 > sampRate*60 | idxFallingEdge1 < length(sync1)-sampRate*60 | idxRisingEdge2 > sampRate*60 | idxFallingEdge2 < length(sync2)-sampRate*60), %SJ: used to be sampRates(1), sampRates(1), sampRates(2), sampRates(2)
            tempEdgeThresh  = tempEdgeThresh+500;
            idxRisingEdge1  = checkEdge(sync1,tempEdgeThresh,'rising');
            idxRisingEdge2  = checkEdge(sync2,tempEdgeThresh,'rising');
            idxFallingEdge1 = checkEdge(sync1,tempEdgeThresh,'falling');
            idxFallingEdge2 = checkEdge(sync2,tempEdgeThresh,'falling');
            
            if edgeTurnedNan==0 & any(isnan([idxRisingEdge1 idxRisingEdge2 idxFallingEdge1 idxFallingEdge2])),
                %- its possible that the code above can actually step over the max pulse amplitude...
                %     if so, checkEdge outputs nan... so detect and backstep if that happens
                edgeTurnedNan = 1;
                tempEdgeThresh  = tempEdgeThresh-500;
                idxRisingEdge1  = checkEdge(sync1,tempEdgeThresh,'rising');
                idxRisingEdge2  = checkEdge(sync2,tempEdgeThresh,'rising');
                idxFallingEdge1 = checkEdge(sync1,tempEdgeThresh,'falling');
                idxFallingEdge2 = checkEdge(sync2,tempEdgeThresh,'falling');
            end
        end
        if any(isnan([idxRisingEdge1 idxRisingEdge2 idxFallingEdge1 idxFallingEdge2]))
            fprintf('\n idx rising for falling edge is nan... that will mess stuff up below.  what up?');
            keyboard
        end
        if tempEdgeThresh~=edgeThresh,
            if (idxRisingEdge1 > sampRate*60 | idxFallingEdge1 < length(sync1)-sampRate*60 | idxRisingEdge2 > sampRate*60 | idxFallingEdge2 < length(sync2)-sampRate*60), %SJ: used to be sampRates(1), sampRates(1), sampRates(2), sampRates(2)
                fprintf('\n Uh Oh, DC12 is expected throughout the recording, but it is not detected during either the first or last minute (or neither)');
                keyboard;
                %- not good to get here... something is wrong with DC12.  Last ditch effort is switching to DC09
                DC09or12=9;
                sync1 = sync{3}; sync2 = sync{4};
                [edgeThresh]    = startingThresh(sync1, sync2);
                idxRisingEdge1  = checkEdge(sync1,edgeThresh,'rising');
                idxRisingEdge2  = checkEdge(sync2,edgeThresh,'rising');
                idxFallingEdge1 = checkEdge(sync1,edgeThresh,'falling');
                idxFallingEdge2 = checkEdge(sync2,edgeThresh,'falling');
            else
                fprintf('\n Good news, DC12 wasnt detected in the first or last minute, but changing the edgeThresh fixed it (does this actually happen?)');
                edgeThresh = tempEdgeThresh;
                keyboard
            end
            
        end
    end
end

%- and now a check that rising and falling edges are in the same ballpark... if not try stepping over different thresholds to correct
met = 0;
deltaRiseEdge = (idxRisingEdge1-idxRisingEdge2);
deltaFallEdge = (idxFallingEdge1-idxFallingEdge2);
deltaDuration = (idxFallingEdge1-idxRisingEdge1) - (idxFallingEdge2 - idxRisingEdge2); %- this should be within 100ms... easier to catch with xcorr though, so set to 1000ms here

if abs(deltaRiseEdge) > sampRate*5 | abs(deltaFallEdge) > sampRate*5 | abs(deltaDuration)>sampRate, %SJ: used to be sampRates(1), sampRates(1), sampRates(1)
    % are pulses even in the same ballpark? sometimes, there are aberrant DC09 pulses that happen before the start of the actual train.
    % if this condition is met, maybe it detected those. lets try with higher threshold and see if that gets rid of noise and keeps good pulses
    tempEdgeThresh   = edgeThresh;
    bestEdgeRiseFall = nan;
    bestEdgeDur      = nan;
    %- jw trying a new version... search all thresholds for the best (rather than finding first to match criterion)
    while tempEdgeThresh < MAX_THRESH_PULSES,
        tempEdgeThresh  = tempEdgeThresh + 500;
        idxRisingEdge1  = checkEdge(sync1,tempEdgeThresh,'rising');
        idxRisingEdge2  = checkEdge(sync2,tempEdgeThresh,'rising');
        idxFallingEdge1 = checkEdge(sync1,tempEdgeThresh,'falling');
        idxFallingEdge2 = checkEdge(sync2,tempEdgeThresh,'falling');
        deltaRiseEdgeTry = (idxRisingEdge1-idxRisingEdge2);
        deltaFallEdgeTry = (idxFallingEdge1-idxFallingEdge2);
        deltaDurationTry = (idxFallingEdge1-idxRisingEdge1) - (idxFallingEdge2 - idxRisingEdge2);
        
        %- edges are closer (not clear this is critical, as long as they are within a second (5 sec thresh)
        if abs(deltaRiseEdgeTry)<abs(deltaRiseEdge) | abs(deltaFallEdgeTry)<abs(deltaFallEdge),
            bestEdgeRiseFall = tempEdgeThresh;
            deltaRiseEdge    = deltaRiseEdgeTry;
            deltaFallEdge    = deltaFallEdgeTry;
        end
        %- duration should be less than 100ms... 0 is not uncommon. only check if edges are wihtin acceptable bounds so we can take it if it exists
        %if abs(deltaRiseEdgeTry)<sampRates(1)*5 & abs(deltaFallEdgeTry)<sampRates(1)*5, %#ok<*AND2>
        if abs(deltaDurationTry) < abs(deltaDuration), %#ok<*NOCOL>
            bestEdgeDur   = tempEdgeThresh;
            deltaDuration = deltaDurationTry;
        end
        % end
        
    end
    if differentFileSystems~=1 %only run the following block if the files are from the same systems, and so are expected to start and end
        % at similar times
        if isnan(bestEdgeRiseFall) || abs(deltaRiseEdge) > sampRate*5 || abs(deltaFallEdge) > sampRate*5 , %#ok<*OR2> %- remove duration requirement... too strict here %SJ: used to be sampRates(1), sampRates(1)
            met=2;
            fprintf('\n Uh Oh, pulses not matching up in length even after trying different edgeThreshs (delta: rise=%d, fall=%d, duration=%d', deltaRiseEdge, deltaFallEdge, deltaDuration);
            keyboard;
            %- does this happen (yes, NIH076 session 7/16
        else
            met=1;
            fprintf('\n Good news, sync signals had very different start/end pulses, or different lengths, but changing edgeTresh fixed it');
            if ~isnan(bestEdgeDur) edgeThresh = bestEdgeDur;   %- improving duration is the best
            else                   edgeThresh = bestEdgeRiseFall; %- improving start/end pulse match is good
            end
        end
    else
        % cant make any assumptions about how distant we expect the file
        % starts to be, since this could be NK/BR match
        
        if ~isnan(bestEdgeDur); edgeThresh = bestEdgeDur;   %- improving duration is the best
        end
    end
end


%- make sure the indices are up to date with edgeThresh
idxRisingEdge1  = checkEdge(sync1,edgeThresh,'rising');
idxRisingEdge2  = checkEdge(sync2,edgeThresh,'rising');
idxFallingEdge1 = checkEdge(sync1,edgeThresh,'falling');
idxFallingEdge2 = checkEdge(sync2,edgeThresh,'falling');
deltaRiseEdge = (idxRisingEdge1-idxRisingEdge2);
deltaFallEdge = (idxFallingEdge1-idxFallingEdge2);
deltaDuration = (idxFallingEdge1-idxRisingEdge1) - (idxFallingEdge2 - idxRisingEdge2); %- this should be within 100ms... easier to catch with xcorr though, so set to 1000ms here


%% validate that the segments of data corresponding to the first pulse and last pulse in the two nsps is actually the same, update lags from XCORR if not
%
% do an xcorr and make sure that the lag on which the first pulses have been identified is really the best.
% previously there was a cutoff on the correlation coefficient, but sometimes the two NSPs
% (in TRC files especially) seem to digitize incoming pulses differently. the first barrage of pulses appears as one pulse in
% NSP2. so the corrcoef will tell us its not the same segment. but it is. so instead, use this xcorr to determine thats the ideal lag

%- taking 1 minute of data, b/c this is enough to encompass multiple pulses of either DC09 or DC12
segLength_xcorr = round(60*5*sampRate); %- 5min %SJ: used to be sampRates(1)
lagXlim = [-1 1]*20*sampRate; %- +/- 20 sec %SJ: used to be sampRates(1)

if length(sync1) < segLength_xcorr
    segLength_xcorr = round(segLength_xcorr/2);
    if length(sync1) < segLength_xcorr;
        fprintf('\nfile is super short - are you sure we should align this?\n')
        keyboard
    end
end

%SJ: Adding in case where the rising edge is within 5 min (segLength_xcorr) of the falling edge (inly for micro now)
% So we are going to redefine segLength_xcorr as the length from idkRisingEdge2 to the end of the data
% Could probably add in the same for sync1 but I haven't discovered a case of this yet
if isMicro && idxRisingEdge2+segLength_xcorr > numel(sync2)
    keyboard %SJ: verify this works, so far have only seen it in align micro file 160618_1528 with noreref 160618_1530
    segLength_xcorr = numel(sync2) - idxRisingEdge2;
end

%- check first pulse
[acor,lag] = xcorr(double(sync2(idxRisingEdge2:idxRisingEdge2+segLength_xcorr)),double(sync1(idxRisingEdge1:idxRisingEdge1+segLength_xcorr)));
[~,I] = max(abs(acor));
lagDiff = lag(I);
if idxRisingEdge2+lagDiff>0,                   %- assume the fix applies to the 2nd nsp, unless that violates length
    idxRisingEdge1mod = idxRisingEdge1;
    idxRisingEdge2mod = idxRisingEdge2+lagDiff;
else
    idxRisingEdge1mod = idxRisingEdge1-lagDiff;
    idxRisingEdge2mod = idxRisingEdge2;
end

%- check last pulse (Falling)
[acorF,lagF] = xcorr(double(sync2(idxFallingEdge2-segLength_xcorr:idxFallingEdge2)),double(sync1(idxFallingEdge1-segLength_xcorr:idxFallingEdge1)));
[~,IF] = max(abs(acorF));
lagDiffF = lagF(IF);
if idxFallingEdge2+lagDiffF < originalSig2Length, %- assume the fix applies to the 2nd nsp, unless that violates length
    idxFallingEdge1mod = idxFallingEdge1;
    idxFallingEdge2mod = idxFallingEdge2+lagDiffF;
else
    idxFallingEdge1mod = idxFallingEdge1-lagDiffF;
    idxFallingEdge2mod = idxFallingEdge2;
end

%- any better than the pulse detect method?
deltaRiseEdgeMod = (idxRisingEdge1mod -idxRisingEdge2mod);
deltaFallEdgeMod = (idxFallingEdge1mod-idxFallingEdge2mod);
deltaDurationMod = (idxFallingEdge1mod-idxRisingEdge1mod) - (idxFallingEdge2mod - idxRisingEdge2mod);

usedXcorrCorrection = 0;
%if deltaDurationMod < deltaDuration, %- if duration better take it.. but maybe that doesn't make sense
if lagDiffF~=0 | lagDiff~=0, %- always use the updated estaimte from xcorr (if anything changes)
    % its better, so use it
    fprintf('\n XCORR FIX: xcorr tweaked rising edge or falling edge estimate. \n new vs old delta: rising edge (%d vs %d), falling edge (%d vs %d), duration (%d vs %d)',deltaDurationMod,deltaDuration,deltaRiseEdgeMod,deltaRiseEdge,deltaFallEdgeMod,deltaFallEdge);
    
    if 0 & deltaDurationMod < deltaDuration & (deltaRiseEdgeMod > deltaRiseEdge | deltaFallEdgeMod > deltaFallEdge),
        fprintf('\n XCORR FIX: duration improves (%d vs %d) but rising edge (%d vs %d) or falling edge (%d vs %d) gets worse? \nprobably ok but check...',deltaDurationMod,deltaDuration,deltaRiseEdgeMod,deltaRiseEdge,deltaFallEdgeMod,deltaFallEdge);
        %keyboard;
    end
    
    fprintf('\n xcorr tweaks on pulse indices improves duration match (%d vs %d), so take it.', deltaDurationMod, deltaDuration);
    idxRisingEdge1  = idxRisingEdge1mod;
    idxRisingEdge2  = idxRisingEdge2mod;
    idxFallingEdge1 = idxFallingEdge1mod;
    idxFallingEdge2 = idxFallingEdge2mod;
    deltaRiseEdge = (idxRisingEdge1-idxRisingEdge2);
    deltaFallEdge = (idxFallingEdge1-idxFallingEdge2);
    deltaDuration = (idxFallingEdge1-idxRisingEdge1) - (idxFallingEdge2 - idxRisingEdge2); %- this should be within 100ms... easier to catch with xcorr though, so set to 1000ms here
    
    usedXcorrCorrection = 1;
    
    %- update met flag if this improves something above
    if abs(deltaRiseEdgeMod) < sampRate*5  &  abs(deltaFallEdgeMod) < sampRate*5  &  abs(deltaDuration)<sampRate, %SJ: used to be sampRates(1), sampRates(1), sampRates(1)
        if met==2,
            fprintf('\n Update: using xcorr to update rise/fall indices now satisfies delta criteria');
        end
        met=1;
    end
end

%- save the final values
transforms{1,4}(2:5) = [idxRisingEdge1 idxFallingEdge1 idxFallingEdge1-idxRisingEdge1 size(signal1,2)]; % {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLength
transforms{1,5}(2:5) = [idxRisingEdge2 idxFallingEdge2 idxFallingEdge2-idxRisingEdge2 size(signal2,2)]; % {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLength
transforms{1,6}      = [lagDiff lagDiffF usedXcorrCorrection]; % {6} = XCORR correction on indices: rising diff, falling diff, used correction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate figure summarizing alignment
fS = 12;
xTimeMin_nsp1 = [1:max([length(sync{1}) length(sync{3})])]/(sampRate*60);  %- time in minutes %SJ: used to be sampRates(1)
xTimeMin_nsp2 = [1:max([length(sync{2}) length(sync{4})])]/(sampRate*60);  %- time in minutes %SJ: used to be sampRates(2)
yLimits = [-750 MAX_THRESH_PULSES+100];
%yLimits = [min([0 cell2mat(sync)])-500 max([maxThresh cell2mat(sync)])+500];



%- visualize sync signals
figure(999); clf; set(gcf,'color','w','name','NSP ALIGMENT RESULTS','units','normalized','outerposition',[0 0 1 1]);
subplot(421)
plot(xTimeMin_nsp1(1:length(sync{1})),sync{1});
set(gca,'fontsize',fS,'box','off','tickdir','out')
ylim([yLimits(1) yLimits(2)]);
if DC09or12==12; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
    set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    hold on; plot(idxRisingEdge1/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(1)
    hold on; plot(idxFallingEdge1/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(1)
end
if exist('optionalOverwriteLabels')
    title(sprintf('<%s DC 12>',topType))
else
    title(['<NSP1 DC 12>'])
end
ylabel('DC Signal (mV)')

subplot(423)
plot(xTimeMin_nsp1(1:length(sync{3})),sync{3})
set(gca,'fontsize',fS,'box','off','tickdir','out')
ylim([yLimits(1) yLimits(2)]);
if DC09or12==9; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
    set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    hold on; plot(idxRisingEdge1/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(1)
    hold on; plot(idxFallingEdge1/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(1)
end
xlabel('Time (min)');
if exist('optionalOverwriteLabels')
    title(sprintf('<%s DC 09>',topType))
else
    title(['<NSP1 DC 09>'])
end
ylabel('DC Signal (mV)')


%- NSP2 sync signals
subplot(422); cla
plot(xTimeMin_nsp2(1:length(sync{2})),sync{2})
set(gca,'fontsize',fS,'box','off','tickdir','out')
ylim([yLimits(1) yLimits(2)]);
if DC09or12==12; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
    set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    hold on; plot(idxRisingEdge2/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(2)
    hold on; plot(idxFallingEdge2/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(2)
end
title(sprintf('<%s DC 12>',bottomType))

subplot(424)
plot(xTimeMin_nsp2(1:length(sync{4})),sync{4})
set(gca,'fontsize',fS,'box','off','tickdir','out')
ylim([yLimits(1) yLimits(2)]);
if DC09or12==9; hold on; limits = xlim; horizontal = line([limits(1) limits(2)],[edgeThresh edgeThresh]);
    set(horizontal,'Color','r','linewidth',0.8,'LineStyle','--');
    hold on; plot(idxRisingEdge2/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(2)
    hold on; plot(idxFallingEdge2/(sampRate*60), edgeThresh, 'co','MarkerSize',18,'linewidth',4); %SJ: used to be sampRates(2)
end
xlabel('Time (min)');
title(sprintf('<%s DC 09>',bottomType))

[ax, h] = suplabel(sprintf('%s\nUsed DC%02d for alignment',sessTitle,DC09or12),'t');
set(h,'FontSize',fS+4)


%- first pulse xcorr
subplot(425)
plot(lag,acor); axis tight;
set(gca,'fontsize',fS,'box','off','tickdir','out', 'xlim',lagXlim)
xlabel('lag (samples)'); ylabel('corr'); h = title(sprintf('xcorr of timeseries around first pulses. max lag= %d samp',lagDiff));

subplot(426);
if  DC09or12==12
    preSamp = min([idxRisingEdge1 idxRisingEdge2 201])-1;
    xTime = -preSamp:1200; %- show up to 200 ms before and 500 ms after the first rising edge
else
    preSamp = min([idxRisingEdge1 idxRisingEdge2 1501])-1;
    xTime = -preSamp:11000; %- show up to 200 ms before and 500 ms after the first rising edge
end
plot(xTime*1000/sampRate,double(sync1(idxRisingEdge1+xTime)),'r:','linewidth',6) %SJ: used to be sampRates(1)
hold on
plot(xTime*1000/sampRate,double(sync2(idxRisingEdge2+xTime)),'k','linewidth',2) %SJ: used to be sampRates(1)
axis tight
ylim([yLimits(1) yLimits(2)]);
set(gca,'fontsize',fS,'box','off','tickdir','out')
xlabel('ms'); h = title('zoomed in timeseries around first pulse');

legend(sprintf('Red=%s',topType),sprintf('Black=%s',bottomType));

%- last pulse xcorr
subplot(427)
plot(lagF,acorF); axis tight;
set(gca,'fontsize',fS,'box','off','tickdir','out', 'xlim',lagXlim)
xlabel('lag (samples)'); ylabel('corr'); h = title(sprintf('xcorr of timeseries around last pulses. max lag= %d samp',lagDiffF));

subplot(428);
if  DC09or12==12
    postSamp = min([length(sync1)-idxFallingEdge1 length(sync2)-idxFallingEdge2 201])-1;
    xTime = -1200:postSamp; %- show up to 500 ms before and up to 200 ms after the last rising edge
else
    postSamp = min([length(sync1)-idxFallingEdge1 length(sync2)-idxFallingEdge2 1501])-1;
    xTime = -11000:postSamp; %- show up to 500 ms before and up to 200 ms after the last rising edge
end
plot(xTime*1000/sampRate,double(sync1(idxFallingEdge1+xTime)),'r:','linewidth',6) %SJ: used to be sampRates(1)
hold on
plot(xTime*1000/sampRate,double(sync2(idxFallingEdge2+xTime)),'k','linewidth',2) %SJ: used to be sampRates(1)
axis tight
ylim([yLimits(1) yLimits(2)]);
set(gca,'fontsize',fS,'box','off','tickdir','out')
xlabel('ms'); h = title('zoomed in timeseries around last pulse');

legend(sprintf('Red=%s',topType),sprintf('Black=%s',bottomType));

%- save the fig
if ~isempty(overwriteLabels)
    fig2pngSimple(gcf,fullfile(writeFigsHere,sprintf('NSP_alignment_plot_%s.png',suffixChans)));
else
    fig2pngSimple(gcf,fullfile(writeFigsHere,'NSP_alignment_plot.png'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%------------ BACK TO CHECKING STUFF -------------------%%%

%- final check (comes after figure is made so user can see whats up)... is duration super tight? (<100ms)
if met==2  |   abs(deltaDuration)>sampRate/10, %SJ: used to be sampRates(1)
    fprintf('\n Index difference for first pulse detected = %0.2f',abs(idxRisingEdge1-idxRisingEdge2))
    fprintf('\n Index difference for firstlast pulse detected = %0.2f',abs(idxFallingEdge1-idxFallingEdge2))
    fprintf('\n Different in NSP segment duration = %0.2f',abs(deltaDuration))
    
    [resp] = (input(sprintf('\nShould we go ahead and align using these pulses? I.e.:\n1) Are they valid/intact pulses? (row used for alignment) \n2) Do they look like an identical timeseries (potentially with a shift)? \n3) Are the x indices of the rise and fall of the first pulse matched (bottom right)? (OK for whats b/w the rise and fall to be messed up)\nY or N? \n If N, data will still be written out but not aligned\n '),'s'));
    %resp = 'N'; %forcing
    if (resp=='Y') || (resp=='y')
    elseif (resp=='N') || (resp=='n')
        
        close gcf

        if strcmp(systemType,'CV')
            if differentFileSystems==1
                signal = signal2; % Why do we do this? SJ
                trimEnd1=0; trimEnd2=0;
            else
                signal = cat(1,signal1,signal2);
                trimEnd1 = 0; trimEnd2 = 0;
            end
        else
            signal = cell(1,2); % force them to be equal sizes anyway (necessary for step 5)
            minSize = min(size(signal2,2),size(signal1,2));
            trimEnd1 = size(signal1,2) - minSize; 
            trimEnd2 = size(signal2,2) - minSize;
            if trimEnd1 < 0
               fprintf('%s\n','ERROR!!! Why is the trim end already negative??? Ask SJ!');
               keyboard
            elseif trimEnd1 > 0 %Going to trim signal1 - is that okay?
                if trimEnd2 ~= 0
                    fprintf('%s\n','ERROR!!! Why are we trying to trim both??? Ask SJ!');
                    keyboard
                end
                if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
                    signal1 = signal1(:,1:minSize); % Trimming end of signal1
                else %Not okay, need to pad signal2
                    trimEnd2 = -1*trimEnd1;
                    trimEnd1 = 0;
                    signal2 = cat(2,signal2,NaN(size(signal2,1),abs(trimEnd2))); %adding NaN values to end
                end
            elseif trimEnd2 > 0 % Trim the end of signal2
                signal2 = signal2(:,1:minSize);
            elseif trimEnd1 == 0 && trimEnd2 == 0
                % Not trimming enything
            else
                fprintf('%s\n','ERROR!!! How do we get here??? Ask SJ!');
                keyboard
            end          
            signal{1,1} = signal1; 
            signal{1,2} = signal2;
        end
        
        transforms{1,1} = [0 trimEnd1 0];  %- NSP1: trim start, trim end, add drift
        transforms{1,2} = [0 trimEnd2 0];  %- NSP1: trim start, trim end, add drift

        %- update final lengths before premature return
        transforms{1,4}(5) = [size(signal1,2)]; % {5} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
        transforms{1,5}(5) = [size(signal2,2)]; % {6} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
        
        syncPulse      = 1;
        syncStartIndex = 1;
        errorStr       = 4;
        
        [done] = writeAlignmentSummary(transforms, errorStr, writeFigsHere, sessTitle, keepFirstSignal, isMicro, topType, bottomType, suffixChans, zeroNeeded);

        close gcf
        return  %- break out early because of error
        
    end
end




%% CHECK TO SEE IF GOOD TO GO  --- JW made this strict... needs to be at most 1 sample drift (physically synched NSPs should always be 0)
% there were pulses. so, check if these NSPs were successfully synced
% conditions: offset of starting pulses, combined with stretch/compression needed, should be less than 5 ms
% i.e., if file 1 started 3 ms later, but file 1s first-last pulse duration was also 3 ms shorter, then total = 0 ms, file is good.
%   but if file 1 started 3 ms later, but file 1 was also longer by 3 ms, then we'll consider that a 6 ms offset. no good

deltaRiseEdge = (idxRisingEdge1-idxRisingEdge2);
deltaDuration = (idxFallingEdge1-idxRisingEdge1) - (idxFallingEdge2 - idxRisingEdge2);
%if abs(deltaRiseEdge + deltaDuration) < minGap  %% OLD WAY, too lenien
if abs(deltaDuration)<=1,  % SYNCHRONIZED, but could have an offset or uneven file lengths, so fix that and return
    
    %- first shift sig1 or sig2 so rising edge is perfectly aligned
    trimStart1 = 0;  trimStart2 = 0; trimEnd1 = 0; trimEnd2 = 0;
    %SJ edits:
    if deltaRiseEdge > 0, %rise2 after rise1, going to trim start of signal1 - is that okay?
        if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
            signal1     = signal1(:,deltaRiseEdge+1:end);
            trimStart1  = deltaRiseEdge;
        else %Not okay, need to pad signal2
            signal2     = cat(2,NaN(size(signal2,1),abs(deltaRiseEdge)),signal2); %adding NaN values to beginning
            trimStart2 = -deltaRiseEdge;
        end
            
    elseif deltaRiseEdge < 0, %rise1 after rise2, going to trim start of signal2
        signal2     = signal2(:,-deltaRiseEdge+1:end);
        trimStart2  = -deltaRiseEdge;
    else
        % Not doing anything
    end
    %- then trim so they are exactly the same length
    minSize  = min(size(signal1,2),size(signal2,2));
    trimEnd1 = size(signal1,2) - minSize;
    trimEnd2 = size(signal2,2) - minSize;
    if trimEnd1 < 0
       fprintf('%s\n','ERROR!!! Why is the trim end already negative??? Ask SJ!');
       keyboard
    elseif trimEnd1 > 0 %Going to trim signal1 - is that okay?
        if trimEnd2 ~= 0
            fprintf('%s\n','ERROR!!! Why are we trying to trim both??? Ask SJ!');
            keyboard
        end
        if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
            signal1 = signal1(:,1:minSize); % Trimming end of signal1
        else %Not okay, need to pad end of signal2
            trimEnd2 = -1*trimEnd1;
            trimEnd1 = 0;
            signal2 = cat(2,signal2,NaN(size(signal2,1),abs(trimEnd2))); %adding NaN values to end
        end
    elseif trimEnd2 > 0 % Trim the end of signal2
        signal2 = signal2(:,1:minSize); 
    elseif trimEnd1 == 0 && trimEnd2 == 0
        % Not trimming enything
    else
        fprintf('%s\n','ERROR!!! How do we get here??? Ask SJ!');
        keyboard
    end       
    
    % trim end samples so that they are exactly the same size
    if strcmp(systemType,'CV')
        signal      = cat(1,signal1,signal2);
    else
        signal      = cell(1,2);
        signal{1,1} = signal1;
        signal{1,2} = signal2;
    end
    transforms{1,1} = [trimStart1 trimEnd1 0];  %- NSP1: trim start, trim end, add drift
    transforms{1,2} = [trimStart2 trimEnd2 0];  %- NSP2: trim start, trim end, add drift
    
    %- update final lengths
    transforms{1,4}(5) = [size(signal1,2)]; % {5} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
    transforms{1,5}(5) = [size(signal2,2)]; % {6} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
    
    syncPulse      = 1; %- use first NSP, but both are the same so doesn't really matter
    syncStartIndex = 1+trimStart1;  %-  in split_BR, the analog channels from NSP 1 should be offset by addSignal before extracting pulses
    %errorStr       = 0; SJ: Not super clear why we need this. In the align summary, this is used to say they are already aligned, which is not true...
    if DC09or12==9
        errorStr = 3;
    elseif DC09or12==12
        errorStr = 2;
    end
    FigAndAbort    = 2; %- make fig, but everything is good
    
    [done] = writeAlignmentSummary(transforms, errorStr, writeFigsHere, sessTitle, keepFirstSignal, isMicro, topType, bottomType, suffixChans, zeroNeeded);

    close gcf
    
    
    
    %SJ: Adding this in (same thing as in next section)
    if strcmpi(systemType,'BR') || (isMicro == 1 && ~strcmpi(systemType,'CV')) % Only do this if aligning BR or a micro session
        triggers1 = get_triggers(sync1',sampRate);
        triggers2 = get_triggers(sync2',sampRate);
        threshMS = 10; mywin = 200;
        [pulses1_use, pulses2_use] = pulsealign(triggers1{1,1}, triggers2{1,1}, sampRate, threshMS, mywin, 0, 0, 1);
        ms_field = 'mstime';
        [alignInfo] = logalign_microVsEcog({pulses1_use},{pulses2_use},ms_field);

        %keyboard %check path below
        writeOutAlignInfo = fullfile(writeFigsHere,'alignmentStats.txt'); % -SJ
        fid  = fopen(writeOutAlignInfo,'w');
        if fid==-1; fprintf('cant open file..\n'); keyboard; end
        fprintf(fid,'NumPoints=%0.4f\n',alignInfo.reg_numPointFit);
        fprintf(fid,'Intercept=%0.4f\n',alignInfo.reg_intercept);
        fprintf(fid,'Slope=%0.4f\n',alignInfo.reg_slope);
        fprintf(fid,'R^2=%0.4f\n',alignInfo.reg_Rsquare);
        fprintf(fid,'MaxDev=%0.4f\n',alignInfo.reg_maxDev);
        fprintf(fid,'MedianDev=%0.4f\n',alignInfo.reg_medianDev);
        fprintf(fid,'PulsesUsed=%d\n',DC09or12);
        fclose(fid);
    end
    
    return %- already done, break out early becaus no drift correction required
end




%-----------------------------------------------------------------------------------------------------------%
%---- DRIFT CORRECTION... only make it here if deltaPulseOnset - deltaPulseDurations > minGap  -------------%
%-----------------------------------------------------------------------------------------------------------%

%% new drift correction... keep the whole time series. function figures out longer/shorter
[signal1Fixed, signal2Fixed, whichNSPunstretched,  unstretchedOffset, transformFixed] = correctDriftFull(signal1, signal2, idxRisingEdge1, idxRisingEdge2, idxFallingEdge1, idxFallingEdge2, keepFirstSignal);


% SJ: Now write out alignmentStats.txt
% So essentially do the same thing, but this time with the sync pulses instead of the signals:
%
if strcmpi(systemType,'BR') || isMicro == 1 % Only do this if aligning BR or a micro session
    [sync1Fixed, sync2Fixed, whichNSPunstretched_sync,  unstretchedOffset_sync, transformFixed_sync] = correctDriftFull(sync1, sync2, idxRisingEdge1, idxRisingEdge2, idxFallingEdge1, idxFallingEdge2, keepFirstSignal);
    
    if ~isequal(transformFixed,transformFixed_sync) || ~isequal(whichNSPunstretched,whichNSPunstretched_sync) || ~isequal(unstretchedOffset,unstretchedOffset_sync)
        fprintf('%s\n','ERROR!!! Something using signal vs sync (pulses) is not the same!!! WHY?? Ask SJ!');
        keyboard
    end
    
    if isMicro ~= 1
        % If ecog NSP-NSP alignment, then we want to output the _transformed files    
        keyboard % SJ: make sure these paths are correct
        if contains(writeFigsHere,'STIM_MAP')
            noreref_path = strrep(writeFigsHere,['raw' filesep 'STIM_MAP'],'eeg,noreref');
        else
            noreref_path = strrep(writeFigsHere,'raw','eeg,noreref');
        end
        split_filename = fullfile(noreref_path,['DC' sprintf('%02d',DC09or12) '_nsp1_transformed']);
        %split_filename = sprintf('%s/DC%02d_nsp1_transformed',noreref_path,9);
        [fchan,msg] = fopen(split_filename,'w','l');
        assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
        fwrite(fchan,sync1Fixed,'int16');
        fclose(fchan);
        
        split_filename = fullfile(noreref_path,['DC' sprintf('%02d',DC09or12) '_nsp2_transformed']);
        [fchan,msg] = fopen(split_filename,'w','l');
        assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
        fwrite(fchan,sync2Fixed,'int16');
        fclose(fchan);
    end
    
    if ~strcmpi(systemType,'CV') %SJ: only doing this is not a CV session (should the ~CV requirement be at the top of this if statement instead, like in the previous pulsealign call...?) I don't think it matters, there's just that extra sync check
        triggers1 = get_triggers(sync1Fixed',sampRate);
        triggers2 = get_triggers(sync2Fixed',sampRate);
        threshMS = 10; mywin = 200;
        [pulses1_use, pulses2_use] = pulsealign(triggers1{1,1}, triggers2{1,1}, sampRate, threshMS, mywin, 0, 0, 1);
        ms_field = 'mstime';
        [alignInfo] = logalign_microVsEcog({pulses1_use},{pulses2_use},ms_field);

        %keyboard %check path below
        writeOutAlignInfo = fullfile(writeFigsHere,'alignmentStats.txt'); % -SJ
        fid  = fopen(writeOutAlignInfo,'w');
        if fid==-1; fprintf('cant open file..\n'); keyboard; end
        fprintf(fid,'NumPoints=%0.4f\n',alignInfo.reg_numPointFit);
        fprintf(fid,'Intercept=%0.4f\n',alignInfo.reg_intercept);
        fprintf(fid,'Slope=%0.4f\n',alignInfo.reg_slope);
        fprintf(fid,'R^2=%0.4f\n',alignInfo.reg_Rsquare);
        fprintf(fid,'MaxDev=%0.4f\n',alignInfo.reg_maxDev);
        fprintf(fid,'MedianDev=%0.4f\n',alignInfo.reg_medianDev);
        fprintf(fid,'PulsesUsed=%d\n',DC09or12);
        fclose(fid);
    end

end


% write out that in the way split_CV and split_BR prefer it
if strcmp(systemType,'CV')
    signal = cat(1,signal1Fixed,signal2Fixed);
else
    signal  = cell(1,2);
    signal{1,1} = signal1Fixed; signal{1,2} = signal2Fixed;
end
syncPulse      = whichNSPunstretched;              %-syncPulse1 means NSP1 will be used for sync pulses, with an offset
syncStartIndex = unstretchedOffset;

transforms{1,1}     = transformFixed{1,1};    % {1} NSP1: trim start, trim end, add drift
transforms{1,2}     = transformFixed{1,2};    % {2} NSP2: trim start, trim end, add drift
transforms{1,3}     = transformFixed{1,3};    % {3  DRIFT: drift rate, drift added, whichNSP
transforms{1,4}(5)  = [size(signal1Fixed,2)]; % {5} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
transforms{1,5}(5)  = [size(signal2Fixed,2)]; % {6} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh


if DC09or12==9
    errorStr = 3;
elseif DC09or12==12
    errorStr = 2;
end

[done] = writeAlignmentSummary(transforms, errorStr, writeFigsHere, sessTitle, keepFirstSignal, isMicro, topType, bottomType, suffixChans, zeroNeeded);

close gcf
return




end %- end of main function align_nsps



%----------------------------------------------------------------------------------------------------%
%% helper functions
%----------------------------------------------------------------------------------------------------%

%%%--------------------------------------------------------------------------------------------------------------%%%
function [signal1Fixed, signal2Fixed, whichNSPnotmod, unstretchedOffset, transformsFixed] = correctDriftFull(signal1,signal2,idxRE1,idxRE2,idxFE1,idxFE2,keepFirstSignal)
% function correctDriftFull returns the entire time series corrected (not just time between pulses)


%- define lengths of segments by difference in rising and falling edges
segLength(1) = idxFE1-idxRE1;
segLength(2) = idxFE2-idxRE2;
driftRate    = segLength(1)/segLength(2);
driftSeg     = abs(segLength(1)-segLength(2));
[~,iShort]   = min(segLength); %- which segment is shorter? the short one will get data added

% %% debugging, 5/14/19 -  trim the end before applying drift correction
% confirmed. drift before trim is the same as trim before drift
% minDur = min([size(signal1,2) size(signal2,2)]);
% trimEnd1 = size(signal1,2)-minDur;  trimEnd2 = size(signal2,2)-minDur; %- trim is always to remove, and is positive
% signal1 = signal1(:,1:minDur);
% signal2= signal2(:,1:minDur);


%- full file duration stuff
sigLength(1)  = size(signal1,2);
sigLength(2)  = size(signal2,2);

%- define the data to be stretched
% SJ: edit to consider keepFirstSignal
removeDrift = false;
if iShort == 1, %Going to try to stretch signal1 - is this okay?
    if keepFirstSignal ~= 1 %Okay to stretch signal1
        iMod = 1;
        whichNSPnotmod = 2;
    else %Not okay- need to shrink signal2 instead
        iMod = 2;
        whichNSPnotmod = 1;
        removeDrift = true;
    end
else
    iMod = 2;
    whichNSPnotmod = 1;
end

% Be extra sure
if keepFirstSignal == 1 && iMod == 1
    fprintf('%s\n','ERROR!!!!!!! keepFirstSignal is flagged, but modifying Signal1. Ask SJ!');
    keyboard
end

driftFull     = round(driftSeg * sigLength(iMod) / segLength(iMod)); %- extrapolate drift to entire time series

%- shouldn't be too much... 30 is most I've seen for 1hour session
if driftFull>100,
    fprintf('\n CORRECTING DRIFT: heads up, extrapolated drift is >50 samples (%d)',driftFull);
    keyboard;
end

[signal1Fixed, signal2Fixed, ~, transformsFixed, unstretchedOffset] = applyAlignCorrection(signal1, signal2, iMod, driftRate, driftFull, removeDrift, keepFirstSignal,'idxRE',[idxRE1,idxRE2]);



%SJ commented out: (now in applyAlignCorrection.m)
%- prepare to stretch the slower signal
% numChans      = size(sigToMod,1);

% if ~removeDrift % If you are already changing signal2 or you don't care which gets changed, proceed normally
%     if iShort ~= iMod
%         fprintf('%s\n','ERROR!!!!! iShort and iMod are not equal?? Ask SJ!');
%         keyboard
%     end
%     modSignal = zeros(numChans,sigLength(iMod) + driftFull);
%     %- where to add in the extra samples
%     correctionInterval = floor(sigLength(iMod)/driftFull);
%     if correctionInterval<1000,
%         fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
%         keyboard;
%     end
%     %- indices that will get a second copy (steps of correction interval, starting 1/2 way in correction interval
%     idx = zeros(1,driftFull);
%     for ii=1:driftFull
%         %-  put it in the middle of the correction interval because that makes sense
%         %       and so an extra sample isn't added to the end on even drifts
%         idx(ii) = floor(correctionInterval/2)+1 + correctionInterval*(ii-1); % Why do we have the +1??
%     end
%     %- setup template for how extra samples are filled.. just a copy of the immediately preceeding sample
%     repN = ones(1,sigLength(iMod));
%     repN(idx) = 2;
%     %- make the drift correction
%     for jj=1:numChans
%         modSignal(jj,:) = repelem(sigToMod(jj,1:sigLength(iMod)),repN); %- fill gap with a copy of the preceeding sample
%     end
% else %Need to delete from signal2 instead
%     modSignal = sigToMod;
%     % Where to remove samples
%     correctionInterval = floor(sigLength(iMod)/driftFull);
%     if correctionInterval<1000,
%         fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
%         keyboard;
%     end    
%     %- indices that will be removed (steps of correction interval, starting 1/2 way in correction interval
%     idx = zeros(length(driftFull));    
%     for ii=1:driftFull
%         %-  remove from the middle of the correction interval because that makes sense
%         %       and so an extra sample isn't added to the end on even drifts
%         idx(ii) = floor(correctionInterval/2)+1 + correctionInterval*(ii-1); % Why do we have the +1??
%     end
%     modSignal(:,idx) = [];
% end

% %- tweak the rising edge index of the streched signal to account for the stretch
% %if iShort==1, %iMod?
% if ~removeDrift  % Go ahead and add Drift to signal1 or 2
%     if iMod == 1
%         idxRE1shift = idxRE1 + sum(idx<idxRE1); %- any correction intervals before idxRE1 will offset it
%         deltaRE = idxRE1shift-idxRE2;
%         signal1Fixed = modSignal;
%         signal2Fixed = signal2;
%     else
%         idxRE2shift = idxRE2 + sum(idx<idxRE2);
%         deltaRE = idxRE1-idxRE2shift;
%         signal1Fixed = signal1;
%         signal2Fixed = modSignal;
%     end
% else %removing Drift from signal2
%     if iMod ~= 2
%         fprintf('%s\n','ERROR!!!!!!! removeDrift flagged, but trying to delete from Signal1 (Should never happen!). Ask SJ!');
%         keyboard
%     else
%         idxRE2shift = idxRE2 - sum(idx<idxRE2); %subtract the number of idxs that were removed before the rising edge, from the rising edge
%         deltaRE = idxRE1-idxRE2shift; %new delta Rising edge is probably a little larger now
%         signal1Fixed = signal1;
%         signal2Fixed = modSignal;   
%     end
% end
% 
% %- trim the time series so there is no delta rising edge
% trimStart(1) = 0; trimStart(2) = 0;
% if deltaRE>0,
%     %- idx1 greater than idx2, so trim 1 (cant add imaginary stuff to 2 to match) -SJ: well now we can
%     if keepFirstSignal ~= 1 %okay, go ahead and trim from start of signal1
%         signal1Fixed = signal1Fixed(:,deltaRE+1:end);
%         trimStart(1) = deltaRE;
%     else %Not okay, need to pad signal2
%         signal2Fixed = cat(2,NaN(size(signal2Fixed,1),abs(deltaRE)),signal2Fixed);
%         trimStart(2) = -deltaRE;
%     end
% elseif deltaRE<0,
%     signal2Fixed = signal2Fixed(:,abs(deltaRE)+1:end);
%     trimStart(2) = abs(deltaRE);
% end
% 
% %- and last trim the end so they are the same length
% minDur = min([size(signal1Fixed,2) size(signal2Fixed,2)]);
% trimEnd1 = size(signal1Fixed,2)-minDur;  
% trimEnd2 = size(signal2Fixed,2)-minDur;
% 
% if trimEnd1 < 0
%    fprintf('%s\n','ERROR!!! Why is the trim end already negative??? Ask SJ!');
%    keyboard
% elseif trimEnd1 > 0 %Going to trim signal1 - is that okay?
%     if trimEnd2 ~= 0
%         fprintf('%s\n','ERROR!!! Why are we trying to trim both??? Ask SJ!');
%         keyboard
%     end
%     if keepFirstSignal ~= 1 %okay, go ahead and trim from end of signal1
%         signal1Fixed = signal1Fixed(:,1:minDur); % Trimming end of signal1
%     else %Not okay, need to pad end of signal2
%         trimEnd2 = -1*trimEnd1;
%         trimEnd1 = 0;
%         signal2Fixed = cat(2,signal2Fixed,NaN(size(signal2Fixed,1),abs(trimEnd2))); %adding NaN values to end
%     end
% elseif trimEnd2 > 0 % Trim the end of signal2
%     signal2Fixed = signal2Fixed(:,1:minDur);
% elseif trimEnd1 == 0 && trimEnd2 == 0
%     % Not trimming enything
% else
%     fprintf('%s\n','ERROR!!! How do we get here??? Ask SJ!');
%     keyboard
% end
% 
% %- update NSP sync selection (non-stretched NSP should be used for getting sync data, with the offset here)
% % If you have keepFirstSignal, then you will always be modifying only signal2 and whichNSPnotmod will always
% % be 1. Because Sig1 was not modified ever, trimStart(1) will always be 0
% % If you don't have keepFirstSignal, then sig1 or sig2 will be modified, but you will never end up with a
% % negative trimStart(whichNSPnotmod). Do a check just to make sure:
% % (I don't think this matters for micro but make sure it always works for the other cases!)
% 
% if keepFirstSignal == 1 && whichNSPnotmod==1 && trimStart(whichNSPnotmod)==0
%     % Good
% elseif keepFirstSignal ~= 1 && trimStart(whichNSPnotmod) >= 0
%     % Good
% else
%     fprintf('%s\n','ERROR!!!!!! You should never get here!! Ask SJ!');
%     keyboard;
% end
% 
% unstretchedOffset = 1+trimStart(whichNSPnotmod); %- this will be used to offset the DC channels from the unstreched NSP in split_BR
% 
% %Assign output values
% if removeDrift
%     driftAdded = -driftFull;
% else
%     driftAdded = driftFull;
% end
% if iMod == 1
%     driftNSP1 = driftAdded;
%     driftNSP2 = 0;
% else
%     driftNSP1 = 0;
%     driftNSP2 = driftAdded;
% end
% 
% %- update the transform outputs
% transformsFixed{1,1} = [trimStart(1) trimEnd1  driftNSP1]; %- NSP1: trim start, trim end, add drift
% transformsFixed{1,2} = [trimStart(2) trimEnd2  driftNSP2]; %- NSP2: trim start, trim end, add drift
% transformsFixed{1,3} = [ driftRate driftAdded iMod];                %- DRIFT: drift rate, drift added, whichNSP

end



% %%%--------------------------------------------------------------------------------------------------------------%%%
% function [corrected] = correctDriftTrim(data,idxRE1,idxRE2,idxFE1,idxFE2)
% % function correctDriftTrim returns only data between first and last pulse... use correctDriftFull for entire time series
%
%
% %- define lengths of segments by difference in rising and falling edges
% seg1Length = idxFE1-idxRE1;
% seg2Length = idxFE2-idxRE2;
% drift      = abs(seg1Length-seg2Length);
%
% numChans   = size(data,1);
% correctedTrim = zeros(numChans,min(seg1Length,seg2Length)+drift+1); %- why +1?
%
%
% %- where to add in the extra samples
% correctionInterval = floor(min(seg1Length,seg2Length)/drift);
% idx = zeros(length(drift));
% for i=1:drift
%     idx(i) = (correctionInterval+1)+correctionInterval*(i-1);
% end
%
% %- setup template for how extra samples are filled.. just a copy of the immediately preceeding sample
% repN = ones(1,min(seg1Length,seg2Length)+1); %- why +1?
% repN(idx) = 2;
%
%
% %- which file is shorter? that is the one that will get the added samples
% if seg1Length < seg2Length
%     slowFile=1;
%     segLength = seg1Length;
% else
%     slowFile=2;
%     segLength = seg2Length;
% end
%
% %- make the correction
% if slowFile == 1
%     for i=1:numChans
%         correctedTrim(i,:) = repelem(data(i,idxRE1:idxRE1+segLength),repN); %- fill gap with a copy of the preceeding sample
%     end
% elseif slowFile ==2
%     for i=1:numChans
%         correctedTrim(i,:) = repelem(data(i,idxRE2:idxRE2+segLength),repN); %- fill gap with a copy of the preceeding sample
%     end
% end
% end
%



%----------------------------------------------------------------------------------------------------%
%check for rising/falling edges
%----------------------------------------------------------------------------------------------------%
function [index] = checkEdge(data, thresh, edgeType)

if strcmp(edgeType, 'falling')
    for i=length(data)-1:-1:1  % start at length-1 so returned value is within bounds
        %if data(i) > thresh
        if data(i) > thresh  & data(i+1) <= thresh, % jw added the second part so actually looking for transition edge, not just high on last sample
            index = i+1;
            return;
        end
    end
end
if strcmp(edgeType, 'rising')
    for i=2:1:length(data)  % start at 2 so returned value is within bounds
        %if data(i) > thresh
        if data(i) > thresh & data(i-1) <= thresh,    % jw added the second part so actually looking for transition edge, not just high on first sample
            index = i-1;
            return;
        end
    end
end

if exist('index')==0
    index = NaN;
end

end



%- find starting pulse edge
function [threshUse] = startingThresh(sync1, sync2)

%- find the range that doesnt include outliers, maybe 5 to 95 percentile, (looks like need 99th percentile to detect pulses)
%    if that is >1V, split put the thresh in the middle;  if <1V, then probabably just say all zeros?

pct1 = prctile([sync1],[5 10 99 100],2);
pct2 = prctile([sync2],[5 10 99 100],2);
rngInner = pct1(3)-pct1(2); %- try to avoid outliers
rngOuter = pct1(4)-pct1(1); %- if dont see pulses with 99th percentile, maybe here?


%-
if rngInner(1)>1000,
    threshUse = pct1(2)+500; %- start 500 mV above bottom
    if threshUse<pct2(2) | threshUse>pct2(3),
        fprintf('\n threshold for sync1 does not work for sync2');
        %keyboard;
    end
    
elseif rngOuter(1)>1000,
    threshUse = pct1(1)+500;
    if threshUse<pct2(1) | threshUse>pct2(4),
        fprintf('\n threshold for sync1 does not work for sync2');
        %keyboard;
    end
    
else
    threshUse = 6000; %- actual range is 5000, so this should result in no pulses;
end

end %-startingThresh


%----------------------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------------------%
function [signalRE1, signalRE2, signalFE1, signalFE2] = convertFs(RE1, RE2, FE1, FE2, syncFs, signalFs)
signalRE1 = RE1/syncFs * signalFs;
signalRE2 = RE2/syncFs * signalFs;
signalFE1 = FE1/syncFs * signalFs;
signalFE2 = FE2/syncFs * signalFs;
end


%----------------------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------------------%
function [done] = writeAlignmentSummary(transforms, errorStr,writeFigsHere, sessTitle, keepFirstSignal, isMicro, topType, bottomType, suffixChans, zeroNeeded)
%- write the alignment summary file into raw

if ~isempty(suffixChans)
    fid = fopen(fullfileEEG(writeFigsHere,sprintf('NSP_alignment_summary_%s.txt',suffixChans)),'w+');
else
    fid = fopen(fullfileEEG(writeFigsHere,'NSP_alignment_summary.txt'),'w+');
end

fprintf(fid,'Alignment Summary: %s\n\n',sessTitle);

if     errorStr==0
    fprintf(fid,'NSPs Already Aligned\n');
    keyboard % SJ: Obsolete now, should never get here, right?
elseif errorStr==1
    fprintf(fid,'No DC09 or DC12 pulses, No Alignment Done\n');
elseif errorStr==2
    fprintf(fid, 'Aligned Successfully Using DC12\n');
elseif errorStr==3
    fprintf(fid, 'Aligned Successfully Using DC09\n');
elseif errorStr==4
    fprintf(fid, 'Pulses But Not Matching\n');
end

if keepFirstSignal
    keep_msg = ' (Signal1 unchanged)';
else
    keep_msg = '';
end
if isMicro == 1
    micro_msg = ' (MICRO ALIGNMENT)';
else
    micro_msg = '';
end

fprintf(fid,['Signal1=' topType ', Signal2=' bottomType keep_msg micro_msg '\n\n']);

% {1} = NSP1 add/remove start, end, drift,  {2} = NSP2 add/remove start, end, drift
% {3} = DRIFT: drift rate, drift added, whichNSP
% {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
% {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh

fprintf(fid,'\n      [Initial duration,  RisingPulseEdge,  FallingPulseEdge,   DeltaRiseFall,    FinalLength]');
fprintf(fid,'\n%s    %d,            %d,               %d,           %d,            %d', topType, transforms{1,4});
fprintf(fid,'\n%s    %d,            %d,               %d,           %d,            %d', bottomType, transforms{1,5});
fprintf(fid,'\nDiff     %d,              %d,                 %d,             %d,               %d', transforms{1,4}-transforms{1,5});

fprintf(fid,'\n\n      [ZeroAmountBeginning,  ZeroAmountEnd]');
fprintf(fid,'\n%s        %d,               %d', topType, zeroNeeded(1), zeroNeeded(2)); %SJ 
fprintf(fid,'\n%s        %d,               %d', bottomType, zeroNeeded(3), zeroNeeded(4));

fprintf(fid,'\n\n      [xCorrLag Rising,  xCorrLag Falling,      AdjustedEdgesFromXCorr]');
fprintf(fid,'\n            %d,                   %d,                   %d', transforms{1,6});
fprintf(fid,'\n\n      [trim from start,   trim from End,         Add for drift]');
fprintf(fid,'\n%s        %d,                   %d,                   %d', topType, transforms{1,1});
fprintf(fid,'\n%s        %d,                   %d,                   %d', bottomType, transforms{1,2});
fprintf(fid,'\n\n      [Drift Rate,       NumSamplesInserted      WhichSignalGotInsertion]');
fprintf(fid,'\n         %.6f                 %d,                   %d', transforms{1,3});
fclose(fid);

done = 1;

end
