function [br_timeStamp, IPIViolationFlag] = correctSplitNEV(br_timeStamp, whichDC, postProc)

    if isempty(br_timeStamp)
       error('br_timeStamp is empty'); 
    end
    
    if isempty(postProc)
       error('postProc is empty'); 
    end
    
    IPIViolationFlag = 0;

    verbose = 0;
    
    nonClockResetTimeStampThresh = 200;

    samplingRate = postProc.samplingFreq; 
    original_timestamps = postProc.nonJunkTimeStamps;

    new_br_timeStamps = br_timeStamp;
    negative_diffs_indices = find(diff(new_br_timeStamps)<0);

    clockResetIncrement = 1;

    cumulative_timeStampUpToBeginningOfSeg = original_timestamps(1);

    for seg=1:length(original_timestamps)-1

        % first check for a packet loss, in which both timestamps
        % are large (first segment was not after clock reset)
        check1 = (original_timestamps(seg) > nonClockResetTimeStampThresh) && (original_timestamps(seg+1) > nonClockResetTimeStampThresh);
        % check for a packet loss, in which first segment was
        % result of clock reset
        check2 = (original_timestamps(seg) < nonClockResetTimeStampThresh) && (original_timestamps(seg+1) > nonClockResetTimeStampThresh);

        if check1 || check2 % checks for non-clock reset split

           if verbose
               fprintf('cumulative_timeStampUpToBeginningOfSeg = %d\n', cumulative_timeStampUpToBeginningOfSeg);
               fprintf('comparing seg %d to seg %d\n', seg, seg + 1);
               fprintf('seg %d has timestamp %d\n', seg, original_timestamps(seg));
               fprintf('seg %d has %d samples\n', seg, postProc.nonJunkCellLengths(seg));
           end
           
           if original_timestamps(seg) ~= original_timestamps(seg+1) % sometimes the universe gives you a free sample
               cumulative_timeStampUpToBeginningOfSeg = cumulative_timeStampUpToBeginningOfSeg + postProc.nonJunkCellLengths(seg)*(30000/samplingRate) + postProc.samplesAdded(seg);
           end

           continue
            % packet loss, add to timestamps in second segment


        else


           %********************************************************************
           %********************************************************************
           %********************************************************************
%            figure();
%            subplot(1, 2, 1);
% 
%            plot(1:length(new_br_timeStamps), new_br_timeStamps, '.');
% 
%            hold on;
% 
%            title(sprintf('before comparing seg %d and seg %d', seg, seg +1));

%           plot(1:length(new_br_timeStamps), repmat(cumulative_timeStampUpToBeginningOfSeg, length(new_br_timeStamps), 1), 'k-');

           %********************************************************************
           %********************************************************************
           %********************************************************************





           if verbose
               fprintf('clock reset!!\n');
               fprintf('cumulative_timeStampUpToBeginningOfSeg = %d\n', cumulative_timeStampUpToBeginningOfSeg);
               fprintf('comparing seg %d to seg %d\n', seg, seg + 1);
               fprintf('seg %d has timestamp %d\n', seg, original_timestamps(seg));
           end
           
           % where is the pulse right after the clock reset
           firstPulse_afterCurrentClockReset = negative_diffs_indices(clockResetIncrement)+1;

           % we only want to apply the timestamp adjustment to the current clock-reset-segment's pulses
           % so where do this segment's pulses end
            if length(negative_diffs_indices) ~= clockResetIncrement
                % more clock resets ahead, use the next index
                lastPulse_afterCurrentClockReset = negative_diffs_indices(clockResetIncrement + 1);
            else
                % reached the last clock reset, add adjustment to rest of the pulses
                lastPulse_afterCurrentClockReset = length(new_br_timeStamps);
            end






           %********************************************************************
           %********************************************************************
           %********************************************************************  

%           plot(firstPulse_afterCurrentClockReset, new_br_timeStamps(firstPulse_afterCurrentClockReset), 'r*', 'MarkerSize', 20);     

           %********************************************************************
           %********************************************************************
           %********************************************************************






            % best estimate of timeStamp corresponding to end of the current segment
            time_upToNegativeDiff_in30 = cumulative_timeStampUpToBeginningOfSeg + postProc.nonJunkCellLengths(seg)*(30000/samplingRate);

            if verbose
                
                fprintf('timestamp before switch point %d\n', new_br_timeStamps(firstPulse_afterCurrentClockReset-1));
                fprintf('timestamp at switch point %d\n', new_br_timeStamps(firstPulse_afterCurrentClockReset));

                fprintf('time_upToNegativeDiff_in30: %d\n', time_upToNegativeDiff_in30);
                fprintf('difference between ts @ switch_point-1 and time_upToNegativeDiff_in30: %d\n', new_br_timeStamps(firstPulse_afterCurrentClockReset-1) - time_upToNegativeDiff_in30);
            end
            
            if whichDC == 12 || whichDC == 9
                
                secondsOfSamples_pastLastPulse = (time_upToNegativeDiff_in30 - new_br_timeStamps(firstPulse_afterCurrentClockReset-1))/(samplingRate);
                
                if verbose
                    fprintf('samples go on for %0.2f seconds after last pulse in this segment\n', secondsOfSamples_pastLastPulse);
                end
                
                if whichDC == 12
                    warning_cutoff = 16;
                elseif whichDC == 9
                    warning_cutoff = 3;
                end
                
                if secondsOfSamples_pastLastPulse > warning_cutoff
                   fprintf('samples continued on for more than %d seconds after last pulse. Pulses should occur every ~%d seconds. Make sure this is a classic clock reset split before proceeding -- setting IPIViolationFlag = 1\n', secondsOfSamples_pastLastPulse, warning_cutoff); 
                   IPIViolationFlag = 1;
                end
                
            end

            % add the estimated timestamp to this segment pulses
            new_br_timeStamps(firstPulse_afterCurrentClockReset:lastPulse_afterCurrentClockReset) = new_br_timeStamps(firstPulse_afterCurrentClockReset:lastPulse_afterCurrentClockReset)+time_upToNegativeDiff_in30;






           %********************************************************************
           %********************************************************************
           %********************************************************************

%             plot(1:length(new_br_timeStamps), repmat(time_upToNegativeDiff_in30, length(new_br_timeStamps), 1), 'b-');
% 
%             subplot(1, 2, 2);
% 
%             plot(1:length(new_br_timeStamps), new_br_timeStamps, '.');
% 
%             title(sprintf('after comparing seg %d and seg %d', seg, seg +1));
% 
%             hold on;
% 
%             plot(1:length(new_br_timeStamps),  repmat(time_upToNegativeDiff_in30, length(new_br_timeStamps), 1), 'b-');

           %********************************************************************
           %********************************************************************
           %********************************************************************





            %keyboard;        

            cumulative_timeStampUpToBeginningOfSeg = cumulative_timeStampUpToBeginningOfSeg + postProc.nonJunkCellLengths(seg)*(30000/samplingRate);
            clockResetIncrement = clockResetIncrement + 1;

        end

    end

    if any(diff(new_br_timeStamps)<0)
       error('finished the adjustment loop, but negative diffs remain. This should not happen'); 
    end

    br_timeStamp = new_br_timeStamps;
    
end