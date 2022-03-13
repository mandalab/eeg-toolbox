
function [alignStats] = logalign_microVsEcog(beh_ms,eeg_offset,ms_field)

% allocate for stuph
b = zeros(2,length(eeg_offset));
eeg_start_ms = zeros(length(eeg_offset),1);
eeg_stop_ms  = zeros(length(eeg_offset),1);

strWarning = '';    %this will hold the warning messages if any are generated
strStats = ''; %this will hold a copy of the stats dumped to the screen

% loop over beh and eeg sync and get slopes
for f = 1:length(eeg_offset)
    % get slope and offset for each eeg file
    bfix = beh_ms{f}(1);
    [b(:,f),bint,r,rin,stats] = regress(eeg_offset{f}, [ones(length(beh_ms{f}),1) beh_ms{f}-bfix]);
    b(1,f) = b(1,f) - bfix*b(2,f);
    
    % calc max deviation
    act=[ones(length(beh_ms{f}),1) beh_ms{f}]*b(:,f);
    maxdev    = max(abs(act - eeg_offset{f}));
    maxDevThresh = 5;
    numAboveDevThr   = length(find(abs(act - eeg_offset{f})>maxDevThresh));
    pcntAboveDevThr  = 100*numAboveDevThr/length(act);
    numAbove2DevThr  = length(find(abs(act - eeg_offset{f})>maxDevThresh*2));
    pcntAbove2DevThr = 100*numAbove2DevThr/length(act);
    meddev{f} = median(abs(act - eeg_offset{f}));
    rawDevs{f}=act - eeg_offset{f};
    
    %     % save regression stats in string that is dumped at the end of the function
    %     %fprintf('%s:\n', eeg_files{f});
    %     fprintf('\tRegression Slope  = %f\n', b(2,f));
    %     fprintf('\tR^2               = %f\n', stats(1));
    %     fprintf('\tMedian Deviation  = %f ms\n', meddev{f});
    %     fprintf('\tMax Deviation     = %f ms\n', maxdev);
    %     fprintf('\tDeviations>%.1fms  = %.1f%%\n', maxDevThresh,pcntAboveDevThr);
    %     fprintf('\tDeviations>%.1fms = %.1f%%\n', maxDevThresh*2,pcntAbove2DevThr);
    
    
    alignStats = struct;
    alignStats.reg_numPointFit = length(eeg_offset{f});
    alignStats.reg_intercept   = b(1,f);
    alignStats.reg_slope       = b(2,f);
    alignStats.reg_Rsquare     = stats(1);
    alignStats.reg_maxDev      = maxdev;
    alignStats.reg_medianDev   = meddev{f};
    alignStats.deviations = sprintf('\tDeviations>%.1fms  = %.1f%%\n', maxDevThresh,pcntAboveDevThr);
end
end
