function signalLag = alignPulses(pulses_1,pulses_2,figSaveDir,whichDC,toSave)
%extractUtahEphys
%Written by: Anthony I. Jang
%
%Align DC12 pulses that are recorded in the blackrock and nk machines using
%cross-correlation
% Inputs:
% pulses_1 & pulses_2: DC12 time stamps from the two machines
%
% Output:
% signalLag: lag (in samples) of signal 2 relative to signal 1
% (e.g. if signalLag < 0, s2 came before s1, and vice versa.)
%
% Need to add: time dilation metric so we can stretch/compress one signal
% relative to the other if there's any time dilation

if ~exist(toSave);
    toSave = 0;
end

% Binarize time stamps
pulses_1(pulses_1==0) = []; pulses_2(pulses_2==0) = [];
s1_bin = zeros(max(pulses_1),1); s1_bin(pulses_1) = 1;
s2_bin = zeros(max(pulses_2),1); s2_bin(pulses_2) = 1;
[r,lags] = xcorr(s1_bin,s2_bin);
signalLag = mean(lags(r==max(r)));


% Embed time stamps onto binary time series, so the actual time stamps are
% preserved after trimming the binary time series
s1_ts = s1_bin; s1_ts(s1_ts==1) = find(s1_ts==1);
s2_ts = s2_bin; s2_ts(s2_ts==1) = find(s2_ts==1);

% if signalLag < 0, s2 came before s1, and vice versa.
if signalLag < 0
    % Align time stamps
    s1_al_ts = s1_ts;
    s2_al_ts = s2_ts(abs(signalLag)+1:end);
    % Align binary time series
    s1_al_bin = s1_bin;
    s2_al_bin = s2_bin(abs(signalLag)+1:end);
    % Unaltered time series (for plotting later)
    ts_unaltered = s1_al_ts;
else
    % Align time stamps
    s1_al_ts = s1_ts(abs(signalLag)+1:end);
    s2_al_ts = s2_ts;
    % Align binary time series
    s1_al_bin = s1_bin(abs(signalLag)+1:end);
    s2_al_bin = s2_bin;
    % Unaltered time series (for plotting later)
    ts_unaltered = s2_al_ts;
end

% Recover the actual time stamps
s1_al_ts = s1_al_ts(s1_al_ts~=0);
s2_al_ts = s2_al_ts(s2_al_ts~=0);
ts_unaltered = ts_unaltered(1:min([length(s1_al_bin),length(s2_al_bin)]));
ts_unaltered = ts_unaltered(ts_unaltered~=0);

% AIJ: don't do the time stretch/compression for now.. Want to finish
% pipeline first.
% Note: resample function doesn't work bc there are too many number of
% samples
% Run a regression to see if time was compressed for one signal
% compared to the other. Apply a stretch/compression to signal 2 using the
% regression slope
% % % numAlignedTS = min([length(s1_al_ts),length(s2_al_ts)]);
% % % b = regress(s1_al_ts(1:numAlignedTS), [ones(numAlignedTS,1) s2_al_ts(1:numAlignedTS)]);
% % % % If slope <1, s2 is compressed compared to s1 --> need to stretch (1/slope > 1)
% % % % If slope >1, s2 is stretched compared to s1 --> need to compress (1/slope < 1)
% % % 
% % % s2_adjustedSampleNum = round(length(s2_al_bin) * (1/b(2)));
% % % s2_al_bin_resampled = resample(s2_al_bin,s2_adjustedSampleNum,length(s2_al_bin));



% Save a figure with the alignment information
showBlocks = 10;    % Show how many divisions of the alignment?

if whichDC==12
    pulseTrainStart_i = find(diff(ts_unaltered)>1000)+1;
    show_pulseTrainStart_i = pulseTrainStart_i(floor(pulseTrainStart_i  :  length(pulseTrainStart_i)/showBlocks  :  length(pulseTrainStart_i)-1))';
    show_pulseTrainEnd_i = pulseTrainStart_i(floor(pulseTrainStart_i  :  length(pulseTrainStart_i)/showBlocks  :  length(pulseTrainStart_i)-1)+1)' - 1;
elseif whichDC==9
    numPulses = 5;
    pulseTrainStart_i = floor(100 : length(ts_unaltered)/showBlocks : length(ts_unaltered));
    show_pulseTrainStart_i = pulseTrainStart_i;
    show_pulseTrainEnd_i = pulseTrainStart_i + numPulses;
end


showBuffer = 200;
sampleRate = 1000;
customColorMap = [1,1,1  ;  0,0,1  ;  1,0,0];
masterFig = figure('units','normalized','outerposition',[0,0,1,1],'visible','off');
for i = 1:length(show_pulseTrainStart_i)
    x_i = ts_unaltered(show_pulseTrainStart_i(i)):ts_unaltered(show_pulseTrainEnd_i(i));
    x_i = min(x_i)-showBuffer : max(x_i)+showBuffer;
    
    % Make lines thicker by padding each binary pulse
    toPlot_s1 = s1_al_bin(x_i); s1_i = find(toPlot_s1==1);
    toPlot_s2 = s2_al_bin(x_i); s2_i = find(toPlot_s2==1);
    padding = 5;
    for pp = 1:padding
        toPlot_s1(s1_i+pp) = 1; toPlot_s2(s2_i+pp) = 1;
        toPlot_s1(s1_i-pp) = 1; toPlot_s2(s2_i-pp) = 1;
    end
    
    figureHandle2 = subplot(showBlocks,2,i*2);
    imagesc([toPlot_s1,toPlot_s2.*2]'); hold on; colormap(customColorMap);
    set(gca,'Xtick',[],'YTick',[1,2],'YTickLabel',{'Signal 1','Signal 2'});
    title(sprintf('%.01f minutes',min(x_i)/sampleRate/60));
    if i == length(show_pulseTrainStart_i), xlabel('Time'); end
    setFigFontSize(16,figureHandle2);
end
figureHandle2 = subplot(showBlocks,2,1:2:showBlocks*2); hold on;
title(sprintf('Max r lag = %d',signalLag));
lagPlot = plot(lags,r,'k'); xlabel('lag (samples)'); ylabel('Cross correlation r');
ar = area([signalLag-floor(length(lags)/100),signalLag+floor(length(lags)/100)],[max(get(gca,'yLim')),max(get(gca,'yLim'))],min(get(gca,'yLim')));
ar.FaceAlpha = 0.1; ar.LineStyle = 'none'; ar.FaceColor = 'b';
setFigFontSize(16,figureHandle2);
uistack(lagPlot);
set(masterFig, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');

% Save the figure
if toSave==1
print(masterFig,figSaveDir,'-dpng','-r300');
end


end



function setFigFontSize(fontSize,figureHandle)
% Sets the font size of all text in a currently open figure
set(gca,'FontSize',fontSize);
%     figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontSize',fontSize);
end
