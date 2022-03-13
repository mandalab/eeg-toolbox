function [isoMetrics,classResults,isoScore_stability,spkRate_stability,invalidTimeRange] = getChanUnitQuality(sorted_allUnitNum_all,sorted_timeStamp_all,sorted_waveFormValues_all,figSavePath,binFileDir,samplerate)
%% Get the unit quality metrics and save a figure with the channel's unit information
%
% 2/19/2020 - SJ added: if segTime_1/2 == 0, segTime_1/2 = 1; end 
%             because can't 0 index. Was this the right thing to do?

%- figure plotting relies on histogram, a function added in Matlab R2014a
if verLessThan('matlab','8.4'),
    fprintf('\n Plotting code below requires matlab R2014 or later; \n   tweak the code (histogram function) or upgrade matlab to run');
    keyboard
end


%- make some decisions here about what the standard output is.
CALC_ISOSCORE_fromPlexon = 1; %- calculate isolation metrics from the snippet time series saved out by plexon.  This was AJ's original hack of the Joshua method
CALC_ISOSCORE_fromBin    = 1; %- calculate isolation metrics by grabbing a fresh snippet from the .bin file and fresh noise snippets, exactly the way Joshua did it.
CALC_CLASSIFY_btwnUnits  = 0; %- use LDA to compute classification accuracy between population of snippets between all pairs of units.  Classifier accuracy was always high, >95%, even if isolation score is crap.  Dont use it.
CALC_ISOSCORE_stability  = 0; %- calc the isolation score in windows across the recording.  very slow computationally, so skip it.

snippetSecs = 3;
maxWaves = 5000;

% Get the save name for this figure
slashes = [strfind(figSavePath,'/') strfind(figSavePath,'\')];
figSaveName = figSavePath(slashes(end)+1:end);
for k = 1:length(figSaveName), if strcmp(figSaveName(k),'_'), figSaveName(k)='-'; end, end

%- initialize output variables
classResults       = [];
isoScore_stability = [];
spkRate_stability  = [];
invalidTimeRange   = {};


isoMetrics         = struct;
tempStruct         = struct;
tempStruct.sigNoiseR = []; tempStruct.isoScore = []; tempStruct.fnScore = []; tempStruct.fpScore = [];
if CALC_ISOSCORE_fromPlexon,
    isoMetrics.unitVSnoise   = tempStruct; %- plexon spike vs plexon noise, as output by manual sort.  This is the original metric used by AJ and JW for isolationScore
    %isoMetrics.allNoise  = tempStruct; %- dont use this anymore. this was plexon spike vs noise+all other spikes
end
if CALC_ISOSCORE_fromBin,
    isoMetrics.joshua    = tempStruct; %- upsampled spike from .bin vs freshly grabbed noise crossings, where noise is defined based on spike amplitude.
end
tempStruct2 = struct;
tempStruct2.sigNoiseR = {}; tempStruct2.isoScore = {}; tempStruct2.fnScore = {}; tempStruct2.fpScore = {};
isoMetrics.btwnUnits = tempStruct2; %- plexon spike vs spike

makeTemplatePlot = 0;
classifyTrialNum = 10; % ten-fold validation



%%% START THE PROCESS %%%%
% Get rid of invalidated waveforms
sorted_allUnitNum     = sorted_allUnitNum_all(sorted_allUnitNum_all>=0);
sorted_timeStamp      = sorted_timeStamp_all(sorted_allUnitNum_all>=0);
sorted_waveFormValues = sorted_waveFormValues_all(sorted_allUnitNum_all>=0,:);
unitCount             = length(unique(sorted_allUnitNum(sorted_allUnitNum > 0)));

sorted_timeStamp_min  = sorted_timeStamp/1000/60;


if unitCount >= 1,   % Unit exists
    
    % Get spike quality measures
    
    % Need to load the binary for (1) plotting raw time series in spike plot and (2) calculating Joshua's isoscore the "real" way, where snippets are grabbed from bin
    tLoadBin = tic;
    fid     = fopen(binFileDir,'r');
    binData = fread(fid, 'int16');
    fclose(fid);
    fprintf('loadBin[%.1f]',toc(tLoadBin));
    
    
    % All unit numbers (non-noise)
    allUnitNums = unique(sorted_allUnitNum)';
    allUnitNums(allUnitNums <= 0) = [];
    
    % Joshua Isolation metrics
    tStartQuality=tic;
    for ui = 1:length(allUnitNums)
        u = allUnitNums(ui);
        
        
        
        if CALC_ISOSCORE_fromPlexon,
            tStartIso1 = tic;
            
            % Original way in the Zaghloul lab: using 34-sample waveforms where noise is just noise
            spkWaves     = sorted_waveFormValues_all(sorted_allUnitNum_all==u,:);
            nzWaves_orig = sorted_waveFormValues_all(sorted_allUnitNum_all==0,:);
            
            % new way in the Zaghloul lab: using 34-sample waveforms where noise is everything except current unit (closer to Joshua's actual method, but in practice doesn't really make a difference so skip it)
            %nzWaves_allOther = sorted_waveFormValues_all(sorted_allUnitNum_all~=u,:);
            
            % Subsample to a maximum of 'maxWaves'
            [spkWaves_orig,nzWaves_orig]         = subsampleWaves(spkWaves,nzWaves_orig,maxWaves);
            %[spkWaves_allOther,nzWaves_allOther] = subsampleWaves(spkWaves,nzWaves_allOther,maxWaves);
            
            % compute the metrics from PLEXON-saved waveforms
            [isoMetrics.unitVSnoise.sigNoiseR(ui,1),     isoMetrics.unitVSnoise.isoScore(ui,1),     isoMetrics.unitVSnoise.fnScore(ui,1),     isoMetrics.unitVSnoise.fpScore(ui,1),     outClustInfo] = klUnitIsolation_AIJ(spkWaves_orig,    nzWaves_orig);
            %[isoMetrics.allNoise.sigNoiseR(ui,1), isoMetrics.allNoise.isoScore(ui,1), isoMetrics.allNoise.fnScore(ui,1), isoMetrics.allNoise.fpScore(ui,1), outClustInfo] = klUnitIsolation_AIJ(spkWaves_allOther,nzWaves_allOther);
            
            fprintf('isoN[%.1fs]',toc(tStartIso1));
        end
        
        
        if CALC_ISOSCORE_fromBin,
            tStartIso2 = tic;
            
            % Joshua's actual method using 144-sample interpolated waveform
            preSpikeMS  = 0.5;  preSpikeSamples  = round(preSpikeMS*(samplerate/1000));
            postSpikeMS = 1;    postSpikeSamples = round(postSpikeMS*(samplerate/1000));
            bufferMS    = 0.5;  bufferSamples    = round(bufferMS*(samplerate/1000));
            
            % Get the spk waves
            u_timeStamp_samples = floor(sorted_timeStamp(sorted_allUnitNum==u)*(samplerate/1000));
            % Subsample
            if length(u_timeStamp_samples) > maxWaves
                i_rand = randperm(length(u_timeStamp_samples));
                u_timeStamp_samples = sort(u_timeStamp_samples(i_rand(1:maxWaves)));
            end
            % Get rid of timestamps that are out of bounds with the buffer
            if any(u_timeStamp_samples<preSpikeSamples | u_timeStamp_samples > length(binData)-postSpikeSamples)
                fprintf('\nGot rid of some timestamps...(out of bounds)\n');
                u_timeStamp_samples(u_timeStamp_samples<(preSpikeSamples+bufferSamples) | u_timeStamp_samples > length(binData)-(postSpikeSamples+bufferSamples)) = [];
            end
            
            spkWaves_joshua = getCubicSplineWaves(u_timeStamp_samples,binData,samplerate,preSpikeMS,postSpikeMS,bufferMS);
            
            % Set up the noise waves using the spike times
            troughValues = spkWaves_joshua(:,mean(spkWaves_joshua,1)==min(mean(spkWaves_joshua,1)));
            q = quantile(troughValues,0.98);
            nzThresh = mean(troughValues(troughValues>q))/2;
            if isnan(nzThresh), nzThresh = min(troughValues)/2; end %- happens when very few spikes, so 98th pcntile is the online sample
            timeStamp_threshCrossing = find(binData<=nzThresh);
            n_timeStamp_samples = timeStamp_threshCrossing([500;diff(timeStamp_threshCrossing)]>1);
            % Remove noise events within the spike events Noise should not be within 1.5ms from time zero for a spike event
            n_samples = zeros(1,max([n_timeStamp_samples;u_timeStamp_samples])+preSpikeSamples+postSpikeSamples);
            n_samples(n_timeStamp_samples) = 1;
            % Remove u_timeStamp_samples that is out of bounds when there are buffers (pre & postSpikeSamples)
            validCross = ones(1,length(n_samples));
            for ii=-preSpikeSamples:postSpikeSamples, validCross(u_timeStamp_samples+ii)=0; end %- smear the spike crossing in time
            crossKeep = n_samples.*validCross;
            n_timeStamp_samples = find(crossKeep); %- timing of valid crosses, in samples
            % Subsample
            if length(n_timeStamp_samples) > maxWaves
                i_rand = randperm(length(n_timeStamp_samples));
                n_timeStamp_samples = sort(n_timeStamp_samples(i_rand(1:maxWaves)));
            end
            
            nzWaves_joshua = getCubicSplineWaves(n_timeStamp_samples,binData,samplerate,preSpikeMS,postSpikeMS,bufferMS);
            
            [isoMetrics.joshua.sigNoiseR(ui,1),   isoMetrics.joshua.isoScore(ui,1),   isoMetrics.joshua.fnScore(ui,1),   isoMetrics.joshua.fpScore(ui,1),   outClustInfo] = klUnitIsolation_AIJ(spkWaves_joshua,  nzWaves_joshua);
            
            fprintf('isoJ[%.1fs]',toc(tStartIso2));
        end
    end
    
    
    % Get between-unit isolation score
    unitPairs = [];
    if unitCount > 1
        tStartIso3 = tic;
        
        unitPairs = cat(1,unitPairs,nchoosek(allUnitNums,2));
        btwnUnitSNR = []; btwnUnitIso = []; btwnUnitFN = []; btwnUnitFP = [];
        for uPair = 1:size(unitPairs,1)
            u1 = unitPairs(uPair,1);
            u2 = unitPairs(uPair,2);
            % Original way using 34-sample waveforms where noise is just noise
            spkWaves_1 = sorted_waveFormValues_all(sorted_allUnitNum_all==u1,:);
            spkWaves_2 = sorted_waveFormValues_all(sorted_allUnitNum_all==u2,:);
            % Subsample to a maximum of 'maxWaves'
            [spkWaves_1,spkWaves_2] = subsampleWaves(spkWaves_1,spkWaves_2,maxWaves);
            [sigNoiseR_temp,isoScore_temp,fnScore_temp,fpScore_temp,outClustInfo] = klUnitIsolation_AIJ(spkWaves_1,spkWaves_2);
            btwnUnitSNR = cat(1,btwnUnitSNR,[unitPairs(uPair,:),sigNoiseR_temp]);
            btwnUnitIso = cat(1,btwnUnitIso,[unitPairs(uPair,:),isoScore_temp]);
            btwnUnitFN  = cat(1,btwnUnitFN,[unitPairs(uPair,:),fnScore_temp]);
            btwnUnitFP  = cat(1,btwnUnitFP,[unitPairs(uPair,:),fpScore_temp]);
        end
        for uuu = 1:length(allUnitNums)
            %- for this unit's cell entry, only include the btwn comparisons to this unit (used to output all comparisons for same channel)
            iRows = find(unitPairs(:,1)==uuu | unitPairs(:,2)==uuu);  %- rows of the comparison matrix with this unit
            isoMetrics.btwnUnits.sigNoiseR = cat(1,isoMetrics.btwnUnits.sigNoiseR,{btwnUnitSNR(iRows,:)});
            isoMetrics.btwnUnits.isoScore  = cat(1,isoMetrics.btwnUnits.isoScore,{btwnUnitIso(iRows,:)});
            isoMetrics.btwnUnits.fnScore   = cat(1,isoMetrics.btwnUnits.fnScore,{btwnUnitFN(iRows,:)});
            isoMetrics.btwnUnits.fpScore   = cat(1,isoMetrics.btwnUnits.fpScore,{btwnUnitFP(iRows,:)});
        end
        
        fprintf('isoBTW[%.1fs]',toc(tStartIso3));
    else
        isoMetrics.btwnUnits.sigNoiseR = cat(1,isoMetrics.btwnUnits.sigNoiseR,{[nan,nan,nan]});
        isoMetrics.btwnUnits.isoScore  = cat(1,isoMetrics.btwnUnits.isoScore,{[nan,nan,nan]});
        isoMetrics.btwnUnits.fnScore   = cat(1,isoMetrics.btwnUnits.fnScore,{[nan,nan,nan]});
        isoMetrics.btwnUnits.fpScore   = cat(1,isoMetrics.btwnUnits.fpScore,{[nan,nan,nan]});
    end
    fprintf('<total%.0fs>',toc(tStartQuality)); 
    
    
    if CALC_CLASSIFY_btwnUnits,
        % Use LDA to classify different clusters
        % 1. Between spike and noise clusters
        % 2. All pairwise between units
        % after doing this for many units we found that classification accuracy was always >95%, even for crappy isolation.  seems like not a great metric
        fprintf('; Classify '); tStartClassify=tic;
        clusterPairs = [];
        for u = allUnitNums, clusterPairs = cat(1,clusterPairs,[0,u]); end
        if unitCount > 1,    clusterPairs = cat(1,clusterPairs,nchoosek(allUnitNums,2)); end
        %keyboard
        accuracy = []; accuracy_all = {};
        for comparison = 1:size(clusterPairs,1)
            acc_thisComp = [];
            
            thisComp = clusterPairs(comparison,:);
            cluster_1 = sorted_waveFormValues(sorted_allUnitNum == thisComp(1),:); size_1 = size(cluster_1,1);
            cluster_2 = sorted_waveFormValues(sorted_allUnitNum == thisComp(2),:); size_2 = size(cluster_2,1);
            minNumPts = min([size_1,size_2]);
            subSize = floor(minNumPts/10);
            minNumPts = subSize*10;
            for classifyTrial = 1:classifyTrialNum
                tempClust_1 = cluster_1(randperm(size_1,minNumPts),:);
                tempClust_2 = cluster_2(randperm(size_2,minNumPts),:);
                
                testing = [tempClust_1(1:subSize,:);tempClust_2(1:subSize,:)];
                testing_group = double([ones(subSize,1);2*ones(subSize,1)]);
                training = []; group = [];
                training = cat(1,training,tempClust_1(subSize+1:end,:)); group = cat(1,group,ones(minNumPts-subSize,1));
                training = cat(1,training,tempClust_2(subSize+1:end,:)); group = cat(1,group,2*ones(minNumPts-subSize,1));
                classified = classify(testing,training,group);
                accuracy(comparison,classifyTrial) = sum(abs(classified-testing_group)==0)/(subSize*2);
                acc_thisComp(classifyTrial,1) = sum(abs(classified(1:subSize)-testing_group(1:subSize))==0)/subSize;
                acc_thisComp(classifyTrial,2) = sum(abs(classified(subSize+1:end)-testing_group(subSize+1:end))==0)/subSize;
            end
            accuracy_all{comparison} = acc_thisComp;
        end
        classResults = [clusterPairs,mean(accuracy,2)];
        fprintf('%.0fs',toc(tStartClassify)); 
    end
    
    
    % Get spike rate and isolation quality vs. time
    fprintf('; vsTime '); tStartVtime=tic;
    
    %- initialize variables
    
    % Spike rate time series
    windowMinutes_1 = 1;
    wSize_1 = 1000*60*windowMinutes_1;
    wStep_1 = wSize_1/2;
    allWindowStartPts_1 = 1:wStep_1:max(sorted_timeStamp)-wSize_1;
    spkRate_stability = nan(unitCount,length(allWindowStartPts_1),1);
    windowMid_1       = nan(length(allWindowStartPts_1),1);
    
    if CALC_ISOSCORE_stability,
        % Isolation quality time series
        windowMinutes_2 = 10;
        wSize_2 = 1000*60*windowMinutes_2;
        wStep_2 = wSize_2/2;
        allWindowStartPts_2 = 1:wStep_2:max(sorted_timeStamp)-wSize_2;
        windowMid_2         = nan(length(allWindowStartPts_2),1);
        
        isoScore_stability = nan(unitCount,length(allWindowStartPts_2),1);
        accuracy_stability = nan(unitCount,length(allWindowStartPts_2),1);
    end
    
    
    % Spike rate vs. time (used for plot)
    for ui = 1:length(allUnitNums)
        u = allUnitNums(ui);
        u_timeStamp = sorted_timeStamp(sorted_allUnitNum==u);
        for w = 1:length(allWindowStartPts_1)
            thisWindowStart = allWindowStartPts_1(w);
            windowMid_1(w)  = thisWindowStart+wSize_1/2;
            spkRate_stability(ui,w) = sum(u_timeStamp>=thisWindowStart & u_timeStamp < thisWindowStart+wSize_1);
        end
        spkRate_stability(ui,:) = (spkRate_stability(ui,:)/wSize_1) * 1000;  % Spikes per sec
        windowMid_1 = windowMid_1/1000/60;   % Times in minutes
    end
    
    
    % Isolation quality for each time window
    if CALC_ISOSCORE_stability,
        for ui = 1:length(allUnitNums)
            u = allUnitNums(ui);
            u_timeStamp = sorted_timeStamp(sorted_allUnitNum==u);
           
            nzWaves      = sorted_waveFormValues(sorted_allUnitNum==0,:);
            nz_timeStamp = sorted_timeStamp(sorted_allUnitNum==0);
            u_spkWaves   = sorted_waveFormValues(sorted_allUnitNum==u,:);
            
            for w = 1:length(allWindowStartPts_2)
                thisWindowStart = allWindowStartPts_2(w);
                windowMid_2(w)  = thisWindowStart+wSize_2/2;
                windowSpkWaves  = u_spkWaves(u_timeStamp>=thisWindowStart & u_timeStamp < thisWindowStart+wSize_2,:);
                windowNzWaves   = nzWaves(nz_timeStamp>=thisWindowStart & nz_timeStamp < thisWindowStart+wSize_2,:);
                
                % If the spike count is low, don't run these
                if size(windowSpkWaves,1) >= 20
                    
                    % Joshua
                    quickQuality = 1; %- only use max 1000 waveforms when doing time series
                    [SNR2,isoScore_stability(ui,w),fnScore2,fpScore2,outClustInfo2] = klUnitIsolation_AIJ(windowSpkWaves,windowNzWaves,quickQuality);
                    
                    % Hill
                    minNumPts = floor(min(size(windowSpkWaves,1),size(windowNzWaves,1))/10)*10;
                    subSize = minNumPts/10;
                    tempClust_1 = windowSpkWaves(randperm(size(windowSpkWaves,1),minNumPts),:);
                    tempClust_2 = windowNzWaves(randperm(size(windowNzWaves,1),minNumPts),:);
                    
                    testing = [tempClust_1(1:subSize,:);tempClust_2(1:subSize,:)];
                    testing_group = double([ones(subSize,1);2*ones(subSize,1)]);
                    training = []; group = [];
                    training = cat(1,training,tempClust_1(subSize+1:end,:)); group = cat(1,group,ones(minNumPts-subSize,1));
                    training = cat(1,training,tempClust_2(subSize+1:end,:)); group = cat(1,group,2*ones(minNumPts-subSize,1));
                    classified = classify(testing,training,group);
                    accuracy_stability(ui,w) = sum(abs(classified-testing_group)==0)/(subSize*2);
                else
                    isoScore_stability(ui,w) = nan;
                    accuracy_stability(ui,w) = nan;
                end
            end
            windowMid_2 = windowMid_2/1000/60;   % Times in minutes
            
            badTimes = [];
            invalidTimeRange{ui,1} = badTimes;
        end % for u = 1:unitCount
    end
    fprintf('%.0fs; SaveFigs ',toc(tStartVtime)); tSaveFigs=tic;
    
    
    
    % Get plots!
    subplot_rows = 6; subplot_cols = 9;
    colors = {'w','y','g','c','r','b','m','k','m','m','m'}; %- add a few extra magentas in the end in the odd case there are more than 6 sorts
    fontSize = 16;
    
    
    % Make template plot
    if makeTemplatePlot == 1
        figure;
        for ii = 1:subplot_rows*subplot_cols
            subplot(subplot_rows,subplot_cols,ii); hold on;
            text(0.5,0.5,num2str(ii),'fontSize',20);
        end
    end
    % Make one big plot with all the small plots
    masterFig = figure('units','normalized','outerposition',[0,0,1,1],'visible','off');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get mean waveforms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fontSize = 12;
    figureHandle = subplot(subplot_rows,subplot_cols,[1:3,10:12]); hold on;
    set(gca,'color','black'); plots = [];
    ui = 1;
    for unitNum = [0,allUnitNums]
        waveForm = sorted_waveFormValues(sorted_allUnitNum == unitNum,:);
        
        % Calculate mean & standard error
        wf_mean = mean(waveForm,1);
        wf_std  = std(waveForm,1);
        plots(ui) = plot(1:length(wf_mean),wf_mean,'Color',colors{unitNum+1},'LineWidth',3);
        
        xStep = 1:length(wf_mean);
        ii    = length(xStep):-1:1;
        fillX = [xStep xStep(ii)];
        fillY = [wf_mean+wf_std wf_mean(ii)-wf_std(ii)];
        hFill = fill(fillX,fillY,colors{unitNum+1});
        set(hFill,'LineStyle','none','FaceAlpha',0.4);
        
        ui = ui+1;
    end
    set(gca,'XLim',[1,size(waveForm,2)],'XTick',0:15:size(waveForm,2),'XTickLabel',[0:15:size(waveForm,2)]/(samplerate/1000));
    xlabel('Time(ms)'); ylabel('Voltage (uV)');
    legendLabel = {};
    for ui = 1:length(allUnitNums), legendLabel{ui} = sprintf('Unit %d',allUnitNums(ui)); end
    legendLabel = ['Noise',legendLabel];
    legendHandle = legend(plots,legendLabel); legendHandle.TextColor = 'w';
    setFigFontSize(fontSize,figureHandle);
    xRange = get(gca,'XLim'); yRange = get(gca,'YLim');
    text(xRange(2)/4,yRange(2)-yRange(2)/9,figSaveName,'Color','w','FontSize',20);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time vs. peak-valley plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fontSize = 12;
    peakValley = max(sorted_waveFormValues,[],2) - min(sorted_waveFormValues,[],2);
    figureHandle = subplot(subplot_rows,subplot_cols,[4:6,13:15]); hold on; set(gca,'color','black');
    plot(sorted_timeStamp_min,peakValley,'.','Color','w','MarkerSize',3);
    ylabel('Voltage (Peak - Trough)'); xlabel('Time(m)');
    setFigFontSize(fontSize,figureHandle);
    figureHandle = subplot(subplot_rows,subplot_cols,[7:9,16:18]); hold on; set(gca,'color','black');
    setFigFontSize(fontSize,figureHandle);
    for unitNum = [0,allUnitNums]
        x = sorted_timeStamp_min(sorted_allUnitNum == unitNum);
        y = peakValley(sorted_allUnitNum == unitNum);
        plot(x,y,'.','Color',colors{unitNum+1},'MarkerSize',3);
    end
    ylabel('Voltage (Peak - Trough)'); xlabel('Time(m)');
    % Get PC1 vs. PC2 plot
    [COEFF,SCORE] = pca(sorted_waveFormValues);
    figureHandle = subplot(subplot_rows,subplot_cols,[22:24,31:33]); hold on; set(gca,'color','black');
    plot(SCORE(:,1),SCORE(:,2),'.','Color','w','MarkerSize',3);
    xlabel('PC1'); ylabel('PC2');
    set(gca,'XTick',[],'YTick',[]);
    setFigFontSize(fontSize,figureHandle);
    figureHandle = subplot(subplot_rows,subplot_cols,[25:27,34:36]); hold on; set(gca,'color','black');
    for unitNum = [0,allUnitNums]
        plot(SCORE(sorted_allUnitNum == unitNum,1),SCORE(sorted_allUnitNum == unitNum,2),'.','Color',colors{unitNum+1},'MarkerSize',3);
    end
    xlabel('PC1'); ylabel('PC2');
    set(gca,'XTick',[],'YTick',[]);
    setFigFontSize(fontSize,figureHandle);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get pdf and cdf of interspike intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxIsiToShow  = 500;
    isiThresh     = 2; %- 2ms abs refract
    cdfSubplotNum = 1;
    fontSize = 12;
    for unitNum = allUnitNums
        if unitNum ~= 0 && cdfSubplotNum <=3
            isi = diff(sorted_timeStamp(sorted_allUnitNum == unitNum,1));
            %isi(isi > maxIsiToShow) = []; %- JW commented out: blanking out long ISI will mess up the calculation of "fraction <2ms"
            isiBins = 0:1:maxIsiToShow;
            [oc, isiBins] = hist(isi,isiBins);
            isiCDF = [];
            for ii=1:length(oc),
                isiCDF(ii) = sum(oc(1:ii))./sum(oc);
            end
            %isiCDF = isiCDF(2:end); %- AJ
            %isiCDF = isiCDF(2:end-1); %- cut out endpoint so plot doesn't have weird uptick for ISI>500ms (JW)
            isiCDF = isiCDF(1:end-1); %- cut out endpoint so plot doesn't have weird uptick for ISI>500ms (JW). and include 1st point so thresh is applied correctly
            figureHandle = subplot(subplot_rows,subplot_cols,18+cdfSubplotNum); hold on; set(gca,'color','black');
            histHandle = histogram(isi,isiBins); xlabel('ISI (ms)');
            histHandle.FaceColor = colors{unitNum}; histHandle.EdgeColor = colors{unitNum+1};
            set(gca,'XLim',[-100,maxIsiToShow]);
            setFigFontSize(fontSize,figureHandle);
            figureHandle = subplot(subplot_rows,subplot_cols,27+cdfSubplotNum); hold on;
            plot(isiCDF,'LineWidth',2);
            plot([1 1]*isiThresh,[0,1],'r-');
            text(50,0.95,sprintf('%.3f (%dms)',isiCDF(isiThresh-1),isiThresh));
            set(gca,'XLim',[-100, maxIsiToShow]);
            setFigFontSize(fontSize,figureHandle);
            cdfSubplotNum = cdfSubplotNum+1;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the quality measures - Joshua
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fontSize = 14;
    figureHandle = subplot(subplot_rows,subplot_cols,[37:39,46:48]); hold on;
    title('Quality (Joshua)')
    toPlot = [isoMetrics.unitVSnoise.sigNoiseR(:,1)';isoMetrics.unitVSnoise.isoScore(:,1)';isoMetrics.unitVSnoise.fnScore(:,1)';isoMetrics.unitVSnoise.fpScore(:,1)'];
    barHandle = bar(1:4,toPlot,'LineWidth',2); set(gca,'XTick',1:4,'XTickLabel',{'SNR','isoScore','fnScore','fpScore'},'XTickLabelRotation',45);
    for ui = 1:length(allUnitNums), set(barHandle(ui),'FaceColor',colors{allUnitNums(ui)+1}); end
    plot([1.5,4.5],[1,1],'k--','LineWidth',1);
    legendLabel = {};
    for ui = 1:length(allUnitNums), legendLabel{ui} = sprintf('Unit %d',allUnitNums(ui)); end
    legend(barHandle,legendLabel);
    % Get the x-location of the bars
    bar_xLoc = bsxfun(@plus, barHandle(1).XData, [barHandle.XOffset]')';
    for barNum = 1:size(bar_xLoc,1)
        for subBarNum = 1:size(bar_xLoc,2)
            text(bar_xLoc(barNum,subBarNum),toPlot(barNum,subBarNum),sprintf('%.2f',toPlot(barNum,subBarNum)),'FontSize',16,'Rotation',45);
        end
    end
    setFigFontSize(fontSize,figureHandle);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot snippets of time series
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take snippets spanning 'snippetSecs' long from first 1/5 and last 1/5
    % of the time series
    snippetSamples = snippetSecs * samplerate;
    temp_ts = sort(sorted_timeStamp);
    segTime_1 = floor(temp_ts(floor(length(temp_ts)*(1/5)))*30);
    segTime_2 = floor(temp_ts(floor(length(temp_ts)*(4/5)))*30);
    if segTime_2 + snippetSamples >= length(binData)
        segTime_2 = length(binData) - snippetSamples - 1;
    end
    
    fontSize = 14;
    figureHandle = subplot(subplot_rows,subplot_cols,40:42); hold on;
    if segTime_1 == 0, segTime_1 = 1; end % SJ added (can't 0 index)
    plot(binData(segTime_1:segTime_1+snippetSamples),'k');
    set(gca,'XTick',0:samplerate:snippetSamples,'XTickLabel',0:snippetSecs,'XLim',[0,snippetSamples]);
    ax(1) = gca; yLim(1,:) = get(gca,'YLim');
    xlabel('Time(s)');
    setFigFontSize(fontSize,figureHandle);
    figureHandle = subplot(subplot_rows,subplot_cols,49:51); hold on;
    if segTime_2 == 0, segTime_2 = 1; end % SJ added (can't 0 index)
    plot(binData(segTime_2:segTime_2+snippetSamples),'k');
    set(gca,'XTick',0:samplerate:snippetSamples,'XTickLabel',0:snippetSecs,'XLim',[0,snippetSamples]);
    ax(2) = gca; yLim(2,:) = get(gca,'YLim');
    xlabel('Time(s)');
    for a = 1:2
        set(ax(a),'YLim',[min(yLim(:,1)),max(yLim(:,2))]);
    end
    setFigFontSize(fontSize,figureHandle);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot spkRate, quality measures vs. time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fontSize = 12;
    figureHandle = subplot(subplot_rows,subplot_cols,[43:45,52:54]); hold on; set(gca,'color','black');
    for ui = 1:length(allUnitNums)
        plot(windowMid_1,spkRate_stability(ui,:),'Color',colors{allUnitNums(ui)+1},'LineWidth',2);
    end
    ylabel('Spikes/sec');
    xlabel('Time(m)');
    setFigFontSize(fontSize,figureHandle);
    
    masterFig = tightfig(masterFig);
    set(masterFig, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    if isempty(strfind(figSavePath,'.png')), figSavePath=strcat(figSavePath,'.png'); end
    print(masterFig,figSavePath,'-dpng','-r300');
    
    fprintf('%.0fs || ',toc(tSaveFigs));
    
else % if unitCount >= 1
    
    masterFig = figure('units','normalized','outerposition',[0,0,1,1],'visible','off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time vs. peak-valley plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fontSize = 12;
    peakValley = max(sorted_waveFormValues,[],2) - min(sorted_waveFormValues,[],2);
    % wholeFigHandle = figure('units','normalized','outerposition',[0,0,1,1],'visible','off');
    figureHandle = subplot(2,2,1); hold on; set(gca,'color','black');
    plot(sorted_timeStamp_min,peakValley,'.','Color','w','MarkerSize',3);
    ylabel('Voltage (Peak - Trough)'); xlabel('Time(m)');
    setFigFontSize(fontSize,figureHandle);
    
    % Get PC1 vs. PC2 plot
    [COEFF,SCORE] = pca(sorted_waveFormValues);
    figureHandle = subplot(2,2,3); hold on; set(gca,'color','black');
    plot(SCORE(:,1),SCORE(:,2),'.','Color','w','MarkerSize',3);
    xlabel('PC1'); ylabel('PC2');
    set(gca,'XTick',[],'YTick',[]);
    setFigFontSize(fontSize,figureHandle);
    
    masterFig = tightfig(masterFig);
    set(masterFig, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    print(masterFig,figSavePath,'-dpng','-r300');
end









%%%%%%%%%%%%%%%%%%%%%%%%%
% Other utility functions:  tightfig and setFigFontSize
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hfig = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
%
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
%
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

if nargin == 0
    hfig = gcf;
end

% There can be an issue with tightfig when the user has been modifying
% the contnts manually, the code below is an attempt to resolve this,
% but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
set(hfig, 'WindowStyle', 'normal');

% 1 point is 0.3528 mm for future use

% get all the axes handles note this will also fetch legends and colorbars as well
hax = findall(hfig, 'type', 'axes');

% get the original axes units, so we can change and reset these again later
origaxunits = get(hax, 'Units');

% change the axes units to cm
set(hax, 'Units', 'centimeters');

% get various position parameters of the axes
if numel(hax) > 1
    %         fsize = cell2mat(get(hax, 'FontSize'));
    ti = cell2mat(get(hax,'TightInset'));
    pos = cell2mat(get(hax, 'Position'));
else
    %         fsize = get(hax, 'FontSize');
    ti = get(hax,'TightInset');
    pos = get(hax, 'Position');
end

% ensure very tiny border so outer box always appears
ti(ti < 0.1) = 0.15;

% we will check if any 3d axes are zoomed, to do this we will check if they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax)
    
    set(hax(i), 'LooseInset', ti(i,:));
    %         set(hax(i), 'LooseInset', [0,0,0,0]);
    
    % get the current viewing angle of the axes
    [az,el] = view(hax(i));
    
    % determine if the axes are zoomed
    iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
    
    % test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.')
    end
    
end

% we will move all the axes down and to the left by the amount necessary to just show the bottom and leftmost axes and labels etc.
moveleft = min(pos(:,1) - ti(:,1));
movedown = min(pos(:,2) - ti(:,2));

% we will also alter the height and width of the figure to just encompass the topmost and rightmost axes and lables
figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);

% move all the axes
for i = 1:numel(hax)
    set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
end

origfigunits = get(hfig, 'Units');
set(hfig, 'Units', 'centimeters');

% change the size of the figure
figpos = get(hfig, 'Position');
set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);

% change the size of the paper
set(hfig, 'PaperUnits','centimeters');
set(hfig, 'PaperSize', [figwidth, figheight]);
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition',[0 0 figwidth figheight]);

% reset to original units for axes and figure
if ~iscell(origaxunits)
    origaxunits = {origaxunits};
end

for i = 1:numel(hax)
    set(hax(i), 'Units', origaxunits{i});
end

set(hfig, 'Units', origfigunits);

%      set(hfig, 'WindowStyle', origwindowstyle);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets the font size of all text in a currently open figure
function setFigFontSize(fontSize,figureHandle)
set(gca,'FontSize',fontSize);
%     figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontSize',fontSize);



% Joshua's actual method using 144-sample interpolated waveform
function waves_joshua = getCubicSplineWaves(timeStamp_samples,binData,samplerate,preSpikeMS,postSpikeMS,bufferMS)
samplePts        = 144;
preSpikeSamples  = round(preSpikeMS*(samplerate/1000));
postSpikeSamples = round(postSpikeMS*(samplerate/1000));
bufferSamples    = round(bufferMS*(samplerate/1000));

% Get the waveform snippets from the bin file and upsample using cubic spline
% Since we need to re-align to the trough after upsampling, get
% more than 144 samples, re-align, then cut off the edges.
waves_joshua  = nan(length(timeStamp_samples),samplePts);
samplePtsBuff = (preSpikeSamples+postSpikeSamples+2*bufferSamples) * (samplePts/(preSpikeSamples+postSpikeSamples));
for ss = 1:length(timeStamp_samples)
    if timeStamp_samples(ss)>(preSpikeSamples+bufferSamples) && (timeStamp_samples(ss) < (length(binData)-(postSpikeSamples+bufferSamples))-1)
        snippet = binData(timeStamp_samples(ss)-(preSpikeSamples+bufferSamples):timeStamp_samples(ss)+(postSpikeSamples+bufferSamples)-1)';
        xx      = 1:(length(snippet)-1)/(samplePtsBuff-1):length(snippet);
        snippet = spline(1:length(snippet),snippet,xx);
        % Get the index for the trough(minimum) that's relevant (not something from the buffer region)
        preMinThresh  = bufferSamples +      (preSpikeMS-0.2)*(samplerate/1000);  % 0.2ms before where the trough should be
        postMinThresh = bufferSamples + preSpikeSamples + 0.2*(samplerate/1000);  % 0.2ms after where the trough should be
        [~, i_min]    = min(snippet(xx>=preMinThresh & xx<=postMinThresh));
        i_min         = i_min+find(xx>=preMinThresh,1,'first')-1; %- convert preMinThresh back into index 
        %i_min = find(snippet==min(snippet(xx>=preMinThresh & xx<=postMinThresh))); %- AJ's way.  fine but zerod out stim sections can cause an issue
        waves_joshua(ss,:) = snippet(i_min-samplePts*(1/3):i_min+samplePts*(2/3)-1);
    end
end

