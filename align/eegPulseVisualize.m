function [hFig_out EKG1] = eegPulseVisualize(rawPath, ekgTagName, channelNdx)
% function ekgPulseVisualize(eegDir)
%   Description: visualizes the .sync pulses
%
%   --Input:
%       rawPath: path to .21E file or directory containing raw files.
%                 all subdirectories labeled "raw" or [day]_[time] will be searched for raw files.
%
%       e.g.  /Volumes/Kareem/data/eeg/NIH018/raw/160306_1619/DA8661N4.21E  will extract and plot just DA8661N4.21E and .EEG
%       e.g.  /Volumes/Kareem/data/eeg/NIH018/raw/                    will  extract & plot all days, all sessions
%
%       ekgTagName:  tag of electrodes with pulse signal, usually 'EKG' or 'DC09'; can pass in a two-element cell array to exactly specify tags (e.g. {'EKG1', 'EKG2'} or {'DC09'}
%       channelNdx [optional] - if passed AND file contains duplicate channel names (usally an error), plot the channelNdx'th duplicate channel
%                               e.g. file contains 2 DC10 channels, plot the 2nd DC10 channel if channelNdx=2. Default 1
%
%   --Outputs:
%             --plot(s) that shows pulses (return handles to the figs as well)
%
%   UPDATED 10/2015 so it can handle original EEG-1100 and "new" EEG-1200 extended file formats JW
%

FIX_FIGURE_NUM_root = 102;

FIGURE_FOR_EACH_RAW  = 0; %each raw in its own figure
FIGURE_WITH_ALL_RAWS = 1; %single figure with subplots

if ~exist('channelNdx', 'var'), channelNdx = 1; end;
FIX_FIGURE_NUM       = FIX_FIGURE_NUM_root+(channelNdx-1)*100; %-  hack to make mulitple figures when diagnising which repeated channel is legit
hFig_out = []; %- list of figure handles created by this function

fprintf('\n\nSearching for .header and .signal files in %s:\n', rawPath);
eegRootList = {};
eegRootList = findRaws(rawPath, eegRootList);
eegRootList = unique(eegRootList);

subjStr = '';
if ~isempty(strfind(rawPath,'eeg/NIH'))
    iSubjStr = strfind(rawPath,'eeg/NIH')+[4:9];
    subjStr = sprintf(' -- %s -- ', rawPath(iSubjStr) );
end

RAW_ROOT = 0;
if strcmp('raw',rawPath(end-2:end)) || strcmp('raw',rawPath(end-3:end-1))
    RAW_ROOT = 1;
end


STIM_ROOT = 0;
if strcmp('STIM_MAP',rawPath(end-7:end)) || strcmp('STIM_MAP',rawPath(end-8:end-1))
    STIM_ROOT = 1;
end


fprintf('\n')
if isempty(eegRootList), fprintf('> No files found!!\n\n');
else
    % dont plot data from NSP2, if multiple NSPs
    indexNSP2 = contains(eegRootList, 'ieeg2') | contains(eegRootList,'INST1');
    indexNSP2 = find(indexNSP2);
    eegRootListNSP2 = eegRootList(indexNSP2);
    eegRootList(indexNSP2) = [];
    
    
    numRow = min( [5 length(eegRootList)] );
    numCol = 0 + (ceil(length(eegRootList)/numRow));
    
    
    if (length(eegRootList)==1),  FIGURE_FOR_EACH_RAW=1; end;
    
    isDC09DC10 = 0; %- flag that will make DC10 negative and red
    if length(ekgTagName)==2,   ekgTag2out = ekgTagName{2}; chanStr=sprintf('%s-%s',ekgTagName{1},ekgTagName{2});
        if strcmp(ekgTagName{1}(1:2),'DC')&strcmp(ekgTagName{2}(1:2),'DC') isDC09DC10 = 1; ekgTagName{3}='DC11'; end
    else                        ekgTag2out = 'gnd';         chanStr = ekgTagName{1}; end
    
    %- grab all the data (and the embedded dates for sorting the output plots)
    for iList=1:length(eegRootList)
        
        %-create a clean version of the path for the figure title (and commandline output)
        clnTitle = eegRootList{iList};
        %if ~isempty(strfind(clnTitle,'raw')), clnTitle = clnTitle(strfind(clnTitle,'raw')+4:end); end
        
        %- raw file names are crazy now... so just use the session name
        if ~isempty(strfind(clnTitle,'raw')),
            clnTitle = clnTitle(strfind(clnTitle,'raw')+4:end);
            parts    = strsplit(clnTitle,'/');
            clnTitle = parts{1};
            if contains(eegRootList{iList},'STIM_MAP') & length(parts)>1,
                clnTitle = ['STIM_MAP/' parts{2}];
            end
            
            %- and then add back in some deets about the recording device
            [~,eegRootFile,~] = fileparts( eegRootList{iList} );
            if     contains(eegRootFile,'INST')|contains(eegRootFile,'ieeg'),     clnTitle = [clnTitle '(BR)'];
            elseif contains(eegRootFile,'EEG') |contains(eegRootFile,'cervello'), clnTitle = [clnTitle '(CV)'];
            else clnTitle = [clnTitle '(NK)']; end
        end
        
        clnTitle(find(clnTitle=='_' | clnTitle=='\'))=' ';
        fprintf('... loading %d of %d: %s ...',iList,length(eegRootList),clnTitle)
        
        [EKG1, EKG2, sampRate, strTime, numTime, annOut, EKG3] = grabPulses(eegRootList{iList}, ekgTagName, channelNdx);
        
        fprintf(' [%s]\n', strTime);
        
        list_clnTitle{iList} = clnTitle;
        list_EKG1{iList}     = EKG1;
        list_EKG2{iList}     = EKG2;
        list_EKG3{iList}     = EKG3;
        list_sampRate{iList} = sampRate;
        list_strTime{iList}  = strTime;
        list_numTime(iList)  = numTime;
        list_annOut{iList}   = annOut; %-structure with annotation information
    end
    
    
    %- open the figure and make full screen so axis are sized correctly before automatic save
    if FIGURE_WITH_ALL_RAWS
        if FIX_FIGURE_NUM==0, hFig = figure;
        else                  hFig = figure(FIX_FIGURE_NUM); set(hFig,'name',sprintf('%s All Raw Files <%s>',subjStr, ekgTagName{1})); set(hFig,'Color','w','units','normalized','outerposition',[0 0 1 1]); end
        clf;
        hFig_out(end+1) = hFig;
    end
    
    
    %- sort the list by embedded raw date/time
    [sorted, iSort] = sort(list_numTime);
    
    % extract pulse data from NSP2 files as well, if they exist, and
    % cross-correlate with matching NSP1 pulses
    lags = NaN(1,15);
    for iList=1:length(eegRootListNSP2)
        [EKG1, EKG2, sampRate, strTime, numTime, annOut] = grabPulses(eegRootListNSP2{iList}, ekgTagName, channelNdx);
        list_EKG1_NSP2{iList}     = EKG1;
        list_EKG2_NSP2{iList}     = EKG2;
        
        if contains(eegRootListNSP2{iList},'ieeg'),
            temp_name_match = strrep(eegRootListNSP2{iList},'ieeg2','ieeg1');
        else
            temp_name_match = strrep(eegRootListNSP2{iList},'INST1','INST0');
        end
        nsp1_file_match = find(strcmp(temp_name_match,eegRootList));
        
        [seq,lag] = xcorr(EKG1,list_EKG1{nsp1_file_match},50);
        
        [~,I] = max(abs(seq));
        lags(iList) = lag(I);         % sensor 2 leads sensor 1 by 350 samples
        if abs(lag(I)) > 5;
            fprintf('NSPs are not synced within 5 ms of one another!! To be corrected when files are split...')
        end
    end
    
    
    % External counter for lag label
    iLag = 1;
    
    %- now plot all the loaded data!
    for iList = iSort
        
        clnTitle = list_clnTitle{iList};
        %clnTitle = strrep(clnTitle,'ieeg','blackrock'); %- old way... but now we only show the session folder name anyway
        %clnTitle = strrep(clnTitle,'INST','blackrock');
        EKG1     = list_EKG1{iList};
        EKG2     = list_EKG2{iList};
        EKG3     = list_EKG3{iList};
        sampRate = list_sampRate{iList};
        strTime  = list_strTime{iList};
        numTime  = list_numTime(iList);
        annOut   = list_annOut{iList};
        
        plotNum  = find(iSort==iList);
        
        if FIGURE_WITH_ALL_RAWS
            figure(hFig);
            
            sp = subplot(numRow,numCol,plotNum);
            t = [1:length(EKG1)]/sampRate/60;  %convert to minutes
            if isDC09DC10,
                if isempty(EKG2); EKG2 = zeros(1,size(EKG1,2)); end
                if isempty(EKG3); EKG3 = zeros(1,size(EKG1,2)); end  %- EKG3 is DC11
                plot(t, EKG1, 'b', t, -EKG2,'r', t, -EKG3,'m');
                set(gca,'ylim',[-5500 5500]);
                
                [gapX,gapY] = find(EKG1==-999);
                if size(gapX ~=0)
                    vertical1 = line([t(gapY) t(gapY)],[max(EKG1) -max(EKG2)]);
                    set(vertical1,'Color','k','linewidth',6)
                    hold on
                    vertical2 = line([t(gapY) t(gapY)],[max(EKG1) -max(EKG2)]);
                    set(vertical2,'Color','y','linewidth',4,'LineStyle','--')
                end
            else
                plot(t, EKG1-EKG2);
            end
            
            set(gca,'box','off','tickdir','out','fontsize',15,'xlim',[t([1 end])]);
            %title( sprintf('%s :: [%s]',clnTitle,strTime),'fontsize', 15);
            title( sprintf('%s',clnTitle));
            if size(eegRootListNSP2,2) > 0 && contains(eegRootList(iList),'INST1'),
                str = sprintf('NSP1-NSP2 Lag = %0.2f ms',lags(iLag));
                pos = get(sp,'position');
                dim = [pos(1)+(pos(4)*0.1) pos(2)+(pos(4)*0.3) pos(3)*0.1 pos(4)*0.1];
                t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
                t.BackgroundColor = 'w';
                if abs(lags(iLag)) > 5;
                    t.Color = 'k';
                    t.BackgroundColor = 'y';
                end
                iLag = iLag + 1;
            end
            
            if sum(~cellfun('isempty',strfind(ekgTagName,'DC'))), strUnits = '(mV)'; else strUnits = '(uV)'; end
            if plotNum==1, ylabel(sprintf('%s-%s %s',ekgTagName{1},ekgTag2out, strUnits));  xlabel('time (min)'); end
            %hT=text(mean(get(gca,'xlim')),max(get(gca,'ylim'))*.99,sprintf('[%s]',strTime));
            %set(hT,'fontsize',18,'HorizontalAlignment','center','VerticalAlignment','top');
            %legend(strTime)
            
        end
        
        if FIGURE_FOR_EACH_RAW,
            figure(FIX_FIGURE_NUM+plotNum); clf
            set(gcf,'name',clnTitle, 'color', 'w')
            hFig_out(end+1) = gcf;
            
            t = [1:length(EKG1)]/sampRate;  %convert to seconds for the individual figures
            
            if length(find(EKG2~=0))>0, OK_EKG2 = 1; else OK_EKG2 = 0; end
            if sum(~cellfun('isempty',strfind(ekgTagName,'DC'))), strUnits = '(mV)'; else strUnits = '(uV)'; end
            
            annFont = 15;
            
            
            %%- top plot (or only plot) is the difference between EKG1 and 2
            subplot(1+OK_EKG2*2,1,1); thisAx(1) = gca;
            plot(t, EKG1-EKG2, '-'); hold on
            yAnn = 0;
            for iAnn = 1:length(annOut.samp),
                hPt = plot(t(annOut.samp(iAnn)),yAnn,'*b','MarkerSize',15);
                hTx = text(t(annOut.samp(iAnn)),yAnn,annOut.text{iAnn},'Rotation',90,'FontSize',annFont);
            end
            axis tight
            grid on
            box off
            
            set(gca,'fontsize',13);
            ylabel(sprintf('%s-%s %s',ekgTagName{1},ekgTag2out, strUnits));
            xlabel('time (s)', 'fontsize',13)
            title(clnTitle,'fontsize', 15, 'Interpreter','none');
            legend(strTime)
            
            %%- if EKG2 is non-zero, then plot EGK1 and EGK2 separately
            if OK_EKG2,
                subplot(2+OK_EKG2,1,2); thisAx(2) = gca;
                plot(t, EKG1,'r'); hold on;
                grid on; box off; set(gca,'fontsize',13);
                if ~isempty(strfind(ekgTagName{1},'DC')), strUnits = '(mV)'; else, strUnits = '(uV)'; end
                ylabel(sprintf('%s %s',ekgTagName{1},strUnits));
                %title(clnTitle,'fontsize', 15);
                yAnn = 0;
                for iAnn = 1:length(annOut.samp),
                    hPt = plot(t(annOut.samp(iAnn)),yAnn,'*b','MarkerSize',15);
                    hTx = text(t(annOut.samp(iAnn)),yAnn,annOut.text{iAnn},'Rotation',90,'FontSize',annFont);
                end
                
                subplot(3,1,3); thisAx(3) = gca;
                plot(t, EKG2,'b'); hold on;
                grid on; box off; set(gca,'fontsize',13);
                if ~isempty(strfind(ekgTagName{2},'DC')), strUnits = '(mV)'; else, strUnits = '(uV)'; end
                ylabel(sprintf('%s %s',ekgTagName{2},strUnits));
                xlabel('time (s)', 'fontsize',13)
                yAnn = 0;
                for iAnn = 1:length(annOut.samp),
                    hPt = plot(t(annOut.samp(iAnn)),yAnn,'*m','MarkerSize',15);
                    hTx = text(t(annOut.samp(iAnn)),yAnn,annOut.text{iAnn},'Rotation',90,'FontSize',annFont);
                end
            end
            
            linkaxes(thisAx,'x')
            
            
            %%% ---  If annotation present, make the same plot in a vertical axis --- %%%
            if ~isempty(annOut.samp),
                figure(FIX_FIGURE_NUM+50+plotNum); clf
                set(gcf,'name',clnTitle, 'color', 'w')
                
                
                subplot(1,1+OK_EKG2,1);  vertAx(1) = gca;
                plot(EKG1,t,'r'); hold on; grid on;
                if sum(~cellfun('isempty',strfind(ekgTagName,'DC')))
                    xlabel(sprintf('%s (mV)',ekgTagName{1}));
                else
                    xlabel(sprintf('%s (uV)',ekgTagName{1}));
                end
                yAnn = 0;
                for iAnn = 1:length(annOut.samp),
                    hPt = plot(yAnn,t(annOut.samp(iAnn)),'*b','MarkerSize',15);
                    hTx = text(yAnn,t(annOut.samp(iAnn)),annOut.text{iAnn},'Rotation',0,'FontSize',annFont);
                end
                set(gca,'fontsize',13);
                set(vertAx(1),'YDir','reverse','Ylim',[0 max(t)]);
                
                if OK_EKG2,
                    subplot(1,2,2);  vertAx(2) = gca; set(gca,'YDir','reverse');
                    plot(EKG2,t,'b'); hold on;  grid on;
                    xlabel(ekgTag2out);
                    yAnn = 0;
                    for iAnn = 1:length(annOut.samp),
                        hPt = plot(yAnn,t(annOut.samp(iAnn)),'*m','MarkerSize',15);
                        hTx = text(yAnn,t(annOut.samp(iAnn)),annOut.text{iAnn},'Rotation',0,'FontSize',annFont);
                    end
                    set(vertAx(2),'YDir','reverse','Ylim',[0 max(t)]);
                end
                set(gca,'fontsize',13);
                
                linkaxes(vertAx,'y')
                
            end
            
        end
    end
    
    %- save a copy of the master figure to subjet's raw directory
    if FIGURE_WITH_ALL_RAWS && (RAW_ROOT || STIM_ROOT),
        if STIM_ROOT, fileSaveName = 'align_PlotPulseChannels_STIM';
        else          fileSaveName = 'align_PlotPulseChannels';
        end
        
        while exist(fullfile(rawPath, sprintf('%s.png',fileSaveName)),'file'), fileSaveName = [fileSaveName '_x']; end %- prevent overwrite
        
        try
            fig2pngSimple(hFig,fullfile(rawPath, sprintf('%s.png',fileSaveName))); %- new version... just save the fig (no screen shot)
        catch
            fprintf('\n uh oh... reverting to screenshot method');
            if ismac,
                reply = input('\n\n PRESS RETURN TO TAKE SCREENSHOT OF THE MASTER FIG AND DROP IT IN SUBJECT/RAW!!!!\n','s');
                figure(hFig);
                pause(1);
                unix(sprintf('screencapture -T 1 "%s"', fullfile(rawPath, sprintf('%s.png',fileSaveName))));  % quick view on mac... this makes the image match the screen, otherwise JW uses a bunch of functions outside of toolbox...
            else
                %- this has become a very big file...
                saveas(hFig, fullfile(rawPath, sprintf('%s.fig',fileSaveName)), 'fig');
                %saveas(hFig, 'align_PulseChannels.pdf', 'pdf');
            end
        end
        
    end
    
end
fprintf('\n')

end %eegPulseVisualize




%%%-----------------------------------------------------------------------%%%
%%% Extract the EKG channels and Raw file creation date
function [EKG1, EKG2, sampRate, strTime, numTime, annOut, EKG3] = grabPulses(eegRoot, ekgTagName, channelNdx)
% Parameter Notes
%   channelNdx [optional] - if passed AND file contains duplicate channel
%                       names (an error), plot the channelNdx'th channel
%
%
% 12/2017 melkalliny - will now handle raw blackrock files
% 3/2018 melkalliny - will now handle raw cervello files
% 10/2018 JW -- add a 3rd trace for DC11 so we can see when microStim was applied



[eegRootPath,eegRootFile,~] = fileparts(eegRoot);
if  contains(eegRootFile,'ieeg') | contains(eegRootFile,'INST'), %if blackrock
    
    %% OLD WAY OF DOING IT: assume its on ns3, and there is nothing else there
    % get analog DC
    %EEG_file = strcat(eegRoot,'.ns3');
    %NS3FileDir = EEG_file; NS3FileDir = strrep(NS3FileDir,'ns2','ns3');  %- assume channels on NS3, that is the most common
    
    %     if exist(NS3FileDir,'file'),
    %         NS3data = concatOpenNSx(NS3FileDir); % For analog DC channel
    %         rawData_NS3 = NS3data.Data;
    %     else
    %         fprintf('\n cant find the expected pulse file... return?');
    %         keyboard
    %     end
    
    %% NEW WAY: use jacksheetBR to figure out which NSx file to look in and which channels within that NSx contain DC9-12
    %- get ananlog PULSE channels.  Following code will identify the appropriate NSx file and index within for just DC9-12
    jackBR = dir(fullfile(eegRootPath,'jacksheetBR_local*.csv'));  %- would like to use complete here, but for some subjects the NSx files were renamed
    if length(jackBR)==0,
        jackBR = makeJacksheetBR(eegRootPath,'',0,1,0);
    else
        jackBR = readtable(fullfile(eegRootPath,jackBR.name));
    end
    
    
    
    %- initialize return variables in case need to break out
    EKG1=zeros(1,100); EKG2=nan(1,100); EKG3=nan(1,100);
    sampRate=1000; strTime=jackBR.RawDir{1}; numTime=0;
    annOut.samp = [];
    annOut.text = {};

    
    %- for now get ALL analog channels that are recorded on this NSP
    iDC = find(jackBR.PhysicalChan>=1001 & jackBR.PhysicalChan<=1032 & ~contains(jackBR.FileName,'.nev') & contains(jackBR.FileName,eegRootFile)); %- DC channels from THIS nsp (not NEV)
    if ~isempty(iDC),
        NS3file     = jackBR.FileName{iDC(1)};
        NS3FileDir  = fullfile(eegRootPath,NS3file);
        NS3data     = concatOpenNSx(NS3FileDir); %- For analog DC channel
        jackNSX     = jackBR(strcmp(jackBR.FileName, NS3file),:); %- jacksheet just for this NSx file; this is ordered to match the outputs of openNSx
        iChanNSP    = find(ismember(jackNSX,jackBR(iDC,:))); %- find rows in jackNSX that match the DC channels identified above
        %keyboard
        rawData_NS3 = NS3data.Data(iChanNSP,:);
        NS3data.ElectrodesInfo = NS3data.ElectrodesInfo(iChanNSP);
    else
        fprintf('\n cant find the expected pulse file... return?');
        EKG1=nan(100,1); EKG2=nan(100,1); EKG3=nan(100,1);
        sampRate=1000; strTime=jackBR.RawDir{1}; numTime=0;
        annOut.samp = [];
        annOut.text = {};
        %keyboard
        return;
    end
    
    
    
    %convert requested tag name into blackrock name... first make sure its "DC"
    if sum(contains(ekgTagName,'DC'))~=length(ekgTagName),
        fprintf('\n Uh Oh... blackrock files only setup to visualize DC channels for now (not chan %s)',ekgTagName{1});
        %keyboard;
        return;
    end
    dcNum = [];
    for ii=1:length(ekgTagName),
        dcNum = [dcNum str2num(ekgTagName{ii}(3:end))];
    end
    dcNum = dcNum-8; %- conversion from DC to ain...  DC starts at 9, ain starts at 1, so subtract 8
    
    %rawData_NS3 = getNSxData(NS3data);
    DCtimeseries = {};
    sampRate = 1000;
    for DC=dcNum;
        if contains(eegRootFile,'ieeg2') | contains(eegRootFile,'INST1')
            analogChanName = sprintf('ain%d',DC+16);
        else
            analogChanName = sprintf('ain%d',DC);
        end
        allChanNames = {NS3data.ElectrodesInfo.Label}';
        
        if iscell(rawData_NS3)
            allChanNames = allChanNames(1:size(rawData_NS3{1},1));
        else
            allChanNames = allChanNames(1:size(rawData_NS3,1));
        end
        analogDC = 0;
        if DC <= size(allChanNames,1)
            %analogDC = rawData_NS3(not(cellfun('isempty', strfind(allChanNames, analogChanName))),:);
            if iscell(rawData_NS3);
                analogDC = rawData_NS3{1}(DC,:); %assumption that blackrock recorded 1-4 in order
            else
                analogDC = rawData_NS3(DC,:); %assumption that blackrock recorded 1-4 in order
            end
            [fsorig, fsres] = rat(NS3data.MetaTags.SamplingFreq/sampRate); % resample
            splitIndex = find(analogDC==999.999);
            analogDC=resample(double(analogDC),fsres,fsorig);
            if size(splitIndex,2) ~=0
                analogDC(round(splitIndex(1) / 2)) = -999;
            end
            
            tempIndex = find(analogDC == -999);
            
            % multiply by gain
            maxDigital = NS3data.ElectrodesInfo(1).MaxDigiValue;
            maxAnalog = NS3data.ElectrodesInfo(1).MaxAnalogValue;
            gain = 1/(double(maxDigital)/double(maxAnalog));
            analogDC = analogDC .* gain;
            
            analogDC(tempIndex) = deal(-999);
            
            %resamplerate = NS3data.MetaTags.SamplingFreq/sampRate;
            %analogDC = analogDC(floor(linspace(1,length(analogDC),length(analogDC)/resamplerate)));
            DCtimeseries{find(DC==dcNum)} = analogDC;
        end
        % get timestamp of analog DC09 pulses
    end;
    EKG1 = double(DCtimeseries{1});
    EKG2 = double(DCtimeseries{2});  EKG3 = double(DCtimeseries{3});
    
    % getting output variables in order
    %EEG_file = strrep(EEG_file,'ns3','ns2');
    %NS3FileDir = EEG_file;
    temp = strfind(NS3FileDir,'/'); dateTime = NS3FileDir(temp(end-1)+1:temp(end)-1);
    numTime=0; % isnt critical. for debugging purposes, lets make it the starting timestamp
    strTime = dateTime;
    
    annOut.samp = [];
    annOut.text = {};
    
    
elseif contains(eegRootFile,'cervello') | contains(eegRootFile,'EEG'); %if cervello
    EEG_file = strcat(eegRoot,'.trc');
    chansGrab = [];
    EEG = readCervelloTRC(EEG_file,chansGrab);  % in uV
    actualName_ALL = {EEG.electrodes.tagName}'; % dont use this for writing out names
    GAIN           = [EEG.electrodes.measurement_unit]';  %- units or gain in volts... 1e-6 means microvolts reported
    fileStemDate   = EEG.fileStemDate;
    if numel(unique(GAIN)) > 1; fprintf('error - not all chans have same gain');
        keyboard; end
    if unique(GAIN) ~= 1*10^-6; fprintf('unexpected- raw data from Cervello not in uV');
        keyboard; end;
    sampRate=EEG.srate; % find sample rate
    
    
    RESAMP=0; %- no need to resample here... and its a slow step
    if RESAMP,
        samprate_orig=EEG.srate; % find sample rate
        samprate_new=1000;
        [fsorig, fsres] = rat(samprate_orig/samprate_new); % resample
        % initialize matrix of resampled data
        [temp] = resample((EEG.data(1,:)),fsres,fsorig); % fx resamples to 2nd entry.incorrectly labeled as 'orig'!!
        chan=cell(size(temp,2),1);
        EEG_data_resamp_temp = zeros(size(EEG.data,1),size(temp,2));
        % resample and assign
        for channel=1:size(EEG.data,1);
            EEG_data_resamp_temp(channel,:) = resample((EEG.data(channel,:)),fsres,fsorig);
            chan{channel}=resample((EEG.data(channel,:)),fsres,fsorig);
        end
        EEG.data = EEG_data_resamp_temp;
        clear EEG_data_resamp_temp
        sampRate = 1000;
    end
    
    
    DCtimeseries = {};
    for DC=[1:4];
        
        % load correct ain data
        analogChanName = sprintf('ain%d',DC);
        temp  = strfind({EEG.electrodes.tagName},analogChanName);
        match = find(not(cellfun('isempty', temp)));
        if ~isempty(match),
            analogDC = EEG.data(match(1),:);
            
            %- optional... fix the pulses here they way they will be fixed for alignment... else commment this line out for the "raw" version of the pulses
            analogDC = fix_CV_analog(analogDC, analogChanName);
            
        else
            analogDC = nan(1,EEG.pnts);
        end
        % get timestamp of analog DC09 pulses
        %[syncPulses]=get_triggers(double(analogDC'),samprate_new);
        
        
        DCtimeseries{DC} = analogDC;
    end
    EKG1 = double(DCtimeseries{1}); EKG2 = double(DCtimeseries{2});
    EKG3 = double(DCtimeseries{3});
    
    % getting output variables in order
    NS3FileDir = EEG_file;
    temp = strfind(NS3FileDir,'/'); dateTime = NS3FileDir(temp(end-1)+1:temp(end)-1);
    numTime=0; % isnt critical. strTime is the important one
    strTime = dateTime;
    annOut.samp = [];
    
    
else; %if nihon kohden
    
    VERBOSE = 0; % 1=some info, 2=LOTS of INFO
    
    % gets filepath for .EEG file
    EEG_file=[eegRoot '.EEG'];
    
    % Same as above for .21E file
    ELEC_file=[eegRoot, '.21E'];
    
    %% gets all the codes and names from the .21E file
    [allCodes,allNames] = textread2(ELEC_file,'%s%s','delimiter','='); %textread reads data from a txt file and write to multiple outputs
    endRange  = find(strcmp(allCodes,'[SD_DEF]')); %finds the range of the electrodes
    allCodes  = allCodes(1:endRange-1);  %stores the codes
    allNames  = allNames(1:endRange-1);  %stores the names
    %disp([allCodes,allNames])
    %goodCodes = [0:36 74 75 100:253];
    goodCodes = [0:36 42:73 74:77 100:253];  %include DC channels, 42-73, plus mark channels 76-77
    badNames  = {'E'};
    actualCode_ALL = {};
    actualName_ALL = {};
    actualNameWbad_ALL = {};  %jw added this... makes it easier to track the file offset for the target channel
    
    
    %% Gets to data in the waveform block
    fid = fopen(EEG_file);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skipping EEG device block: 128 bytes, same for new (EEG-1200) and old (EEG-1100) filetypes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    deviceBlockLen=128; %skips the first 128 bytes
    fseek(fid,deviceBlockLen,'bof');  %fseek(fileID, offset, origin) moves to specified position in file. bof=beginning of file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reading EEG1 control Block (contains names and addresses for EEG2 blocks) -- 896 bytes for new and old filetypes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=fread(fid,1,'*uint8');  if VERBOSE==2, fprintf('block ID: %d\n',x); end
    x=fread(fid,16,'*char');  if VERBOSE==2, fprintf('device type: %s\n',x); end;  if strcmp(x(1:9)','EEG-1200A'),NEW_FORMAT=1;else NEW_FORMAT=0; end;
    x=fread(fid,1,'*uint8');  if VERBOSE==2, fprintf('number of EEG2 control blocks: %d\n',x); end
    numberOfBlocks=x;
    if numberOfBlocks > 1
        % we think we will never have this
        % throw an error for now and re-write code if necessary
        fprintf('ERROR: %d EEG2 control blocks detected (only expecting 1).\n');
        return
    end
    % if numberOfBlocks is ever > 1, the following should be a for loop
    blockAddress=fread(fid,1,'*int32');  if VERBOSE==2, fprintf('address of block %d: %d\n',i,blockAddress); end
    x=fread(fid,16,'*char');             if VERBOSE==2, fprintf('name of EEG2 block: %s\n',x); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reading EEG2m control block (contains names and addresses for waveform blocks)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fseek(fid,blockAddress,'bof');          if VERBOSE==2, fprintf('\nin EEG21 block!\n');  end
    x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('block ID: %d\n',x); end
    x=fread(fid,16,'*char');                if VERBOSE==2, fprintf('data format: %s\n',x); end
    numberOfBlocks=fread(fid,1,'*uint8');   if VERBOSE==2, fprintf('number of waveform blocks: %d\n',numberOfBlocks); end
    if numberOfBlocks > 2
        % we think we will never have this
        % throw an error for now and re-write code if necessary
        fprintf('ERROR: %d waveform blocks detected (only expecting 1).\n');
        return
    end
    % if numberOfBlocks is ever > 1, the following should be a for loop
    blockAddress=fread(fid,1,'*int32'); if VERBOSE==2, fprintf('address of block %d: %d\n',i,blockAddress); end
    x=fread(fid,16,'*char');            if VERBOSE==2, fprintf('name of waveform block: %s\n',x); end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %Reading waveform block -- if New format the original waveform block will contain a single channel (channel 1) for exactly 1 second..
    %%%%%%%%%%%%%%%%%%%%%%%
    fseek(fid,blockAddress,'bof');      if VERBOSE==2, fprintf('\nin EEG waveform block!\n'); end
    x=fread(fid,1,'*uint8');            if VERBOSE==2, fprintf('block ID: %d\n',x); end
    x=fread(fid,16,'*char');            if VERBOSE==2, fprintf('data format: %s\n',x); end
    x=fread(fid,1,'*uint8');            if VERBOSE==2, fprintf('data type: %d\n',x); end
    L=fread(fid,1,'*uint8');            if VERBOSE==2, fprintf('byte length of one data: %d\n',L); end
    M=fread(fid,1,'*uint8');            if VERBOSE==2, fprintf('mark/event flag: %d\n',M); end
    
    %%- annonomous function to convert binary to decimal.  input is binary string created with dec2bin
    bcdConverter2 = @(strDec2bin)  10*bin2dec(strDec2bin(1:4)) + bin2dec(strDec2bin(5:8));
    
    % get the start time
    T_year   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    T_month  = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    T_day    = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    T_hour   = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    T_minute = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    T_second = bcdConverter2(dec2bin(fread(fid,1,'*uint8'),8));
    strTime  = sprintf('%d/%d/%d %02d:%02d:%02d',T_month,T_day,T_year,T_hour,T_minute,T_second); if NEW_FORMAT, strTime=[strTime '(N)']; end
    numTime  = datenum(T_year,T_month,T_day,T_hour,T_minute,T_second);
    
    % get the sampling rate
    x=fread(fid,1,'*uint16');  %fprintf('sample rate (coded): %d\n',x);
    switch(x)
        case hex2dec('C064'),
            actSamplerate=100;
        case hex2dec('C0C8'),
            actSamplerate=200;
        case hex2dec('C1F4'),
            actSamplerate=500;
        case hex2dec('C3E8'),
            actSamplerate=1000;
        case hex2dec('C7D0'),
            actSamplerate=2000;
        case hex2dec('D388'),
            actSamplerate=5000;
        case hex2dec('E710'),
            actSamplerate=10000;
        otherwise
            fprintf('UNKNOWN SAMPLING RATE\n');
    end
    if VERBOSE==2, fprintf('Sampling rate: %d Hz\n',actSamplerate); end
    
    
    % get the number of 100 msec block
    num100msBlocks=fread(fid,1,'*uint32');         if VERBOSE==2, fprintf('Length of Session: %2.2f hours\n',double(num100msBlocks)/10/3600); end
    numSamples  = actSamplerate*num100msBlocks/10; if VERBOSE==2, fprintf('number of samples: %d\n',numSamples); end
    AD_off      = fread(fid,1,'*int16');           if VERBOSE==2, fprintf('AD offset at 0 volt: %d\n',AD_off); end
    AD_val      = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('AD val for 1 division: %d\n',AD_val); end
    bitLen      = fread(fid,1,'*uint8');           if VERBOSE==2, fprintf('bit length of one sample: %d\n',bitLen); end
    comFlag     = fread(fid,1,'*uint8');           if VERBOSE==2, fprintf('data compression: %d\n',comFlag); end
    numChannels = fread(fid,1,'*uint8');           if VERBOSE==2, fprintf('number of RAW recordings: %d\n',numChannels); end
    
    
    if (numChannels==1 & numSamples-actSamplerate==0) && NEW_FORMAT==0, fprintf('\n expecting old format, but 1 channel for 1 second'); keyboard; end
    if (numChannels>1)                                && NEW_FORMAT==1, fprintf('\n expecting new format, but >1 channel ');           keyboard; end
    
    if NEW_FORMAT || (numChannels==1 && numSamples-actSamplerate==0)
        
        if VERBOSE, fprintf('** New File Format **'); end
        
        %- seek the file location of the new wave data... need to make a pit stop in the new EEG2 header, which will provide the direct address
        
        waveformBlockOldFormat = 39 + 10 + 2*actSamplerate + double(M)*actSamplerate; %- with new format the initial waveform block contains 1 channel (10 bytes of info) for 1 second (2bytes x 1000)
        controlBlockEEG1new    = 1072;
        %controlBlockEEG2new    = 20+24*double(numberOfBlocks); % this isn't working, so read the wave start location from controlBlockEEG2
        blockAddressEEG2       = blockAddress + waveformBlockOldFormat + controlBlockEEG1new;% + controlBlockEEG1new + controlBlockEEG2new;
        
        
        
        %- EEG2' format
        addTry = blockAddressEEG2;
        fseek(fid,addTry,'bof');                if VERBOSE, fprintf('--EEG2-prime format--\n'); end
        x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('block ID: %d\n',x); end
        x=fread(fid,16,'*char');                if VERBOSE==2, fprintf('data format: %s\n',x); end
        x=fread(fid,1,'*uint16');               if VERBOSE==2, fprintf('number of waveform blocks: %d\n',x); end;
        x=fread(fid,1,'*char');                 if VERBOSE==2, fprintf('reserved: %s\n',x); end
        x=fread(fid,1,'*int64');ii=1;           if VERBOSE==2, fprintf('address of block %d: %d\n',ii,x);    end; waveBlockNew = x;
        
        
        %- EEG2' waveform format
        fseek(fid,waveBlockNew,'bof');          if VERBOSE, fprintf('--EEG2-prime WAVE format--\n'); end
        x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('block ID: %d\n',x); end
        x=fread(fid,16,'*char');                if VERBOSE==2, fprintf('data format: %s\n',x); end
        x=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('data type: %d\n',x); end
        L=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('byte length of one data: %d\n',L); end
        M=fread(fid,1,'*uint8');                if VERBOSE==2, fprintf('mark/event flag: %d\n',M); end
        
        %- now things get a little different with the new header
        x=fread(fid,20,'*char');                if VERBOSE==2, fprintf('start time string: %s\n',x); end
        x=fread(fid,1,'*uint32');               if VERBOSE==2, fprintf('data interval (sample rate): %d\n',x);                end; actSamplerate  = double(x);
        x=fread(fid,1,'*uint64');               if VERBOSE==2, fprintf('Length of Session: %2.2f hours\n',double(x)/10/3600); end; num100msBlocks = double(x);
        
        numSamples  = actSamplerate*num100msBlocks/10; if VERBOSE==2, fprintf('number of samples: %d\n',numSamples); end
        AD_off      = fread(fid,1,'*int16');           if VERBOSE==2, fprintf('AD offset at 0 volt: %d\n',AD_off); end
        AD_val      = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('AD val for 1 division: %d\n',AD_val); end
        bitLen      = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('bit length of one sample: %d\n',bitLen); end
        comFlag     = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('data compression: %d\n',comFlag); end
        reserveL    = fread(fid,1,'*uint16');          if VERBOSE==2, fprintf('reserve length: %d\n',reserveL); end
        x           = fread(fid,reserveL,'*char');     if VERBOSE==2, fprintf('reserve data: %s\n',x); end
        
        numChannels = fread(fid,1,'*uint32');          if VERBOSE==2, fprintf('number of RAW recordings: %d\n',numChannels); end
        
    end
    
    
    %- parse the channel names -- connect .21E information with .EEG information
    listChanStringCode = {};
    listActualName     = {};
    for k=1:numChannels
        x=fread(fid,1,'*int16');  %reads in 1 byte every time you iterate the loop
        chanCode(k)=x; %and stores it in chanCode(k)
        if (VERBOSE) fprintf(' Index %d ''name'': Channel %d\n',k,x); end
        chanCodeString=sprintf('%04d',x); %format data into string. Same as chanCode except the format is string.
        matchingRow=find(strcmp(chanCodeString,allCodes)); %looks for this particular string in allCodes and stores its locations in matchingRow
        actualName=allNames{matchingRow};
        
        listChanStringCode{end+1} = chanCodeString;
        listActualName{end+1} = actualName;
        
        if ~ismember(chanCode(k),goodCodes) %if not a member of goodCodes
            if (VERBOSE) fprintf(' chan %d (%s) is a bad channel code and excluded\n',chanCode(k),actualName); end;
            goodElec(k)=false;
        elseif any(strcmp(actualName,badNames)) %or if it's part of badNames
            if (VERBOSE) fprintf(' chan %d (%s) is a bad address\n',chanCode(k),actualName); end;
            goodElec(k)=false;
        else
            if (VERBOSE) fprintf(' chan %d (%s) is good!\n',chanCode(k),actualName); end
            goodElec(k)=true;
        end
        
        % save out the names for the jacksheet
        if goodElec(k); actualName_ALL(end+1)=allNames(matchingRow); actualCode_ALL{end+1}=chanCodeString; end %if it is a good electrode, append it to the jacksheet
        actualNameWbad_ALL(end+1)=allNames(matchingRow);
        
        fseek(fid,6,'cof'); %skipping the six most sig. bits of 'name'
        
        %finds the channel sensitivity
        chan_sensitivity = fread(fid,1,'*uint8');   if VERBOSE==2, fprintf('channel sensitivity: %d\n',chan_sensitivity); end
        chan_unit        = fread(fid,1,'*uint8');   if VERBOSE==2, fprintf('         unit: %d\n',chan_unit); end
        switch chan_unit,
            case 0; CAL=1000;%microvolt
            case 1; CAL=2;%microvolt
            case 2; CAL=5;%microvolt
            case 3; CAL=10;%microvolt
            case 4; CAL=20;%microvolt
            case 5; CAL=50;%microvolt
            case 6; CAL=100;%microvolt
            case 7; CAL=200;%microvolt
            case 8; CAL=500;%microvolt
            case 9; CAL=1000;%microvolt
        end
        GAIN(k)=CAL/double(AD_val);%OK TO ASSUME THIS IS CONSTANT FOR ALL ELECS?
    end
    %disp([listChanStringCode' listActualName']);
    
    %- starting point of filepointer for reading the data
    fStart = ftell(fid);
    
    tReadAll = tic;
    %fprintf('\nReading Data...') %\n means new line
    d=fread(fid,[double(numChannels+1) double(numSamples)],'*uint16'); %reads the content into an array
    %fprintf('done reading in %.3fs\n', toc(tReadAll))
    
    dEvntMrk = d((numChannels+1),:);  %additional element in time series (chan+1) is 16-bit event/mark data, where bits 7-14 encode DC09-DC13 triggers
    trigDC09 = bitget(dEvntMrk,7);
    trigDC10 = bitget(dEvntMrk,8);
    trigDC11 = bitget(dEvntMrk,9);
    trigDC12 = bitget(dEvntMrk,10);
    trigDC13 = bitget(dEvntMrk,11);
    trigDC14 = bitget(dEvntMrk,12);
    trigDC15 = bitget(dEvntMrk,13);
    trigDC16 = bitget(dEvntMrk,14);
    
    d=d([goodElec],:);
    %fprintf('Removing offsets... total time %.3f s\n', toc(tReadAll))
    mark1  = int16(d(find(~cellfun('isempty',strfind(actualCode_ALL,'76'))),:));  %mark channels 76 and 77 are signed
    mark2  = int16(d(find(~cellfun('isempty',strfind(actualCode_ALL,'77'))),:));
    
    d_int16=int16(int32(d)+int32(AD_off)); %convert to int16
    
    
    %GAIN_DC = 500 / 2730 ; % E11FFmt.pdf says "Ox0555 corresponds to 500 mV"; Ox555 = 1365;  looks like actually OxAAA-->500mV (2730)
    GAIN_DC = 500 / 1365 ; % E11FFmt.pdf says "Ox0555 corresponds to 500 mV"; Ox555 = 1365;  looks like actually OxAAA-->500mV (2730)
    
    %%- Now it finds the EKG channels and returns the waveforms
    if iscell(ekgTagName) && length(ekgTagName)==2,
        iTemp1 = find(~cellfun('isempty',strfind(actualName_ALL,ekgTagName{1})));
        iTemp2 = find(~cellfun('isempty',strfind(actualName_ALL,ekgTagName{2})));
        if isempty(iTemp1)+isempty(iTemp2)==2,
            fprintf('\nERROR: cant find sync pulse tag %s or %s in the following list:\n', ekgTagName{1},ekgTagName{2});
            disp([listChanStringCode' listActualName']);
            error('revise sync pulse tag name(s) and call eegPulseVisualize again');
        elseif length(iTemp1)>1 & length(iTemp2)>1,
            %designed to catch error with NIH017, where both EKG lines were labeled 'EKG' instead of 'EKG1' and 'EKG2'
            fprintf('\nWARNING: sync pulse tag %s and %s both produce multiple hits... assuming chaCodeString missed suffix:\n', ekgTagName{1},ekgTagName{2});
            iTemp1 = iTemp1(1);
            iTemp2 = iTemp2(2);
        else
            if length(iTemp1)>1, fprintf('\nWARNING: syncpulse tag#1 %s had two hits, taking the first one', ekgTagName{1}); iTemp1=iTemp1(1); end
            if length(iTemp2)>1, fprintf('\nWARNING: syncpulse tag#2 %s had two hits, taking the first one', ekgTagName{2}); iTemp2=iTemp2(1); end
        end
        %- put in a condition where just one of the tags can be missing from this RAW... important when visualizing entire raw folder and some are missing DC9 or 10
        if ~isempty(iTemp1), idx(1) = iTemp1;  EKG1 = double(d_int16(idx(1),:))*GAIN_DC; %convert to volts
        else EKG1 = zeros(size(EKG2)); end
        if ~isempty(iTemp2), idx(2) = iTemp2;  EKG2 = double(d_int16(idx(2),:))*GAIN_DC;
        else EKG2 = zeros(size(EKG1)); end
        
        
    else
        iTemp1 = find(~cellfun('isempty',strfind(actualName_ALL,ekgTagName{1})));
        if isempty(iTemp1),
            fprintf('\nERROR: cant find sync pulse tag %s in the following list:\n', ekgTagName{1});
            disp([listChanStringCode' listActualName']);
            error('revise sync pulse tag name(s) and call eegPulseVisualize again');
        elseif exist('channelNdx', 'var') && channelNdx <= length(iTemp1)
            idx = iTemp1(channelNdx);
        else
            idx = iTemp1(1);
        end
        
        EKG1 = double(d_int16(idx,:)) * GAIN_DC; %convert to volts
        EKG2 = zeros(size(EKG1));
    end
    sampRate = actSamplerate;
    EKG3 = nan(size(EKG1));
    
    
    %%- Grab info from the LOG file if present
    annOut.samp = [];
    annOut.text = {};
    if ~isempty((trigDC10 > 0) | (trigDC11 > 0) | (trigDC12 >0))
        annStruct = nk_parseAnnotation(eegRoot);
        annTimes = [];
        annStr   = {};
        for iAnn=1:length(annStruct),
            annTimes(iAnn) = (annStruct(iAnn).timeSec+1)*actSamplerate;  %-add 1 sample to avoid indexing 0 (could cause problem for end of file?)
            annStr{iAnn}   = sprintf('  %s',annStruct(iAnn).str);
        end
        annOut.samp = annTimes;
        annOut.text = annStr;
    end
    
    
    %%- use the following code to test digital pulses from DC09
    TESTING_DIGITAL_PULSES = 0;
    
    if TESTING_DIGITAL_PULSES,
        thisAx = [];
        
        figure(90); clf
        subplot(5,1,1);  thisAx(end+1) = gca;
        plot(mark1);      title('mark1');
        subplot(5,1,2);   thisAx(end+1) = gca;
        plot(mark2);       title('mark2');
        subplot(5,1,3);   thisAx(end+1) = gca;
        plot(trigDC09,'r');  title('trig DC09');
        subplot(5,1,4);   thisAx(end+1) = gca;
        plot(EKG1,'b'); hold on; plot(trigDC09*200,'r');     title(ekgTagName{1});
        subplot(5,1,5);   thisAx(end+1) = gca;
        plot(EKG2);      title(ekgTagName{2});
        linkaxes(thisAx,'x')
        
        figure(91); clf
        plot(EKG1,'b'); hold on; plot(trigDC09*200,'r');     title(ekgTagName{1});
        
        figure(92); clf
        plot(dEvntMrk,'b'); hold on; plot(trigDC09,'r');     title('all mark events');
        
        filestem = eegRoot;
        if find(trigDC09>0)
            updowns = [trigDC09(1) diff(trigDC09)];
            uptimes = find(updowns==1);
            fprintf('Digital pulse trigger found... extracting sync file:');
            chanfile = sprintf('%s.trigDC09.sync.txt', filestem);
            fchan = fopen(chanfile,'w','l');
            fprintf(fchan,'%d\n',uptimes);
            fclose(fchan);
            fileattrib(chanfile, '+x', 'a'); %JW - change files to "executable"... helps for sorting in mac finder
        end
    end
end
end %grabPulses


function [varargout] = textread2(filename, format_str, varargin)
% quick update of out-dated textread to textscan
fd = fopen(filename, 'r');
c = textscan(fd, format_str, varargin{:});
fclose(fd);
for i = 1:length(c)
    varargout(i) = c(i);
end
if ~iscell(varargout), varargout = {varargout}; end
end

%#ok<*UNRCH>
%#ok<*NOCOL>

