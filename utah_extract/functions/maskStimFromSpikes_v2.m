function maskStimFromSpikes(subj,overrideDataPath56, msGrowDC10pulse, testOnSomeChans, sessionDirPick)
%  function  maskSpikeTimeSeries
%
%   will loop over processed sessions, and will create MASKED versions of any stim session spike_reref or spike_noreref that is found (for sorting stimulation data)
%
%   Inputs:
%      subj -- subject string, "NIH069"
%      overrideDataPath56 -- (standard=''), use this to pass in location of a local copy of the spike data for processing
%      msGrowDC10pulse    -- (standard=5), pads the DC10 & 11 pulse with this value before and after.
%      testOnSomeChans    -- (standard=0), use a positive number to indicate the number of channels to split.  helpful for testing
%      sessionDirPick     == (standard=''), pass a string of a session to process, otherwise code will automatically search all
%
%   -pass in a subject directory on FRNU56
%   -a list of sessions to create mask of... will apply to spikes_reref and spikes_noreref
%   msGrowDC10pulse default 5ms, amout of time to add to beginnging and end of DC10(and/or DC11) pulse
%   testOnFewChan... set to 1,2..5 however many minimum number of channels you want if you are just testing
%
%  JW 2/2019
%

if nargin<5,
    sessionDirPick  = '';
end
if nargin<4,
    testOnSomeChans = 0;
end
if nargin<3,
    msGrowDC10pulse = 5;
end
if nargin<2,
    overrideDataPath56 = '';
end

%- defulat, fprintf status of every channel
VERBOSE = 1;


numMSgrow = msGrowDC10pulse; %- how much to expand the DC10, or DC11 stimulation pulses when creating a mask


% Get subject info
subjectInfo = getMicroSubjInfo_v11(subj,overrideDataPath56);

% convert subject string into subject number for conditional below (DC11 only because a standard stim pulse with NIH062')
subjNum = str2num(subj(4:6));


% Directories for the critical input csv files
fileDir_bookkeeping = fullfile(subjectInfo.dataPath56 ,subj,'_extraction_notes',sprintf('micro_PickSess2Extract.xlsx'));



% Load utahInfo CSV files:  get the "SESSIONS" and the "CHAN STRINGS"
tSessList = readtable(fileDir_bookkeeping);
tSessList = tSessList.folderName(strcmpi(tSessList.toExtract,'y') & strcmpi(tSessList.isStim,'y')); %- only look at the "to extract" sessions that are "isStimSess".

if ~isempty(sessionDirPick),
    tSessList = {sessionDirPick};
    VERBOSE   = 0;  %- usually this call will come from extractMicroPhys... only provide summary outputs if so
    fprintf('\n Checking for masked files in %d session folders', length(tSessList));
else
    fprintf('\n Found %d sessions to extract and isStim... checking for masked files now',length(tSessList));
end


%- Loop over each "to extract" session
for iS=1:length(tSessList),
    
    thisSess    = tSessList{iS};
    
    sessionsDir = fullfile(subjectInfo.dataPath56,subj,'processed_SPK',thisSess);
    fprintf('\n    session: %s',thisSess);
    
    
    
    %- a sessionDir likely contains multiple splits (noreref, reref).
    %    There is a *chance* that these were created at different times with different pulse structs,
    %    where the pulse struct changed because we figured out some other fix/error in the blackrock splitting and concatination functions
    %  to protect against this, save the pulseStruct SEPARATELY for each split, and create a mask separately for each split
    
    
    
    
    %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
    %%- Next, loop over every channel in spikes_reref.  Load, apply mask, and save masked version
    spkDirList = {};
    spkDir = fullfile(sessionsDir,'sort_reref');   if exist(spkDir,'dir'), spkDirList{end+1} = spkDir; end
    spkDir = fullfile(sessionsDir,'sort_noreref'); if exist(spkDir,'dir'), spkDirList{end+1} = spkDir; end
    
    extraFilesDir = '_extractionInfo'; %- sub-directory in each of the spike_XXref folders that has all the extra extraction/alignment info. Underscore so at top of directory in finder
    
    
    
    
    %- loop over target directories
    for iDir=1:length(spkDirList),
        
        %- first get a list of channels to mask
        sourceDir   = spkDirList{iDir};
        chanListGet = dir(fullfile(sourceDir,'*.bin'));  %- will not get anything from pruned...
        if length(chanListGet)==0, continue; end
        
        
        
        %- now grab the pulseStruct for this spkDir and create the mask
        [maskList NO_DC10or11 maskStr] = createStimMask(sourceDir,extraFilesDir,numMSgrow,subjNum,thisSess);
        clear maskData; %- trying to trigger an error below if jackTable doesn't ID proper pulseStruct...
        
        
        %- break out here if there are no pulses for mask-making (but make a dummy target dir so we know what is up)
        if NO_DC10or11,
            sourceDir = spkDirList{1};
            destDir = sprintf('%s(stimMask%d)_SKIPnoDC10or11',sourceDir,numMSgrow);
            if ~exist(destDir,'dir')
                mkdir(destDir);
                saveMaskFigs(destDir,maskList);
            end
            continue; %- jump to the next spkDirList directory
        end
        
        
        %- check for existance of destination dir and channels
        destDir = sprintf('%s(stimMask%d)',sourceDir,numMSgrow);
        if exist(destDir,'dir'),
            chanListPut = dir(fullfile(destDir,'*.bin'));
            if length(chanListPut)==length(chanListGet),
                fprintf('\n All %d channels already masked and saved in %s... skipping this directory',length(chanListPut), destDir);
                continue;
            elseif testOnSomeChans>0 && length(chanListPut)>=testOnSomeChans,
                fprintf('\n  %d of %d test channels already masked and saved in %s... \n  skipping this directory (rerun with testSomeChans==0 later)',length(chanListPut),testOnSomeChans,destDir);
                continue;
            end
        else
            %- create the destination directory and populate the subfolder with info
            mkdir(destDir);
            
            %- create the extra files folder and copy everything over from the source
            [SUCCESS,MSG,MGG2] = copyfile(fullfile(sourceDir, extraFilesDir),fullfile(destDir, extraFilesDir)); %- copy the entire extras folder in one shot
            
            %- add a line to readme about the masking done here
            fidTXT = fopen(fullfile(destDir,extraFilesDir,'_extractMicro_readme.txt'),'a+');
            fprintf(fidTXT, '\n  STIMULATION artifact masked (set to zero) during DC10 or 11 pulses; %d ms margin added before and after each pulse',numMSgrow);
            fclose(fidTXT);
            
            %- save the mask figures in the extras folder
            mkdir(fullfile(destDir,extraFilesDir,'stimMaskFigs'));
            saveMaskFigs(fullfile(destDir,extraFilesDir, 'stimMaskFigs'),maskList);
        end
        
        
        %- binaries exist, now grab the spike-channel jackTable so we can match the currect pulse file if more than one NSP
        jackOut = fullfile(destDir,extraFilesDir,'jacksheetBR_justSplitChan.csv');
        tJackSplit = readtable(jackOut); %- just the split channels contained within this folder..
        
        %- test on a few channels... j
        if testOnSomeChans~=0, chanListGet = chanListGet(1:min([testOnSomeChans length(chanListGet)])); end
        
        fprintf('\n Creating %d masked channel files in %s:', length(chanListGet),destDir); tStart=tic;
        %- loop over channels in this directory.  load; mask; save
        for iChan=1:length(chanListGet),
            if VERBOSE, fprintf('\n %s ',chanListGet(iChan).name); tStartChan=tic; end
            
            %- load the raw time series
            rawSpikeBin = fullfile(sourceDir, chanListGet(iChan).name);
            fid         = fopen(rawSpikeBin,'r');
            spkData     = fread(fid,'int16');
            fclose(fid);
            if ~iscolumn(spkData), spkData=spkData';  end
            
            
            %- for the first channel file (or for any channel if two NSPs), grab the correct mask and make sure the mask duration matches up
            if iChan==1 | size(maskList,1)>1,
                
                %- find this channel in the jacktable and match up the appropriate mask (matters if two NSPs, which gives two masks)
                thisChanStr = chanListGet(iChan).name(1:end-4);
                iJack       = find(strcmp(tJackSplit.SortChanName,thisChanStr));  %- new way, channel names are [NSPoffset]ChanName.bin
                thisNSxFile = tJackSplit.FileName{iJack};
                iMask       = find(contains(thisNSxFile,maskList(:,1)));
                maskData    = maskList{iMask,2};
                
                
                %- now make sure the mask duration matches up, else element-by-element multiplcation will fail below
                nSampMask  = length(maskData);
                nSampSpike = length(spkData);
                if nSampMask<nSampSpike,
                    maskData(end+1:nSampSpike) = 1;
                elseif nSampMask>nSampSpike,
                    maskData = maskData(1:nSampSpike);
                end
                if ~iscolumn(maskData), maskData=maskData'; end
            end
            
            
            %- now mask and saveout result
            rawSpikeBin = fullfile(destDir, chanListGet(iChan).name);
            fid         = fopen(rawSpikeBin,'w+');
            count       = fwrite(fid,spkData.*maskData,'int16');
            fclose(fid);
            if count~=nSampSpike, fprintf('\ Error... fwrite issue?'); keyboard; end
            
            
            %- save a figure showing the result for the first channel
            if iChan==1,
                fprintf(' (making example fig '); tFig = tic;
                %- create a variable representing time in minutes
                dtMS     = 1/30;
                tSampMS  = [0:nSampSpike]*dtMS;
                tSampMin = tSampMS/1000/60;
                tPlot    = tSampMin;  xStr = 'min';
                iTuse    = 1:nSampSpike;
                
                %-  Now plot the data
                hFig = figure(133); clf;
                set(gcf,'color','w','name','spike mask alignemnt','units','norma','position',[.25 .25 .5 .75]);
                
                axList = [];
                axList=[axList subplot(311)];
                plot(tPlot(iTuse),(1-maskData(iTuse))*1,'r'); hold on;
                set(gca,'box','off','tickdir','out','fontsize',15);
                ylabel('mask off/on');
                title(sprintf('%s: DC10andDC11 widened by %d ms to create mask',subj,numMSgrow))
                
                axList=[axList subplot(312)];
                plot(tPlot(iTuse),spkData(iTuse),'k'); hold on
                set(gca,'box','off','tickdir','out','fontsize',15);
                ylabel('voltage (uV)')
                tStr = sprintf('Spike Time Series: %s: %s, %s',subj,thisSess,chanListGet(iChan).name);  tStr(find(tStr=='_'))='-';
                title(tStr);
                
                axList=[axList subplot(313)];
                plot(tPlot(iTuse),spkData(iTuse).*maskData(iTuse),'k'); hold on
                set(gca,'box','off','tickdir','out','fontsize',15);
                xlabel(sprintf('time (%s)',xStr));
                ylabel('voltage (uV)')
                title(sprintf('Masked Time Series: (%s)',maskStr))
                
                linkaxes(axList,'x');
                set(gca,'xlim', tPlot(iTuse([1 end])));
                set(gca,'ylim',[-350 200])
                
                %- save the figure into the "others" folder
                figOutPath = fullfile(destDir, extraFilesDir, 'stimMaskFigs', sprintf('maskExample_%dms.png',numMSgrow));
                fig2pngSimple(hFig,figOutPath)
                pause(2);
                close(hFig);
                
                fprintf(' [%.1fs])',toc(tFig));
            end
            if VERBOSE, fprintf(' [%.1fs]',toc(tStartChan)); end
            
        end %- for iChan
        if VERBOSE, fprintf('\n '); end
        fprintf(' [total duration %.1fs]',toc(tStart));
        
        
    end %- for iDir (noreref, reref)
    
end %- for tSessList

end %- main function



%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%
%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%%%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%% %%%%%%%



%%%%%  KEY HELPER FUNCTION --- create the mask %%%%%%
% create the stim Mask. Will have to make 2 if two NSPs recording spikes
function [maskList NO_DC10or11 maskStr] = createStimMask(targetDir,extraFilesDir,numMSgrow,subjNum,thisSess)
%
%  pass in targetDir, extraFilesDir, numMSgrow -- to pull the pulse files and create the mask
%  pass in subjNum and thisSess                      -- to properly label output figure
%

%- default return: assume nothing worked
maskList = {};  NO_DC10or11 = 1; 


%- code below was looking for one pulse struct at a time... need to pull those out of the pulseStructCell now....
pulseStructPath = fullfile(targetDir,extraFilesDir,'pulseStructCell.mat');
if ~exist(pulseStructPath,'file'),
    fprintf('\n ERROR: unexpected error... cant find pulseStructCell');
    keyboard;
    return;
end
x = load(pulseStructPath);               %- load the matfile. Contains a cell array: e.g., {'INST0'}    {'YkGiDmVfJGdfG-20190712-091113-INST0.ns5'}    {1×1 struct}
pulsesStructCell = x.pulsesStructCell;   % size of cell array = numNSP x 3


%- loop over pulse structs and create a different mask for each
for iNSP=1:size(pulsesStructCell,1),
    
    %- isolate pulses and nsp suffix for this NSP
    pulses        = pulsesStructCell{iNSP,3};
    pulseNSPsuff  = pulsesStructCell{iNSP,1}; %- INST0 or INST1... just used for naming the output figures
    
    if     isfield(pulses,'ain2_ts')  & isfield(pulses,'ain3_ts'),
        dc10 = pulses.ain2_ts;
        dc11 = pulses.ain3_ts;
        NO_DC10or11 = 0;
    elseif isfield(pulses,'ain18_ts') & isfield(pulses,'ain19_ts'),
        dc10 = pulses.ain18_ts;
        dc11 = pulses.ain19_ts;
        NO_DC10or11 = 0;
    elseif isfield(pulses,'din2_ts') & isfield(pulses,'din3_ts'),
        fprintf('\n heads up... using din_ts for the first time. Does it work? Check on transition detection below');
        keyboard; 
        dc10 = pulses.din2_ts; %- this could be from NSP1 or NSP2... for din we dont have a din17 din18, etc in the pulseStrruct
        dc11 = pulses.din3_ts;
        NO_DC10or11 = 0;
    else
        fprintf('\n ERROR: maskStimFromSpikes cant find pulse channel.');
        dc10 = [];
        dc11 = [];
    end
    
    %- could add a check here that IF the ain pulses existed BUT seemed to be missing the recording, then use din.  Not sure if that is a condition though.
    
    if subjNum<62,
        dc11 = dc10; %- DC11 started being used to indicate microstimulation with NIH062
    end
    
    maskStr = sprintf('stimPulse grown +/-%d ms',numMSgrow); %- add ~5ms before and after DC10/11 to completely mask the stim
    
    
    %- Pulse threshold detector from makeJacksheetBR
    %- find the range that doesnt include outliers, maybe 5 to 95 percentile, (looks like need 99th percentile to detect pulses)
    %    if that is >1V, split put the thresh in the middle;  if <1V, then probabably just say all zeros?
    bncData     = double([dc10; dc11]);  iBNC=[1 2];
    bncData_mV  = bncData*5000/32764;      %- convert to mV so threshold can be set in terms of mV
    Y = prctile(bncData_mV,[5 10 99 100],2);
    rngInner = Y(:,3)-Y(:,2); %- try to avoid outliers
    rngOuter = Y(:,4)-Y(:,1); %- if dont see pulses with 99th percentile, maybe here?
    threshUse = nan(size(iBNC));
    for iiBNC=1:length(iBNC),
        if ~isempty(rngInner)
            if rngInner(iiBNC)>1000,
                threshUse(iiBNC)  = Y(iiBNC,2)+rngInner(iiBNC)/2;
            elseif rngOuter(iiBNC)>1000,
                threshUse(iiBNC)  = Y(iiBNC,1)+rngOuter(iiBNC)/2;
            else
                threshUse(iiBNC) = 6000; %- range 5000, so this should result in no pulses;
            end
            bncPulseCnt(iBNC(iiBNC)) = length(find( bncData_mV(iiBNC,2:end)>=threshUse(iiBNC) & bncData_mV(iiBNC,1:end-1)<threshUse(iiBNC)));
        else
            bncPulseCnt(iBNC(iiBNC)) = 0 ;
        end
    end
    
    
    
    %- loop over stim DC channels and create binary version. attempt to detect and remove long runs of "high" due to stimulator being unplugged
    for whichDC=[10 11],
        
        %- bookkeeping
        iBNC = whichDC-9;
        if whichDC==10, thisDC = dc10;
        else            thisDC = dc11;
        end
        
        dcBin = zeros(size(thisDC));  iUp = [];
        if bncPulseCnt(iBNC)>10,%- sometimes disconnected value is >500 the whole time...
            %if range(thisDC)>5000, %- sometimes disconnected value is >500 the whole time...
            %dcBin(thisDC>max([0.9*max(thisDC) 500])) = 1;
            dcBin(thisDC>threshUse(1)) = 1;
            iUp   = find(dcBin(1:end-1)==0 & dcBin(2:end)==1);
            iDown = find(dcBin(1:end-1)==1 & dcBin(2:end)==0);
            if length(iUp)>1,       dcBin(1:iUp(1))=0;       end  %- set section before first up to 0
            if iDown(end)<iUp(end), dcBin(iUp(end):end) = 0; end  %- set section after last up to 0
            dcBin(iUp-[0:numMSgrow-1]')=1;  dcBin(iDown+[1:numMSgrow]')=1;
            if length(iUp)+length(iDown)<10,
                dcBin = zeros(size(thisDC)); %- real stim session should have MANY more pulses than 5.. trying to catch unplugged bnc
            end
        end
        
        %- bookkeeping
        if whichDC==10, dc10bin = dcBin;  dc10up = iUp;
        else            dc11bin = dcBin;  dc11up = iUp;
        end
        
    end
    
    
    
    %-upsample binary mask to 30kS so its matched with the spike time series
    maskData=zeros(length(dc11bin)*30,1);
    iHi = find(dc10bin>0 | dc11bin>0);
    if length(iHi)>0,
        maskData(iHi*30+[0:29]')=1;
    end
    maskData = 1-maskData; %-invert for subsequent steps... so now high means
    
    
    %- little sanity check
    pcntSessMasked = 100.0 * length(iHi)/length(dc10bin);
    maskSummaryStr = sprintf('\n MASK applied to %.1f%% of the %.1f min session (%d total stim pulses detected)',pcntSessMasked,length(dc10)/1000/60, length(iHi));
    disp(maskStr);
    if pcntSessMasked > 30,
        fprintf('\n WARNING: mask shouldnt be on that often. Should be way less than <30%%.  About to break to figure it out');
        %- issue... if stimulation stop and possibly if stimulator is disconnected, DC10/11 can transition high and stay there
        MASK_WARNING = 1;
    else
        MASK_WARNING = 0;
    end
    
    
    %- make a plot to check it
    if 1 | MASK_WARNING,
        tMin = [1:length(dc10)]/1000/60;
        hFigMask = figure(10+iNSP); clf
        set(gcf,'color','w','name','Mask Construction');
        
        subplot(311)
        plot(tMin,dc10/max(dc10),'b'); hold on;
        plot(tMin,.1+dc10bin*.8,'r');
        plot(tMin(dc10up),ones(size(dc10up)),'*')
        ylabel('dc10')
        set(gca,'box','off','tickdir','out')
        strCln = sprintf('NIH%03d :: %s',subjNum,thisSess);  strCln(find(strCln=='_'))='-';
        title(strCln);
        
        subplot(312)
        plot(tMin,dc11/max(dc11),'b'); hold on;
        plot(tMin,.1+dc11bin*.8,'r');
        plot(tMin(dc11up),ones(size(dc11up)),'*')
        ylabel('dc11')
        set(gca,'box','off','tickdir','out')
        
        subplot(313)
        tMinMask =  [1:length(maskData)]/30/1000/60; %- mask is at 30kS, so use a different time var here
        plot(tMinMask(1:30:end),maskData(1:30:end)/max(maskData),'b'); hold on;
        plot(tMinMask(1:30:end),.1+maskData(1:30:end)*.8,'r');
        ylabel('MASK on/off')
        set(gca,'box','off','tickdir','out')
        title(maskSummaryStr);
        
        if MASK_WARNING, keyboard; end
        
        
        %- figure will get saved outside of this function
        
        %- save the figure to "other" folder, where the pulseStruct lives
        %figOutPath = strrep(pulseMat,'pulseStruct','maskCreation');
        %figOutPath = strrep(figOutPath,'.mat','.png');
        %fig2pngSimple(hFigMask,figOutPath)
        %pause(2);
        %close(hFigMask);
    end
    
    %- create cell array output struct with mask for each NSP
    maskList{iNSP,1} = pulseNSPsuff;   %- NSP suffix "INST0"... example pulse file name is "pulseStruct_INST0.mat"
    maskList{iNSP,2} = maskData;       %- mask data
    maskList{iNSP,3} = hFigMask;       %- figure handle for saving figure outside of this function
end
end %- function createStimMask



%%%- NEXT HELPER... saving the figures
function saveMaskFigs(destDir,maskList)

for iFig=1:size(maskList,1),
    
    hFigOut = maskList{iFig,3};
    if ~isempty(hFigOut) && isgraphics(hFigOut),
        figOutPath = fullfile(destDir,sprintf('maskCreation_%s.png',maskList{iFig,1}));
        fig2pngSimple(hFigOut,figOutPath);
        pause(2);
        close(hFigOut);
    end
end
end % function saveMaskFigs


