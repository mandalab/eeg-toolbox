function fixSPKextractionInfo(subj, tJackSplit, thisSess, isStimSess, processedDirUse, rawDir, subjectInfo, fidOut, compMemoryGB)

processedDir = processedDirUse;

devNums = unique(tJackSplit.MicroDevNum);
extractOutputs = {'sort_reref' 'sort_unfilt'};
GET_COMMON_AVE = 1;
extraFilesDir   = '_extractionInfo';
commonAveDir = fullfile(processedDir, 'sort_reref',extraFilesDir,'rerefStuff'); 



% Get Filter parameters
samprate_br  = 30000;
samprate_LFP = 1000;
samprate_nk  = 1000;
CONVERT_2_uV = 0;


butterOrder         = 2;

% Spikes
samprate_spikes     = samprate_br;
nyquist_freq_spikes = (samprate_spikes/2);
spikeBand           = [300 3000];
wo_spike=spikeBand/nyquist_freq_spikes; % bandpass 300Hz to 3000HZ %mountainsort: 600 to 6000 Hz
[b_spike,a_spike]   = butter(butterOrder,wo_spike);
str_spike           = sprintf('Spike Band: butterworth bandpass %d-%d Hz; filtfilt for zero phase shift (effective order %d)',spikeBand,butterOrder*2); %- filtfilt doubles order

% Spikes - unfiltered
unfiltBand          = [1,5000];
wo_spike_unfilt     = unfiltBand/nyquist_freq_spikes; % "unfiltered": 1Hz - 5000Hz
[b_spike_unfilt,a_spike_unfilt] = butter(butterOrder,wo_spike_unfilt);
str_spike_unfilt    = sprintf('"Unfilt" Spike: butterworth bandpass %d-%d Hz; filtfilt for zero phase shift (effective order %d)',unfiltBand,butterOrder*2); %- filtfilt doubles order





%- Identify which NSPs the sorted channels come from, and include a pulse structure (for subsequent alignemnt) for any NSP included
fileNames = unique(tJackSplit.FileName,'stable'); %- stable to keep in order of jacksheetBR.  probably doesn't matter tho

%- add a bit of info to the command line output:  #channels and recording duration
fprintf(1,      '\n          -- Extracting %d micro chan; %d devNums; %d ns5/6 files; %.1fs min recording: ', height(tJackSplit), length(devNums), length(fileNames), tJackSplit.DurationMin(1));
fprintf(fidOut, '\n          -- Extracting %d micro chan; %d devNums; %d ns5/6 files; %.1fs min recording: ', height(tJackSplit), length(devNums), length(fileNames), tJackSplit.DurationMin(1));


%- create a single date/time string for this meta extraction... this can be used to link together all the extractions (LFP, noreref, reref)... 
%    in the past I grabbed a fresh datestr when each extraction type was created, but this could lead to confusion about whether the same extraction code was used for all
sessExtractDateStr = datestr(now,21); 


%- Loop over ns5/6 files and extract
pulsesStructCell   = {};
postProcStructCell = {};
fprintf(1, '(making pulse struct:'); ticPulse = tic;
for iF = 1:length(fileNames),

    %- NSP suffix from jackTableBR... from a row with this fileName
    nspSuffix = tJackSplit.NSPsuffix{find(strcmp(tJackSplit.FileName,fileNames{iF}),1,'first')};

    if contains(fileNames{iF},'ns5')
        fileType = 'ns5';
    elseif contains(fileNames{iF},'ns6')
        fileType = 'ns6';
    end
    ns3File = strrep(fileNames{iF},fileType,'ns3');  ns3File = fullfile(rawDir,ns3File);
    nevFile = strrep(fileNames{iF},fileType,'nev');  nevFile = fullfile(rawDir,nevFile);
    if ~exist(ns3File); ns3File = strrep(ns3File,'ns3','ns4');  end


    %- extract the pulse up times and pulse time series, to append to LFP.mat for subsequent alignment
    [pulses] = makePulsesStruct('ns3_fpath',ns3File,'nev_fpath',nevFile);


    %- no longer saving the individual pulse structs... instead save the cell array of pulse structs in each extraction folder (e.g., noreref, reref, etc)

    %- this cell array will be attached to the LFP.mat below; and also save to processed_SPK/_extaction_notes/pulseStructCell.mat
    pulsesStructCell{iF,1} = nspSuffix;
    pulsesStructCell{iF,2} = fileNames{iF};
    pulsesStructCell{iF,3} = pulses;
end
fprintf(1,'%.1fs)',toc(ticPulse));


%- now the big loop over nsX files.  will be >1 if using 2 NSPs, or for some reason mixed ns6 and ns5 data (did that with NIH050 microwires)
for iF = 1:length(fileNames)

    %- these are the channels that we plan to split from this file
    tJackFileSplit = tJackSplit(strcmp(tJackSplit.FileName,fileNames{iF}),:);

    thisNSPsuffix  = tJackFileSplit.NSPsuffix{1};


    %- list of device number for the channels to split
    devNums  = unique(tJackFileSplit.MicroDevNum);
    devNames = {};
    for iiDev=1:length(devNums)
        iDev = find(subjectInfo.deviceNum==devNums(iiDev));
        devNames{end+1} = subjectInfo.newChanName{iDev,1}; %- create cell array of device names matching the device nums on this NSPs file
    end

    %- pull a "raw" jacksheet from the file of interest... best way to be 100% sure that channel names are matching channel numbers
    NSxFileDir     = fullfile(rawDir,fileNames{iF});
    dataFromFile   = openNSx(NSxFileDir,'noread');
    labelsFromFile = {dataFromFile.ElectrodesInfo.Label};
    labelsFromFile = cellfun(@deblank,labelsFromFile,'UniformOutput',false); % remove whitespace

    %- index of channels to pull from nsX
    iChanRead      = find(ismember(labelsFromFile,tJackFileSplit.ChanName));

    % determine how many channels to load at once based on computer memory
    [chanGroups, fileSizeGB] = getChanGroupsToLoad(NSxFileDir,iChanRead,compMemoryGB);
    if isempty(chanGroups), fprintf('\n uh oh... not enough memory to load a single channel'); keyboard; end

    %- extract group by group and compute the commonAve for each unique device in the file
    [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
    fprintf(fidOut,'\n LOADING ALL DATA TO COMPUTE COMMON AVERAGE BEFORE SPLITTING FILES  <current Memory: used %.1f GB, free %.1f GB, compressed %.1f GB>',memUsedGB,memFreeGB,memCompressGB);
    fprintf(fidOut,'\n [filesize %.0f GB; %d groups of <=%d chan to load]', fileSizeGB, length(chanGroups), length(chanGroups{1}));
    for grp = length(chanGroups):-1:1, %- start with last group because it can have fewer channels and we are gonna skip loading below (so want to end with most channels)

        %- load the data (slow step)
        fprintf(fidOut,'\n   loading group %d (%d chan)',grp,length(chanGroups{grp})); tLoad = tic;
        fprintf(1,'.');
        clear NSxData; %- make room because openNSx is going to use a lot of temp space?
        [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
        NSxData = concatOpenNSx(NSxFileDir,CONVERT_2_uV,1,chanGroups{grp}); % The zero should actually be CONVERT_2_uV


        %- grab the concat NSx info... that way if things change in the future (i.e., gain; or whether we throw away "junk", etc... we'll know what was in this one)
        physio_nsx_postProc = NSxData.postProc;  %- get the "post process" string returned from concatOpenNSx and save it to LFP and spike structures
        postProcStructCell{iF,1} = thisNSPsuffix;
        postProcStructCell{iF,2} = fileNames{iF};
        postProcStructCell{iF,3} = NSxData.postProc; %- will be filled out below

        NSxData = NSxData.Data;
        [memUsedGB2, memFreeGB2, memTotalGB2, memCompressGB2] = getMemoryUsage();
        fprintf(fidOut,' [%.1fs; %.1fGB --> %.1f GB used;  %.1f -> %.1f GB free]',toc(tLoad),memUsedGB, memUsedGB2,memFreeGB,memFreeGB2);
        fprintf(1,'.');

        %- now that we know how long the time series is we can initailize the commonAve time series (1 per device)
        if grp==length(chanGroups),
            numSampBR = size(NSxData,2);
            commonAve = zeros(length(devNums),numSampBR);
            cAveCount = zeros(length(devNums),1);
            exampChan = commonAve; %- single trace from each device for the output figure
        end

        %- create a list of the channel names and device numbers that were just loaded
        chanStrLoad = labelsFromFile(chanGroups{grp});
        chanDevLoad = []; chanIsStim = []; %chanNotStim = []; - SJ changed to IsStim, chanNotStim not used anywhere else
        for iC=1:length(chanStrLoad), 
            chanDevLoad(iC) = tJackFileSplit.MicroDevNum(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC})); 
            chanIsStim(iC)  = tJackFileSplit.microStimPigtail(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC}));  %- 0 for non-stim, >0 for stim
        end;

        %- chanGroup might contain multiple devices, so have each channel contribute to the ave of its own device
        if GET_COMMON_AVE,
            for iD=1:length(devNums),
                iChanD = find(chanDevLoad==devNums(iD) & chanIsStim==0); %- all channels from this device EXCLUSING microStim channels, which have weird noise even if not stim session
                if length(iChanD)>0,
                    if cAveCount(iD)==0,
                        exampChan(iD,:) = double(NSxData(iChanD(1),:)); %- save the first trace from each device to output in a figure below
                    end
                    commonAve(iD,:) = commonAve(iD,:)+sum(NSxData(iChanD,:),1); %- int16 can sum, dont need to convert to double here
                    cAveCount(iD)   = cAveCount(iD)+length(iChanD);
                end
            end % for iD
        end
    end % for grp


    %- save out 30kS common average if outputing spikes or lfp_reref
    if GET_COMMON_AVE,

        %- output the results
        fprintf(fidOut,'\n saving common average time series(s) '); tStart = tic;
        for iD=1:length(devNums),
            %- convert to a mean
            commonAve(iD,:) = commonAve(iD,:) / cAveCount(iD);

            % Save common average as .bin that could be loaded to plexon (only getting saved in processed, not processedLFP, because we only save noreref for LFP
            commonAvFileName_raw = sprintf('globalAvg_30kRaw_(%s)_[dev%d].bin', devNames{iD}, devNums(iD));
            fid = fopen(fullfile(commonAveDir,commonAvFileName_raw), 'w');
            fwrite(fid,commonAve(iD,:),'int16');  %Save for plexon... confirm this is saving in the correct format
            fclose(fid);

            % Save the reref figure for one channel... probably should do this separately for each device
            figSaveDir = fullfile(commonAveDir,sprintf('rerefPlot_(%s)_[dev%d]', devNames{iD}, devNums(iD)));
            getRerefFigs(subj,exampChan(iD,:),commonAve(iD,:),figSaveDir,'');
        end
        clear exampChan;
        fprintf(fidOut,' DONE [%.1fs] !\n', toc(tStart));
        fprintf(1,'.');

        %- sanity check figure...
        MAKE_GLOBAL_FIG = 0;
        if MAKE_GLOBAL_FIG,
            figure(1);clf
            for i=1:size(commonAve,1)
                hold on
                subplot(1,size(commonAve,1),i)
                plot(NSxData(:,1:1000:end)'); hold on;
                plot(commonAve(1,1:1000:end),'k','linewidth',4); hold on;
                tStr = sprintf('%s%s%s',processedDir,'-array ', subjectInfo.oldChanName{i,1});
                tStr(tStr=='_')='-';
                title(tStr);
            end
            keyboard
            % from the examplie I just looked it (NIH066), seems like it would make more sense to regress out the common signal rather than subtract from everything...
            %  most follow common but some do not.
        end

    end  %- if GET_COMMON_AVE



    %- Now loop over each channel, and output the requested data
    %

    %- create matricies to hold the LFP data... it will all be written at once after processing all channels
    if any(contains(extractOutputs,'lfp')),
        fprintf(fidOut,'\n Initalizing LFP output matricies.  [%.1f GB --> ', getMemoryUsage()); tInitLFP = tic;
        lfpDS = [1:(samprate_br/samprate_LFP):numSampBR];
        numSampLFP = length(lfpDS);
        if any(strcmp(extractOutputs,'lfp_noreref')), lfp_noreref = zeros(length(iChanRead),numSampLFP,'int16'); end
        clear tLFPdata; %- to become a table of lfp channels
        fprintf(fidOut,' %.1f GB;  %.1fs]', getMemoryUsage(), toc(tInitLFP));
    end

    %- time series of zeros used in place of common average for "noreref" outputs
    zeroCommon = zeros(1,numSampBR);

    [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
    fprintf(fidOut,'\n RE-LOADING DATA TO SPLIT FILE  <current Memory: used %.1f GB, free %.1f GB, compressed %.1f GB>',memUsedGB,memFreeGB,memCompressGB);
    fprintf(fidOut,'\n [filesize %.0f GB; %d groups of <=%d chan to load]', fileSizeGB, length(chanGroups), length(chanGroups{1}));
    for grp = 1:1 % ONLY DOING THIS ONCE

        %- first group currently loaded and in memory, so start reload on group 2 (avoiding a reload of the same data)
        if grp > 1,
            fprintf(fidOut,'\n   loading group %d (%d chan)',grp,length(chanGroups{grp})); tLoad = tic;
            clear NSxData; %- make room because openNSx is going to use a lot of temp space?
            [memUsedGB, memFreeGB, memTotalGB, memCompressGB] = getMemoryUsage();
            NSxData = concatOpenNSx(NSxFileDir,CONVERT_2_uV,1,chanGroups{grp});
            NSxData = NSxData.Data;
            [memUsedGB2, memFreeGB2, memTotalGB2, memCompressGB2] = getMemoryUsage();
            fprintf(fidOut,' [%.1fs; %.1fGB --> %.1f GB used;  %.1f -> %.1f GB free]',toc(tLoad),memUsedGB, memUsedGB2,memFreeGB,memFreeGB2);
            fprintf(1,'.');

            %- create a list of the channel names and device numbers that were just loaded
            chanStrLoad = labelsFromFile(chanGroups{grp});
            chanDevLoad = []; for iC=1:length(chanStrLoad), chanDevLoad(iC) = tJackFileSplit.MicroDevNum(strcmp(tJackFileSplit.ChanName,chanStrLoad{iC})); end;
        end


        %- now loop over each channel separately and create all the split offs
        fprintf(fidOut,'\n splitting channels:');
        for iC = 1:1 % ONLY DOING THIS ONCE

            tStartChan     = tic;
            chanDataRaw    = double(NSxData(iC,:));
            thisChanIndex  = chanGroups{grp}(iC);
            thisChanStrOld = chanStrLoad{iC};            %-channel name
            thisChanStrOut = tJackFileSplit.SortChanName{strcmp(tJackFileSplit.ChanName,thisChanStrOld)};  %- spike-sorting.bin filename [nspSuffix]ChanName.bin (use the new table column to make sure single source of this name)
            thisChanDevNum = chanDevLoad(iC);            %-actual device number, used for file name outputs
            thisChanDevInd = find(devNums==thisChanDevNum); %-inded into devNums, which is also index into commonAve

            fprintf(fidOut,'\n   %s (old string=%s; devNum=%d; indexThisNSP=%d)',thisChanStrOut,thisChanStrOld,thisChanDevNum,thisChanIndex);


            %- NSxData is 30kS raw data.
            for iOut=1:length(extractOutputs),
                switch extractOutputs{iOut},

                    case 'sort_unfilt',
                        %- sort_unfilt.  subtract the global mean, filter 1-5000, output the channel
                        numeratorCoeff   = b_spike_unfilt;
                        denominatorCoeff = a_spike_unfilt;
                        subCommon        = commonAve(thisChanDevInd,:);
                        strFilt          = str_spike_unfilt;
                        processedDirUse  = processedDir;

                    case 'sort_reref',
                        %- sort_reref.  subtract the global mean, filter 300-3000, output the channel
                        numeratorCoeff   = b_spike;
                        denominatorCoeff = a_spike;
                        subCommon        = commonAve(thisChanDevInd,:);
                        strFilt          = str_spike;
                        processedDirUse  = processedDir;

                    case 'sort_noreref',
                        %- sort_reref.  subtract the global mean, filter 300-3000, output the channel
                        numeratorCoeff   = b_spike;
                        denominatorCoeff = a_spike;
                        subCommon        = zeroCommon;
                        strFilt          = str_spike;
                        processedDirUse  = processedDir;

                    case 'lfp_noreref',
                        %- lfp_noreref.   low pass filter and downsample
                        numeratorCoeff   = b_LFP;
                        denominatorCoeff = a_LFP;
                        subCommon        = zeroCommon;
                        strFilt          = str_LFP;
                        processedDirUse  = processedDirLFP;

                    case 'sort_raw',
                        %- dont filter. just extract
                        numeratorCoeff   = nan;
                        denominatorCoeff = nan;
                        subCommon        = zeroCommon;
                        strFilt          = 'raw spike time series split from ns5 or ns6 without reref or filtering';
                        processedDirUse  = processedDir;
                end

                fprintf(1,'.'); %- let the user know something is happening...

                %- possible to filter all channels at once?  one at a time takes 10x longer than just saving out the channel (30 vs 3sec per channel)
                %keyboard
                %fprintf('\n filering all channels at once (super memory intensive)'); tFiltStart = tic;
                %filteredAll = [filtfilt(numeratorCoeff,denominatorCoeff,[double(NSxData)-(ones(length(chanGroups{grp}),1)*subCommon)]')]'; %- filtfilt operates along first non-zero dimension, so make it time x channels then convert back
                %fprintf(' [%.1f]',toc(tFiltStart));

                %- apply zero-phase-shift filter...
                if strcmp(extractOutputs{iOut}, 'sort_raw')
                    filtered = chanDataRaw; %- just split the channels without reref or filter
                else
                    filtered = filtfilt(numeratorCoeff,denominatorCoeff,chanDataRaw-subCommon);
                end

                %- different steps for LFP and SPIKES
                if contains(extractOutputs{iOut},'lfp'),
                    %- downsample the LFPs
                    tLFPdata(thisChanIndex,:) = tJackFileSplit(strcmp(tJackFileSplit.ChanName,thisChanStrOld),:); %- make a table matched to the LFP data
                    if     strcmp(extractOutputs{iOut},'lfp_noreref'), lfp_noreref(thisChanIndex,:) = int16(filtered(lfpDS));
                    elseif strcmp(extractOutputs{iOut},'lfp_reref'),   lfp_reref(thisChanIndex,:)   = int16(filtered(lfpDS)); end

                else
                    %- spike output .bin to disk
                    spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},sprintf('%s.bin',thisChanStrOut));
                    fchan = fopen(spkBinFile,'w','l');
                    fwrite(fchan,filtered,'int16'); %- convert to int16 before writing
                    fclose(fchan);
                end

                %- save a copy of the jacksheet table in each output directory.  LFPs also get a copy embedded in .mat
                if iF ==1 & grp ==1 & iC == 1,    %- first channel dumped... only need to save the table once per outputDir
                    if contains(extractOutputs{iOut},'sort'),
                        extraFilesDirUse = extraFilesDir; %- for spikes put all extraction/alignment stuff in a subfolder
                    else
                        extraFilesDirUse = ''; %- for LFPs there is not much in the folder so put them there
                    end
                    
                    if ~exist(fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse),'dir')
                        mkdir(fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse))
                    end
                    
                    jackOut = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'jacksheetBR_justSplitChan.csv');
                    writetable(tJackSplit,jackOut); %- just the split channels contained within this folder..

                    %- also make a little helper readme file with extraction date/time, computer name, filter settings
                    fidTXT = fopen( fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'_extractMicro_readme.txt'), 'a+');
                    fprintf(fidTXT, '\n\n\n**** extract MicroPhys on %s ***',  sessExtractDateStr);
                    fprintf(fidTXT, '\n  full directory path: %s ',       fullfile(processedDirUse,extractOutputs{iOut}) );
                    fprintf(fidTXT, '\n  run on computer %s, by user %s ',char(java.net.InetAddress.getLocalHost.getHostName),char(java.lang.System.getProperty('user.name')) );
                    fprintf(fidTXT, '\n  filtering:  %s ',          strFilt);
                    if all(subCommon==0),
                        fprintf(fidTXT, '\n  referencing: NO common average rerefrence subtracted during extraction');
                    else
                        fprintf(fidTXT, '\n  referencing: Common average rereference (by device) subtracted before filtering');
                        fprintf(fidTXT, '\n               %d stimulation micro channels removed from reference groups', sum(tJackFileSplit.microStimPigtail>0));
                        fprintf(fidTXT, '\n               total num channels per reference (device) group: ');
                        for iii=1:length(cAveCount), fprintf(fidTXT, '%d, ',cAveCount(iii)); end
                    end
                    fclose(fidTXT);

                    %- if spike folder, save the pulseStruct and postProcStruct so it can be appended from this folder instead of requiring "other"
                    if contains(extractOutputs{iOut},'sort'),

                        %- output copy of filtered first channel
                        spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('%s.bin',thisChanStrOut));
                        fchan = fopen(spkBinFile,'w','l');
                        fwrite(fchan,filtered,'int16'); %- convert to int16 before writing
                        fclose(fchan);

                        %- output copy of raw first channel
                        spkBinFile = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('%s_raw.bin',thisChanStrOut));
                        fchan = fopen(spkBinFile,'w','l');
                        fwrite(fchan,chanDataRaw,'int16'); %- convert to int16 before writing
                        fclose(fchan);

                        %- and save a full copy of the pulseStructCell.mat, to be used with alignment
                        pulseCellOut  = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('pulseStructCell.mat',nspSuffix));
                        save( pulseCellOut, 'pulsesStructCell',  '-v7.3');

                        %- and the postProc cell... THIS GIVES DEEPER INFO ABOUT THE NS5/6 SPLIT (openNSx)
                        postProcCellOut  = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,sprintf('postProcStructCell.mat',nspSuffix));
                        save( postProcCellOut, 'postProcStructCell',  '-v7.3');
                    end
                    
                    % Get rid of previous ones
                    good_jack = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'jacksheetBR_justSplitChan.csv');
                    bad_jack = fullfile(processedDirUse,extractOutputs{iOut},'jacksheetBR_justSplitChan.csv');
                    good_readme = fullfile(processedDirUse,extractOutputs{iOut},extraFilesDirUse,'_extractMicro_readme.txt');
                    bad_readme = fullfile(processedDirUse,extractOutputs{iOut},'_extractMicro_readme.txt');

                    if exist(good_jack,'file') && exist(bad_jack,'file')
                        delete(bad_jack)
                    end
                    if exist(good_readme,'file') && exist(bad_readme,'file')
                        delete(bad_readme)
                    end

                end %- copy files

            end %- for extractionOut

            %- sanity check
            if 0 & iC==1,
                figure(1); clf;
                subplot(211);  plot(double(lfp_noreref(2,1:1000)),'r'); title('noreref');
                subplot(212);  plot(double(lfp_reref(2,1:1000)),'r');   title('reref');
            end

            fprintf(fidOut,' [%.1fs]', toc(tStartChan));

        end %- channels in this chanGroup

        clear NSxData;

    end %- for all chanGroups

end %- for iF = 1:length(fileNames)   loop over file stems (i.e., INST0 and INST1)

%fprintf(1     ,'  Extraction finished for %s [%.1f min]\n',thisSess,toc(tStartSess)/60);
%fprintf(fidOut,'  Extraction finished for %s [%.1f min]\n',thisSess,toc(tStartSess)/60);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the stim session masks using default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- if a freshly extracted stim session, create a masked version of the spike time series for sorting.
if isStimSess & any(contains(extractOutputs,'sort')) %2/24/20 - SJ changed 'spike' to 'sort'
    %- now loop over extracted sessions and created masked version of splits for any stim session
    testOnSomeChans = 0;  pickSessDir = thisSess; %- JW found that 5ms margin worked well with NIH069 and 66...
    maskStimFromSpikes_v2(subj,overrideDataPath56, maskStimMsGrowDC10pulse, testOnSomeChans, pickSessDir);
end



end

