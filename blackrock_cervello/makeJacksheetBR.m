function  jackTable = makeJacksheetBR(rawDirPath,rootEEGdir,overwrite,GET_RANGE_and_PULSES,SKIP_KEYBOARD)
%
%  goal of this function is to create a table that list ALL recorded channels, sample rates, and physical electrodes
%   this will be created for each raw directory, and can be used to feed a masterJack sheet or other stuff
%   intent is to run it on FRNU56, but can be run in an EEG/raw directory as well (but will not have info about ns5/6 data not in the dir
%
%  So far seems to work on Cervello-controlled data.  Try on some old patients.
%      And then try pulling actual range and pulses from analog files... if too slow than skip that
%      then run on all subjects, in FRNU and FRNU56.  note if created from complete set
%
%
%--inputs:         --
%   rawDirPath can be local or on FRNU56
%      ex '/Volumes/56C/UTAH_C/NIH061/data_raw/180710_1443'
%      ex '/Volumes/JW24TB/data24TB/eeg_new/NIH066/raw/180712_1559'
%
%--optional inputs:--
%   rootEEGdir, path to local copy of subject, only used if rawDirPath is on FRNU56, in which case a copy is made locally
%      ex '/Volumes/JW24TB/data24TB/eeg_new/'
%      ex ''  % pass in empty string to bypass this functionality
%
%   overwrite: (0), 1, or -1:  if 0, skips processing if jacksheet already exists
%                              if -1, creates new jacksheet even if it already exists, but still not overwrite it
%
%   GET_RANGE_and_PULSES: (1) or 0, use 0 for a faster first-pass of the subject, 1 calculated range and count for all pulses
%
%   SKIP_KEYBOARD: (0).  keyboards mostly happen in microNameConvert, where it is iffy if things are working corectly.  Pass 1 if you know all those keyboards are OK.  
%
%
%-- EXAMPLE USAGE
%
%   (1) to quickly create a local jacksheetBR for a specific raw directory of a local subject (for phys2chan mapping).
%    jt = makeJacksheetBR('/Volumes/JW24TB/data24TB/eeg_new/NIH066/raw/181103_1507','',0,0);
%
%
%   (2) to create a "complete" jacksheetBR on FRNU56 that measures voltage ranges and pulse counts on ALL channels (including ns5)
%    jt = makeJacksheetBR('/Volumes/56C/UTAH_C/NIH059/data_raw/180405_1838','',0,1);
%
%
%   (3) to copy the "complete" jacksheets to a local copy (will copy even if already made)
%    jt = makeJacksheetBR('/Volumes/56C/UTAH_C/NIH059/data_raw/180405_1838','/Volumes/JW24TB/data24TB/eeg_new/NIH066/raw/181103_1507',0,1);
%
%
%%%%%%%%%%%%%
%
% JHW 12/2018  created it
% JHW 3/2019  updated... added some fields (e.g. MicroDevNum), changd the range calculation (now 5-95% for physio, .001-99.999 for digital
% SNJ 5/11/2020: New input for overwrite. Previously, overwrite specified as "0" would just grab the
%                jacksheetBR_local that already in the folder if it existed. By specifying overwrite as 
%                -1, now it will create a new jacksheetBR_local even if it already exists, but still NOT 
%                overwrite it. This is to perform a check in split_BR that ensures the "ChanName" column 
%                has not been tampered with.
% SNJ 7/16/2020: - Changed isspace to deblank for electrode labels (isspace does not work because the extra
%                  characters are actually cntrl, not white space)
%                - Added in a few more checks for overwrite
%                - changed instances of readtable to readtableSafe
%

%- rootEEGdir is optional... it is used to push a "complete" copy to a subject folder if the subject folder exists
if nargin<5
    SKIP_KEYBOARD=0;  %- (0).. catches some potential errors in microNameConvert that require a user check. set to 1 if confident all issues are worked out for a subject
end
if nargin<4
    GET_RANGE_and_PULSES=1;  %- (1).. use 0 when trying to debugrange of analog data for each channel, pulse counts for BNC inputs
end
if nargin<3
    overwrite = 0; %- default is to skip the folder if a complete table already exists
end
if nargin<2
    rootEEGdir = '';
elseif contains(rawDirPath,'data_raw')
    %- check to see if rootEEGdir/subj contains the same raw dir, if so, put a copy there as well
    if length(rootEEGdir)>0 && ~exist(rootEEGdir,'dir')
        fprintf('\n error, rootEEGdir specified but doesnt exist: %s',rootEEGdir);
        keyboard;
        return;
    end
end
%SJ added this in to be safe
if overwrite ~= 1 && overwrite ~= 0 && overwrite ~= -1
    fprintf('%s\n',['ERROR!! overwrite = ' num2str(overwrite) ' (overwrite variable should be equal to 1, 0, or -1).']);
    keyboard
end

jackTable = [];

%- define the target file
if GET_RANGE_and_PULSES==0, pulseSuffix = '_noPulses'; else pulseSuffix = ''; end
if contains(rawDirPath,'data_raw')
    targetFile      = sprintf('jacksheetBR%s_complete.csv',pulseSuffix); %- what will be written with current settings
    targetFileIdeal = 'jacksheetBR_complete.csv';  %- this version takes more time, but if exists use it
    targetFiles2cut = {'jacksheetBR_noPulses_complete.csv'};
else
    targetFile = sprintf('jacksheetBR%s_local.csv',pulseSuffix);
    targetFileIdeal = 'jacksheetBR_local.csv';
    targetFiles2cut = {'jacksheetBR_noPulses_local.csv'};
end


%- try to figure out the subject string
subjStr = '';
iNIH = strfind(rawDirPath,'NIH');
if isempty(iNIH), iNIH = strfind(rawDirPath,'TRE'); end
if length(iNIH)>0
    subjStr = rawDirPath(iNIH(end)+[0:5]); %- want NIHXYZ
end
%- now check it.
if length(subjStr)~=6 | ~any(contains(subjStr,{'NIH' 'TRE'}))
    fprintf('\n subject string not correct: %s\n should only affect population of micro chanName and deviceNum');
    keyboard;
end


%- split the full path into parts
[rawRootDir,thisRawDir,~] = fileparts(rawDirPath);


%- get a directory listing for thisRawDir
fprintf('\n PROCESSING %s in makeJacksheetBR', fullfileEEG(rawRootDir,thisRawDir));
if ~exist(fullfileEEG(rawRootDir,thisRawDir),'dir')
    fprintf('\n Uh Oh... directory doesnt exist! returning');
    %keyboard;%
    return;
end


%- Check for and fix an error with the ChanName strings...
fixChanNameStr = 0;  %- seems that all fixes are in place
if fixChanNameStr & exist(fullfileEEG(rawRootDir,thisRawDir,targetFileIdeal),'file')
    jackTable = readtableSafe( fullfileEEG(rawRootDir,thisRawDir,targetFileIdeal) ); %SJ changed from readtable
    fixCount = 0;
    for iC=1:height(jackTable)
        thisChan    = jackTable.ChanName{iC};
        thisChanFix = deblank(thisChan); %(1:find(~isspace(thisChan),1,'last'));
        if length(thisChanFix)<length(thisChan)
            fixCount = fixCount+1;
            jackTable.ChanName{iC} = thisChanFix;
        end
    end
    if fixCount>0
        writetable(jackTable,fullfileEEG(rawRootDir,thisRawDir,targetFileIdeal));
        fprintf('\n updated table with %d fixes to chanName',fixCount);
    else
        fprintf('\n no fixes to chan name required');
    end
    return;
end


%- overwrite existing file or not?
if overwrite==0
    %- if ideal or current file exists skip processing
    if exist(fullfileEEG(rawRootDir,thisRawDir,targetFileIdeal),'file')
        %fprintf('\n already made %s and overwrite is zero, skipping this folder',targetFileIdeal);
        copy2localSubj(rawRootDir,subjStr,thisRawDir,rootEEGdir,targetFileIdeal);
        jackTable = readtableSafe( fullfileEEG(rawRootDir,thisRawDir,targetFileIdeal) ); %SJ changed from readtable
        fprintf(' --> %s already created',targetFileIdeal);
        return;
    elseif exist(fullfileEEG(rawRootDir,thisRawDir,targetFile),'file')
        %fprintf('\n already made %s and overwrite is zero, skipping this folder',targetFile);
        copy2localSubj(rawRootDir,subjStr,thisRawDir,rootEEGdir,targetFile);
        jackTable = readtableSafe( fullfileEEG(rawRootDir,thisRawDir,targetFile) ); %SJ changed from readtable
        fprintf(' --> %s already created',targetFile);
        return;
    end
end
fileList = dir(fullfileEEG(rawRootDir,thisRawDir));
fileList = {fileList.name}';
fileList = fileList(contains(fileList,'.nev') | contains(fileList,'.ns'));
fprintf('---> %d nev or nsx files found', length(fileList));
if length(fileList)==0
    fprintf('\n Uh oh, directory empty. no files found... returning');
    %keyboard;
    return;
end


%- automatically move nev and nsX files into a sub-directory if they can't be read... will save headaches down the line
corruptFileDir = 'cant_read_file';

%- determine whether there are 1 or 2 NSPs and decide which is which
%-   numbering of the two NSPs matters for mapping physical electrode to channel name (particuarly when eCog spans two NSPs)
fileListStem = {};
for iF=1:length(fileList)
    [~,fileListStem{iF},~] = fileparts(fileList{iF});
end
uniqueStems = unique(fileListStem);
numNSPs = length(uniqueStems);
nspList = ones(size(fileList));
if numNSPs<0 | numNSPs>2, fprintf('\n error: more than one unique fileListStem... split the directory'); keyboard; return;
elseif numNSPs==2
    iNSP2 = contains(fileList,{'ieeg2','INST1','utah'});
    if all(iNSP2)
        %- all files from NSP2?  Probably means matlab controlled recording with one utah array on each
        iNSP2 = contains(fileList,{'utah2','utah_u','utahF'});
    end
    if all(iNSP2) | sum(iNSP2)==0
        fprintf('\n');
        disp(fileList);
        fprintf('didnt automatically detect which was the second NSP from list above.  Why?');
        
        keyboard;
        nspList(strcmp(fileListStem,uniqueStems{2}))=2; %- just set it to the second unique string in this case;
    else
        %- a subset of files from NSP2, that is what we expect
        nspList(iNSP2)=2;
    end
end




%- output columns (define here to make sure NEV and NSx creates the same table elemnts)
%strColumns = {'ChanName','NSP','NSPchan','PhysicalChan','SplitterBox','ChanMilliVolt','SampFreq','RawDir','FileName','DurationMin','RangeMilliV','PulseCount'}; % jw's original version
strColumns = {'ChanName','ChanNameNew','NSPsuffix','NSP','NSPchan','PhysicalChan','MicroDevNum','SampFreq','RecFilterHz','ChanMilliVolt','RawDir','FileName','DurationMin','RangeMilliV','PulseCount','SamplesAdded'};

%-
REPORTED_UNUSUAL_FILTER = 0 ; %- just a flag of whether this has been reported... if so, dont need to report separately for each channel

%- loop over files and extract channel information into a table
Tall = [];
for iF=1:length(fileList)
    
    thisFile     = fileList{iF};
    thisFilePath = fullfileEEG(rawRootDir,thisRawDir,thisFile);
    fprintf('\n %s',thisFilePath);
    
    %- calculate nsp Suffix here and add a column for it
    [~,NSPsuffix,~] = fileparts(thisFile);
    if contains(NSPsuffix,'_beh'),    NSPsuffix = NSPsuffix(1:strfind(NSPsuffix,'_beh')-1);    end
    if contains(NSPsuffix,'_premie'), NSPsuffix = NSPsuffix(1:strfind(NSPsuffix,'_premie')-1); end
    NSPsuffixTry = NSPsuffix(find(NSPsuffix=='-' | NSPsuffix=='_',1,'last')+1:end); %- trying to pull out "_INST0", or "_micro", etc
    if ~any(contains(NSPsuffixTry,{'INST','ieeg','utah','micro'}))
        NSPsuffixTry = NSPsuffix(find(NSPsuffix=='-' | NSPsuffix=='_',2,'last')+1:end); %- use the 2nd to last underscore or dash.. this catches "utah_m"
    end
        
    if ~any(contains(NSPsuffixTry,{'INST','ieeg','utah','micro'}))
        fprintf('\n unexpected nsp suffix: "%s"',NSPsuffixTry);
        keyboard;
        thisNSPsuffix = 'error';
    else
        thisNSPsuffix = NSPsuffixTry;
    end
    
    
    thisChanNameNew = '-';
    thisMicroDevNum = nan; %-
    
    
    %- is this NSP1 or NSP2?  numbering matters for mapping physical electrode to channel name
    thisNSP = nspList(iF);
    
    if contains(thisFile,'.nev')
        
        %- NEV as data source has pros and cons.
        %-   pro -- get all the channels at once
        %-   con -- not sure if its possible to tell which were recorded; also sample rate, max analog value, etc are not reported
        try
            if GET_RANGE_and_PULSES
                fprintf(' --> reading all'); tReadStart = tic;
                NEVdata = openNEV(thisFilePath,'nosave','noread','nomat'); % load the data
            else
                fprintf(' --> reading fast'); tReadStart = tic;
                NEVdata = openNEV(thisFilePath,'nosave','noread','nomat','t:0:20'); % load the data
            end
            fprintf(' [%.1fs]',toc(tReadStart));
            fileInfoNEV = NEVdata.MetaTags;
            chanInfoNEV = NEVdata.ElectrodesInfo; %-
        catch
            fprintf('\n error reading NEV file, moving to "cant_read" subdirectory and continuing to NSx');
            if ~exist(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir),'dir')
                mkdir(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir));
            end
            movefile(thisFilePath,fullfileEEG(rawRootDir,thisRawDir,corruptFileDir,thisFile));
            continue;
        end
        
        %- confirm whether 128 or 256 channel nsp  [NEV lists everything]
        if     length(chanInfoNEV)==272, is256nsp = 1;  nspOffset=256;
        elseif length(chanInfoNEV)==144, is256nsp = 0;  nspOffset=128;
        else keyboard; end
        
        
        %- loop over channels and pull some info.
        T = struct2table(chanInfoNEV);
        %- only channels we care about are the serial inputs (not spikes)
        T = T(nspOffset+[1:16],:); %- what is the name of the digital channels?
        
        for iChan=1:height(T)
            thisChan = T.ElectrodeLabel{iChan}';
            thisChan  = deblank(thisChan);
            if sum(isstrprop(thisChan,'alphanum'))~=length(thisChan), fprintf('\n deblank missed something? check the char string %s',thisChan); keyboard; end
            T.ChanName{iChan}  = thisChan;
            T.ChanNameNew{iChan}  = thisChanNameNew;
            T.MicroDevNum(iChan) = thisMicroDevNum;
            T.NSP(iChan)       = thisNSP;
            
            if T.NSP(iChan) ~= 1 && T.NSP(iChan) ~= 2
                keyboard
            end
            
            T.SampFreq(iChan)  = fileInfoNEV.SampleRes;  %- cant get this from NEV.... fileInfoNEV.SampleRes;
            T.FileName{iChan}  = thisFile;
            T.NSPsuffix{iChan} = thisNSPsuffix;
            T.RawDir{iChan}    = thisRawDir;
            
            %- new field JW is adding... recorded filter band.  Look at a few of these and decide whether to cut/keep
            hiFreqHz = double(T.HighFreqCorner(iChan))/1e3; %- data saved as milliHz, divide by 1e3 to make it Hz.  
            loFreqHz = double(T.LowFreqCorner(iChan)) /1e3; %- data saved as milliHz, divide by 1e3 to make it Hz
            if T.HighFreqOrder(iChan)==0 & T.HighFilterType(iChan)==0, hiFreqHz = 0; end %- Confirm order and filter-type are non-zero.  If either of those are zero, then no filter applied.
            if T.LowFreqOrder(iChan)==0  & T.LowFilterType(iChan)==0,  loFreqHz = 0; end %- Confirm order and filter-type are non-zero.  If either of those are zero, then no filter applied.
            if     hiFreqHz==0 & loFreqHz==0,  filterStr = 'none';
            elseif hiFreqHz-floor(hiFreqHz)>0, filterStr = sprintf('%.1f-%.0f', hiFreqHz, loFreqHz);
            else                               filterStr = sprintf('%.0f-%.0f', hiFreqHz, loFreqHz); end
            T.RecFilterHz{iChan} = filterStr;
            
            
            %- convert max analog value to milivolts so BNC and physio values can be combined into a single table
            %thisMaxMiliVolt = double(T.MaxAnalogValue(iChan));
            %if contains(T.AnalogUnits{iChan},'uV'),thisMaxMiliVolt = thisMaxMiliVolt/1e3; end
            T.ChanMilliVolt(iChan) = 5000; %TTL
            
            
            %- NINDS mapping with 512 channel system:
            %      NSP1,  chan 001-128 --> physical channel 1-128    <cervello accessible>
            %      NSP2,  chan 001-128 --> physical channel 129-256  <cervello accessible>
            %      NSP1,  chan 129-256 --> physical channel 257-384  <research only, no cervello access>
            %      NSP2,  chan 129-256 --> physical channel 385-512  <research only, no cervello access>
            %      NSP1,  chan 257-272 --> digital input ain 01-16 (call it 1001 to 1016)
            %      NSP2,  chan 257-272 --> digital input ain 17-32 (call it 1017 to 1032)
            
            %- convert elecID + nspNum to physical electrode number
            elecID = double(T.ElectrodeID(iChan));
            if     thisNSP==1 & elecID<=128,      physElec = elecID;
            elseif thisNSP==1 & elecID>nspOffset, physElec = 1000+elecID-nspOffset;
            elseif thisNSP==1,                    physElec = elecID+128; %- should only happen for 256nsp
            elseif thisNSP==2 & elecID<=128,      physElec = elecID+128;
            elseif thisNSP==2 & elecID>nspOffset, physElec = 1000+elecID-nspOffset+16;
            elseif thisNSP==2,                    physElec = elecID+256; %- should only happen for 256nsp
            end
            
            splitBoxRng  = [1:64:511 1001 1017; 64:64:512 1016 1032];
            splitBoxStr  = {'A(1-64)','B(65-128)','C(129-192)','D(193-256)','E(257-320)','F(321-384)','G(385-448)','H(449-512)','BNC(NSP1)','BNC(NSP2)'};
            thisSplitBox = splitBoxStr{find(physElec>=splitBoxRng(1,:) & physElec<=splitBoxRng(2,:))};
            
            T.NSPchan(iChan)      = T.ElectrodeID(iChan);
            T.PhysicalChan(iChan) = physElec;
            T.SplitterBox{iChan}  = thisSplitBox;
            
            T.DurationMin(iChan)  = fileInfoNEV.DataDurationSec/60;
            
            T.RangeMilliV(iChan)  = -1;
            
            if GET_RANGE_and_PULSES
                %T.PulseCount(iChan)   = length(getBlackrockPulses(NEVdata,9-1+iChan));  %- original way... DC version is basically the same but allows for postProc struct to fix 
                T.PulseCount(iChan)   = length(getBlackRockPulsesDC(NEVdata,9-1+iChan)); %  but allows for postProc struct to fix time gaps in NEV, but that shouldn't change our counts here
                T.SamplesAdded(iChan) = nan; %- for the jacksheet only compute this with NSx files, as that is the number that will feed into the NEV file
            else
                T.PulseCount(iChan)   = nan;
                T.SamplesAdded(iChan) = nan; %- for the jacksheet only compute this with NSx files, as that is the number that will feed into the NEV file
            end
            
        end
        T = T(:,strColumns);
        Tall = vertcat(Tall, T);
        
        
    else
        fprintf(' --> reading header'); tReadStart = tic;
        NSxdata = openNSx(thisFilePath,'noread'); %- 'noread' --> just get the header info
        fprintf(' [%.1fs]',toc(tReadStart));
        if ~isstruct(NSxdata) && (NSxdata==-1 | NSxdata==-99)
            fprintf('\n error reading NSxdata.... moving to "cant_read" subdirectory and continuing');
            keyboard
            if ~exist(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir),'dir')
                mkdir(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir));
            end
            movefile(thisFilePath,fullfileEEG(rawRootDir,thisRawDir,corruptFileDir));
            continue;
        end
        fileInfoNSx = NSxdata.MetaTags;
        chanInfoNSx = NSxdata.ElectrodesInfo; %-
        
        %- confirm whether 128 or 256 nsp... NEV header would tell us, else look to see if BNC input occurs below 256
        %- assume we'll shift any physio channels above 128,
        % if they dont exist, only question is how much to shift BNC channels
        is256nsp = 1;  nspOffset=256;
        if chanInfoNSx(end).ElectrodeID<=144 & chanInfoNSx(end).MaxAnalogValue==5000, is256nsp = 0;  nspOffset=128; end
        
        
        %- now grab extra data if requested
        if GET_RANGE_and_PULSES
            fprintf(' --> reading rng'); tReadStart = tic;
            NSxdata2 = concatOpenNSx(thisFilePath,0,fileInfoNSx.SamplingFreq);  %- use "skipfactor" to grab a sample per second for all channels (used for Range)
            
            if ~isstruct(NSxdata2) && NSxdata2==-99
                fprintf('\n caught openNSx error (infinite loop, with CZs new escape');
                keyboard;
                fprintf('\n error reading NSxdata.... moving to "cant_read" subdirectory and continuing');
                keyboard
                if ~exist(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir),'dir')
                    mkdir(fullfileEEG(rawRootDir,thisRawDir,corruptFileDir));
                end
                movefile(thisFilePath,fullfileEEG(rawRootDir,thisRawDir,corruptFileDir));
                continue;
            end
            
            fprintf(' [%.1fs]',toc(tReadStart));
            if size(NSxdata2.Data,2)>1,
                chanRange = diff(prctile(double(NSxdata2.Data),[5 95],2),1,2); %- use pcntile to protect against outlier points; must convert to double so range is not limited to int16 values
            else
                chanRange = nan(size(chanInfoNSx));
            end
            %--physio:-32764          32764            -8191              8191         'uV
            %--BNC   :-32767          32767            -5000              5000         'mV              '
            bncPulseCnt = nan(size(chanInfoNSx));
            iBNC = find(ismember(double([chanInfoNSx.ElectrodeID]),nspOffset+[1:16]));
            if length(iBNC)>0
                
                fprintf(' --> reading pulses'); tReadStart = tic;
                NSxdata3  = concatOpenNSx(thisFilePath, 0, 1, iBNC); %- need to do a full resolution grab here... slower step
                fprintf(' [%.1fs]',toc(tReadStart));
                bncData = double(NSxdata3.Data);
                
                %- compute the range of the data
                chanRange(iBNC) = diff(prctile(bncData,[0.001 99.999],2),1,2); %- for BNC make the range VERY wide because pulses are rare
                
                %- find the range that doesnt include outliers, maybe 5 to 95 percentile, (looks like need 99th percentile to detect pulses)
                %    if that is >1V, split put the thresh in the middle;  if <1V, then probabably just say all zeros?
                bncData_mV  = bncData*5000/32764;      %- convert to mV so threshold can be set in terms of mV
                Y = prctile(bncData_mV,[5 10 99 100],2);
                rngInner = Y(:,3)-Y(:,2); %- try to avoid outliers
                rngOuter = Y(:,4)-Y(:,1); %- if dont see pulses with 99th percentile, maybe here?
                threshUse = nan(size(iBNC));
                for iiBNC=1:length(iBNC)
                    if rngInner(iiBNC)>1000
                        threshUse(iiBNC)  = Y(iiBNC,2)+rngInner(iiBNC)/2;
                    elseif rngOuter(iiBNC)>1000
                        threshUse(iiBNC)  = Y(iiBNC,1)+rngOuter(iiBNC)/2;
                    else
                        threshUse(iiBNC) = 6000; %- range 5000, so this should result in no pulses;
                    end
                    bncPulseCnt(iBNC(iiBNC)) = length(find( bncData_mV(iiBNC,2:end)>=threshUse(iiBNC) & bncData_mV(iiBNC,1:end-1)<threshUse(iiBNC)));
                end
            else
                fprintf(' --> reading SamplesAdded'); tReadStart = tic;
                NSxdata3 = concatOpenNSx(thisFilePath,0,1,1);    %- grab the first channel in the list at high resolution for SamplesAdded calculation. Unfortunately whole file is still read, but less memory used so can be faster
                fprintf(' [%.1fs]',toc(tReadStart));
            end
            
            %- now look at the high resolution data (single channel or DC channels)... where there any samples added?
            if isfield(NSxdata3.postProc,'samplesAdded')
                samplesAdded = sum(NSxdata3.postProc.samplesAdded);
                if samplesAdded>0, fprintf('<%d added>',samplesAdded); end
            else
                samplesAdded = 0;
            end
            
            
        else
            chanRange    = nan(size(chanInfoNSx));
            bncPulseCnt  = nan(size(chanInfoNSx));
            samplesAdded = nan(size(chanInfoNSx));
        end
        
        
        
        %- loop over the channels and pull info
        numRow = length(chanInfoNSx); %- assume more than
        if numRow==1
            chanInfoNSx(2) = chanInfoNSx; %- hack so table columns are cell arrays even if just one entry
        end
        T = struct2table(chanInfoNSx);
        T = T(ismember(T.ElectrodeID, fileInfoNSx.ChannelID),:); %- only take the recorded channels
        if numRow==1
            T=T(1,:); %- cut back out the fake row used to trick struct2table
        end
        %- take the table craeted from infoNSX, and manipulate/add fields then trim out the unused fields in the end
        for iChan=1:height(T)
            thisChan = T.Label{iChan};
            %thisChan = thisChan(1:find(~isspace(thisChan),1,'last'));
            thisChan = deblank(thisChan); %SJ changed because the 'spaces' are actually cntrl characters, could also use thisChan = thisChan(1:find((~isspace(thisChan)&~isstrprop(thisChan,'cntrl')),1,'last')); 
            T.ChanName{iChan}  = thisChan;
            T.NSP(iChan)       = thisNSP;
            T.ChanNameNew{iChan}  = thisChanNameNew;
            T.MicroDevNum(iChan) = thisMicroDevNum;
            T.SampFreq(iChan)  = fileInfoNSx.SamplingFreq;
            T.FileName{iChan}  = thisFile;
            T.NSPsuffix{iChan} = thisNSPsuffix;
            T.RawDir{iChan}    = thisRawDir;
            
            %- new field JW is adding... recorded filter band.  Look at a few of these and decide whether to cut/keep
            hiFreqHz = double(T.HighFreqCorner(iChan))/1e3; %- data saved as milliHz, divide by 1e3 to make it Hz.  
            loFreqHz = double(T.LowFreqCorner(iChan)) /1e3; %- data saved as milliHz, divide by 1e3 to make it Hz
            if T.HighFreqOrder(iChan)==0 & T.HighFilterType(iChan)==0, hiFreqHz = 0; end %- Confirm order and filter-type are non-zero.  If either of those are zero, then no filter applied.
            if T.LowFreqOrder(iChan)==0  & T.LowFilterType(iChan)==0,  loFreqHz = 0; end %- Confirm order and filter-type are non-zero.  If either of those are zero, then no filter applied.
            if     hiFreqHz==0 & loFreqHz==0,  filterStr = 'none';
            elseif hiFreqHz-floor(hiFreqHz)>0, filterStr = sprintf('%.1f-%.0f', hiFreqHz, loFreqHz);
            else                               filterStr = sprintf('%.0f-%.0f', hiFreqHz, loFreqHz); end
            %- something weird happens with ns6 files sometimes having a filter.  JW thinks its from also recording an ns4 at the same time?
            if ~strcmp(filterStr,'0.3-7500') & contains(thisFilePath,'.ns6')
                if REPORTED_UNUSUAL_FILTER==0
                    fprintf('\n unusual filter string for ns6... should be none or 0.3-7500. recording value from file but heads up');
                    REPORTED_UNUSUAL_FILTER = 1;
                end
                filterStr = '0.3-7500';  %- JW check with blackrock. All our old files look good.  Was a known bug in old central that wrong filter was reported.  So OK to overwrite in jacksheet.
            end
            T.RecFilterHz{iChan} = filterStr;
            
            %- convert max analog value to milivolts so BNC and physio values can be combined into a single table
            thisMaxMiliVolt = double(T.MaxAnalogValue(iChan));
            if contains(T.AnalogUnits{iChan},'uV'),thisMaxMiliVolt = thisMaxMiliVolt/1e3; end
            T.ChanMilliVolt(iChan) = thisMaxMiliVolt;
            
            if thisMaxMiliVolt==8.191, scaleFactor = 8.191/32764; else scaleFactor = 5000/32764; end
            
            
            %- NINDS mapping with 512 channel system:
            %      NSP1,  chan 001-128 --> physical channel 1-128    <cervello accessible>
            %      NSP2,  chan 001-128 --> physical channel 129-256  <cervello accessible>
            %      NSP1,  chan 129-256 --> physical channel 257-384  <research only, no cervello access>
            %      NSP2,  chan 129-256 --> physical channel 385-512  <research only, no cervello access>
            %      NSP1,  chan 257-272 --> digital input ain 01-16 (call it 1001 to 1016)
            %      NSP2,  chan 257-272 --> digital input ain 17-32
            
            %- convert elecID + nspNum to physical electrode number
            elecID = double(T.ElectrodeID(iChan));
            if     thisNSP==1 & elecID<=128,      physElec = elecID;
            elseif thisNSP==1 & elecID>nspOffset, physElec = 1000+elecID-nspOffset;
            elseif thisNSP==1,                    physElec = elecID+128; %- should only happen for 256nsp
            elseif thisNSP==2 & elecID<=128,      physElec = elecID+128;
            elseif thisNSP==2 & elecID>nspOffset, physElec = 1000+elecID-nspOffset+16;
            elseif thisNSP==2,                    physElec = elecID+256; %- should only happen for 256nsp
            end
            
            splitBoxRng = [1:64:511 1001 1017; 64:64:512 1016 1032];
            splitBoxStr = {'A(1-64)','B(65-128)','C(129-192)','D(193-256)','E(257-320)','F(321-384)','G(385-448)','H(449-512)','BNC(NSP1)','BNC(NSP2)'};
            thisSplitBox = splitBoxStr{find(physElec>=splitBoxRng(1,:) & physElec<=splitBoxRng(2,:))};
            
            T.NSPchan(iChan)      = T.ElectrodeID(iChan);
            T.PhysicalChan(iChan) = physElec;
            T.SplitterBox{iChan}  = thisSplitBox;
            
            %- grab duration. can be a cell array, so grab the last one (or sum them?)
            durationMin = sum(fileInfoNSx.DataDurationSec);
            T.DurationMin(iChan)  = durationMin/60;
            T.RangeMilliV(iChan)  = chanRange(iChan)*scaleFactor;
            T.PulseCount(iChan)   = bncPulseCnt(iChan);
            T.SamplesAdded(iChan) = samplesAdded;
        end
        T = T(:,strColumns);
        Tall = vertcat(Tall, T);
    end
    
end
%- last check... did we make a table?
if isempty(Tall)
    fprintf('\n no channel data pulled from any files in the session folder, returning');
    return;
end


%- sort the final table by physical channel numbers, this will put all BNC channels together at the bottom
%jackTable = sortrows(Tall,{'NSP','PhysicalChan'}); %- sort by NSP, and then channel number
jackTable = sortrows(Tall,{'PhysicalChan'});        %- allow BNC channels to all group at the bottom


%- clean up PulsCount.  Set to zero for all channels where not recorded
jackTable.PulseCount(jackTable.PhysicalChan<1000) = -1;  %- this will leave nans in BNCs if there was an issue with detection


%- truncate numerical entries... makes for a cleaner output
trunkFields = {'DurationMin' 'RangeMilliV'};
decPlaces   = 100; %- keep X decimal places
for iF=1:length(trunkFields)
    rawD = jackTable.(trunkFields{iF});
    jackTable.(trunkFields{iF}) = floor(rawD*decPlaces)/decPlaces; 
end



%- but first, are there any micro-recordings?  if so, call microNameConvert to populate MicroDevNum field
i30k = find(jackTable.SampFreq==30000 & jackTable.PhysicalChan<1000 & ~contains(jackTable.FileName,'.nev')); %- MOSTLY micros... in early subjects (<~NIH040) occasionally used on eCog and/or stim DC signals
subjNoMicroIEEG30k = {'NIH038' 'NIH068' 'NIH073' 'TRE006' 'NIH075'};
if length(i30k)>0 & ~contains(subjStr,subjNoMicroIEEG30k)
    
    subjInfo = getMicroSubjInfo_v11(subjStr);
    if isempty(subjInfo)
        fprintf('\n ERROR: jacksheet contains 30kHz chanenls, but micro details not described in getMicroSubjInfo.  Make temporary jacksheet, populate getMicroSubj, and rerun\n');
        disp(jackTable(i30k,:))
        keyboard;
    else
        %- sometimes the microNameConvert can throw and error, maybe make a temp jacksheet first?
        targetFileTemp = fullfileEEG(rawRootDir,thisRawDir,sprintf('%s_tempPreMicroName.csv',targetFile(1:end-4)));
        writetable(jackTable,targetFileTemp);
        
        [~,jackTable] = microNameConvert(subjInfo,jackTable,SKIP_KEYBOARD); %- this will populate ChanNameNew for 30kHz microChans, and MicroDevNum for all channels
        
        delete(targetFileTemp); %- microNameConvert can hang up or crash... best to make a test jacksheet and delete in case that happens
    end
else
    %- no 30k channels, so show that we evaluated that by setting all MicroDev nums to -1
    jackTable.MicroDevNum = -ones(height(jackTable),1);
end


%- prepare to write the table
tableOutFile  = fullfileEEG(rawRootDir,thisRawDir,targetFile);  %- on FRNU56, so have all the files

%- write the table
if overwrite == -1 %SJ
    % Not overwriting!!
else
    writetable(jackTable,tableOutFile);
    %- if file was written successfully, clean up and copy
    if exist(tableOutFile,'file')
        fprintf('\n   successfully wrote %s', tableOutFile);
        
        if strcmp(targetFile,targetFileIdeal)
            for iCheck=1:length(targetFiles2cut)
                thisCut = fullfileEEG(rawRootDir,thisRawDir,targetFiles2cut{iCheck});
                if exist(thisCut,'file')
                    delete(thisCut);
                    fprintf('\n deleted %s',thisCut);
                end
            end
        end
        
        %- copy to local subject folder
        copy2localSubj(rawRootDir,subjStr,thisRawDir,rootEEGdir,targetFile);
    else
        fprintf('\n   Uh oh... table didnt write correctly');
        keyboard
    end
end





%%%%%%%%% LOCAL FUNCTION DEFINITION %%%%%%%%%%%%%%%%
%- if data_raw, then this will be a "complete" file, and it should get copied to the processed subject folder
function copy2localSubj(rawRootDir,subjStr,thisRawDir,rootEEGdir,targetFile)

if contains(rawRootDir,'data_raw')
    
    %- check to see if rootEEGdir/subj contains the same raw dir, if so, put a copy there as well
    tryTarget = {};
    tryTarget{1} = fullfileEEG(rootEEGdir,subjStr,'raw',thisRawDir);
    tryTarget{2} = fullfileEEG(rootEEGdir,subjStr,'raw','STIM_MAP',thisRawDir);
    if strcmp(thisRawDir(end-3:end),'_beh')
        %-sometimes FRNU56 has _beh at the end of the session folder, but FRNU does not
        tryTarget{3} = fullfileEEG(rootEEGdir,subjStr,'raw',thisRawDir(1:end-4));
        tryTarget{4} = fullfileEEG(rootEEGdir,subjStr,'raw','STIM_MAP',thisRawDir(1:end-4));
    end
    
    
    %- test out the different possibilities
    rootEEGtarget = '';
    for iTry=1:length(tryTarget)
        if exist(tryTarget{iTry},'dir')
            rootEEGtarget = fullfileEEG(tryTarget{iTry},targetFile);
        end
    end
    
    
    %- copy to subject folder if subject folder exists
    tableOutFile = fullfileEEG(rawRootDir,thisRawDir,targetFile);
    if length(rootEEGtarget)>0
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile(tableOutFile,rootEEGtarget);
        if SUCCESS
            fprintf('\n--> successful copy to %s  <------\n\n',rootEEGtarget);
        else
            fprintf('\n--> Uh oh, copy to %s failed <----\n\n',rootEEGtarget);
            keyboard
        end
    end
end


