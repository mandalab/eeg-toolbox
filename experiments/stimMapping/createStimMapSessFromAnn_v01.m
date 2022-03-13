function createStimMapSessFromAnn_v01(rootEEGdir, subj)
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
%
%   --Outputs:
%             --plot(s) that shows pulses
%
%   UPDATED 10/2015 so it can handle original EEG-1100 and "new" EEG-1200 extended file formats JW
%


%- first find the STIM_MAP directory
%rawPath       = fullfileEEG(rootEEGdir,subj,'raw');
rawStimFolder = dir(fullfileEEG(rootEEGdir,subj,'raw','*STIM*'));
if isempty(rawStimFolder),
    fprintf('\n WARNING: no raw/*STIM* (ideally raw/STIM_MAP) folder found for %s \n  cant create behavioral/stimMapAnn from annotation file!', subj);
    return;
end
rawStimFolder = fullfileEEG(rootEEGdir,subj,'raw',rawStimFolder(find([rawStimFolder.isdir])).name);
rawPath       = rawStimFolder;
if ~strcmp(rawStimFolder(end-7:end),'STIM_MAP'),
    error(sprintf('\n ERROR: found directory %s, \n         but should be named raw/STIM_MAP; \n       Rename and rerun prepAndAlign\n',rawStimFolder));
end




FIX_FIGURE_NUM       = 101;

FIGURE_FOR_EACH_RAW  = 1; %each raw in its own figure


fprintf('\n\nSearching for .21E and .EEG files in %s:\n', rawPath);
eegRootList = {};
eegRootList = findRaws(rawPath, eegRootList);
eegRootList= unique(eegRootList);




fprintf('\n')
if isempty(eegRootList), fprintf('> No files found!!\n\n');
else
    numRow = min( [5 length(eegRootList)] );
    numCol = 0 + (ceil(length(eegRootList)/numRow));
    
    
    
    %- grab all the data (and the embedded dates for sorting the output plots)
    for iList=1:length(eegRootList)
        
        %-create a clean version of the path for the figure title (and commandline output)
        clnTitle = eegRootList{iList};
        if ~isempty(strfind(clnTitle,'raw')), clnTitle = clnTitle(strfind(clnTitle,'raw'):end); end
        clnTitle(find(clnTitle=='_' | clnTitle=='\'))=' ';
        
        fprintf('... loading %d of %d: %s ...',iList,length(eegRootList),clnTitle)
        
        %disp(eegRootList{iList});
        whichRoot = iList-1;  %- this will become the sessionNumber
        [DC09, DC10, DC11, DC12, sampRate, strTime, numTime, annOut, stimMappingDir] = grabPulsesDC(eegRootList{iList}, whichRoot);
        
        
        
        fprintf(' [%s]\n', strTime);
        
        list_clnTitle{iList} = clnTitle;
        list_DC09{iList}     = DC09;
        list_DC10{iList}     = DC10;
        list_DC11{iList}     = DC11;
        list_DC12{iList}     = DC12;
        list_sampRate{iList} = sampRate;
        list_strTime{iList}  = strTime;
        list_numTime(iList)  = numTime;
        list_annOut{iList}   = annOut; %-structure with annotation information
        list_stimMappingDir{iList} = stimMappingDir;
    end
    
    
    
    %- sort the list by embedded raw date/time
    [sorted, iSort] = sort(list_numTime);
    
    
    
    %- now plot all the loaded data!
    for iList = iSort
        
        clnTitle = list_clnTitle{iList};
        DC09     = list_DC09{iList};
        DC10     = list_DC10{iList};
        DC11     = list_DC11{iList};
        DC12     = list_DC12{iList};
        strTime  = list_strTime{iList};
        numTime  = list_numTime(iList);
        annOut   = list_annOut{iList};
        stimMappingDir = list_stimMappingDir{iList} ;
        
        plotNum  = find(iSort==iList);
        
        
        
        figure(FIX_FIGURE_NUM+plotNum); clf
        set(gcf,'name',clnTitle, 'color', 'w');
        set(gcf,'units','pixels','position',[50    50   1200   1200]);
        
        t = [1:length(DC10)]/sampRate;  %convert to seconds for the individual figures
        strUnits = '(mV)';
        annFont = 15;
        thisAx = [];
        
        % keyboard
        %%- top plot (or only plot) is the difference between
        for iSub=1:4,
            if     iSub==1, DCout = DC09; DCstr = 'DC09';
            elseif iSub==2, DCout = DC10; DCstr = 'DC10';
            elseif iSub==3, DCout = DC11; DCstr = 'DC11';
            elseif iSub==4, DCout = DC12; DCstr = 'DC12';
            end
            
            subplot(4,1,iSub); thisAx(iSub) = gca;
            if ~isempty(DCout), plot(t, DCout, 'r-'); hold on; end
            yAnn = 0;
            if iSub>1,
                for iAnn = 1:length(annOut.samp),
                    hPt = plot(annOut.samp(iAnn)/sampRate,yAnn,'*b','MarkerSize',15); hold on;
                    hTx = text(annOut.samp(iAnn)/sampRate,yAnn,annOut.text{iAnn},'Rotation',90,'FontSize',annFont);
                end
            end
            axis tight;
            set(gca,'ylim',[-500 5500]);
            grid on;
            box off;
            
            set(gca,'fontsize',18);
            ylabel(sprintf('%s %s', DCstr, strUnits));
            xlabel('time (s)', 'fontsize',13);
            if iSub==1,
                title(clnTitle,'fontsize', 15, 'Interpreter','none');
                legend(strTime);
            end
        end
        linkaxes(thisAx,'x');
        
        
        
        figFile = fullfileEEG(stimMappingDir,'DCchannels.png');
        if exist(figFile,'file'),
            figFile = fullfileEEG(stimMappingDir,'DCchannels_x.png');
        end
        
        set(gcf,'Color','w','units','normalized','outerposition',[0 0 1 1]);
        fig2pngSimple(gcf,figFile);
        
        %%
        %set(gcf,'PaperOrientation','portrait','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
        %saveas(gcf,fullfileEEG(stimMappingDir,'DCchannels.png'));
        %export_fig(gcf,figFile,'-png');
        %print(gcf, fullfileEEG(stimMappingDir,'DCchannels.png'), '-dpng', '-r0'); %- saveas looks a little better.. both are shitty
        
    end
    
    
end
fprintf('\n');

end %eegPulseVisualize




%%%-----------------------------------------------------------------------%%%
%%% Extract the EKG channels and Raw file creation date
function [DC09, DC10, DC11, DC12, sampRate, strTime, numTime, annOut, stimMappingDir] = grabPulsesDC(eegRoot, whichRoot)
% Parameter Notes

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
fileStemDate = sprintf('%02d%02d%02d_%02d%02d',T_year,T_month,T_day,T_hour,T_minute);     % new version: file stem of extracted channels (new version.... YYMMDD_HHMM -- JHW 11/2013

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

fclose(fid);

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
iTemp09 = find(~cellfun('isempty',strfind(actualName_ALL,'DC09')));
iTemp10 = find(~cellfun('isempty',strfind(actualName_ALL,'DC10')));
iTemp11 = find(~cellfun('isempty',strfind(actualName_ALL,'DC11')));
iTemp12 = find(~cellfun('isempty',strfind(actualName_ALL,'DC12')));

if length(iTemp09)>1, fprintf('\n WARNING: more than one DC09 found, taking the first'); iTemp09=iTemp09(1); end
if length(iTemp10)>1, fprintf('\n WARNING: more than one DC10 found, taking the first'); iTemp10=iTemp10(1); end
if length(iTemp11)>1, fprintf('\n WARNING: more than one DC11 found, taking the first'); iTemp11=iTemp11(1); end
if length(iTemp12)>1, fprintf('\n WARNING: more than one DC12 found, taking the first'); iTemp12=iTemp12(1); end

DC09 = double(d_int16(iTemp09,:))*GAIN_DC; %convert to volts
DC10 = double(d_int16(iTemp10,:))*GAIN_DC; %convert to volts
DC11 = double(d_int16(iTemp11,:))*GAIN_DC;
DC12 = double(d_int16(iTemp12,:))*GAIN_DC;


sampRate = actSamplerate;

%%- Grab info from the LOG file if present
annOut.samp = [];
annOut.text = {};
%if ~isempty((trigDC10 > 0) | (trigDC11 > 0) | (trigDC12 >0))
if 1,  %- always make the ann output
    annStruct = nk_parseAnnotation(eegRoot);
    annTimes = [];
    annStr   = {};
    annStrPlot = {}; %- just for the graphical outpu
    for iAnn=1:length(annStruct),
        annTimes(iAnn)   = (annStruct(iAnn).timeSec+1)*actSamplerate;  %-add 1 sample to avoid indexing 0 (could cause problem for end of file?)
        annStrPlot{iAnn} = sprintf('  %s',annStruct(iAnn).str);
        annStr{iAnn}     = sprintf('ANNOTATE \t %s',annStruct(iAnn).str);
    end
    annOut.samp = annTimes;
    annOut.text = annStrPlot;
end


%- output a sync for DC10 as well, this will be useful for newer stimMapping sessions
numTrigDC09 = length(find([trigDC09(1) diff(trigDC09)]==1));
numTrigDC10 = length(find([trigDC10(1) diff(trigDC10)]==1));
numTrigDC11 = length(find([trigDC11(1) diff(trigDC11)]==1));
numTrigDC12 = length(find([trigDC12(1) diff(trigDC12)]==1));
numAnn      = length(annTimes);


%- missing trigger... need to manually extract peaks or skip this subject
if numTrigDC10+numTrigDC11+numTrigDC12==0 && numAnn>0,
    fprintf('\n Uh oh... found annotations but no triggered DC outputs.  Probably will have to split and manually extract DC10 times?');
    fprintf('\n    or perhaps NOT extract stimMapping Session?');
    keyboard;
end


%- if stimulation, then DC10, 11, or 12 should have pulses, and possibly annotation file contains channel info
if numTrigDC10>0 ||  numTrigDC11>0 || numTrigDC12>0 || numAnn>0,
    
    
    rawDir         = eegRoot;
    subjDir        = eegRoot(1:strfind(eegRoot,'/raw'));
    stimMappingDir = fullfileEEG(subjDir,'behavioral','stimMapAnn',sprintf('session_%d_(%s)',whichRoot,fileStemDate)); if ~exist(stimMappingDir,'dir'), mkdir(stimMappingDir); end
    
    
    
    %-output annotation file (ONLY TO BEHAVIORAL SESSION FOLDER, NOT TO SPLIT HERE!!!!5
    %... it will also be merged with the updown files below
    if numAnn>0
        chanfile    = fullfile(stimMappingDir, 'annotation.txt');
        fchan       = fopen(chanfile,'w','l');
        fprintf(fchan,'1 \t FILENAME \t %s \n1 \t EEGSTEM \t %s \n', chanfile(min(strfind(chanfile, fileStemDate)):end),fileStemDate); %- make first entry "FILENAME"; trim off the path info
        for iOut = 1:length(annTimes),
            fprintf(fchan,'%d \t %s\n',annTimes(iOut), annStr{iOut});
        end
        fclose(fchan);
        fprintf('\nAnnotation events found (%d total) and extracted to: %s', length(annTimes), chanfile);
    end
    
    
    
    
    
    %- loop through trig arrays and see whether they have any non-zero entries... if so, export a up-down file
    for trigOut = [10 11 12],
        trigStr  = sprintf('DC%02d',trigOut); %-create variable representing trigDC09, 10, 11, etc
        thisTrig = eval(sprintf('trig%s',trigStr)); %-create variable representing trigDC09, 10, 11, etc
        
        %- additional outputs if DC10, 11, or 12 had trigger events: this is used to determine stimulation timining 2/2014
        if sum(thisTrig)>0,
            trigDCout   = double(thisTrig);
            
            %-create list of pulse start (up) and stop (down) times (in units of sample), and a string for each event
            %strUpDown   = {sprintf('%s \t PULSE_HI',trigStr),sprintf('%s \t PULSE_LO',trigStr)};
            strUpDown   = {sprintf('PULSE_HI \t %s',trigStr),sprintf('PULSE_LO \t %s',trigStr)};
            updowns     = [trigDCout(1) diff(trigDCout)];  %diff requires double input for proper functionality
            updownTimes = find(updowns==1  |  updowns==-1);
            updownStr   = {strUpDown{ ((updowns(updownTimes)-1)*-.5)+1 }} ; %-convert -1-->2 and 1->1
            
            %-merge annotation and pulse times.
            timesMerge  = [updownTimes annTimes];
            textMerge   = [updownStr   annStr  ];
            [y,iSort]   = sort(timesMerge);
            timesMerge  = timesMerge(iSort);
            textMerge   = textMerge(iSort);
            
            %-output updown file, with annotation if present
            chanfile    = fullfile(stimMappingDir, sprintf('trig%s.updown.txt', trigStr));
            fchan       = fopen(chanfile,'w','l');
            fprintf(fchan,'1 \t FILENAME \t %s \n1 \t EEGSTEM \t %s \n', chanfile(min(strfind(chanfile, fileStemDate)):end),fileStemDate); %- make first entry "FILENAME"; trim off the path info
            for iOut = 1:length(timesMerge)
                fprintf(fchan,'%d \t %s\n',timesMerge(iOut), textMerge{iOut});
            end
            fprintf(fchan,'%d \t SESS_END \n', size(d_int16,2)); %- make first entry "FILENAME"; trim off the path info
            fclose(fchan);
            fprintf('\nDigital trigger events found on %s... extracted %d pulse up and down times to:\n   %s', trigStr, length(updownTimes), chanfile);
            
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