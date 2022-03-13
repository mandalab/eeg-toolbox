function allTags = split_CV(subj, rootEEGdir, EEG_file, BATCH)

% split_CV - Splits a cervello .TRC datafile into separate channels
% the specified directory.
%
% FUNCTION:
%    split_cv(subj,br_dir,EEG_file,...)
%
% INPUT ARGs:
% subj = 'TJ022'
% rootEEGdir = '/Users/damerasr/damerasr/data/eeg/'
% EEG_file = '/Users/damerasr/damerasr/data/eeg/160614_1008/x.EEG'
% BATCH = 0 or 1, if 1, skip keyboards that ask for user input
%
% Output
%	output is saved out to: output_dir = '/Users/damerasr/damerasr/data/eeg/[SUBJ]/eeg.noreref/' set in line 49
%
%
% Edited from previous version so that manual input of the ouput_dir would not be necessary
%
% 12/2013... now uses jacksheetMaster.txt to guide channel number outputs
% 10/2015... now can handel "new" EEG-1200 extended file formats
% 06/2016... convert nihon khoden to cervello split format
% 1/2018 melkalliny - integrated into NihonKohden/Blackrock-compatible prepAndAlign
% 12/2018... ME added and JW tweaked nsp aligment check
% 8/28/2020. SJ: removed pulseFs from call to align_nsps (redundant b/c same as signalFs, also align_nsps only accepts 1 input now)
% 9/2020     SJ: Implementatoin for jacksheetCV_local, copy jacksheetCV_localUsed, remove ChanNames.txt
%                   - if it already exists, check results of current process against it
%                   - Special case for NIH063
%                   - **IMPORTANT: Fix at end where nontrigchans includes DC channels (previously only
%                   looked for ain channels so DC channels were overwritten by versions without
%                   fix_CV_analog applied...
% 9/18/2020  SJ: (from Jiali's account): Adding sort to comparison between names_txt (chanNames.txt) and jschans jacksheetCV_local.csv b/c order doesn't matter
% 9/29/2020  SJ: Centralizing place where channel names are extracted from jacksheetBR/CV- now calling new function getChansFromJack
% 10/2/2020  SJ: fprintf missing a "t"
%

%%old fx inputs
% INPUT ARGs:
% subj = 'TJ022'
% cv_dir = '.../eeg/TJ022/rawCervello/160521_1431'% must contain a single .EEG file
% tagNameOrder={'RFA'; 'RFB'; 'ROFA'; 'ROFB'; 'RAT'; 'RST'; 'RPT'; 'LOF'; 'E';'LF'; 'LAT'; 'LPT'; 'RAH'; 'EKG'; 'RMH'; 'RPH'}

% Edited from previous version so that manual input of the ouput_dir would not be necessary
% need to change namesFromBRCell into namesFromCVCell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VERBOSE = 0 ; %[0,1] = output info about each channel's remapping


tagNameOrder = getTagNames(subj, rootEEGdir);


% load the jacksheetMaster, or create it, or give a warning that it can't be created...
subjDir             = fullfile(rootEEGdir, subj);
rawDir              = fullfile(subjDir, 'raw');
[rawSubDir,~,~]     = fileparts(EEG_file); %- path to the raw file being processed, so FRNU/data/eeg/NIHXXX/raw/121401_1415/
jackMaster_file_new = fullfile(subjDir, 'docs/jacksheetMaster.csv');
sessTitle           = EEG_file(strfind(EEG_file,subj):end);

% create jacksheetMaster if it doesn't exist or hasn't looked at raws
if ~exist(jackMaster_file_new, 'file')
    createMasterJack(subj, rootEEGdir);
end

jacktable = getJackTable(subj, rootEEGdir);
jackMaster_names = jacktable.chanName;
jackMaster_chans = [1:length(jackMaster_names)];


% read the raw_info file to identify channels intentionally EXCLUDED from jacksheetmaster... no reason to give a warning if those are found below
raw_info_file = fullfile(subjDir, 'docs/raw_info.csv');
raw_info      = readtableSafe(raw_info_file);
chanDontSplit = raw_info.chanName(find(raw_info.in_jackSheetMaster==0));
chanDontSplit{end+1} = ''; %- dont report error or try to split channels with no name


% 1. define and create directories
outputDir = fullfile(rootEEGdir,'/',subj,'/eeg.noreref/'); %define noreref directory
if ~exist(outputDir);
    mkdir(outputDir);
    disp('eeg.noreref directory has been created');
end


folderOutDir = strrep(EEG_file,'raw','eeg.noreref'); % for split files
folderOutDir = strrep(folderOutDir,'STIM_MAP',''); % for raws from stim_map
temp = strfind(folderOutDir,'/'); folderOutDir(temp(end):end) = [];
if ~exist(folderOutDir,'file'); mkdir(folderOutDir);  end



% extract data from raw files & resample


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%- open cervello datafile, read all data and channel info
% extract data from raw files & resample
% should be no resampling necessary
chansGrab = []; %- set to string '0' to ONLY get header info, or string with channel list




%% Open files, check validity
EEG = readCervelloTRC(EEG_file,chansGrab);
rawPath = fileparts(EEG_file);

%- with cervello we routinely record ain17,18,19,20,din2 on the 2nd NSP even if there were no physio channels on second nsp, so check for that
lastChanNSP1 = find(contains({EEG.electrodes(1:min([length(EEG.electrodes) 128])).tagName},{'ain' 'din'},'IgnoreCase',true),1,'last');
if length(EEG.electrodes)>lastChanNSP1
    physChanNSP2 = sum(~contains({EEG.electrodes(lastChanNSP1+1:end).tagName},{'ain' 'din'},'IgnoreCase',true));
    if lastChanNSP1~=128
        fprintf('\n heads up... two nsps, but cervello only recording %d from the first nsp\n',lastChanNSP1);
    end
else
    physChanNSP2 = 0;
end

% For NIH063, for some unknown reason CRV can see the second half of NSP1
% and reads it in as NSP2 and Check to see if the session is before October 1
if strcmp('NIH063',char(regexp(rawSubDir,'NIH\d{3}','match'))) && str2double(char(regexp(rawSubDir,'(?<=\d{2})\d{4}(?=_\d{4})','match'))) < 1000
    lastChanNSP1 = 188;
    physChanNSP2 = 0; %Don't want it to align anything since it's only 1 NSP
end

%% old version, assumed first nsp had 128 channels if two nsps... but what if deactivated? that happens in NIH063 (and probably others)
%physChanAbove128 = sum(~contains({EEG.electrodes(129:end).tagName},{'ain' 'din'},'IgnoreCase',true));
%if size(EEG.data,1) <= 128 | physChanAbove128==0, % all physio on a single NSP, no alignment needed

if physChanNSP2 == 0 || strcmp('NIH091',char(regexp(rawSubDir,'NIH\d{3}','match'))) % all physio on a single NSP, no alignment needed
    signal = EEG;
    syncPulse = 1;
    syncStartIndex = 1;
    pulses_to_use  = 1; % use sync pulses from the first nsp
    
else
    
    DC09or12 = 12;
    if (size(find(strcmpi({EEG.electrodes.tagName},'ain4')),2) < 1) || (size(find(strcmpi({EEG.electrodes.tagName},'ain20')),2) < 1)
        fprintf('\nDC12 signal not present on both NSPs- will try aligning with DC09\n')
        DC09or12 = 9;
        %keyboard; %- jw wants to see when this happens
        if (size(find(strcmpi({EEG.electrodes.tagName},'ain1')),2) < 1) && (size(find(strcmpi({EEG.electrodes.tagName},'ain17')),2) < 1)
            DC09or12 = NaN;
            fprintf('\n ERROR: should never get here unless naming was funny (i.e., ainp1 and ainp17) or channels not recorded. check into it (call JW or ME');
            keyboard;
            error('nsp aligment error','split_CV error on line 132, no paired ain channels found when trying to align NSPs');
        end
    end
    
    signalFs = EEG.srate;
    %pulseFs  = signalFs; %SJ: redundant. Why do we even have this??
    signal1  = EEG.data(1:lastChanNSP1,:);  %- use "lastChanNSP1" because if channels are disabled on first NSP will be less than 128 channels
    signal2  = EEG.data(lastChanNSP1+1:end,:);
    
    
    %- define the preferred sync channels
    sync = cell(1,4);
    if DC09or12 == 12
        syncStr = {'ain4', 'ain20','ain1', 'ain17'};
        sync{1,1} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain4')),:);
        sync{1,2} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain20')),:);
        %- grab 3 and 4 for visualization and debugging in case 1 and 2 dont work as expected
        sync{1,3} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain1')),:);
        sync{1,4} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain17')),:);
    elseif DC09or12 == 9
        syncStr = {'ain1', 'ain17'};
        sync{1,3} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain1')),:);
        sync{1,4} = EEG.data(find(strcmpi({EEG.electrodes.tagName},'ain17')),:);
        % there is no DC12, so dont bother
    end
    
    
    %-
    %[signal, pulses_to_use, ~, errorStr, transforms] = align_nsps(sync,signal1,signal2,[pulseFs signalFs],'CV',DC09or12,rawSubDir,sessTitle,0);
    [signal, pulses_to_use, ~, errorStr, transforms] = align_nsps(sync,signal1,signal2,signalFs,'CV',DC09or12,rawSubDir,sessTitle,0); %SJ 8/28/2020: removed pulseFs (redundant b/c same as signalFs, also align_nsps only accepts 1 input now)
    
    EEG.data = signal;
    
    %- used to write alignment summary here... now happens inside of align_nsps
    
end



actualName_ALL = {EEG.electrodes.tagName}'; % dont use this for writing out names
GAIN           = [EEG.electrodes.measurement_unit]';  %- units or gain in volts... 1e-6 means microvolts reported
fileStemDate   = EEG.fileStemDate;

if numel(unique(GAIN)) > 1; fprintf('error - not all chans have same gain');
    keyboard; end
if unique(GAIN) ~= 1*10^-6; fprintf('unexpected- raw data from Cervello not in uV');
    keyboard; end

%- check on sample rate and give warning if not 1024.
samprate_orig=EEG.srate; % find sample rate
if samprate_orig < 1000
    fprintf('\ncervello sample rate less than expected (1000 Hz). use another filesystem if possible\n otherwise, make sure data analyzers know!\n')
    % need to add a way for batch to skip the keyboard here
    if ~BATCH
        keyboard;
    end
end
%samprate_new=1000;
%[fsorig, fsres] = rat(samprate_orig/samprate_new); % resample


%% do not resample anymore - 11/2018 - melkalliny
% initialize matrix of resampled data
% [temp] = resample((EEG.data(1,:)),fsres,fsorig); % fx resamples to 2nd entry.incorrectly labeled as 'orig'!!
% chan=cell(size(temp,2),1);
% EEG_data_resamp_temp = zeros(size(EEG.data,1),size(temp,2));
% for channel=1:size(EEG.data,1)
%     EEG_data_resamp_temp(channel,:) = resample((EEG.data(channel,:)),fsres,fsorig);
%     chan{channel}=EEG_data_resamp_temp(channel,:);
% end


% % instead, write out at native sampleing rate
%
EEG_data_resamp_temp = zeros(size(EEG.data,1),size(EEG.data,2));
chan = cell(size(EEG.data,2),1);

for channel=1:size(EEG.data,1)
    EEG_data_resamp_temp(channel,:) = EEG.data(channel,:);
    chan{channel}=EEG.data(channel,:);
end

EEG.data = EEG_data_resamp_temp;
clear EEG_data_resamp_temp

if ~exist(fullfile(rawSubDir,'jacksheetCV_local.csv'),'file')
    fprintf('\n%s\n','ERROR! jacksheetCV_local.csv does not exist!'); %SJ: fprintf missing a "t"
    keyboard
end



jacksheetCV = readtableSafe(fullfile(rawSubDir,'jacksheetCV_local.csv'));

%SJ replacing the following commented code with this (Centralized pulling of chanNames from
%jacksheet_local):

jschans = getChansFromJack(jacksheetCV);

% jschans = cell(numel(jacksheetCV.ChanName),1);
% for ii = 1:numel(jacksheetCV.ChanName)
%     if strcmp(jacksheetCV.ChanNameNew{ii},'-')
%         jschans{ii} = jacksheetCV.ChanName{ii};
%     else
%         jschans{ii} = jacksheetCV.ChanNameNew{ii};
%     end
% end
% 
% % Now we need to remove the 'exclude' channels (din*, ain17-20) -> from
% % getChanNames_CV (even though these are written to jacksheetCV in the next step after)
% exclude_chans = ~ismember(jschans,regexpCell(jschans,'^(din.*|ain1[7-9]|ain20)$'));
% jschans = jschans(exclude_chans);


%- also make a copy of jacksheetCV_local so eeg.noreref has low-level info about electrode numbers
jackCSV = fullfile(folderOutDir,'jacksheetCV_localUsed.csv'); %- call eeg.noreref copy "localUsed" to deliniate the source file in raw from the copy in eeg.noreref
writetable(jacksheetCV,jackCSV);


channellabels = jschans';

%SJ: got rid of ChanNames.txt in pipeline, but if it still exists, check to
%make sure it's the same as when we use jacksheetCV instead!
% cut raw channels that were not specified in master
% to do so, pull chanNames from .txt made while aligning w/ jacksheet
temp = strfind(EEG_file,'/');
channelTxtPath = strcat(EEG_file(1:temp(end)),'ChanNames.txt');
if exist(channelTxtPath,'file')
    fprintf('\n%s','ChanNames.txt still exists. Checking against jacksheetCV to make sure we get the same result ...');
    fid = fopen(channelTxtPath,'r');
    if fid==-1
        fprintf('Should not happen... how does it  exist then? \n')
        %fprintf('Need to run createMasterJack, so that this .txt file can be created \n')
        keyboard
    end
    names_txt = textscan(fid,'%s');
    fclose(fid);
    names_txt = names_txt{1};
    
    if ~isequal(sort(names_txt),sort(jschans)) %SJ: Adding sort because they are not always in the same order
        fprintf('\n%s\n','ERROR!! Not the same! Get SJ');
        keyboard
    else
        fprintf('%s\n',' we do! Deleting chanNames.txt...');
        delete(channelTxtPath)
    end
    
end

% 
% namesFromBRCell = {};
% %namesFromBRCell{1} = names{1};
% namesFromBRCell{1} = jschans;
%namesFromBRCell = jschans;

% SJ: This wasn't even used for anything...
% actualCode_ALL = {};
% for name = 1:length(names{1})
%     actualCode_ALL{name,1} = sprintf('%0.4d',name);
% end

% make jacksheet.txt using chans from both NSPs
%channellabels = namesFromBRCell{1}';
%channellabels = namesFromBRCell';

%SJ: In the future, I think we can get rid of this and when prep&align
%calls in jacksheet.txt, just use jacksheetCV_local and get rid of
%ain17-20, din, etc like usual

jackFile = fullfile(folderOutDir,'/','jacksheet.txt');
fileOut = fopen(jackFile,'w','l');
if fileOut==-1; error('Jacksheet output directory is not found.'); end
for channel=1:length(channellabels)
    if isempty(strfind(channellabels{channel},'DC')) %fix needed. ain1-4?
        fprintf(fileOut,'%d %s\n',channel,channellabels{channel});
    end
end
for DC=1:4
    %  assuming there were 4 DC chans rec'd
    if DC==1; fprintf(fileOut,'%d DC0%s\n',length(channellabels)-4+DC,num2str(DC+8)); %SJ: used to be DC==1 || DC==2, but that doesn't make sense or else you get DC010
    else; fprintf(fileOut,'%d DC%s\n',length(channellabels)-4+DC,num2str(DC+8)); end
end
fclose(fileOut);

% make params.txt
paramsFile=fullfile(folderOutDir,'/','params.txt');
fileOut = fopen(paramsFile,'w','l');
if fileOut==-1; error('params output directory is not found.'); end
fprintf(fileOut,'samplerate %0.11f\n',samprate_orig);
fprintf(fileOut,'dataformat ''int16''\n');
fprintf(fileOut,'gain %d\n',(GAIN(1)/(1*10^-6))); % default, if violated, will be error above
fclose(fileOut);

% make sourcetype.txt (which system recorded raw)
sourceFile=fullfile(folderOutDir,'/','sourcetype.txt');
fileOut = fopen(sourceFile,'w','l');
fprintf(fileOut,'cervello raw\n');
fclose(fileOut);


% %%- cervello TRC is always a single file... never splits with 2 NSPs
%SJ: No reason to be assigning these again...
%namesFromBRCell = namesFromBRCell{1};
%channellabels   = namesFromBRCell';


% exclude raw chans that were recDontUse in jacksheetMaster
%[~,temp] = intersect(namesFromBRCell,chanDontSplit);
[~,temp] = intersect(channellabels',chanDontSplit);
chan(temp(:)) = []; channellabels(temp(:)) = [];
%namesFromBRCell(temp(:)) = [];
% exclude raw chans that were otherwise not present in jackMaster
%[mismatches, index] = setdiff(namesFromBRCell,jackMaster_names,'stable');
[mismatches, index] = setdiff(channellabels',jackMaster_names,'stable');
if length(mismatches) ~=0
    fprintf(['\n  WARNING: raw channel list contains chans that jacksheetMaster.csv does not.'...
        '\n Channels will not be extracted.  Correct the header of the raws, or add it to element_info.csv to split, '...
        '\n or: suppress in element_info with recNotUsed column, delete jacksheetMaster, and rerun' ...
        '\n if just added a new raw: delete jacksheetMaster and re-run prepAndAlign so the channel can be fixed' ...
        '\n Channels in raws but not jacksheetMaster are:\n']) ;
    disp(mismatches);
    fprintf('\n NOTE: probably shouldnt be possible to get here. so contact wittig if you get here');
    keyboard %SJ: I haven't checked this for new cases anyway
    chan(index(:)) = []; channellabels(index(:)) = [];
    %namesFromBRCell(index(:)) = [];
end


%SJ: the following wasn't used anywhere...
% exclude = {'din' 'ain17' 'ain18' 'ain19' 'ain20'};
% updatedChans = actualName_ALL;
% validChanIndex = 1:length(actualName_ALL);
% for y= 1:length(exclude)
%     temp = strfind(updatedChans(validChanIndex), exclude{y});
%     validChanIndex = find(~not(cellfun('isempty', temp)));
% end

% get analog DC, resample, and save analog files + trigDC09, trigDC12
%if ~contains(EEG_file,'cervello2');
% already have ain data loaded
for DC=1:4
    if (strcmp(subj, 'NIH087') || strcmp(subj, 'NIH091')) && DC > 1
        continue
    end
    
    %- loop over one channel at a time and pull out the up transitions
    addDCs = [0 16];
    for iter=1:2
        % load correct ain data
        analogChanName = sprintf('ain%d',DC+addDCs(iter));
        temp = strfind({EEG.electrodes.tagName},analogChanName);
        
        if DC == 1 && strcmp(subj, 'NIH087') && iter == 1 %Temporary fix for first Micromed Patient NIH087
            temp = strfind({EEG.electrodes.tagName}, 'MKR1+');
        elseif DC > 1 && strcmp(subj, 'NIH087') && iter > 1
            continue
        end
        
        if DC == 1 && strcmp(subj, 'NIH091') && iter == 1 %Temporary fix for first Micromed Patient NIH091
            temp = strfind({EEG.electrodes.tagName}, 'MKR2+');
        elseif DC > 1 && strcmp(subj, 'NIH091') && iter > 1
            continue
        end        

        match = find(not(cellfun('isempty', temp)));
        

        
        
        if ~isempty(match)
            analogDC = EEG.data(match(1),:);
        else
            analogDC = nan(1,EEG.pnts);
        end
        
        
        %% get timestamp of analog DC09 pulses
        
        % this is a janky but so far successful way of dealing with unusual
        % pulse structure - which we initially thought was related to
        % stimulation being applied, but now know can be present even on
        % record-only sessions
        %% JW commenting this out on 12/2018 after fixing the voltage scaling for DC channels in readCervelloTRC.. does that break anything? YES... so JW made "fix_CV_analog" below
        if 0
            %- fill up temp with a count of samples where analogDC is greater than some thresh
            sampAboveThresh = [];  counter = 1;
            for thresh=6000:-100:-500;
                sampAboveThresh(counter, [1:2]) = [thresh length(find(analogDC(1:end-1) > thresh & analogDC(2:end) <= thresh))]; %- now count transitions 
                %sampAboveThresh(counter, [1:2]) = [thresh length(find(analogDC(:) > thresh))]; %- was counting samples above
                counter=counter+1;
            end
            
            %- if temp found ANY samples > ANY thresh
            if sum(sampAboveThresh(:,2))>0
                %- find the biggest jump in super-thresh samples... this should roughly mark the peak of the pulses (i.e., ~5V)
                diffs      = diff(sampAboveThresh(1:end-1,2)); %- end-1 so the thresh can equal the thresh below max
                threshJump = sampAboveThresh(find(diffs==max(diffs),1,'first')+1,1);
                th
                
                %- in some weird cases ain1 is hovering at 1000mV for 1/4 of the session... that screws up this approach
                newHiVal =  mean(analogDC(analogDC > threshJump));
                analogDC(analogDC < 0)          = newHiVal;           %- convert negatives to mean positive.
                analogDC(analogDC > threshJump) = newHiVal;           %- convert positives to mean positive.
            end
            
        else
            analogDC = fix_CV_analog(analogDC,analogChanName );
        end
        
        
        %- this normalizes the traces and takes 25% as threshold crossing.. So if weird Cervello negative blips are in there this will not get the meat of the pulses9
        [syncPulses]=get_triggers(double(analogDC'),samprate_orig);
        
        
        % write out analog files
        if iter==1
            split_filename = sprintf('%s/DC%02d',folderOutDir,DC+8);
            [fchan,msg] = fopen(split_filename,'w','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
            fwrite(fchan,analogDC,'int16');
            fclose(fchan);
        end
        
        
        % write out trig files for NSP1 *AND* NSP2
        % not writing out trigDC12.updown
        if (DC==1 || DC==2 || DC==4)
            if iter==1
                if     DC==1; DCsuffix='trigDC09_nsp1.sync.txt';
                elseif DC==2; DCsuffix='trigDC10_nsp1.syncStim.txt';
                elseif DC==4; DCsuffix='trigDC12_nsp1.syncNSP.txt'; end
                fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
                fprintf(fchan,'%d\n',syncPulses{1}');
                fclose(fchan);
                
            elseif iter==2
                if     DC==1; DCsuffix='trigDC09_nsp2.sync.txt';
                elseif DC==2; DCsuffix='trigDC10_nsp2.syncStim.txt';
                elseif DC==4; DCsuffix='trigDC12_nsp2.syncNSP.txt'; end
                fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
                fprintf(fchan,'%d\n',syncPulses{1}');
                fclose(fchan);
            end
        end
        
        
        % now write out the "winner" of the nsp1 nsp2 alignment
        %   note: sync signal has already been cut according to startIndex - dont cut here
        if (pulses_to_use==2) && (iter==2)
            
            if     DC==1; DCsuffix='trigDC09.sync.txt';
            elseif DC==2; DCsuffix='trigDC10.syncStim.txt';
            elseif DC==4; DCsuffix='trigDC12.syncNSP.txt'; end
            fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
            fprintf(fchan,'%d\n',syncPulses{1}');
            fclose(fchan);
            
        elseif (pulses_to_use==1) && (iter==1)
            
            if     DC==1; DCsuffix='trigDC09.sync.txt';
            elseif DC==2; DCsuffix='trigDC10.syncStim.txt';
            elseif DC==4; DCsuffix='trigDC12.syncNSP.txt'; end
            fchan = fopen(fullfile(folderOutDir,'/',DCsuffix),'w','l');
            fprintf(fchan,'%d\n',syncPulses{1}');
            fclose(fchan);
        end
        
    end % for iter = 1:2
    
end  %- for DC=[1:4]


debuggingPlot=0;
if debuggingPlot==1
    subplot(2,1,1)
    t = [1:length(analogDC)]/1000/60;
    plot(t,analogDC,'b')
    hold on
    tempPulses = syncPulses{1}/1000/60;
    scatter(tempPulses,5000.*ones(1,length(tempPulses)),50,'r')
    xlim([12 13])
    keyboard
    subplot(2,1,2)
    t = [1:length(analogDC)]/1000/60;
    plot(t,analogDC,'b')
end


% split and write all channels into noreref
fprintf('\nwriting files:')
ticker=0;
tick_inc=10;
nontrigchans=find(cellfun('isempty', regexp(channellabels, '(ain)|(DC)'))); % only for non-DC chans
for thisChan=nontrigchans
    if thisChan/length(nontrigchans)*100>=ticker
        fprintf(' %2.0f%%',ticker)
        ticker=ticker+tick_inc;
    end
    split_filename = sprintf('%s/%s',folderOutDir,channellabels{thisChan});
    [fchan,msg] = fopen(split_filename,'w','l');
    assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
    fwrite(fchan,EEG.data(thisChan,:),'int16');
    fclose(fchan);
    if ismac
        try
            % JHW - change files to "executable"... helps for sorting in mac finder
            % MST - avoid permission error (only owner can change permission)
            fileattrib(chanfile, '+x', 'a');
        catch
        end
    end
end

allTags = {};

end
