function [sortFiles_struct sortFiles_cell] = getBR_rawFileList(subj56, rootEEGdir, rootFRNU56dir, shortTest, onlyReturnUtahInfo)
% Create/update the bookkeeping CSV and XLSX files for processing iEEG data using blackrock acquisition system
%
%  input: subj56             = subject string on FRNU56/UTAH_X, e.g. 'NIH050' or 'NIH053_micro"
%         rootEEGdir         = directory with subject EEG folder where rawFileList will get dumped (e.g., '/Volumnes/JW24TB/eeg')
%         rootFRNU56dir      = e.g., '/Volumes/56C/UTAH_C'     %- frnu56 drive with this subject's data
%         shortTest          = optional, assumed (0);  if 1 then execute code only looking at first 10 nev files to make sure everything works
%         onlyReturnUtahInfo = optional, assumed (0);  if 1, then DONT write out the rawFileList, instead format the columns for utahInfo_bookkeeping
%
%  output: rawFileList_NIHXXX.xlsx written to disk on FRNU56 and/or local copies of FRNU/EEG
%          sortFiles_struct and sortFiles_cell:  struct and cell of same info, both packaged for creating CSVs used for manual sorting
%
%
% 2016?   - Written by: Anthony I. Jang
% 05/2018 - Modified by Mostafa El-Kalliny for simpler file-referencing purpose
% 05/2018 - tweaked by JW; add input for rootFRNU56dir
% 11/2018 - another JW overhaul with more columns/info
% 03/2019 - rawFileList will now be saved as xlsx intead of csv.  this protects agains comma's in the notes section messing things up

%-
shortStr = '';
if nargin<5,
    onlyReturnUtahInfo = 0; % Changed from forUtahInfo, SNJ
end
if nargin<4,
    shortTest = 0;
elseif shortTest == 1,
    fprintf('\n\n >> HEADS UP, about to run a short test... rerun with shortTest=0 after things checkout. <<\n');
    pause(2);
    shortStr='_shortTest';
end


%- subj56 is the name of the subject on FRNU56... could have a suffix, like "NIH053_micro".  Make "subjCln" be ONLY the NIH or TRE part of the name (as would be expected on FRNU/eeg)
subjCln = subj56;
if length(subj56)>6,
    if strcmp(subj56(1:3),'NIH')|strcmp(subj56(1:3),'TRE'),
        subjCln = subj56(1:6);
    else
        fprintf('\n unexpected subject name (%s)... check this',subj56);
        keyboard;
    end
end

%- confirm rootEEGdir/subjCln/raw exists, if not, wont be able to write the result
WRITE_TO_LOCAL_EEG = 1;
if ~isempty(rootEEGdir) && ~exist( fullfile(rootEEGdir,subjCln,'raw'),'dir'),
    fprintf('\n WARNING: cant find local EEG dir %s. Break and correct, passing in empty string to bypass making a local copy of jacksheet.', rootEEGdir);
    return;
elseif isempty(rootEEGdir)
    WRITE_TO_LOCAL_EEG = 0;
end

%- confirm rootFRNU56dir exists
if ~exist(rootFRNU56dir,'dir'),
    fprintf('\n ERROR: specified frnu56 path (%s) not found. \nconfirm its mounted (control-k) and you can access it through finder, then rerun',rootFRNU56dir);
    keyboard;
    return;
end


%---------------------------------------------------------------------------------------------------------%
%- find subject subfolder within rootFRNU56dir
dir56 = dir(rootFRNU56dir);
dir56 = {dir56([dir56.isdir]).name};               %- convert to cell array of dir names
iDir  = find(contains(dir56,subj56));  %- check all directory names for subject string. use contains instead of strcmp incase suffix added on frnu56 (e.g., NIH055_iEEGonly)
if length(iDir)~=1,
    fprintf('\n ERROR: expected exactly 1 subj dir in frnu56, found %d\n resolve issue and rerun',length(iDir));
    keyboard;
    return;
end
frnu56subjPath = fullfile(rootFRNU56dir,dir56{iDir});


%- Check for manual sorts subfolder... this is what we expect on 56PUB or local versions of sorting.
notesDir = '_extraction_notes';  %- usually notes, but if a sorting subdirectory use sort_notes
frnu56subjPathManSort = fullfile(frnu56subjPath,'manual_sorts');
if exist(frnu56subjPathManSort,'dir'),
    if exist(fullfile(frnu56subjPath,'data_raw')),
        fprintf('\n ERROR: found a subj/maunal_sorts folder AND a subj/data_raw folder. Should only be one or the other');
        keyboard;
        error('\n fit it');
    else
        frnu56subjPath = frnu56subjPathManSort;
        notesDir = 'sort_notes';
    end
end


if ~exist(fullfile(frnu56subjPath,notesDir),'dir'),
    mkdir(fullfile(frnu56subjPath,notesDir));
    fprintf('\n HEADS UP: created %s directory for storing result',fullfile(frnu56subjPath,notesDir));
end


%---------------------------------------------------------------------------------------------------------%
%- now find data_rawXXX folder within subject folder.  just take the first one found, could be iEEG, or no suffix, or even utah
dir56 = dir(frnu56subjPath);
dir56 = {dir56([dir56.isdir]).name};  %- convert to cell array of dir names
iDir  = find(contains(dir56,'data_raw'));  %- check all directory names for subject string. use contains instead of strcmp incase suffix added on frnu56 (e.g., NIH055_iEEGonly)
if length(iDir)==0,
    fprintf('\n ERROR: expected 1 or more data directories named %s \ncheck dir, potentially rename or generate dir, and rerun', frnu56subjPath);
    keyboard;
    return;
end
frnu56rawPath = fullfile(frnu56subjPath,dir56{iDir(1)});


%---------------------------------------------------------------------------------------------------------%
%- now look for subfolders within, will loop over those to get the duration/pulse info
dir56 = dir(frnu56rawPath);
dir56 = {dir56([dir56.isdir]).name};  %- convert to cell array of dir names
dir56 = dir56(3:end);                 %- cut out . and ..
dir56 = dir56(~strncmp(dir56,'_',1)); %- if directory starts with an underscore dont include it
rawDirList = dir56;

%- initalize the containers
folderNames = {};
brSessData = struct;
brSessData.DC09count  = [];
brSessData.DC10count  = [];
brSessData.DC11count  = [];
brSessData.DC12count  = [];
brSessData.whereDC09  = {}; %- 'ns3' 'nev', 'n/a'
brSessData.brDuration = [];
brSessData.samplesAdded = [];
brSessData.nevCount   = [];
brSessData.ns2Count   = [];
brSessData.ns3Count   = [];
brSessData.ns5Count   = [];
brSessData.ns6Count   = [];
brSessData.fileNameMicro = {};
brSessData.microChanCount = [];

%- loop over every folder and pull the info from  jacksheetBR
num2check = length(rawDirList);
if shortTest, num2check = min([50 num2check]); end
tStart = tic;
for iRaw=1:num2check,
    
    thisRawDir = fullfile(frnu56rawPath,rawDirList{iRaw});
    folderNames{iRaw} = rawDirList{iRaw};
    fprintf('\n %d of %d) reading %s', iRaw,length(rawDirList),thisRawDir);
    
    nevFile    = dir(fullfile(thisRawDir,'*.nev'));
    brSessData.nevCount(iRaw)=length(nevFile);
    brSessData.ns2Count(iRaw)=length(dir(fullfile(thisRawDir,'*.ns2')));
    brSessData.ns3Count(iRaw)=length(dir(fullfile(thisRawDir,'*.ns3')));
    brSessData.ns5Count(iRaw)=length(dir(fullfile(thisRawDir,'*.ns5')));
    brSessData.ns6Count(iRaw)=length(dir(fullfile(thisRawDir,'*.ns6')));
    
    
    %- initilize the info that would come from the jacksheet... then if no jacksheet possible there are defaults in place
    brSessData.brDuration(iRaw)     = 0;
    brSessData.samplesAdded(iRaw)   = nan;
    for whichDC = 9:12
        fieldName = sprintf('DC%02dcount',whichDC);
        brSessData.(fieldName)(iRaw) = 0;
    end
    brSessData.whereDC09{iRaw}      = '???';
    brSessData.microChanCount(iRaw) = 0;
    brSessData.maxMicroRng(iRaw)    = 0;
    brSessData.fileNameMicro{iRaw}  = '???';
    
  
    
    %- initialization complete; now load the jackTable or create one
    jackSheetComplete = fullfile(thisRawDir,'jacksheetBR_complete.csv');
    overwrite=1; GET_RANGE_and_PULSES=1; SKIP_KEYBOARD=0;
    if ~exist(jackSheetComplete,'file'),
        fprintf('\n   --->  heads up, missing jacksheetBR_complete.csv, so attempting to create for this session');
        jackTable = makeJacksheetBR(thisRawDir,rootEEGdir,overwrite,GET_RANGE_and_PULSES,SKIP_KEYBOARD);
    else
        %- jacktable found.  makeJacksheet will push a copy to rootEEGdir and return a table read
        jackTable = makeJacksheetBR(thisRawDir,rootEEGdir,0,1);
        %jackTable = readtable(jackSheetComplete); %- just read table to skip some checks... 
    end
    
    
    %- 
    if isempty(jackTable),
        fprintf('\n  WARNING: empty jacktable for %s', thisRawDir);
        jackTable.DurationMin(1) = nan;
        jackTable.FileName    = {''};
        jackTable.ChanName    = {''};
        continue;
    else
        %- Jacktable found... a few data checks to see if it should be updated:
        %-     (1) raw_folder name should be accurate (sometimes beh suffix added after jacksheet made)
        %-     (2) any ns6 file should show filter string .3 to 7500 Hz (bug in old central version put the wrong filter in the header)
        if ~strcmp(jackTable.RawDir{1},rawDirList{iRaw}),
            fprintf('\n Found inaccurate RawDir listing in jacksheet.  Probably modified folder name with _beh or other.  Recreating jack now');
            keyboard;
            overwrite=1; 
            jackTable = makeJacksheetBR(thisRawDir,[],overwrite,1);
        end
        
        iNS6 = find(contains(jackTable.FileName,'.ns6'));
        if ~all(strcmp(jackTable.RecFilterHz(iNS6),'0.3-7500')|strcmp(jackTable.RecFilterHz(iNS6),'none')),
            fprintf('\n  Found incorrect filter setting from NSx (blackrock bug confirmed)... re-creating jackTable with correct setting');
            %keyboard;
            overwrite=1; 
            jackTable = makeJacksheetBR(thisRawDir,[],overwrite,1);
        end
        
    end
   
    

    
    %- pull the info from the jackTable
    brSessData.brDuration(iRaw)   = jackTable.DurationMin(1);
    brSessData.samplesAdded(iRaw) = max(jackTable.SamplesAdded); %- this should pickup the 30kS channels. point is whether non-zero
    
    whereDC09 = 'n/a';
    iNS6 = find(contains(jackTable.FileName,'.ns6'));
    iNS4 = find(contains(jackTable.FileName,'.ns4'));
    iNS3 = find(contains(jackTable.FileName,'.ns3'));
    iNS2 = find(contains(jackTable.FileName,'.ns2'));
    iNEV = find(contains(jackTable.FileName,'.nev'));
    for whichDC = 9:12
        fieldName    = sprintf('DC%02dcount',whichDC);
        
        %- try 5 possible names for DC channel:  DC09 == ain1 == ainp1 == ain17 == ainp17
        iDC =       find(strcmp(jackTable.ChanName,sprintf('DC%02d',whichDC)));
        iDC = [iDC; find(strcmp(jackTable.ChanName,sprintf('ain%d', whichDC-8)))];
        iDC = [iDC; find(strcmp(jackTable.ChanName,sprintf('ainp%d',whichDC-8)))];
        iDC = [iDC; find(strcmp(jackTable.ChanName,sprintf('ain%d', whichDC-8+16)))];
        iDC = [iDC; find(strcmp(jackTable.ChanName,sprintf('ainp%d',whichDC-8+16)))];
        if whichDC == 10,
            %- for NIH034, ainp2 was called "stimWave"...
            iDC = [iDC; find(strcmp(jackTable.ChanName,'stimWave'))];
        end
        
        thisDCPulses = nan;
        if ~isempty(iDC),
            thisDCPulses = max(jackTable.PulseCount(iDC)); %- iDC could be NEV and NSX, take the max of anything we find
            
            if whichDC==9,
                %- see if NS3 can be used for alignment. thisDCPulses is max pulse count from any channel (ns3, BNC, etc) that records pulses
                pulseNS3 = max(jackTable.PulseCount(intersect(iNS3,iDC)));
                pulseNEV = max(jackTable.PulseCount(intersect(iNEV,iDC)));
                pulseNS4 = max(jackTable.PulseCount(intersect(iNS4,iDC))); %- early subj used NS4 occasionally
                pulseNS2 = max(jackTable.PulseCount(intersect(iNS2,iDC))); %- NIH042 used NS2 once
                pulseNS6 = max(jackTable.PulseCount(intersect(iNS6,iDC))); %-
                if     ~isempty(pulseNS3) & pulseNS3 >= 0.9 * thisDCPulses,
                    whereDC09    = 'ns3';
                    thisDCPulses = pulseNS3;
                elseif ~isempty(pulseNS4) & pulseNS4 >= 0.9 * thisDCPulses,
                    whereDC09    = 'ns4';
                    thisDCPulses = pulseNS4;
                elseif ~isempty(pulseNS2) & pulseNS2 >= 0.9 * thisDCPulses,
                    whereDC09    = 'ns2';
                    thisDCPulses = pulseNS2;
                elseif ~isempty(pulseNEV) & pulseNEV >= 0.9 * thisDCPulses,
                    whereDC09    = 'nev';
                    thisDCPulses = pulseNEV;
                elseif ~isempty(pulseNS6) & pulseNS6 >= 0.9 * thisDCPulses,
                    whereDC09    = 'ns6';
                    thisDCPulses = pulseNS6;
                else
                    fprintf('\n is it possible to get here?  if so, grab the suffix of the winner');
                    keyboard;
                end
            end
        elseif length(jackTable.ChanName)>1,
            fprintf('\n %d channel names, but didnt find DC%02d... that is unusual',length(jackTable.ChanName),whichDC);
            if contains(jackTable.RawDir{1},'NoDC'),
                fprintf(' --> but RawDir named "NoDC", so its acknowledged already');
            else
                keyboard
            end
        end
        brSessData.(fieldName)(iRaw) = thisDCPulses;
    end
    brSessData.whereDC09{iRaw} = whereDC09;
    
    
    %- how many micro channels exist in this session, what was the max range (i.e., were all of them turned off?)
    iMicroChan = find(jackTable.MicroDevNum>0 & ~isnan(jackTable.MicroDevNum) & jackTable.SampFreq==30000);
    brSessData.microChanCount(iRaw) = length(iMicroChan);
    if length(iMicroChan)>0,
        maxRng = max(jackTable.RangeMilliV(iMicroChan));
        if isnan(maxRng) | maxRng<=0, maxRng=-1; end %- condition it here so later we catch that is should be considered zero
        brSessData.maxMicroRng(iRaw) = maxRng;
    else
        brSessData.maxMicroRng(iRaw) = 0;
    end
    
    %- identify all of the nsX files with micro-electrode data... should have channel names utah and micro, but maybe not
    iNs5or6 =           find(contains(jackTable.FileName,'.ns5'));
    iNs5or6 = [iNs5or6; find(contains(jackTable.FileName,'.ns6'))];
    ns5and6 = unique(jackTable.FileName(iNs5or6));
    
    brSessData.fileNameMicro{iRaw} = ns5and6;
    
    
end
fprintf('\n ----  ALL JACKSHEETs READ in %.1f min  ----\n\n',toc(tStart)/60);


%- clean the duration for output readability
decPlaces   = 10; %- keep X decimal places
for iRaw = 1:length(folderNames)
    rawD = brSessData.brDuration(iRaw);
    brSessData.brDuration(iRaw) = floor(rawD*decPlaces)/decPlaces;
end



%---------------------------------------------------------------------------------------------------------%
% Now package and save the csv:   One for utahInfo and one for rawFileList
%---------------------------------------------------------------------------------------------------------%



%%- (1) UtahInfo -%%

%- this column list matches what updateUtahInfoCSV expects
%bk_columns = {'fileName','folderName','ecog_sessName','toExtract','alignEEG','task','pulses_beh','pulses_stim','br_duration','extracted','sorted','notes'}; %- cut br_pulses...
%bk_columns = {'folderName','ecog_sessName','task','toExtract','num_NS5or6','pulses_beh','pulses_stim','br_duration','extracted','sorted','notes'}; %- cut br_pulses...
%bk_columns = {'folderName','task','toExtract','isStim','num_NS5or6','num_microChan','pulses_beh','pulses_stim','br_duration','extracted','sorted','notes'}; %- cut br_pulses...
bk_columns = {'folderName','task','toExtract','isStim','num_NS5or6','num_microChan','pulses_beh','pulses_stim','br_duration','samplesAdded','maxMicroRng','extracted','sorted','notes'}; %- cut br_pulses...

sortFiles_struct = struct;
for ff = bk_columns, sortFiles_struct.(char(ff)) = {}; end

iOut = 1;
for iRaw = 1:length(folderNames)
    fileNames = brSessData.fileNameMicro{iRaw}; %- could be more than one file with micro data if recorded on two NSPs (i.e., two Utahs, or utah + microwire)
    
    %- single row per folder, vs row per ns5/ns6.   Subsequent code is now setup for single row per session, so lets make that happen
    iName = 1;
    %for iName = 1:length(fileNames),
    %sortFiles_struct(iOut).fileName    = fileNames{iName};
    sortFiles_struct(iOut).folderName   = folderNames{iRaw};
    
    sortFiles_struct(iOut).num_NS5or6   = length(fileNames);
    sortFiles_struct(iOut).num_microChan = brSessData.microChanCount(iRaw);
    
    sortFiles_struct(iOut).pulses_beh   = brSessData.DC09count(iRaw);
    sortFiles_struct(iOut).pulses_stim  = brSessData.DC10count(iRaw) + brSessData.DC11count(iRaw); %- micro stim on DC11
    sortFiles_struct(iOut).br_duration  = brSessData.brDuration(iRaw);
    sortFiles_struct(iOut).samplesAdded = brSessData.samplesAdded(iRaw);
    sortFiles_struct(iOut).maxMicroRng  = brSessData.maxMicroRng(iRaw);

    
    toExtractDefault = '-';
    if sortFiles_struct(iOut).pulses_beh>100,   toExtractDefault = 'y';   end  %-  taking an educated guess on what sessions to extract... only applies to freshly created CSV
    if sortFiles_struct(iOut).num_microChan==0, toExtractDefault = 'n/a'; end  %-  cant extract if no micro chan
    if sortFiles_struct(iOut).maxMicroRng<=0,   toExtractDefault = 'n/a'; end  %-  no reason to extract if signals are all zero... utah was unplugged or turned off
    
    
    isStimDefault = '-';
    if sortFiles_struct(iOut).pulses_stim>100, isStimDefault = 'y';   end  %-  taking an educated guess on what sessions are stim ... only applies to freshly created CSV
    
    %- initailze empty fields with string so writetable/readtable properly formats
    sortFiles_struct(iOut).task      = '-';
    sortFiles_struct(iOut).toExtract = toExtractDefault;
    sortFiles_struct(iOut).isStim    = isStimDefault;
    sortFiles_struct(iOut).extracted = '-';
    sortFiles_struct(iOut).sorted    = '-';
    sortFiles_struct(iOut).notes     = '-';
    
    iOut = iOut+1;
    %end
end

sortFiles_cell = cat(1,bk_columns,squeeze(struct2cell(sortFiles_struct))');

%- dont write anything... just return the struct and cell array
if onlyReturnUtahInfo,
    fprintf('\n data packaged for utahInfo');
    return;
end




%---------------------------------------------------------------------------------------------------------%
%- Now create the csv raw_fileList

%%- (2) rawFileList -%%

bk_columns = {'folderName','num_NEV','num_NS2','num_NS3','num_NS5or6','num_microChan','pulses_DC09','pulses_DC10','pulses_DC11','pulses_DC12','whereDC09','dur_minutes','samplesAdded','task','notes'};

bkStruct = struct;
for ff = bk_columns, bkStruct.(char(ff)) = ''; end

for ff = 1:length(folderNames)
    %bkStruct(ff).fileName     = fileNames{ff};
    bkStruct(ff).folderName    = folderNames{ff};
    bkStruct(ff).num_NEV       = brSessData.nevCount(ff);
    bkStruct(ff).num_NS2       = brSessData.ns2Count(ff);
    bkStruct(ff).num_NS3       = brSessData.ns3Count(ff);
    bkStruct(ff).num_NS5or6    = brSessData.ns5Count(ff)+brSessData.ns6Count(ff);
    bkStruct(ff).num_microChan = brSessData.microChanCount(ff);
    bkStruct(ff).pulses_DC09   = brSessData.DC09count(ff);
    bkStruct(ff).pulses_DC10   = brSessData.DC10count(ff);
    bkStruct(ff).pulses_DC11   = brSessData.DC11count(ff);
    bkStruct(ff).pulses_DC12   = brSessData.DC12count(ff);
    bkStruct(ff).whereDC09     = brSessData.whereDC09{ff};
    bkStruct(ff).dur_minutes   = brSessData.brDuration(ff);    %range( br_pulseTimeStamps.DC12{ff}) / br_sampleRate / 60;
    bkStruct(ff).samplesAdded  = brSessData.samplesAdded(ff);  
    
    bkStruct(ff).task          = '-';
    bkStruct(ff).notes         = '-'; %- default placeholder for empty
end

%bkCell = cat(1,bk_columns,squeeze(struct2cell(bkStruct))');
bkTable = struct2table(bkStruct);


% %- clean the duration for output readability [done above to the underlying data so it gets both outputs]
% trunkFields = {'dur_minutes'};
% decPlaces   = 10; %- keep X decimal places
% for iF=1:length(trunkFields),
%     rawD = bkTable.(trunkFields{iF});
%     bkTable.(trunkFields{iF}) = floor(rawD*decPlaces)/decPlaces;
% end


%disp(bkCell);

%- output directories for rawFileList... make a copy in frnu56 and rootEEG dir
fileDir_rawFileListFRNU56 = fullfile(frnu56subjPath, notesDir,sprintf('rawFileList_%s%s.xlsx',subjCln,shortStr));
fileDir_rawFileListEEGdir = fullfile(rootEEGdir,subjCln,'raw',sprintf('rawFileList_%s%s.xlsx',subjCln,shortStr));

%- grab the task and notes from the old rawFileList if present, then rename that old file
if exist(fileDir_rawFileListFRNU56,'file'),
    %- load the existing table and copy its notes and task info to the new table
    oldTable = readtable( fileDir_rawFileListFRNU56 );
    for iTold=1:height(oldTable),
        iTnew = find(strcmp(bkTable.folderName,oldTable.folderName{iTold}));
        if length(iTnew)==1,
            bkTable.task{iTnew}  = oldTable.task{iTold};
            bkTable.notes{iTnew} = oldTable.notes{iTold};
            if ~strcmp(oldTable.task{iTold},'-') |  ~strcmp(oldTable.notes{iTold},'-'),
                fprintf('\n transfered info from old to new table, folder=%s: task=%s, notes=%s',oldTable.folderName{iTold},oldTable.task{iTold},oldTable.notes{iTold});
            end
        end
    end
    
    %- check to see if new copy is identicle to old copy, if not, save a copy of the old table
    if height(oldTable)>0 && (numel(oldTable)~=numel(bkTable) || height(intersect(bkTable,oldTable))~=height(oldTable)),
        oldFile = sprintf('%s[replaced %s].xlsx',fileDir_rawFileListFRNU56(1:end-4),datestr(now,'yymmdd_HHMM'));   %- add date string so multiple "old" files can exist
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(fileDir_rawFileListFRNU56,oldFile,'f');
        if ~SUCCESS,
            writetable(oldTable,oldFile);
            if   exist(oldFile,'file'),  SUCCESS=1;
            else fprintf('\n ERROR: cant save copy of old table... tried two ways'); keyboard; end
        end
    end
end



%%- Make a copy of the new table.  On FRNU56 and in local EEG
%cell2csv(bkCell,fileDir_rawFileListFRNU56);
writetable(bkTable,fileDir_rawFileListFRNU56);

if WRITE_TO_LOCAL_EEG, writetable(bkTable,fileDir_rawFileListEEGdir); end
if ~exist(fileDir_rawFileListFRNU56,'file') | (WRITE_TO_LOCAL_EEG & ~exist(fileDir_rawFileListEEGdir,'file')),
    fprintf('\n ERROR: output files not saved. take a closer look');
    keyboard;
else
    fprintf('\n FILES SAVED:\n    %s\n    %s\n',fileDir_rawFileListFRNU56,fileDir_rawFileListEEGdir);
end

if shortTest,
    fprintf('\n\n  *** SHORT-TEST COMPLETE... RERUN with shortTest=0 now to create proper files! ***\n');
else
    %- not a short test... delete short test result if present
    shortStr = '_shortTest';
    fileDir_rawFileListFRNU56 = fullfile(frnu56subjPath,'notes',sprintf('rawFileList_%s%s.xlsx',subjCln,shortStr));
    if exist(fileDir_rawFileListFRNU56,'file') delete(fileDir_rawFileListFRNU56); end
end



