function getBR_grabFileList(subj,rootEEGdir,rootFRNU56dir, sessList, JUST_CHECK )
% Grab raw files from FRNU56 and move them to a local rootEEG directory for subject procesing (whether eCog in "eeg" folder; or Micro in "56PUB" drive
%
%  input: subj          = subject string, e.g. 'NIH050'
%         rootEEGdir    = directory with subject EEG folder where rawFileList will get dumped (e.g., '/Volumnes/JW24TB/eeg')
%         rootFRNU56dir = e.g., '/Volumes/56C/UTAH_C'     %- frnu56 drive with this subject's data
%         shortTest     = optional, assumed (0);  if 1 then execute code only looking at first 10 nev files to make sure everything works
%
% 11/2018 - JW turned this into a function that will grab any session with non-empty notes
%
%


if nargin<5,
    JUST_CHECK = 0;
end


%- confirm rootFRNU56dir exists
if ~exist(rootFRNU56dir,'dir'),
    fprintf('\n ERROR: specified frnu56 path (%s) not found. \nconfirm its mounted (control-k) and you can access it through finder, then rerun',rootFRNU56dir);
    keyboard;
    return;
end


%- is target a local EEG copy, or is it 56PUB for manual sort processing?
COPY_2_LocalEEG   = 1; %- assume destination is local EEG for eCog processing... will be updated in a moment if not
if contains(rootEEGdir,fullfile(subj,'data_raw')),
    targetDir      = rootEEGdir;
    tweakCSV       = sprintf('rawFileList_%s.xlsx',subj); %- go for the untweaked version in 56A/B/C/D/E
    COPY_2_LocalEEG = 0;
    if ~exist(targetDir,'dir'), mkdir(targetDir); end 
else
    targetDir      = fullfile(rootEEGdir,subj,'raw');
    tweakCSV       = fullfile(rootEEGdir,subj,'raw',sprintf('rawFileList_%s_tweak.xlsx',subj));
end


%- confirm rootEEGdir/subj/raw exists, if not, wont be able to write the result
if ~exist( targetDir,'dir'),
    fprintf('\n ERROR: cant find local EEG dir %s\n break and correct, cause cant write files without this access.',targetDir);
    keyboard;
    return;
end



%---------------------------------------------------------------------------------------------------------%
%- find subject subfolder within rootFRNU56dir
dir56 = dir(rootFRNU56dir);
dir56 = {dir56([dir56.isdir]).name};    %- convert to cell array of dir names
iDir  = find(contains(dir56,subj));     %- check all directory names for subject string. use contains instead of strcmp incase suffix added on frnu56 (e.g., NIH055_iEEGonly)
if length(iDir)~=1,
    fprintf('\n ERROR: expected exactly 1 subj dir in frnu56, found %d\n resolve issue and rerun',length(iDir));
    keyboard;
    return;
end
frnu56subjPath = fullfile(rootFRNU56dir,dir56{iDir});
if ~exist(fullfile(frnu56subjPath,'notes'),'dir'),
    mkdir(fullfile(frnu56subjPath,'notes'));
    fprintf('\n HEADS UP: created %s directory for storing result',fullfile(frnu56subjPath,'notes'));
end


%- now that we know the subject directory on FRNU56, check to see if rawFileList is found (if sessions not specified)
if isempty(sessList),
    if COPY_2_LocalEEG==0,
        tweakCSV = fullfile(frnu56subjPath,'notes',tweakCSV);
    end
    
    if ~exist( tweakCSV,'file' ),
        fprintf('\n tweakCSV: no sessions specified and cant find rawFileList tweak: %s\n', tweakCSV);
        keyboard;
        return
    end
end


%---------------------------------------------------------------------------------------------------------%
%- now find data_rawXXX folder within subject folder.  just take the first one found, could be iEEG, or no suffix, or even utah
dir56 = dir(frnu56subjPath);
dir56 = {dir56([dir56.isdir]).name};  %- convert to cell array of dir names
iDir  = find(contains(dir56,'data_raw'));  %- check all directory names for subject string. use contains instead of strcmp incase suffix added on frnu56 (e.g., NIH055_iEEGonly)
if length(iDir)==0,
    fprintf('\n ERROR: expected 1 or more data directories named %s$$\n in %s \ncheck dir, potentially rename or generate dir, and rerun', frnu56subjPath);
    keyboard;
    return;
end
frnu56rawPath = fullfile(frnu56subjPath,dir56{iDir(1)});


%---------------------------------------------------------------------------------------------------------%
%- now look for subfolders within, will loop over those to get the duration/pulse info
dir56 = dir(frnu56rawPath);
dir56 = {dir56([dir56.isdir]).name};  %- convert to cell array of dir names
dir56 = dir56(3:end);                 %- cut out . and ..
rawDirList = dir56;




%-----
%- if session lis is empty, then grab all sessions with non-empty rows in "notes" column of tweak
sessNotes = {};
if isempty(sessList),
    
    T = readtable( tweakCSV );
    
    if COPY_2_LocalEEG==0,
        iHasNotes = find(T.pulses_DC09 + T.pulses_DC10 + T.pulses_DC11 > 120); %- two minutes of pulses means take it
    else
        iHasNotes = find(~strcmp(T.task,'-')&~strcmp(T.task,''));
    end
    sessList  = T.folderName(iHasNotes);
    sessNotes = T.task(iHasNotes);
end

MAKE_COPY = ~JUST_CHECK;

%%%%-  Now loop over "sessList" and only pull data from folders that are missing in local EEG
grabSuff = {'/*.ns2' ,'/*.ns3', '/*.nev' ,'/jacksheetBR*'};
if COPY_2_LocalEEG==0, grabSuff = {grabSuff{:}, '/*.ns4' ,'/*.ns5' ,'/*.ns6'}; end %- add the high-sample-rate files
    
totNumCopied = 0;  totNumFilesToCopy = 0;  totNumAlreadyThere = 0;
for iS=1:length(sessList),
    thisSess = sessList{iS};
    if length(sessNotes)>0, thisNote = sessNotes{iS}; else thisNote=''; end
    fprintf('\n\n Sesssion %d of %d:  %s <%s>',iS,length(sessList),thisSess,thisNote);
    filesSrc={};
    for iG=1:length(grabSuff),
        fList = dir(fullfile(frnu56rawPath,thisSess,grabSuff{iG}));
        if isempty(fList),
            fList = dir(fullfile(frnu56rawPath,[thisSess '_beh'],grabSuff{iG})); %- JW modifies the folder names as he scans through the directory
        end
        for iF=1:length(fList),
            filesSrc{end+1,1} = fullfile(fList(iF).folder,fList(iF).name);
            filesSrc{end  ,2} = fList(iF).name;
        end
    end
    totNumFilesToCopy = totNumFilesToCopy+size(filesSrc,1);
    fprintf('  found %d files to potentially copy', size(filesSrc,1));
    
    
    %- Now check to see if these exist on local copy, and if not, copy them
    numAlreadyThere = 0;
    numCopiedSess = 0;
    if size(filesSrc,1)>0,
        if contains(thisSess,'_beh'), fprintf('\n removing "_beh" from session name at destination'); end
        thisSessDst = regexprep(thisSess,'_beh','');
        
        if COPY_2_LocalEEG,
            destDir     = fullfile(rootEEGdir,subj,'raw',thisSessDst);
            destDirStim = fullfile(rootEEGdir,subj,'raw/STIM_MAP',thisSessDst);
        else
            destDir     = fullfile(rootEEGdir,thisSess); %- copy sesion name to 56PUB (dont trim "beh" or put in stim/nostim)
            destDirStim = '';
        end
        if ~exist(destDir,'dir') & ~exist(destDirStim,'dir') & MAKE_COPY,
            mkdir(destDir);
        end
        
        %- loop over files and make the copy
        for iF=1:size(filesSrc,1),
            fileDst     = fullfile(destDir,filesSrc{iF,2});
            fileDstStim = fullfile(destDirStim,filesSrc{iF,2});
            if ~exist(fileDst,'file') & ~exist(fileDstStim,'file'), %- if doesnt exist in RAW or RAW/STIM_MAP, then grab it!
                fileSrc = filesSrc{iF,1};
                fprintf('\n about to copy %s --> %s ', fileSrc, fileDst);
                if MAKE_COPY,
                    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fileSrc,fileDst);
                    if ~SUCCESS,
                        fprintf('\n uh oh');
                        keyboard;
                    else
                        totNumCopied = totNumCopied+1;
                        numCopiedSess = numCopiedSess+1;
                    end
                end
            else
                %fprintf('\n file already exists at destination %s, so dont copy ', fileDst);
                totNumAlreadyThere = totNumAlreadyThere+1;
                numAlreadyThere = numAlreadyThere+1;
            end
        end %- for iF
        fprintf('\n >>>  %d checked,  %d already there,  %d copied',  size(filesSrc,1), numAlreadyThere, numCopiedSess);
        
    end %- if length(fileSrc)>0
    
end

fprintf('\n\n  ------  complete:  %d sessions checked,  %d potential files found,  %d already on local EEG, %d of %d copied \n', length(sessList), totNumFilesToCopy,totNumAlreadyThere, totNumCopied,totNumFilesToCopy-totNumAlreadyThere);

