function numChangesMade = eegModifyAll21Es(subj, rootEEGdir)
%
%  Point of this function is to make it easy to modify all of the 21Es with a changed channel name
%    for instance, if the clinicians mis-labeld EKG1 as EKG, call this function and help push the update
%
%- do this for a whole line instead of just the tag so it is completly unambiguous which line to modify
%            (for instance if there were duplicates it would be ambiguous)
%
% 12/2017 - melkalliny updated so that modification possible for NK + BR

fprintf('\n\n\n called eegModifyAll21Es to identify and implement possible changes to .21E.  Loading 21Es now...\n');

% get channels from raws
rawDir = fullfile(rootEEGdir, subj, 'raw');
raws = findRaws(rawDir, {});
if isempty(raws)
    fprintf('\n No raws found... cant modify 21Es');
    return;
end


%- single text file at root of raw used to note changes made using this function:
rawChangeFile = fullfile(rootEEGdir, subj, 'raw/raw21E_changeHistory.txt');


%- Find the list of .21Es and make sure all corresponding .EEGs are present for subsequent getChanNames call
c21Es           = strcat(raws, '.21E');
cEEGsNK           = strcat(raws, '.EEG');
cEEGsBR           = strcat(raws, '.ns2'); % not flexible 
cEEGsCV           = strcat(raws, '.TRC');
cEEGsfinal = {};
numRaw          = length(c21Es);
codeChansPerRaw = cell(1, numRaw); %- string of codes in .21E indicating each row
for iRaw = 1:numRaw
    d21E = c21Es{iRaw};
    dEEGNK = cEEGsNK{iRaw};
    dEEGBR = cEEGsBR{iRaw};
    dEEGCV = cEEGsCV{iRaw};
    if ~exist(dEEGNK, 'file') && ~exist(dEEGBR, 'file') && ~exist(dEEGCV, 'file');
        assert(exist(dEEGNK, 'file') > 0, 'No EEG file for %s', d21E);
    end
    
    if     contains(raws{iRaw},'ieeg')|contains(raws{iRaw},'INST'), cEEGsfinal{iRaw} = strcat(raws{iRaw},'.ns2');
    elseif contains(raws{iRaw},'TRC'),                             cEEGsfinal{iRaw} = strcat(raws{iRaw},'.TRC');
    else;                                                           cEEGsfinal{iRaw} = strcat(raws{iRaw},'.EEG'); end
    
    
end


numChangesMade = 0;
notDone = 1;
while notDone,
    
    % Get channel names and codes from EEG and 21E file (updated after every .21E change
    for iRaw = 1:numRaw,
        [cRawChans  cRawCodes] = getChanNames(c21Es{iRaw}, cEEGsfinal{iRaw}); % channels as specified in this .21E (can return duplicates if exist, doesn't return empty strings ''
        cCodeChan = {};
        for ii=1:length(cRawCodes),
            cCodeChan{ii} = sprintf('%s=%s',cRawCodes{ii},cRawChans{ii});
        end
        codeChansPerRaw(iRaw) = {cCodeChan};
    end
    
    %- user prompt of which 21E to look at
    fprintf('\n\n Following is a list of .21E files for %s',subj);
    for iRaw = 1:numRaw,
        fprintf('\n %d) %s',iRaw,c21Es{iRaw});
    end
    [iPick, success] = str2num(input(sprintf('\n SELECT a 21E to open and list the contents of. (Enter instance number 1-%d; 0 to quit): ',  length(c21Es)), 's'));
    if isempty(iPick) | iPick<0 | iPick>length(c21Es),
        fprintf('\n invalid selection, idiot. try again');
        continue;
    elseif iPick==0,
        notDone=0;
        return;
    end
    cCodeChan = codeChansPerRaw{iPick};
    
    
    %- list of contents of selected 21E
    fprintf('\n Copy of .21E contents from %s:', c21Es{iPick});
    for iChan=1:length(cCodeChan),
        fprintf('\n%s',cCodeChan{iChan});
    end
    
    %- prompt the change
    fprintf('\nAbove is the list of active (recorded) channels in %s',c21Es{iPick});
    [iChoice, success] = str2num(input(sprintf('\n Enter [1 (default)] to select a line to modify, [2] to see another 21E, [3] to quit \n'),'s'));
    if iChoice==3,
        notDone=0;
        return;
    elseif iChoice==2,
        continue;
    end
    
    %- proceed with modification code
    iOldLine=[];
    while length(iOldLine)~=1,
        oldLine  = input('\nCopy the line to be replaced now, including code and name (e.g., "0062=RAC16", or if multiple, "0062=RAC16,0063=RPI1"): ', 's');        
        lineInputs = strsplit(oldLine,',');
        
        iOldLine = find( strcmp(lineInputs{1,1},cCodeChan) );
        if length(iOldLine)==0, fprintf('\n your selection doesnt match any rows above, try again'); end
    end
    newLine = input('\nEnter the replacement line now,   including code and name (e.g., "0062=RAC16", or if multiple, "0062=RAC16,0063=RPI1"): ', 's');
    newLineInputs = strsplit(newLine,',');
    
    for inputIndex = 1:size(lineInputs,2)
        iOldLine = find( strcmp(lineInputs{1,inputIndex},cCodeChan) );
    %- list of contents of selected 21E
    fprintf('\n\n Intended update to 21E from %s:', c21Es{iPick});
    for iChan=1:length(cCodeChan),
        if iChan==iOldLine,
            numFound = 0;
            for iRaw = 1:numRaw,
                numFound = sum(strcmp(codeChansPerRaw{iRaw},lineInputs{1,inputIndex}))+numFound;
            end
            fprintf('\n%s   < old version "%s" found in %d (of %d) 21Es',newLineInputs{1,inputIndex},lineInputs{1,inputIndex},numFound,numRaw);
        else
            fprintf('\n%s   ',cCodeChan{iChan});
        end
    end
    resp = input(sprintf('\n Make change to all %d affected 21Es (including this one)? [Y]',numFound),'s');
    if isempty(resp), resp='Y'; end
    if resp=='Y' | resp=='y',
        for iRaw = 1:numRaw,
            if sum(strcmp(codeChansPerRaw{iRaw},lineInputs(1,inputIndex))),
                if rewrite21Eline(c21Es{iRaw}, lineInputs{1,inputIndex}, newLineInputs{1,inputIndex},rawChangeFile),
                    fprintf('\n  updated %s',c21Es{iRaw});
                    numChangesMade = numChangesMade+1;
                end
            end
        end
    else
        fprintf('\n  No change made.');
        continue;
    end
    
    end
    
end %while notDone,
end %function



%---------------------------------------------------------------------------------------------------------------------------%
function success = rewrite21Eline(fileName, oldLine, newLine, rawChangeFile)
%- do this for a whole line instead of just the tag so it is completly unambiguous which line to modify (for instance if there were duplicates it would be ambiguous)
success = false;
temp_name = fullfile(tempdir, 'matlab.nk_split.rewriteChannel');

if  contains(fileName,'ieeg') | contains(fileName,'INST') | contains(fileName,'EEG_'), % blackrock or cervello, the same
    temp = strfind(oldLine,'='); % chans are written differently
    oldLine = oldLine(temp+1:end);
    temp = strfind(newLine,'=');
    newLine = newLine(temp+1:end);
    
    temp = strfind(fileName,'/'); % get file that is equiv to 21E
    fileName = fileName(1:temp(end)-1);
    fileName = fullfile(fileName,'ChanNames.txt');
    
    % same as in 21E from here down
    fid_raw = fopen(fileName, 'r');
    fid_tmp = fopen(temp_name, 'w');
    
    line = fgets(fid_raw); %- grab the next line
    anyFound = 0;
    while ischar(line)
        % check line has repeated name with the right channel number and associated renaming
        if strcmp(strtrim(line),oldLine), %- use strtrim to remove white spaces, matching textscan output from getChanNames in oldLine
            anyFound = 1;
            fprintf(fid_tmp, '%s', strrep(line, oldLine, newLine)); %- perserve white spaces and carrage return from .21E
        else
            fprintf(fid_tmp, '%s', line);
        end
        
        line = fgets(fid_raw); %- grab the next line
    end
    success = anyFound;
    
    fclose(fid_tmp);
    fclose(fid_raw);
    if ~success, return; end % didn't find name
    
    backupAffix = '.original';
    [path, name, ext] = fileparts(fileName);
    backupFile = fullfile(path, [name ext backupAffix]);
    if ~exist(backupFile,'file')
        copyfile(fileName, backupFile); % create backup
    end
    copyfile(temp_name, fileName); % create new
    delete(temp_name); % remove temp
    
    %- keep a record of changes in a single file
    fid_rawRecord = fopen(rawChangeFile,'a+');
    fprintf(fid_rawRecord,'\n%s  -->  %s  (%s)',strtrim(oldLine),strtrim(newLine),fileName);
    fclose(fid_rawRecord);
    fprintf('\n NOTE: user-specified change to .21E files:  %s  -->  %s  (%s)',strtrim(oldLine),strtrim(newLine),fileName);
    
    
    
    
else; % correcting channel name extracted from NK
    
    fid_raw = fopen(fileName, 'r');
    fid_tmp = fopen(temp_name, 'w');
    
    line = fgets(fid_raw); %- grab the next line
    anyFound = 0;
    while ischar(line)
        % check line has repeated name with the right channel number and associated renaming
        if strcmp(strtrim(line),oldLine), %- use strtrim to remove white spaces, matching textscan output from getChanNames in oldLine
            anyFound = 1;
            fprintf(fid_tmp, '%s', strrep(line, oldLine, newLine)); %- perserve white spaces and carrage return from .21E
        else
            fprintf(fid_tmp, '%s', line);
        end
        
        line = fgets(fid_raw); %- grab the next line
    end
    success = anyFound;
    
    fclose(fid_tmp);
    fclose(fid_raw);
    if ~success, return; end % didn't find name
    
    backupAffix = '.original';
    [path, name, ext] = fileparts(fileName);
    backupFile = fullfile(path, [name ext backupAffix]);
    if ~exist(backupFile,'file')
        copyfile(fileName, backupFile); % create backup
    end
    copyfile(temp_name, fileName); % create new
    delete(temp_name); % remove temp
    
    %- keep a record of changes in a single file
    fid_rawRecord = fopen(rawChangeFile,'a+');
    fprintf(fid_rawRecord,'\n%s  -->  %s  (%s)',strtrim(oldLine),strtrim(newLine),fileName);
    fclose(fid_rawRecord);
    fprintf('\n NOTE: user-specified change to .21E files:  %s  -->  %s  (%s)',strtrim(oldLine),strtrim(newLine),fileName);
end
end %function success = rewriteChannel