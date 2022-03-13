function jackTable = makeJacksheetCV(rawDirPath,overwrite)
%makeJacksheetCV creates a jacksheetCV_local.csv, which is a 'lite' version
%                of the jacksheetBR but for Cervello data
% This function will return a jacksheetCV that contains the following
% fields: 
%     ChanName: electrode name from .TRC "tagName"
%     ChanNameNew: '-' for now, but can be filled in by getChanNames_CV
%                  (for DC channels) or by user
%     NSP: NSP1 or NSP2
%     PhysicalChan_guess: Physical channel # based on electrode number from
%                         .TRC "chan_record"; Just a guess for now because it doesn't
%                         necessarily map back to the BR channels #s in jacksheetBR
%
% INPUTS:
%     rawDirPath = path to the raw directory ex) '/Volumes/56PROC/eeg_processing/.update/eeg/NIH068/raw/STIM_MAP/181211_1127/
% 
%     overwrite:   (0) If file exists, return that without generating new copy
%                   1  If file exists, generate new copy and return that, overwriting old one
%                  -1  If file exists, generate new copy and return that, WITHOUT overwriting old one
% 
% Called by getChanNames_CV
% 
% Created by Samantha N Jackson 7/17/2020
%
% 7/17/2020 SJ: Created function
%

%CVfile = readCervelloTRC('/Volumes/56PROC/eeg_processing/.update/eeg/NIH068/raw/STIM_MAP/181211_1127/181211_1127_cervello.TRC',[]);

if nargin < 2
    overwrite = 0;
end

if overwrite ~= 1 && overwrite ~= 0 && overwrite ~= -1
    fprintf('%s\n',['ERROR!! overwrite = ' num2str(overwrite) ' (overwrite variable should be equal to 1, 0, or -1).']);
    keyboard
end

%- split the full path into parts
%[rawRootDir,thisRawDir,~] = fileparts(rawDirPath);
targetFile = 'jacksheetCV_local.csv';

if overwrite == 0
    if exist(fullfile(rawDirPath,targetFile),'file')
        jackTable = readtableSafe(fullfile(rawDirPath,targetFile),'file');
        return
    else
        fprintf('\n%s\n',['No jacksheetCV_local.csv ; Creating it now...' ' (' rawDirPath ')'])
    end
end

CVpath = getDirNamesRegexp(rawDirPath,'.*\.TRC$');
if numel(CVpath) > 1 % Char for multiple files will have multiple rows
    fprintf('\n%s\n',['ERROR! Multiple TRC files in this raw folder, ' rawDirPath ' : ']);
    fprintf('\t%s\n',CVpath{:});
    keyboard
else
    CVfile = readCervelloTRC(fullfile(rawDirPath,char(CVpath)),[]);
end

T2 = struct2table(CVfile.electrodes);

% phys_num = regexpCell({CVfile.electrodes.positive_input},'\d{1,3}(?=-)');
% elec_names = {CVfile.electrodes.tagName};
% elec_names_from_pos = regexpCell({CVfile.electrodes.positive_input},'(?<=\d{1,3}-).*');
% if ~isequal(elec_names,elec_names_from_pos)
%     keyboard
% end
% brd_num = {NSxdata.ElectrodesInfo.ElectrodeID};
% brd_names = {NSxdata.ElectrodesInfo.ElectrodeLabel};

% Taken from split_CV:
%- with cervello we routinely record ain17,18,19,20,din2 on the 2nd NSP even if there were no physio channels on second nsp, so check for that
% Find where NSP1 ends -> either 128 or the last ain/din channel before 128
lastChanNSP1 = find(contains(T2.tagName(1:min([size(T2,1) 128])),{'ain' 'din'},'IgnoreCase',true),1,'last');
if size(T2,1)>lastChanNSP1 % Do we have more than 1 NSP?
    physChanNSP2 = sum(~contains(T2.tagName(lastChanNSP1+1:end),{'ain' 'din'},'IgnoreCase',true));
    if lastChanNSP1~=128
        fprintf('\n heads up... two nsps, but cervello only recording %d from the first nsp\n',lastChanNSP1);
    end
else
    physChanNSP2 = 0; %Do we need this for anything...
end

% For NIH063, for some unknown reason CRV can see the second half of NSP1
% and reads it in as NSP2 and Check to see if the session is before October 1
if strcmp('NIH063',char(regexp(rawDirPath,'NIH\d{3}','match'))) && str2double(char(regexp(rawDirPath,'(?<=\d{2})\d{4}(?=_\d{4})','match'))) < 1000
    lastChanNSP1 = 188;
end

jackTable = table('Size',[size(T2,1),4],'VariableTypes',{'string','string','double','double'},'VariableNames',{'ChanName','ChanNameNew','NSP','PhysicalChan_guess'});
jackTable.ChanName = T2.tagName;
jackTable.ChanNameNew = repmat({'-'},numel(jackTable.ChanNameNew),1);


% Find physical channel numbers
% Doesn't seem to be a reliable mapping from CRV to BR, but here are my observations for NIH063:
% Can't just use order/index location because EKG channels are often not placed right after the last
% channel so they will be a higher # in BR but appear right at the end of the list in CV
% NSP1: CV = BR + 31  (start on 32) (even for EKG)
% NSP2: CV = BR + 254 (start on 255)
% *ain/din channnels: ? But for jacksheet purposes, just order them and # 1001+ (no din though)
% *Look at channel naming in makeJacksheerBR ~ line 520+
%


for ii = 1:size(T2,1)
    thisChan = jackTable.ChanName{ii};
    if ii <= lastChanNSP1 % Part of NSP1
        jackTable.NSP(ii) = 1;
        jackTable.PhysicalChan_guess(ii) = T2.chan_record(ii) - 31;
        % Special case for NIH063
        if strcmp('NIH063',char(regexp(rawDirPath,'NIH\d{3}','match'))) && str2double(char(regexp(rawDirPath,'(?<=\d{2})\d{4}(?=_\d{4})','match'))) < 1000
            jackTable.PhysicalChan_guess(ii) = T2.chan_record(ii) + 2; %No clue why this is but this matches jacksheetBR_complete
        end
    else % Part of NSP2
        jackTable.NSP(ii) = 2;
        jackTable.PhysicalChan_guess(ii) = T2.chan_record(ii) - 254 + 128; %To get it to start on 129
    end
    % Check for sync channels and apply different physicalChan
    if contains(thisChan,{'ain','din'}) % should work for ainp too..?
        DC_num = str2double(char(regexp(thisChan,'\d*','match')));
        if contains(thisChan,'ain')
            jackTable.PhysicalChan_guess(ii) = 1000 + DC_num;
            % Populating chanNameNew with DC happens after this function in
            % getChanNames_CV (as with BR)
            
        else %for din
            % Using 2000 here, but no reason why
            jackTable.PhysicalChan_guess(ii) = 2000 + DC_num;
            
        end
    end    
end
% Remove din since this is not part of BR?

% Sort by PhysicalChan_guess (even if we don't know this is correct)
jackTable = sortrows(jackTable,{'PhysicalChan_guess'});

tableOutFile = fullfileEEG(rawDirPath,targetFile);

if overwrite == -1
    % Not overwriting!!
else
    writetable(jackTable,tableOutFile);
    %- if file was written successfully, clean up and copy
    if exist(tableOutFile,'file')
        fprintf('%s\n',['   successfully wrote ' tableOutFile]);
        
%         if strcmp(targetFile,targetFileIdeal)
%             for iCheck=1:length(targetFiles2cut)
%                 thisCut = fullfileEEG(rawRootDir,thisRawDir,targetFiles2cut{iCheck});
%                 if exist(thisCut,'file')
%                     delete(thisCut);
%                     fprintf('\n deleted %s',thisCut);
%                 end
%             end
%         end
        
        %- copy to local subject folder
        %copy2localSubj(rawRootDir,subjStr,thisRawDir,rootEEGdir,targetFile);
    else
        fprintf('%s\n','Uh oh... table didnt write correctly');
        keyboard
    end
end

end

