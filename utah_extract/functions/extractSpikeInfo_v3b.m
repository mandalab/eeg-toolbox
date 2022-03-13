function extractSpikeInfo(subj,rootMicroDir,sessFolders,ONLY_UPDATE_TXT)
% AIJ: This function looks through the trim2eeg folder and extracts spikes from
% all the sessions that have been trimmed.
%
% 1. Creates spikeInfo.mat which contains info on spike times, spike waveforms, quality metrics, etc.
% 2. Saves a PNG figure with each channel's spike quality in a folder named 'sortFigs'
% 3. Creates sorts_INCOMPLETE.txt or sorts_(complete).txt
% 
% Input:
% subj = 'NIH0XX'
% rootMicroPath = directory of the utah folder in external drive (e.g. '/Volumes/56A/UTAH_A')
% sessFolders = [optional parameter] cell array of strings of sessions to process
% ONLY_UPDATE_TXT = 1 (true) or 0 (false)
%
%
%  3/2019:  JW overhauls... now this operates on data_processed/sorted/subDir/sort_txt
%  9/2019:  SJ edits: correct issues with summary .txt file (variable
%          sortNoteSummary) where 0 sorts out of 0 grades was flagged as
%          incomplete and add date, create list of channels that should probably be
%          sorted, more detailed message about channels without grades (blanks)
%  10/2019: SJ edits: fix error with older sorting with only 1 grade
%           column, or cases where sorter only grades 1 column, fixed issue where
%           code would not recognize [NoUnits] spike info, added option to only
%           update sort summary text (ONLY_UPDATE_TXT)
%  1/2020:  SJ edits: adding in case for when someone updates a
%           sortNoteTable for a session that has already been sorted, but
%           has no new sorts -> re-extract the spike info
%  2/13/20: SJ edits: adding edit to above case, compare spikeInfo and sortNoteTable 
%           modification dates instead (previous solution was flagging every time)
%  3-4/2020:SJ edits: Adding in noise units
%  8/17/2020SJ edits: Adding case for dummySort column is all NaN; also add "any"
%  11/30/2020 SJ: Fixed case where dummySort, chanNoise, and globalNoise columns are read with str2double when there are no strings in the column
%
% TO-DO:
%   Update spikeInfo to include any new changes to sortNoteTable, sort
%   summary text, instructions (?), etc
%
%
%
subDir    = fullfile(rootMicroDir,subj);
sortedDir = fullfile(subDir,'sorts_manual');  %- sorts will be in here. 
procDir   = fullfile(subDir,'processed_SPK'); %- binarys used to grab spike waveforms will be in here
samplerate = 30000;

SAVE_ALL_SNIPPET_WAVES = 0; %- optional... include a field with all snippet waveforms. mu and sd is always saved, so this would be on top of that
%- current version of saving snippet waves draws on the plexon output time series. IF this is set to 1 ever, should probably use the longer timeseries from .bin


%- if user doesn't pass in session folder list, then check all of the folders
if nargin<3 || length(sessFolders)==0
    sessFolders = dir(sortedDir);
    sessFolders = {sessFolders([sessFolders.isdir]).name};
    sessFolders = sessFolders(~strcmp(sessFolders,'.') & ~strcmp(sessFolders,'..'));
    fprintf('\n Found %d data_processed sub-folders to check', length(sessFolders));
end
if nargin<4
    ONLY_UPDATE_TXT = 0;
end




tStartAll = tic;

%- loop over processed session folders
for s = 1:length(sessFolders)
    tStartSess = tic;
    
    %- list of already-created sort sub-directories.  these could have been modifed with sorter initials, so
    sessID      = sessFolders{s};
    %sortDirRoot = fullfile(sortedDir,sessID,'sorted'); %- old path... processed/session/sorted/ref_sortedBy??. 
    sortDirRoot = fullfile(sortedDir,sessID);           %- new path... sorts_manual/session/ref_sortedBy??. 
    if exist(sortDirRoot,'dir')
        sortDirList = dir(fullfile(sortDirRoot));
        sortDirList = {sortDirList([sortDirList.isdir]).name};
        sortDirList = sortDirList(~strcmp(sortDirList,'.') & ~strcmp(sortDirList,'..'));
    else
        sortDirList = {};
    end
    
    
    if isempty(sortDirList),
        fprintf('\n No sorted data found: %s\n possibly call grabSortedData here...',sessID);
        %keyboard;
        continue;
    end
    
    
    %- loop over spike directories and look for .txt files, then rename (_grabbed) and copy to sorted/session
    for iDir=1:length(sortDirList)
        tStartSortTxt = tic;
      
        sortDir      = sortDirList{iDir}; %- e.g., spikes_noreref
        sortedByPath = fullfile(sortDirRoot,sortDir);
        sortTxtPath  = fullfile(sortedByPath,'sort_txt');
        sortInfoFile = fullfile(sortedByPath,'_extractionInfo','_extractMicro_readme.txt');
        
        fprintf('\n\n>>>> CHECKING %s/%s/%s to see if spikeInfo is created and contains all available sort.txt <<<<',subj,sessID,sortDir);
        
        
        %- look for text files in the spikes folder
        sortTxtFiles = dir(fullfile(sortTxtPath,'*.txt'));
        sortTxtFiles = {sortTxtFiles.name};
        if length(sortTxtFiles)==0
            fprintf('\n    no sort files in %s, \n skipping that session', sortedByPath);
            continue;
        else
            %- add second column to "sortTxtFiles" that has the "grabbed XYZ" removed.. the clean version will be used to match up with physio channels
            sortTxtFilesCln = {}; %- with "grabbed" removed from the string
            for iF=1:length(sortTxtFiles)
                iGrab = strfind(sortTxtFiles{iF},'_(grabbed');
                if ~isempty(iGrab), sortTxtClean = sprintf('%s.txt',sortTxtFiles{iF}(1:iGrab-1));
                else                sortTxtClean = sortTxtFiles{iF}; end
                sortTxtFilesCln{iF} = sortTxtClean;
            end
        end
        
        
        %- look at the _extractMicro_readme.txt that is linked to these sorts... down below we will confirm this is the same as the extractMicro_readme info linked to the split physio files
        %    hopefully that sanity check alawys comes up negative... but JW saved this file in both places just in case.  If the files are different, next step would be to confirm that the 
        %    resultant splits were identicle, in which case the sort outputs can be used with the physio splits even though they werent created at the same time
        if ~exist(sortInfoFile,'file')
            fprintf('\n expecting "extractMicro_readme.txt" but didnt find it.  This is used to confirm sorts were made from same binaries as are used for pulling spike waveforms');
            sort_extractMicro_str = 'no extractionMicro_readme.txt found with sorted spikes (legacy?)';
            keyboard;
        else
            %- microExtract info file present... read it in.  Will use to compare with physio extraction file down below
            fidInfo = fopen( sortInfoFile, 'r');
            sort_extractMicro_str = '';
            while 1
                tline = fgetl(fidInfo);
                if      ~ischar(tline), break;
                elseif ~isempty(tline), sort_extractMicro_str = sprintf('%s%s\n',sort_extractMicro_str,tline); end
            end
            fclose(fidInfo);
            if isempty(sort_extractMicro_str)
                fprintf('\n sort extraction string is empty.. that shouldnt happen. Figure it out');
                keyboard;
            else
                sort_extractMicro_str = sort_extractMicro_str(1:end-1); %- remove last carriage return
            end
        end
        
        
        %- spikeInfo file
        %spkInfoPath   = fullfile(sortedByPath,sprintf('%s_%s_spikeInfo.mat',subj,sessID)); %- missing NSPsuffix
        spkInfoPath   = fullfile(sortedByPath,sprintf('%s_spikeInfo.mat',sessID)); %- remove the subject from the front
        spkInfoPath_NU= fullfile(sortedByPath,sprintf('%s_spikeInfo[NoUnits].mat',sessID));    % Also need case with [NoUnits] - SJ
        figSaveFolder = fullfile(sortedByPath,'sortFigs');
        
        if exist(spkInfoPath_NU,'file') % SJ
            if exist(spkInfoPath,'file') % SJ
                fprintf('%s\n','Error! We have a [NoUnits] spikeInfo and a normal one! Why is this?');
                keyboard
            else
                spkInfoPath = spkInfoPath_NU;
            end
        end
        
        % Add in flag where you can update (ie replace) the summary .txt
        % file even if it already exists (ie for then the below flag will
        % skip it normally)
        
        %- if spike info exists, confirm it was made using all the sorts in the sort folder
        MAKE_NEW_SPIKEINFO = 1;
        if exist(spkInfoPath,'file')
            temp      = load(spkInfoPath);
            spikeInfo = temp.spikeInfo;
            if isfield(spikeInfo,'sort_filenames') %SJ: What is this field doesn't exist? No else condition!!
                sortsSaved = spikeInfo.sort_filenames;
                if length(sortsSaved)==length(sortTxtFiles)
                    LIA = ismember(sortTxtFiles,sortsSaved); %- for every sort file, is there a corresponding save file?
                    if sum(LIA)==length(LIA)
                        % Compare current sort note table to the one in spikeInfo to see if the sorter has updated any sort grades
                        % SJ
                        spkinfosortnotetable = spikeInfo.sortNoteTable{1,2};
                        oldsortnotetable = removevars(spkinfosortnotetable,{'maxGrade','hasSortTxt'});
                        newsortnotetable_path = [sortedByPath filesep char(getDirNamesRegexp(sortedByPath,'sortNotes\(\d{6}_\d{4}.*\)_sortedBy[^\?].*'))];
                        newsortnotetable = readtable(newsortnotetable_path);
                        if ONLY_UPDATE_TXT == 1 %SJ
                            [~, ~, ~] = updateSortSummaryTxt(subj,sortedByPath,sessID,sortTxtFilesCln);
                            MAKE_NEW_SPIKEINFO = 0; 
                            continue
                        end
                        if strcmp(newerFile(spkInfoPath,newsortnotetable_path),newsortnotetable_path) %If sort table has been modified after spikeInfo
                            %~isequal(oldsortnotetable,newsortnotetable)
                            % If the sortnoteTable in spikeInfo is not the
                            % same as the one in the folder, stop the
                            % program and ask if spikeInfo should be
                            % updated with the new one - I.e., did someone
                            % change this purposely or did they make a
                            % mistake? - You will want to do the extraction
                            % now!
                            fprintf('\n%s\n%s\n','Spike info contains all sorts, BUT it does not have the most current sortNoteTable.','Most likely, the sortNoteTable has since been updated with different grades. Please confirm this is the case and continue to re-extract spike info.');
                            keyboard
                            newsi = input('Would you like to create a new spikeInfo (n) or continue as-is (c)?','s');
                            if strcmp(newsi,'n') || strcmp(newsi,'N')
                                MAKE_NEW_SPIKEINFO = 1;
                            elseif strcmp(newsi,'c') || strcmp(newsi,'C')
                                MAKE_NEW_SPIKEINFO = 0;
                                continue
                            else
                                fprintf('%s\n','Error!!!!! You did not type one of the response options.');
                                keyboard
                            end
                            
                        else
                            fprintf('\n    spike info contains all sorts & spikeInfo has most current sortNoteTable, no action taken: %s', spkInfoPath);
                            MAKE_NEW_SPIKEINFO = 0;  
                            continue

                        end
                    end
                elseif length(sortsSaved)>length(sortTxtFiles)
                    fprintf('\n    kinda a weird place to stop. were some sort files removed?');
                    keyboard;
                end
            end
            
            if MAKE_NEW_SPIKEINFO
                oldSpikeInfo = sprintf('%s[replaced %s].mat',spkInfoPath(1:end-4),datestr(now,'yymmdd_HHMM')); 

                %oldSpikeInfo = sprintf('%s.old', spkInfoPath);
                %while exist(oldSpikeInfo,'file'),
                %    oldSpikeInfo = sprintf('%s.old', oldSpikeInfo);
                %end
                [SUCCESS,MESSAGE,MESSAGEID] = movefile(spkInfoPath,oldSpikeInfo,'f');
                if SUCCESS==0
                    fprintf('\n Uh oh, problem saving the old spike info: %s \n permissions issue?  must fix before moving forward (can manually rename then rerun)',oldSpikeInfo);
                    keyboard
                    continue
                end
                fprintf('\n    existing spike info did NOT contain all sorts and/or new sortNoteTable, creating new spike info, remaing old one %s', oldSpikeInfo);
            end
        else
            fprintf('\n    sorted channels found, and no existing spike info.  Making spikeInfo now');
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %- if code above says there is new info, then make the spike info
        if MAKE_NEW_SPIKEINFO
            %- now that we know we are making a new spike_info, we need to identify the matching sorts_XXX folder to acces the physio binary files
            %physioDirList = dir(fullfile(procDir,sessID,'spikes_*')); %- old path name was spikes_noreref, spikes_reref, etc
            physioDirList = dir(fullfile(procDir,sessID,'sort_*'));
            physioDirList = {physioDirList([physioDirList.isdir]).name};
            if isempty(physioDirList)
                fprintf('\n No spike directory found, so cant process spikeInfo: %s',sessID);
                keyboard
                continue
            else
                %- found the directories, now identify the match
                physioDirListCln = {}; iMatch = [];
                sortDirCln = sortDir(1:find(sortDir=='_',1,'last')-1); if length(sortDirCln)==0, sortDirCln=sortDir; end
                for iD=1:length(physioDirList)

                    physioDirListCln{iD} = physioDirList{iD}(6:end); %- cut off "sorts_"
                    if strcmp(physioDirListCln{iD},sortDirCln)
                        iMatch = [iMatch iD];
                    end
                end
                if length(iMatch)==0
                    fprintf('\n ERROR:  cant find spike dir that matches the sort dir %s', sortDir);
                    keyboard
                    continue
                elseif length(iMatch)>1
                    fprintf('\n ERROR: more than one spike dir matches the sort dir %s', sortDir);
                    disp(physioDirList(iMatch));
                    keyboard
                    continue
                else
                    %- found exactly one match, as expected
                    SpikeBinDir     = physioDirList{iMatch};
                    if contains(SpikeBinDir,'legacy'); SpikeBinDir = sprintf('%s/%s',SpikeBinDir,'v_use/'); end
                    SpikeBinDirPath = fullfile(procDir,sessID,SpikeBinDir); %- where the bin files exist that will be used for channel unit quality
                end
            end


            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            % CONFIRM THAT THE PHYSIO TIME SERIES USED FOR SORTING MATCHES THE TIME SERIES USED FOR PULLING OUT TIMING and VOLTAGE WAVEFORMS
            %    look at the _extractMicro_readme.txt that is linked to these sorts... down below we will confirm this is the same as the extractMicro_readme info linked to the split physio files
            %    hopefully that sanity check alawys comes up negative... but JW saved this file in both places just in case.  If the files are different, next step would be to confirm that the
            %    resultant splits were identicle, in which case the sort outputs can be used with the physio splits even though they werent created at the same time
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            physioInfoFile = fullfile(procDir,sessID,SpikeBinDir,'_extractionInfo','_extractMicro_readme.txt');
            if ~exist(physioInfoFile,'file')
                fprintf('\n expecting "extractMicro_readme.txt" but didnt find it.  This is used to confirm sorts were made from same binaries as are used for pulling spike waveforms');
                physio_extractMicro_str = 'no extractionMicro_readme.txt found with split physio files (seems impossible)';
                keyboard;
            else
                %- microExtract info file present... read it in.  Will use to compare with physio extraction file down below
                fidInfo = fopen( physioInfoFile, 'r');
                physio_extractMicro_str = '';
                while 1
                    tline = fgetl(fidInfo);
                    if      ~ischar(tline), break;
                    elseif ~isempty(tline), physio_extractMicro_str = sprintf('%s%s\n',physio_extractMicro_str,tline); end
                end
                fclose(fidInfo);
                if isempty(physio_extractMicro_str)
                    fprintf('\n sort extraction string is empty.. that shouldnt happen. Figure it out');
                    keyboard
                else
                    physio_extractMicro_str = physio_extractMicro_str(1:end-1); %- remove last carriage return
                end
            end

            %- here is the check.. does sort and physio extraction match.
            if ~strcmp(physio_extractMicro_str,sort_extractMicro_str)
                fprintf('\n  WARNING: process_SPK physio binaries were created at a different time from those used to spike sort.  \n    Should add some additional checks to see split is identicle, else correct for the difference');
                keyboard %SJ: This was commented out before, but do we get here..?
            end    


            %- if the unfilt direcotry exists, then grab those waveforms when making the spike info
            physioDirUnfilt = 'sort_unfilt';
            if exist(fullfile(procDir,sessID,physioDirUnfilt))
                fprintf('\n    found spikes_unfilt folder, so spikeInfo will include unfiltered version of spike snippet');
                getRawWaveform = 1;
            else
                fprintf('\n    did NOT find spikes_unfilt folder, so spikeInfo will include unfiltered version of spike snippet');
                getRawWaveform = 0;
            end



            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            % GRAB THE JACKSHEET SPLIT USED TO MAKE THE binaries that were sorted (get from sortedBy path rather than physio path... could do a confirmation that these two things match)
            %   then identify which rows of the jacksheet correspond to the sort.txt files... this will be used to make correctly-named output channels
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %- get the jacksheet with rename strings from the sorted spike folder
            jacksheetPath = fullfile(sortedByPath,'_extractionInfo','jacksheetBR_justSplitChan.csv'); %- we could get this from the sort folder...

            % if legacy folder, lets copy it in first
            if contains(SpikeBinDir,'legacy')
                jacksheet_toCopy_Path = fullfile(procDir,sessID,'sort_reref','_extractionInfo','jacksheetBR_justSplitChan.csv');
                if exist(jacksheet_toCopy_Path,'file')
                    copyfile(jacksheet_toCopy_Path,jacksheetPath)
                end
            end

            if ~exist(jacksheetPath,'file')

                fprintf('\n ERROR: cant process spikeInfo without %s',jacksheetPath);
                keyboard
                continue

            else
                %- load split jackTable... this contains info on every channel split for spike sorting and will be used to convert the sort_txt filenames back into channels and device numbers
                tJackSplit  = readtable(jacksheetPath);
                %
                %- key fields: tJackSplit.ChanName, tJackSplit.newChanStr, tJackSplit.deviceNum

                iSort2jack = nan(1,length(sortTxtFilesCln));
                for iS=1:length(sortTxtFiles)
                    %- make the comparison with sortTxtFilesCln... this version of the file name has "(grabbedXXX)" cut out
                    iJ = find(strcmp(sortTxtFilesCln{iS}(1:end-4),tJackSplit.SortChanName)); %- cut ".txt" from sortTxtFile
                    if length(iJ)~=1
                        fprintf('\n uh oh... some kind of mismatch between channel string and jacktable: %s',sortTxtFilesCln{iS});
                        keyboard
                    else
                        iSort2jack(iS)=iJ;
                    end
                end
                tJackUse = tJackSplit(iSort2jack,:);
                fprintf('\n    found the following sort channels to process:\n\n');
                %disp(tJackUse);
            end



            %- get the sort-channel information ready here... it will get pushed into the spikeInfo struct further down
            chanNumsPerSort   = tJackUse.PhysicalChan; %- sinlge number per recording that should be stable across recordings unless NSPs are reconfigured (i.e., array moved from one NSP to another)
            deviceNumsPerSort = tJackUse.MicroDevNum;  %-
            newChanStrPerSort = tJackUse.ChanNameNew;  %-
            nsxChanStrPerSort = tJackUse.ChanName;     %-
            nsxFileStrPerSort = tJackUse.FileName;     %-  I think this will all be the same string because I'm looping over nsX files above
            nspSuffixPerSort  = tJackUse.NSPsuffix;


            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %  GRAB THE "sortNotes.xlsx" used to determine which channels to sort.
            %     that table will be included with the spikeInfo,
            %     AND the table will be used to confirm/estimate whether the sorts were actually completed
            %
            %    look for sortNotes(session)_sortedBy??.xlsx.  If found, create a table so those manual notes travel with the spikeInfo
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            % SJ: Moved this to a function (updateSortSummaryTxt.m)
            
            [sortNoteSummary, tSortNotes, flag_sortsIncomplete] = updateSortSummaryTxt(subj,sortedByPath,sessID,sortTxtFilesCln);
            
            %SJ: Look through sortNotes and find which channels have a
            %noise unit - put into tJackUse
            sortNotes = tSortNotes{1,2};
            %sortNotesUse = sortNotes(iSort2jack,:);
            %whichArray = size(sortNotes,1)/2 + 1;
            
            % Append to tJackUse, the columns: arrayNum, chanNoise, globalNoise
            
            [chanArrayMask,arrayList,firstIdxArray] = getArrayFromSortNoteTable(sortNotes);
            
            if any(strcmpi('dummySort',sortNotes.Properties.VariableNames))
                if any(strcmpi('chanNoise',sortNotes.Properties.VariableNames)) || any(strcmpi('globalNoise',sortNotes.Properties.VariableNames))
                    fprintf('%s\n','ERROR!!! Have dummySort column and either/both chanNoise and globalNoise! Ask SJ!');
                    keyboard
                else
                    %You have a dummy sort column... I guess for now let's just treat this as a globalNoise units
                    
                    if isnumeric(sortNotes.dummySort(iSort2jack)) %All numbers or NaN (ie no ???) SJ 11/30/2020
                        col2useDS = sortNotes.dummySort(iSort2jack);
                    else
                        col2useDS = str2double(sortNotes.dummySort(iSort2jack));
                    end
                    if isnan(col2useDS) %This is the same as is all parts are NaN
                        col2useDS = NaN(size(tJackUse,1),1);
                    end
                    tJackUse.globalNoise = col2useDS;
                    keyboard %Verify this works
                end
            else
                if any(strcmpi('globalNoise',sortNotes.Properties.VariableNames))
                    if isnumeric(sortNotes.globalNoise(iSort2jack)) %SJ 11/30/2020
                        col2useGN = sortNotes.globalNoise(iSort2jack);
                    else
                        col2useGN = str2double(sortNotes.globalNoise(iSort2jack));
                    end
                    if isnan(col2useGN)
                        col2useGN = NaN(size(tJackUse,1),1);
                    end
                    tJackUse.globalNoise = col2useGN;
                else
                    tJackUse.globalNoise = NaN(size(tJackUse,1),1);
                end
                if any(strcmpi('chanNoise',sortNotes.Properties.VariableNames))
                    if isnumeric(sortNotes.chanNoise(iSort2jack))  %SJ 11/30/2020
                        col2useCN = sortNotes.chanNoise(iSort2jack);
                    else   
                        col2useCN = str2double(sortNotes.chanNoise(iSort2jack));
                    end
                    if isnan(col2useCN)
                        col2useCN = NaN(size(tJackUse,1),1);
                    end
                    tJackUse.chanNoise = col2useCN;
                else
                    tJackUse.chanNoise = NaN(size(tJackUse,1),1);
                end
            end
            
            % Get rid of dummySort if it exists in tJackUse
            if any(strcmp('dummySort',tJackUse.Properties.VariableNames))
                if iscell(tJackUse.dummySort)
                    if any(~strcmp(tJackUse.dummySort,''))
                        fprintf('%s\n','ERROR!!! tJackUse.dummySort column is not blank!');
                        keyboard
                    end
                else
                    if ~all(isnan(tJackUse.dummySort))
                        fprintf('%s\n','ERROR!!! tJackUse.dummySort column is not blank!');
                        keyboard
                    end
                end
                tJackUse.dummySort = [];
            else
                keyboard
                % Do we ever get here? Probably should only get here if the
                % dummySort column stops getting made in jackSheet
            end
            
            tJackUse.arrayNum = chanArrayMask(iSort2jack);
            
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %  GRAB a copy of the POST-PROC cell array and PULSES cell array so the spikes can be aligned and any future changes in processing can be accounted for
            %-   use the pulseCellStruct that travels with the sorted spikes... that way the pulses and postProc info is about the time series used when sorting the data
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf('\n now the pulseStructCell and postProcStuctCell are saved in the _sorted folder... no need to reach back to the original SPK folder (which *could* change)');
            pulsesStructFile = fullfile(sortedByPath, '_extractionInfo', 'pulseStructCell.mat');
            postProcFile     = fullfile(sortedByPath, '_extractionInfo', 'postProcStructCell.mat');
            if ~exist(pulsesStructFile,'file') | ~exist(postProcFile,'file')
                fprintf('\n ERROR: sorts must be associated with a pulseStruct and postProcCell that can be used for subsequent alignment');
                keyboard
                error('\n fit it!');
            else
                pulsesStructCell = load(pulsesStructFile);
                pulsesStructCell = pulsesStructCell.pulsesStructCell;

                postProcStructCell = load(postProcFile);
                postProcStructCell = postProcStructCell.postProcStructCell;
            end

            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %  NOW CREATE THE SPIKE INFO, populating the static information just gathered (i.e., sortNote table), then loop over sorts and process
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %----------------------------------------------------------------------------------------------------------------------------------------------------------
            %----------------------------------------------------------------------------------------------------------------------------------------------------------

            % Save one structure with all spike info, combining across NSPs and DEVICES.  Each NSP has its own pulseStruct, so channels from eac NSP can be sepearately aligned in a subsequent step
            %- JW is making a single spikeInfo for all devices and NSPs, rather than breakup up by device
            spikeInfo = struct;

            %uniqueUnitIDvar = {'PhysChanNum' 'UnitNum' 'DeviceNum'   'CombinedNum'  'DeviceNumB'   'ChanUnitName' 'NSxChanName' 'NSxFileName' 'NSPsuffix' 'newChanStr'};

            %- this is a copy of the mountainsort info string...
            infoStr = sprintf('this spikeInfo file, generated %s, contains the following fields:\n',datestr(now,21));
            infoStr = sprintf('%s\n createdDate      - a string indicating when this spikeInfo.mat was created',infoStr);
            infoStr = sprintf('%s\n sessUniqueUnitID - a ( #unit x columns ) table of unit identifiers (see table headers or variable names and following description)',infoStr); 
            infoStr = sprintf('%s\n     containing: ( channel_num, unit_num, array_location_num, combined_num, chanUnit_name, nsx_chan_name, nsx_file_name, nsx_file_suffix, newChanName )',infoStr);
            infoStr = sprintf('%s\n          where: channel_num=phys_chan from jacksheet,  unit_num=sort output num,  array_location_num=devNum from jacksheet,  combined_num=devNum.chanNum.unitNum,  array_num=devNum(redundant)',infoStr);
            infoStr = sprintf('%s\n               : chanUnit_name=ChanNameNew from jacksheet with sort suffix,  nsx_chan_name=original ChanName from jacksheet or nsx,  nsx_file_name=full ns5 or ns6 filename',infoStr);  
            infoStr = sprintf('%s\n               : nsx_file_suffix=suffix pulled from nsx filename, ChanNameNew=updated micro name in from getMicroSubjInfo)',infoStr);
            infoStr = sprintf('%s\n timeStamp        - a ( #unit x 1 ) cell array with each cell containing the timestamps (in ms) of the correspondingly indexed unit in sessUniqueUnitID',infoStr);
            infoStr = sprintf('%s\n sessNoiseUnitID - The equivalent of the sessUniqueUnitID, except for the noise units specified by chanNoise or globalNoise, with the following new fields:',infoStr);
            infoStr = sprintf('%s\n               : chanNoise=unit # from sortNoteTable sorted as chanNoise, globalNoise=unit # from sortNoteTable sorted as globalNoise for that array, arrayNum=used to apply globalNoise to that array',infoStr);
            infoStr = sprintf('%s\n timeStampNoise   - a ( #unit x 1 ) cell array equivalent of timeStamp for noise units in sessNoiseUnitID',infoStr);
            infoStr = sprintf('%s\n noiseMask        - a ( #nsp x 3 ) cell array with:',infoStr);
            infoStr = sprintf('%s\n               : chanNoise/globalNoise[array] name',infoStr);
            infoStr = sprintf('%s\n               : ( #samples in LFP x 1 ) binary array where 1s indicate noise during that sample in aligned LFP data. createNoiseMaskSPK.m params: windowDur=5000, numChanReq=2',infoStr);
            infoStr = sprintf('%s\n               : %% of LFP data covered by the mask',infoStr);
            if SAVE_ALL_SNIPPET_WAVES
                infoStr = sprintf('%s\n waveform_plexon     - (optional) a (#unit x 1) cell array containing the waveform from each trial, as output within the plexon.txt sort file',infoStr);
            end
            infoStr = sprintf('%s\n waveForm         - a struct containing the mean and sd waveform (#unit x 210samples) from spike band (300-3000Hz)',infoStr);
            if getRawWaveform
                infoStr = sprintf('%s\n waveForm_raw     - contains waveform mu and sd from unfiltered (1-500Hz bandpass) timeseries [only output if spike_unfilt folder exists]',infoStr);
            end
            infoStr = sprintf('%s\n metrics          - a struct containing various quality scores for each unit (see readme within metrics for complete description)',infoStr);

            infoStr = sprintf('%s\n sessStr          - the session name the spikes were extracted from',infoStr);
            infoStr = sprintf('%s\n sessDurMin       - the duration of session in minutes',infoStr); %[mountainsort way]
            infoStr = sprintf('%s\n startTime_datenum - startTime of the session in matlab datenum format (computed from date string in sessStr)',infoStr); %[mountainsort way]
            infoStr = sprintf('%s\n extractInfoStr   - a string describing how many potential units were found by *manual sort* on each channel',infoStr);
            %infoStr = sprintf('%s\n referencing_info - [mountainsort only] a cell array containing: subject name, nsx filename [mountainsort output has other/more stuff here] ',infoStr);
            infoStr = sprintf('%s\n pulses              - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "pulse" struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries for that NSP',infoStr);
            infoStr = sprintf('%s\n physio_nsx_postProc - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "postProc" struct containing info about openNSx post processing file modifications used to construct the .bin files that were sorted',infoStr);
            infoStr = sprintf('%s\n alignedTo        - a string indicating which file the spikes have been aligned to',infoStr);
            infoStr = sprintf('%s\n alignmentChan    - a string indicating which channel in the pulses struct was used for alignment',infoStr);
            infoStr = sprintf('%s\n sort_dirname     - a string indicating the name of the sorted subdirectory housing the sort.txt files this spikeInfo was built from',infoStr);
            infoStr = sprintf('%s\n sort_filenames   - a cell of filenames of the sort.txt files used to make this spike info. particularly useful for confirming that the spikeInfo is up to date',infoStr);
            infoStr = sprintf('%s\n jackTableSplit   - jacksheetBR table including all "split" micro channels from this session, with device numbers and new channel names, and SortChanName',infoStr);
            infoStr = sprintf('%s\n jackTableUsed    - just the jacksheet for the channels incorporated into this spikeInfo. this table can be used to go back and forth between original and new channel names',infoStr);
            infoStr = sprintf('%s\n sortNoteTable    - (#sort notes x 2) cell array of sortNote.xls filenames and table contents. Table includes letter grades and notes made manually about each channel before and during sorting',infoStr);
            infoStr = sprintf('%s\n sortNoteSummary  - string summarizing #As,Bs,Cs,Ds,Fs reported and #sorted. indicates whether As,Bs,Cs all done or incomplete',infoStr);
            infoStr = sprintf('%s\n flag_sortsComplete - a binary flag indicating whether sorts *appear* to be complete (i.e., all "A" and "B" units + potentially all "C" or "D"',infoStr);


            spikeInfo.readme            = infoStr;

            %%- Closer to mountainsort output
            spikeInfo.createdDate       = datestr(now,21);

            %- data saved for each spike
            spikeInfo.sessUniqueUnitID  = {};    %sessUniqueUnitID - a ( #unit x 10 ) table [init as cell] with each row containing ID info for a unit ( channel_num, unit_num, array_location_num, combined_num, array_num, nsx_channel_name) +ADD NSX file name
            spikeInfo.timeStamp         = {};    %timeStamp - a ( #unit x 1 ) cell array with each cell containing the 30kHz timestamps of the correspondingly indexed unit in sessUniqueUnitID

            %SJ: New for noise units
            spikeInfo.sessNoiseUnitID   = {};
            spikeInfo.timeStampNoise    = {};
            
            %SJ: Adding in noise unit mask (will be empty if there aren't any noise units)
            spikeInfo.noiseMask         = {};
            
            if SAVE_ALL_SNIPPET_WAVES
                spikeInfo.waveForm_plexon   = {};    %waveform: - a ( #unit x 1 ) cell array containing all the snippet waveforms for every trial int16(total spikes x 32 samples).  possibly save this as a separate file or cut
            end

            %- vector value info saved for each unit
            spikeInfo.waveForm          = struct;  %- struct with mean and sd of each units waveform;  filtered and unfiltered
            spikeInfo.waveForm.infoStr  = 'mean and sd in uV, time in ms;  spike-band filtered (300-3000Hz)';
            spikeInfo.waveForm.timeMS   = []; %- updated in the middle of the code when this is computed
            spikeInfo.waveForm.mean     = [];
            spikeInfo.waveForm.sd       = [];
            spikeInfo.waveForm_raw      = spikeInfo.waveForm;
            spikeInfo.waveForm_raw.infoStr  = 'mean and sd in uV, time in ms;  unfiltered (1-5000Hz) after global reref';


            %- metrics contains structs and cells arrays with isolation quality metrics per uniqueUnitID
            spikeInfo.metrics           = struct;  %- struct of struct of struct:   metrics = { unitVSnoise joshua btwnUnits}.  unitVSnoise = {sigNoiseR isoScore fnScore fpScore},  joshua={sigNoiseR isoScore fnScore fpScore}, btwnUnits = cell of cell for all pairwise comparisons on a channel

            infoStr = 'spike Info metrics contains:';
            infoStr = sprintf('%s\n unitVSnoise, Joshua, btwnUnits are each structs containing sigNoiseR,isoScore,fnScore,fpScore, See Joshua et al 2007 for definition. sigNoise is >=0; isoScore between 0 and 1, false neg and positive is ???',infoStr);
            infoStr = sprintf('%s\n     unitVSnoise -- this is metric AJ and JW used in 2016 and 2018 papers. computed using plexon time series from manually identified unit and noise ',infoStr);
            infoStr = sprintf('%s\n     Joshua      -- exact match to Joshua method. This grabs fresh snippes from the .bin file, interpolates to a higher sample rate, and pulls fresh noise samples based on unit amplitude.',infoStr);
            infoStr = sprintf('%s\n     btwnUnits   -- cell array of iso measures (using plexon time series) for all pairwise comparisons of units on one electrode. ',infoStr);;
            infoStr = sprintf('%s\n    Note: JW believes a reasonable and very conservative isolation score for a unit would be minimum of unitVSnoise and btwnUnits',infoStr);
            infoStr = sprintf('%s\n ISIcontamination   -- vector of fraction of ISIs that were less than absolute refractory period of 2ms',infoStr);
            infoStr = sprintf('%s\n spkRate_eachMinute -- cell array with time series of spike rate (spikes/sec) per 1min interval for duration of recording. This is a first order stability measure.',infoStr);

            %- array or cell array saved for each channel
            spikeInfo.metrics.readme  = infoStr;
            spikeInfo.metrics.unitVSnoise = struct;    %- placeholder: populated below with subfields in place
            spikeInfo.metrics.joshua      = struct;
            spikeInfo.metrics.btwnUnits   = struct;


            spikeInfo.metrics.fracISIbelowThr   = [];  % per unit, fraction of spikes with ISI below absolute refractory period of 2ms
            spikeInfo.metrics.spkRate_eachMinute = {};  %- spike rate per minute vs time, in .5min steps

            %- abandoned outputs
            %spikeInfo.metrics.classResults       = {}; %- % Use LDA to classify different clusters % 1. Between spike and noise clusters % 2. All pairwise between units.  Outputs always >95% accuracy... not very useful so dont bother.
            %spikeInfo.metrics.isoScore_stability = {}; %- isoscore vs time.  very expensive to compute, so dont do it
            %spikeInfo.metrics.invalidTimeRange   = {}; %- not actually computed, so dont include in struct

            %- session-level info that is here to match the mountainsort spikeInfo and just because info like this is good to keep with the data
            spikeInfo.sessStr            = sessID;
            spikeInfo.sessDurMin         = tJackUse.DurationMin(1);     %-
            spikeInfo.startTime_datenum  = datenum(sessID, 'yymmdd_HHMM');
            spikeInfo.extractInfoStr     = {};     %- populated line-by-line below with each loaded sort file
            %spikeInfo.referencing_info  = {subj, fileNamesNSX{iF}};  %- Manual version (here) is different from Mountainsort version... just a diagnostic so dont sweat it
            spikeInfo.pulses             = pulsesStructCell; %- will be used for alignment
            spikeInfo.physio_nsx_postProc= postProcStructCell; %- will be used for diagnostic of how time series was split from openNSx

            spikeInfo.alignedTo          = ''; %- placeholder, will be populated during subsequent alignment step
            spikeInfo.alignmentChan      = ''; %- placeholder, will be populated during subsequent alignment step

            %- cell array of sort files
            spikeInfo.sort_dirname       = sortDir;
            spikeInfo.sort_filenames     = sortTxtFiles; %- will contain the "(grabbedXXX)" part of the file name

            %- complete jacktable and subset included in this spikeInfo
            spikeInfo.jackTableSplit     = tJackSplit;
            spikeInfo.jackTableUsed      = tJackUse;

            spikeInfo.sortNoteTable      = tSortNotes;
            spikeInfo.sortNoteSummary    = sortNoteSummary;
            spikeInfo.flag_sortsComplete = ~flag_sortsIncomplete;


            %- grand total of what is being processed... might get split into different NSPs in a minute
            fprintf('\n    processing %d sort files', length(sortTxtFiles));

            %- loop over channels and process them
            for f = 1:length(sortTxtFiles)
                %SJ: Get channel name to be able to compare it to channels
                %with noise units in 'SortChanName' column
                thisSortChanName        = char(regexp(sortTxtFiles{f},'.*(?=_\(grabbed.*)','match'));
                
                thisSpikeFileName       = sortTxtFiles{f};
                thisSpikeFileNameClean  = sortTxtFilesCln{f};                               %- sort.txt filename with "(grabbedXXYYZZ)" removed... need that to associate with the correct raw file
                thisSpikeFilePath       = fullfile(sortTxtPath,thisSpikeFileName);          %- path to actual sort.txt file (which could contain "grabbed)

                %- path to the binary files with physio time series. This path will get passed to getChanUnitQuality AND will be used below to pull out the mean waveform
                thisBinPhysioPath       = fullfile(procDir,sessID,SpikeBinDir,thisSpikeFileNameClean); %- path to matched sort.bin 
                thisBinPhysioPath(strfind(thisBinPhysioPath,'.txt'):end) = '.bin';
                thisBinPhysioPath_raw   = fullfile(procDir,sessID,physioDirUnfilt,thisSpikeFileNameClean); %- this 
                thisBinPhysioPath_raw(strfind(thisBinPhysioPath_raw,'.txt'):end) = '.bin';

                % confirm the physio binary and raw binary are present.. if not code will break below
                if ~exist(thisBinPhysioPath,'file')
                    fprintf('\n ERROR: cant find expected bin file: %s',thisBinPhysioPath);
                    keyboard
                end
                if getRawWaveform && ~exist(thisBinPhysioPath_raw,'file')
                    fprintf('\n ERROR: unfiltered physio folder exists, but cant find expected bin file: %s',thisBinPhysioPath_raw);
                    keyboard
                end


                %- get the sort.txt modification date... this is the date it was grabbed and/or created
                d = dir(thisSpikeFilePath);
                thisSpikeFileDateStr = d.date;


                %- update to command line
                fprintf('\n  %d/%d) Processing channel %s (modified %s; devNum %d; nsxChanNum %d; nsxChanStr %s, nsxFile %s, newChanName %s):',f, length(sortTxtFiles),thisSpikeFileName,thisSpikeFileDateStr, deviceNumsPerSort(f),chanNumsPerSort(f),nsxChanStrPerSort{f},nsxFileStrPerSort{f},newChanStrPerSort{f});
                fprintf('\n    [Load sort.txt: '); tStartChan = tic;



                % Extract spike times
                % Note: importing the text file from plexon causes the values to be reshaped into a N by 1 vector containing both waveform
                % information and the spike time stamps. This is organized so that there is a fixed number of rows for each unit, which starts with
                % the unit number (i.e. 0,1,2), the time stamp, then the waveform values. Code below separates these

                %- first open and look at second line of file ot see what delimeter is being used.  first line could be header, so ignore that
                fid = fopen(thisSpikeFilePath);
                tline = fgetl(fid); tline = fgetl(fid); % jumping past the first line, which might be a header.
                if tline~=-1
                    if contains(tline,',')
                        dlm = ','; %- figure out the delimator for this file
                    else
                        dlm = '';
                    end
                    %- read the sorted text file into a matrix.  Use the delimiter read from line 2 above to adapts to comma or tab based delimiteres.
                    %  always skip the first row to protect against whether or not sorter included headers during export (with the acceptable loss of a single spike)
                    thisSpikeContent = dlmread(thisSpikeFilePath,dlm,1,0); %- jump past header then read all the data
                else
                    thisSpikeContent = [];
                end
                fclose(fid);


                % Get the channel number for this spike....
                chanIDnum    = chanNumsPerSort(f); %- now this is the physical channel number in jacksheetBR (a unique identifier)
                chanIDstr    = thisSpikeFileName(1:strfind(thisSpikeFileName,'.txt')-1);  %- will be used to name the sort figure.  Keep the "(grabedXXX)" nomenclature


                % deal with unusual case of empty file... this is a sort that was opened and obviously didnt have spikes, so empty .txt was saved to indicate complete
                if isempty(thisSpikeContent)
                    fprintf('----No units found for channel (empty file) %s',chanIDstr);

                    %- extraction info saved for each channel (even if no spikes present... allows for easy bookkeeping of number of channels with putative spikes)
                    spikeInfo.extractInfoStr     = cat(1,spikeInfo.extractInfoStr,sprintf('%s (modified %s) --> 0 units, 0 total spikes, 0 noise snippets (empty file)',thisSpikeFileName,thisSpikeFileDateStr));
                    fprintf('\n WARNING: empty file is rare... usually a channel with no units still has unsorted snippets. make sure file was loaded/saved properly');
                    continue;

                elseif size(thisSpikeContent,2) == 2

                    fprintf('\n WARNING: looks like sort file only contains timestamp and unit number.  not setup for this yet. must switch to all .bin based analysis');
                    keyboard;

                end
                
                % Go look for this channel in tJackUse - find it with PhysicalChan == chanIDnum
                idxChan = find(tJackUse.PhysicalChan == chanIDnum);
                arrayNum = tJackUse.arrayNum(idxChan);
                if any(strcmp('chanNoise',tJackUse.Properties.VariableNames)) && tJackUse.chanNoise(idxChan) > 0 %SJ: add "any"
                    % chanUnit - doesn't matter which array it is from.
                    chanUnit = tJackUse.chanNoise(idxChan);
                else
                    chanUnit = 0;
                end
                if  tJackUse.globalNoise(idxChan) > 0
                    globalUnit = tJackUse.globalNoise(idxChan);
                    %arrayNum = tJackUse.arrayNum(idxChan);
                else
                    globalUnit = 0;
                    %arrayNum = 0;
                end
                
                twoNoise = false;
                if (globalUnit ~= 0) && (chanUnit ~= 0)
                    fprintf('%s\n','Heads up!!! There is a case where these is a global noise unit AND a chan noise unit. Verify that this will work.');
                    twoNoise = true;
                    keyboard
                end

                %%%% Get the timestamp, quality, and waveform values of unit(s) in this channel %%%
                % Get unit numbers
                
                % Get rid of manual noise units
                % !!!! Need to also get rid of 0 units (case where a
                % channel only has 0 and a noise unit - do NOT want that to
                % appear in the normal units)
                exc_ManNoiseUnits = thisSpikeContent;
                if globalUnit > 0
                    exc_ManNoiseUnits = exc_ManNoiseUnits(exc_ManNoiseUnits(:,1)~=globalUnit,:);
                end
                if chanUnit > 0
                    %exc_ManNoiseUnits = exc_ManNoiseUnits~=chanUnit;
                    %SJ: I have no idea why the above line was used initially
                    exc_ManNoiseUnits = exc_ManNoiseUnits(exc_ManNoiseUnits(:,1)~=chanUnit,:);
                end
                only_ManNoiseUnits = thisSpikeContent(thisSpikeContent(:,1)>0 & (thisSpikeContent(:,1)==globalUnit | thisSpikeContent(:,1)==chanUnit),:);
                
                %allSortedUnitNum = thisSpikeContent(:,1);
                allSortedUnitNum = exc_ManNoiseUnits(:,1);
                allSortedUnitNum_noise = only_ManNoiseUnits(:,1);
                
                allUnitNum       = unique(allSortedUnitNum(allSortedUnitNum>0))'; %- exclude invalidated (-1) and noise (0)
                rowCount         = length(allUnitNum);
                
                
                
                % Get time stamps
                %timeStamp        = thisSpikeContent(:,2)*1000;   % in ms
                timeStamp        = exc_ManNoiseUnits(:,2)*1000;   % in ms
                timeStamp_noise  = only_ManNoiseUnits(:,2)*1000;
                
                % Get waveForm values
                %waveFormValues   = thisSpikeContent(:,3:end);
                waveFormValues   = exc_ManNoiseUnits(:,3:end);
                waveFormValues_noise = only_ManNoiseUnits(:,3:end);

                %- enforce that waveforms are interger values spanning int16 range.  they will be saved as int16 later.
                if SAVE_ALL_SNIPPET_WAVES
                    if sum(waveFormValues(:)-round(waveFormValues(:)))>0
                        fprintf('\n WARNING: plexon time series saved as float instead of int16... converting');
                        keyboard;
                        if max(abs(waveFormValues(:)))*4 >= 2^15
                            fprintf('\n uh oh... looks like it wasnt scaled by 4, so round and save');
                            waveformValues = round(waveformValues);
                        else
                            waveformValues = round(waveformValues * 4); %- assume its in uV... convert back to int16 (blackrock scales by 4)
                        end
                    end
                end


                %%%% Combine data in order to sort according to units
                combinedMat = cat(2,allSortedUnitNum,timeStamp,waveFormValues);
                combinedMat = sortrows(combinedMat,1);
                
                combinedMat_noise = sortrows(cat(2,allSortedUnitNum_noise,timeStamp_noise,waveFormValues_noise),1);

                sorted_allUnitNum     = combinedMat(:,1);
                sorted_timeStamp      = combinedMat(:,2);
                sorted_waveFormValues = combinedMat(:,3:end);
                
                sorted_allUnitNum_noise     = combinedMat_noise(:,1);
                sorted_timeStamp_noise      = combinedMat_noise(:,2);
                sorted_waveFormValues_noise = combinedMat_noise(:,3:end);                

                unitCount             = length(unique(sorted_allUnitNum(sorted_allUnitNum > 0)));
                unitCount_noise       = length(unique(sorted_allUnitNum_noise(sorted_allUnitNum_noise>0)));
                
                % Get quality measures and save figure for this channel
                fprintf('%.1fs (%d units, %d total spikes);   Quality: ',toc(tStartChan), rowCount, sum(allSortedUnitNum>0)); tStartQual = tic;


                if ~isdir(figSaveFolder), mkdir(figSaveFolder); end;
                figSaveDir = fullfile(figSaveFolder,chanIDstr);
                [isoMetrics,classResults,isoScore_stability,spkRate_stability,invalidTimeRange] = getChanUnitQuality(sorted_allUnitNum,sorted_timeStamp,sorted_waveFormValues,figSaveDir,thisBinPhysioPath,samplerate);
                isoScore_stability_cell = {}; spkRate_stability_cell = {};
                for u = 1:unitCount
                    if ~isempty(isoScore_stability), isoScore_stability_cell{u,1} = isoScore_stability(u,:); end %- usually dont export this
                    spkRate_stability_cell{u,1}  = spkRate_stability(u,:);
                end
                fprintf('%.1fs]',toc(tStartQual));


                % Get the longer versions of waveform mean and sd
                preTroughMS = -2; preTroughSamples  = preTroughMS*(samplerate/1000);
                postTroughMS = 5; postTroughSamples = postTroughMS*(samplerate/1000);
                expectedTroughMS = [-1.0 1.65];  expectedTroughSamples = round(expectedTroughMS*(samplerate/1000))-preTroughSamples+1; %- 1ms before to 1.6 ms after... not symmetric beacause trigger should be before trough
                expectedTroughSamples = expectedTroughSamples(1):expectedTroughSamples(2);
                sampleBuffer = 100;   % Get more samples in order to trough-align
                waveform_mean = []; waveform_sd = [];
                x_30KHz = preTroughSamples:postTroughSamples-1;

                %- open and load the spike-band filtered time series used to sort
                fid     = fopen(thisBinPhysioPath,'r');
                binData = fread(fid, 'int16');
                fclose(fid);

                % optionally open and load the non-filtered (or low-passed (raw)) version of the waveforms (happens automatically if the unfilt folder exists)
                if getRawWaveform==1
                    waveform_raw_mean = []; waveform_raw_sd = [];
                    fid = fopen(thisBinPhysioPath_raw,'r');
                    binDataRaw = fread(fid, 'int16');
                    fclose(fid);
                end

                %- loop over the units
                %SJ: Not even going through here with the noise units... that's okay, right?
                for uu = allUnitNum
                    
                    
                    %- pull the extended time series from the bin file
                    unitTimeStamps  = floor(sorted_timeStamp(sorted_allUnitNum==uu)*(samplerate/1000));
                    outOfBoundCount = sum(unitTimeStamps < abs(preTroughSamples) | unitTimeStamps > length(binData)-postTroughSamples);
                    u_waveforms     = nan(length(unitTimeStamps)-outOfBoundCount,postTroughSamples-preTroughSamples + 2*sampleBuffer);
                    if getRawWaveform==1, u_waveforms_raw = u_waveforms; end

                    %- loop over spikes and grab wider versions of the snippets from the .binaries
                    for uuu = 1:length(unitTimeStamps)
                        if (unitTimeStamps(uuu)+preTroughSamples - sampleBuffer >= 0) && (unitTimeStamps(uuu)+postTroughSamples + sampleBuffer <= length(binData))
                            u_waveforms(uuu,:) = binData(unitTimeStamps(uuu)+preTroughSamples - sampleBuffer : unitTimeStamps(uuu)+postTroughSamples + sampleBuffer-1);
                            if getRawWaveform==1
                                 u_waveforms_raw(uuu,:) = binDataRaw(unitTimeStamps(uuu)+preTroughSamples - sampleBuffer : unitTimeStamps(uuu)+postTroughSamples + sampleBuffer-1);
                            end
                        end
                    end
                    %- might need to convert to double in this step...
                    this_wf_mean    = nanmean(double(u_waveforms),1);
                    this_wf_sd      = std(double(u_waveforms),1,1);
                    trough_region_i = expectedTroughSamples + sampleBuffer;  %- [was hardcoded 40:120, JW changed to ms based] expected region for trough: -1 to +1ms around spike time... if outside of this range then not this spike

                    if getRawWaveform==1
                        this_wf_mean_raw    = nanmean(double(u_waveforms_raw),1);
                        this_wf_sd_raw      = std(double(u_waveforms_raw),1,1);
                    end

                    % Flip waveform if upward
                    outputSign = 1;
                    if this_wf_mean(abs(this_wf_mean)==max(abs(this_wf_mean(trough_region_i)))) > 0
                        this_wf_mean = -this_wf_mean;
                        outputSign = -1;
                    end


                    % Trough-align
                    [~, min_i] = min(this_wf_mean);
                    if min_i < min(trough_region_i) | min_i > max(trough_region_i)
                        %- make the search window more narrow
                        [~, min_ii] = min(this_wf_mean(trough_region_i));
                        min_i = trough_region_i(min_ii); %- convert back to this_wf_mean index
                    end
                    this_wf_mean = this_wf_mean(min_i+preTroughSamples : min_i+postTroughSamples-1);
                    this_wf_sd   = this_wf_sd(  min_i+preTroughSamples : min_i+postTroughSamples-1);


                    %- check alignment
                    if this_wf_mean(x_30KHz==0) ~= min(this_wf_mean(expectedTroughSamples)),
                        fprintf('\n Heads up: Spike Trough align failed... can happen with crazy noise.');
                        figure; plot(x_30KHz,this_wf_mean); hold on; plot([0 0],get(gca,'ylim'),'k--');
                        title(sprintf('%s/%s: Filtered Waveform for %s, unit %d failed trough align',subj, sessID, thisSpikeFileName,uu));
                    end


                    %- save to output structure
                    waveform_mean = cat(1,waveform_mean,this_wf_mean*outputSign); %- unflip so saved data has true spike shape
                    waveform_sd   = cat(1,waveform_sd,this_wf_sd);

                    if getRawWaveform==1
                        %- use the shift identified with non-raw trough alignment.  with raw there can be a big offset in absolute voltage so flip/unflip trick doesn't work anyway
                        this_wf_mean_raw = this_wf_mean_raw(min_i+preTroughSamples : min_i+postTroughSamples-1);
                        this_wf_sd_raw   = this_wf_sd_raw(  min_i+preTroughSamples : min_i+postTroughSamples-1);

                        %- save to the output structure
                        waveform_raw_mean = cat(1,waveform_raw_mean,this_wf_mean_raw); %- for raw we dont need to flip and unflip... use the filtered time series to pick out peak/trough
                        waveform_raw_sd   = cat(1,waveform_raw_sd,this_wf_sd_raw);
                    end
                end

                % Save spikes to output structure
                if ~isempty(isoMetrics.unitVSnoise.sigNoiseR) %- this means there was at least 1 unit from the manual sort

                    % Convert unit number to unique identifier including channel and unit info ###(chan)##(unit)
                    noiseUnits = sorted_allUnitNum<=0;
                    sorted_allUnitNum(noiseUnits)       = [];
                    sorted_timeStamp(noiseUnits)        = [];
                    sorted_waveFormValues(noiseUnits,:) = [];
                    sorted_unitIdentifier = nan(length(sorted_allUnitNum),1);
                    for r = 1:length(sorted_allUnitNum)
                        sorted_unitIdentifier(r,1) = str2double(sprintf('%d%02d',chanIDnum,sorted_allUnitNum(r)));
                    end

                    unitsToWriteOut           = sort(unique(sorted_allUnitNum));
                    sorted_timeStamp_writeOut = cell(length(unitsToWriteOut),1);
                    sorted_waveForm_writeOut  = cell(length(unitsToWriteOut),1);
                    for i=1:length(unitsToWriteOut)
                        sorted_timeStamp_writeOut{i,1} = single(sorted_timeStamp(find(sorted_allUnitNum==unitsToWriteOut(i))))';
                        sorted_waveForm_writeOut{i,1}  = int16(sorted_waveFormValues(find(sorted_allUnitNum==unitsToWriteOut(i)),:))';  %- int16 ok as long as plexon was saved as binary values
                    end
                    uniqueUnitIDvar = {'PhysChanNum' 'UnitNum' 'DeviceNum'  'CombinedNum' 'ChanUnitName' 'NSxChanName' 'NSxFileName' 'NSPsuffix' 'ChanNameNew'};
                    sorted_uniqueUnitID_writeOut = cell(length(unitsToWriteOut),length(uniqueUnitIDvar));
                    for i=1:length(unitsToWriteOut)
                        sorted_uniqueUnitID_writeOut{i,1}  = chanNumsPerSort(f);    % f is index of channel number, taken from filename
                        sorted_uniqueUnitID_writeOut{i,2}  = unitsToWriteOut(i);
                        sorted_uniqueUnitID_writeOut{i,3}  = deviceNumsPerSort(f);
                        sorted_uniqueUnitID_writeOut{i,4}  = deviceNumsPerSort(f)*1e6 + chanNumsPerSort(f)*1e3 + unitsToWriteOut(i); % JW switched to numeric to match Chris's mountainsort output. first digit=array num, next 3 digits=chan num, next 3=unit num
                        sorted_uniqueUnitID_writeOut{i,5}  = sprintf('%s_%02d',newChanStrPerSort{f},unitsToWriteOut(i));
                        sorted_uniqueUnitID_writeOut{i,6}  = nsxChanStrPerSort{f};
                        sorted_uniqueUnitID_writeOut{i,7}  = nsxFileStrPerSort{f};
                        sorted_uniqueUnitID_writeOut{i,8}  = nspSuffixPerSort{f};
                        sorted_uniqueUnitID_writeOut{i,9}  = newChanStrPerSort{f};
                    end

                    %%%  populate the output structure %%%
                    %- data saved for each spike
                    spikeInfo.timeStamp          = cat(1,spikeInfo.timeStamp,sorted_timeStamp_writeOut);
                    spikeInfo.sessUniqueUnitID   = cat(1,spikeInfo.sessUniqueUnitID,sorted_uniqueUnitID_writeOut);

                    %-  mountaionsort saves this as a separate .mat file for... JW agrees not much rationale in saving here. Skip saving altogether for now. if saving to mat file use the wider-time samples grabbed from .bin above for the mu and sd
                    if SAVE_ALL_SNIPPET_WAVES
                        spikeInfo.waveForm_plexon    = cat(1,spikeInfo.waveForm_plexon,sorted_waveForm_writeOut);
                    end

                    %-
                    conv2uV = 0.25; %- for the waveform mean and sd, multiply by blackrock's scalefactor of 0.25 to get microVolts.  This assumes (correctly) that extractMicrophys outputs unscaled int16 values for sorting
                    spikeInfo.waveForm.mean      = cat(1,spikeInfo.waveForm.mean,waveform_mean*conv2uV);
                    spikeInfo.waveForm.sd        = cat(1,spikeInfo.waveForm.sd,  waveform_sd*conv2uV);
                    spikeInfo.waveForm.timeMS    = x_30KHz/30; %- convert from samples to ms. this is always the same, so just a single vector not a matrix of values
                    if getRawWaveform == 1
                        spikeInfo.waveForm_raw.mean   = cat(1,spikeInfo.waveForm_raw.mean,waveform_raw_mean*conv2uV);
                        spikeInfo.waveForm_raw.sd     = cat(1,spikeInfo.waveForm_raw.sd,  waveform_raw_sd*conv2uV);
                        spikeInfo.waveForm_raw.timeMS = x_30KHz/30; %- convert from samples to ms. this is always the same, so just a single vector not a matrix of values
                    end

                    % Isolation metrics
                    allIsoTypes   = fieldnames(isoMetrics);
                    allIsoMetrics = fieldnames(isoMetrics.(char(allIsoTypes(1))));
                    % Initialize isometrics structure if necessary
                    if isempty(fieldnames(spikeInfo.metrics.(allIsoTypes{1}))) %- new way to check, because now initialzing fields above and also adding fields not in isoMetric output
                        for ff = allIsoTypes'
                            for fff = allIsoMetrics'
                                spikeInfo.metrics.(char(ff)).(char(fff)) = [];
                            end
                        end
                    end
                    for ff = allIsoTypes'
                        for fff = allIsoMetrics'
                            spikeInfo.metrics.(char(ff)).(char(fff)) = cat(1,spikeInfo.metrics.(char(ff)).(char(fff)),isoMetrics.(char(ff)).(char(fff)));
                        end
                    end

                    %- cell array info saved for each unit
                    spikeInfo.metrics.fracISIbelowThr    = []; %- placeholder... computed below after all units are evaluated
                    %spikeInfo.metrics.isoScore_stability = cat(1,spikeInfo.metrics.isoScore_stability,isoScore_stability_cell); %- JW decommisioned this one
                    spikeInfo.metrics.spkRate_eachMinute = cat(1,spikeInfo.metrics.spkRate_eachMinute,spkRate_stability_cell);

                    %- extraction info saved for each channel
                    spikeInfo.extractInfoStr     = cat(1,spikeInfo.extractInfoStr,sprintf('%s (modified %s) --> %d units, %d total spikes, %d chanNoise units, %d globalNoise units',thisSpikeFileName,thisSpikeFileDateStr,size(sorted_uniqueUnitID_writeOut,1),size(sorted_unitIdentifier,1),double(chanUnit>0),double(globalUnit>0)));
                    
                    fprintf(' --> %d units, %d total spikes, %d chanNoise units, %d globalNoise units',length(unique(sorted_unitIdentifier)),size(sorted_unitIdentifier,1),double(chanUnit>0),double(globalUnit>0));
                else
                    %- extraction info saved for each channel (even if no spikes present... allows for easy bookkeeping of number of channels with putative spikes)
                    spikeInfo.extractInfoStr     = cat(1,spikeInfo.extractInfoStr,sprintf('%s (modified %s) --> 0 units, 0 total spikes, %d chanNoise units, %d globalNoise units',thisSpikeFileName,thisSpikeFileDateStr, double(chanUnit>0),double(globalUnit>0)));
                    fprintf(' --> 0 units, 0 total spikes, %d chanNoise units, %d globalNoise units', double(chanUnit>0),double(globalUnit>0));
                end
                
                %Now do the exact same thing for noise units, minus all the metrics
                if ~isempty(sorted_allUnitNum_noise)
                    % Convert unit number to unique identifier including channel and unit info ###(chan)##(unit)
                    sorted_unitIdentifier_noise = nan(length(sorted_allUnitNum_noise),1);
                    for rr2 = 1:length(sorted_allUnitNum_noise)
                        sorted_unitIdentifier_noise(rr2,1) = str2double(sprintf('%d%02d',chanIDnum,sorted_allUnitNum_noise(rr2)));
                    end
                    
                    unitsToWriteOut_noise           = sort(unique(sorted_allUnitNum_noise));
                    sorted_timeStamp_writeOut_noise = cell(length(unitsToWriteOut_noise),1);
                    sorted_waveForm_writeOut_noise  = cell(length(unitsToWriteOut_noise),1);
                    for ss2=1:length(unitsToWriteOut_noise)
                        sorted_timeStamp_writeOut_noise{ss2,1} = single(sorted_timeStamp_noise(find(sorted_allUnitNum_noise==unitsToWriteOut_noise(ss2))))';
                        sorted_waveForm_writeOut_noise{ss2,1}  = int16(sorted_waveFormValues_noise(find(sorted_allUnitNum_noise==unitsToWriteOut_noise(ss2)),:))';  %- int16 ok as long as plexon was saved as binary values
                    end
                    uniqueUnitIDvar_noise = {'PhysChanNum' 'UnitNum' 'DeviceNum'  'CombinedNum' 'ChanUnitName' 'NSxChanName' 'NSxFileName' 'NSPsuffix' 'ChanNameNew', 'chanNoise', 'globalNoise', 'arrayNum'};
                    sorted_noiseUnitID_writeOut = cell(length(unitsToWriteOut_noise),length(uniqueUnitIDvar_noise));
                    for tt3=1:length(unitsToWriteOut_noise)
                        sorted_noiseUnitID_writeOut{tt3,1}  = chanNumsPerSort(f);    % f is index of channel number, taken from filename
                        sorted_noiseUnitID_writeOut{tt3,2}  = unitsToWriteOut_noise(tt3);
                        sorted_noiseUnitID_writeOut{tt3,3}  = deviceNumsPerSort(f);
                        sorted_noiseUnitID_writeOut{tt3,4}  = deviceNumsPerSort(f)*1e6 + chanNumsPerSort(f)*1e3 + unitsToWriteOut_noise(tt3); % JW switched to numeric to match Chris's mountainsort output. first digit=array num, next 3 digits=chan num, next 3=unit num
                        sorted_noiseUnitID_writeOut{tt3,5}  = sprintf('%s_%02d',newChanStrPerSort{f},unitsToWriteOut_noise(tt3));
                        sorted_noiseUnitID_writeOut{tt3,6}  = nsxChanStrPerSort{f};
                        sorted_noiseUnitID_writeOut{tt3,7}  = nsxFileStrPerSort{f};
                        sorted_noiseUnitID_writeOut{tt3,8}  = nspSuffixPerSort{f};
                        sorted_noiseUnitID_writeOut{tt3,9}  = newChanStrPerSort{f}; %ChanNameNew
                        if twoNoise %If we have 2 noise units (1 chanNoise and 1 globalNoise), then we do not want to assign both numbers to both rows, only one each!
                            if unitsToWriteOut_noise(tt3) == chanUnit;
                                sorted_noiseUnitID_writeOut{tt3,10} = chanUnit; %chanNoise
                                sorted_noiseUnitID_writeOut{tt3,11} = 0; %globalNoise
                            elseif unitsToWriteOut_noise(tt3) == globalUnit
                                sorted_noiseUnitID_writeOut{tt3,10} = 0; %chanNoise
                                sorted_noiseUnitID_writeOut{tt3,11} = globalUnit; %globalNoise
                            else
                                keyboard % how can this not be equal to the chanNoise unit or GlobalNoise unit??
                            end
                        else
                            sorted_noiseUnitID_writeOut{tt3,10} = chanUnit; %chanNoise
                            sorted_noiseUnitID_writeOut{tt3,11} = globalUnit; %globalNoise
                        end
                        sorted_noiseUnitID_writeOut{tt3,12} = arrayNum; %arrayNum
                    end
                    
                    %%%  populate the output structure %%%
                    %- data saved for each spike
                    spikeInfo.timeStampNoise    = cat(1,spikeInfo.timeStampNoise,sorted_timeStamp_writeOut_noise);
                    spikeInfo.sessNoiseUnitID   = cat(1,spikeInfo.sessNoiseUnitID,sorted_noiseUnitID_writeOut);
                    
                else
                    % Nothing
                end

            end  %- loop over sort files


            % Look at spike times and see how many are within absolute refractory period
            absRef = 2;   % 2ms
            fracISIbelowThr = [];
            for unitID = 1:size(spikeInfo.sessUniqueUnitID,1)
                thisUnitTimeStamp = spikeInfo.timeStamp{unitID};
                contamination     = sum(diff(thisUnitTimeStamp) < absRef)/length(thisUnitTimeStamp);
                fracISIbelowThr  = cat(1,fracISIbelowThr,contamination);
            end
            spikeInfo.metrics.fracISIbelowThr = cat(1,spikeInfo.metrics.fracISIbelowThr,fracISIbelowThr);


            %- convert cell array to table before saving struct
            if exist('uniqueUnitIDvar','var') & ~isempty(spikeInfo.sessUniqueUnitID)
                spikeInfo.sessUniqueUnitID = cell2table(spikeInfo.sessUniqueUnitID, 'VariableNames',uniqueUnitIDvar);
            else
                spikeInfo.sessUniqueUnitID = {'No isolated units in this Session (even though there are sort.txt files)'};
                spkInfoPath = regexprep(spkInfoPath, 'spikeInfo.mat','spikeInfo[NoUnits].mat'); 
            end
            % Do the same for noise units
            if exist('uniqueUnitIDvar_noise','var') & ~isempty(spikeInfo.sessNoiseUnitID)
                spikeInfo.sessNoiseUnitID = cell2table(spikeInfo.sessNoiseUnitID, 'VariableNames',uniqueUnitIDvar_noise);
            else
                spikeInfo.sessNoiseUnitID = {'No isolated noise units in this Session'};
            end
            
            spikeInfo.noiseMask = 'Use createNoiseMaskSPK.m to create mask, or see aligned spikeInfo.';
            
            

            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Now save as a structure
            %%%%%%%%%%%%%%%%%%%%%%%%%
            save(spkInfoPath,'spikeInfo'); %- path defined above
            fprintf('\n All channels sorted and processed, and saved to %s [%.1f min total]',spkInfoPath, toc(tStartSortTxt)/60);

        end  %- if MAKE_NEW_SPIKEINFO
    end %- for iDir sortFolders
    fprintf('\n\n>>> Finished extracting spike information for all %d sort_txt folders in SESSION %s [total time %.1f min]',length(sortDirList), sessID,toc(tStartSess)/60);
end %- for 1:length(sessFolders)
if length(sessFolders)>1,  fprintf('\n\n>>> Finished extracting spike information for %s, %d sessions [total time %.1f min]\n',subj,length(sessFolders),toc(tStartAll)/60); end

end %- end function



%%%%- HELPER FUNCTION
function folders = getDirFileNames(path)

folders = dir(path);

for k = length(folders):-1:1
    
    % remove folders starting with .
    fname = folders(k).name;
    if fname(1) == '.'
        folders(k) = [ ];
    end
    
    
end
% Just get the names
folders = {folders.name};
end

