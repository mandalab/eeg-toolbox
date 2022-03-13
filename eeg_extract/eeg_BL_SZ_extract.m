function events = eeg_BL_SZ_extract(subj, rootEEGdir, timestamps, type)
%EEG_BL_SZ_EXTRACT
%
% This function extracts and pre-processes seizure and/or baseline clips.
% If raws correspond to seizures, a seizure event will be created, all
% other raws will be treated as baselines (i.e. no behavioral processing)
%
%Usage:
%EEG_BL_SZ_EXTRACT(subj,rootEEGdir,timestamps,type)
%
% The type should be 'baseline', 'seizure', or 'all'

%% set up directories, tag names, and create jackMaster_file if it doesn't exist
force_sz_raws = 0; % 1 to try and pull all raws from the seizure table, 0 to process only the input timestamps

%- confirm subject directory exists
subjDir = fullfileEEG(rootEEGdir,subj);
behDir   = fullfileEEG(subjDir,'behavioral','seizure');  % 'behavioral', but could also have 'behavioral_preOp' or 'behavioral_postOp' here
norerefDir = fullfile(subjDir, 'eeg.noreref');
rawDir   = fullfile(subjDir, 'raw');

if ~exist(subjDir,'dir')
    error('Error: root eeg directory: %s \n       does not contain subject: %s',rootEEGdir, subj);
end
if ~exist(behDir, 'dir')
    mkdir(behDir);
end

fprintf('\n\n\n**************************************************************************************'); %%- Command line output
fprintf('\n************************* %s Seizure/Baseline Extraction *************************', subj);     %%- Command line output
fprintf('\n**************************************************************************************\n');     %%- Command line output

%%- check for the jacksheetMaster.txt and create if not found.  This would normally happen during a call to nk_split, but add here for subjects that are already extracted
% import from_eegtoolbox.*; % requires the full eegtoolbox for now, not just the pre-packaged ones with DT

fprintf('\n checking for jacksheetMaster.csv: ');
jackMaster_file = fullfileEEG(subjDir, 'docs/jacksheetMaster.csv');

if ~exist(jackMaster_file, 'file')
    fprintf('creating...\n');
    createMasterJack(subj, rootEEGdir);
end

%% get events from the seizure master table if type= seizure, all, or ignore if type=baseline

if any(strcmpi(type,{'seizure','all'}))
    fprintf('\nLooking for seizure event spreadsheet...\n');
    NKT_DIR_OVERRIDE = []; % empty results in default for eegRawServerGrabFiles
    [sz_info,~,events,isSubjFound, sz_filepath] = load_SZ_info(subj);
    rem_ind=isnan(sz_info.EEG_Change_s); % files which are missing from the server, but present in the table
    sz_info(rem_ind,:)=[];
    if iscell(events) && length(events) == 1
        events = events{1}; % cell --> struct
        events(rem_ind)=[];
    else
        warning('The events structure appears to be in the wrong format for seizures')
    end
    
    if isSubjFound
        % check for raw existance
        all_raws_fullpath   = findRaws(rawDir);
        all_timestamps = cellfun(@fileparts, all_raws_fullpath, 'uniformoutput',0);
        all_timestamps      = fileparts_last(all_timestamps);
        raw_ind = ismember(all_timestamps,timestamps);
        raws   = fileparts_last(all_raws_fullpath(raw_ind));
        [missing_raws, missing_ind] = setdiff(sz_info.rawFile, raws);
        
        if force_sz_raws
            if ~isempty(missing_raws)
                fprintf('Not all raws foundin %s; trying to grab them from /EEG server...\n', rawDir);
                disp(missing_raws);
                
                try
                    eegRawServerGrabFiles(subj, rootEEGdir, [], missing_raws,[], NKT_DIR_OVERRIDE);
                catch e
                    warning('Failed to grab raws with eegRawServerGrabFiles: %s', e.message);
                end
                
                % check for raw existance
                all_raws_fullpath   = findRaws(rawDir);
                all_timestamps = cellfun(@fileparts, all_raws_fullpath, 'uniformoutput',0);
                all_timestamps      = fileparts_last(all_timestamps);
                raw_ind = ismember(all_timestamps,timestamps);
                raws   = fileparts_last(all_raws_fullpath(raw_ind));
                
                [missing_raws, missing_ind] = setdiff(sz_info.rawFile, raws);
%                 events(missing_ind)=[];
            end
        end
        
        events(missing_ind)=[];
        sz_raws_ind = ismember(raws,sz_info.rawFile);
        if ~all(sz_raws_ind) && strcmpi(type,'all')
            fprintf('The following raws (with corresponding timestamp) were found in the ''raw'' directory but do not exist in the seizure table');
            disp([raws(~sz_raws_ind)',timestamps(~sz_raws_ind)'])
            fprintf('Treating these raws as baseline sessions, you will be prompted to create events manually.\n')
        elseif strcmpi(type,'seizure')
            fprintf('Only extracting seizure clips...')
            all_raws_fullpath = all_raws_fullpath(sz_raws_ind); % Should only look at and extract seizure clips
            raws = raws(sz_raws_ind); 
            timestamps = timestamps(sz_raws_ind); 

            sz_raws_ind= ones(size(timestamps));
        else
            error('Unrecognized user input, type must be ''seizure'', ''baseline'', or ''all''')
        end
        
    else
        % check for raw existance        
        all_raws_fullpath   = findRaws(rawDir);
        all_timestamps = cellfun(@fileparts, all_raws_fullpath, 'uniformoutput',0);
        all_timestamps      = fileparts_last(all_timestamps);
        raw_ind = ismember(all_timestamps,timestamps);
        raws   = fileparts_last(all_raws_fullpath(raw_ind));
        
        sz_raws_ind=zeros(1,length(timestamps));
        events=[];
        fprintf('\tSubject %s not found in table %s\n', subj, sz_filepath);
        fprintf('Treating all raws as baseline sessions, you will be prompted to create events manually.\n')
        disp([raws',timestamps'])
    end
    
elseif strcmpi(type,'baseline')
    % check for raw existance
    all_raws_fullpath   = findRaws(rawDir);
    all_timestamps = cellfun(@fileparts, all_raws_fullpath, 'uniformoutput',0);
    all_timestamps      = fileparts_last(all_timestamps);
    raw_ind = ismember(all_timestamps,timestamps);
    raws   = fileparts_last(all_raws_fullpath(raw_ind));
    
    sz_raws_ind=zeros(1,length(timestamps));
    events=[];
    fprintf('\nBaseline mode, treating raws as baseline sessions, you will be prompted to create events manually.\n')
    disp([raws',timestamps'])
else
    error('Unrecognized user input, type must be ''seizure'', ''baseline'', or ''all''')
end

%% split files and update events structure so offsets are in number of samples
% [events, sz_timestamps] = behavioralProcessingSeizure(subj,rootEEGdir);
%fprintf('%d events\n', numel(events));
fprintf('\n All sessions to be processed:\n');
fprintf('\t%s\n', timestamps{:});

createRawInfo(subj, rootEEGdir);
for k=1:length(timestamps)
%     d_timestamp_raw = fullfile(rootEEGdir, subj, 'raw', all_timestamps{i});
    d_timestamp_split = fullfile(rootEEGdir, subj, 'eeg.noreref', timestamps{k});
    
    dirSet=dir(d_timestamp_split);
    dirSet(strncmp({dirSet.name}, '.', 1)) = []; %% remove folders with dots
    dirSet={dirSet.name};
    
    if numel(dirSet) <= 1
% %         fname_eeg = [char(findRaws(d_timestamp_raw)) '.EEG'];
% %         nk_split(subj, rootEEGdir, fname_eeg);
%         nk_split(subj, rootEEGdir, [all_raws_fullpath{i} '.EEG']);
        split_NK(subj, rootEEGdir, [all_raws_fullpath{k} '.EEG']);

    end
end

if ~isempty(events)
    % once the files are split, get the sampling rate and update the events
%     if exist(norerefDir, 'dir')
%         Fs = GetRateAndFormat(struct('eegfile',norerefDir));
%     else
%         fprintf('behavioralProcessing: %s not found. Cannot set eegoffset field until raws are split\n', norerefDir);
%     end
    
    % seizure events will not be aligned to anything, so set eegfile and eegoffset
    for k=1:length(events)
        ndx = find(strcmpi(raws, events(k).rawFile), 1);
        events(k).eegfile = fullfile(norerefDir, timestamps{ndx});
        Fs = GetRateAndFormat(fullfile(norerefDir,timestamps{ndx}));
        events(k).eegoffset = events(k).eegoffset_s * Fs;
        events(k).EEG_Change_ms = events(k).EEG_Change_s * Fs;
        events(k).Ictal_Onset_ms = events(k).Ictal_Onset_s * Fs;
        events(k).Clinical_Onset_ms = events(k).Clinical_Onset_s * Fs;
        events(k).Termination_ms = events(k).Termination_s * Fs;

        % datetime of the start of the seizure
        events(k).date.Hour=str2double(timestamps{ndx}(end-3:end-2));
        events(k).date.Minute=str2double(timestamps{ndx}(end-1:end));
        events(k).date=events(k).date + seconds(events(k).EEG_Change_s);

    end
    events = rmfield(events, 'eegoffset_s');
    events = rmfield(events, 'EEG_Change_s');
    events = rmfield(events, 'Ictal_Onset_s');
    events = rmfield(events, 'Clinical_Onset_s');
    events = rmfield(events, 'Termination_s');
    
    % Change name and save to root events.mat
    rootevntFileStr = sprintf('%s/events.mat',behDir);
    save(rootevntFileStr, 'events', '-v7');  %- version 7 file format is compact... only need version 7.3 if >2GB
    fprintf('\n --- %d events saved in %s ---\n', length(events), rootevntFileStr);
    %
end

% baseline events
if any(~sz_raws_ind)
    bl_timestamps=timestamps(~sz_raws_ind);
    disp('Entering manual event creation for baseline clips.\n')
    disp(bl_timestamps)
    state_vec=input('Please provide a binary state vector with one entry per timestamp, 1==awake, 0==asleep\n');
    create_BL_events(subj,bl_timestamps,state_vec,rootEEGdir)
    
end
fprintf('Final Step: processing seizure data from /eeg.noreref to /eeg.processed and /eeg.noreref...');
processAndReref(subj, rootEEGdir,'timestamps',timestamps);
end

%%
% function [isExtracted, firstExtractedFilename] = findAnyExtractions(subj, rootEEGdir, extractRootName)
% % look for a SINGLE expected extraction file. If we find ONE, assume that
% % extraction is complete
% chans = getLeads(subj, rootEEGdir);
% for i = 1 : length(chans)
%     expectedFilename = fullfileEEG(rootEEGdir,subj,'eeg.noreref',extractRootName, chans{i});
%     if exist(expectedFilename, 'file')
%         isExtracted = true;
%         firstExtractedFilename = expectedFilename;
%         return;
%     end
% end
% isExtracted = false;
% end

function x = fileparts_last(x)
for i=1:length(x), [~,x{i}] = fileparts(x{i}); end
end

