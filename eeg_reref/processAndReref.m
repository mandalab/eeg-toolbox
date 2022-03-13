function [stemsChecked, stemsRerefed] = processAndReref(subj, rootEEGdir, varargin)
% PROCESSANDREREF preprocesses and rereferences
%
% Function processAndReref(subj,rootEEGdir)
%
%   Description: This function implements both bipolar and average
%   re-referencing schemes. It also preprocesses by removing line noise and
%   running a local detrending algorithm and identifies bad channels
%
%   Input:
%         --subj=subject name (ex. 'NIH001')
%         --eegDir=directory where data is held
%                  (ex. '/Users/damerasr/Sri/data/eeg/')
%
%   Key-Value Pair Input
%       timestamps: cell array of timestamps to rereference (default: all in eeg.noreref)
%       disp: 1 to show and 0 to hide figures while they are being saved
%
%   Output:
%     stemsChecked: All timestamps that ARE currently processed and reref'd
%     stemsRerefed: All timestamps that were JUST NOW processed and reref'd
%
% Called by:
%   eegPrepAndAlign
%
% Revision History:
%   05/17 MST - Created. Based on rerefWrapper
%
% Dependencies:
%   eeg_toolbox/basic/
%
%
% See Also: rerefwrapper

% -------------------------------------------------------------------------
Z_THRESH = 3.0; % variance std dev cutoff
DUNGEON_TOOLBOX = 'eeg_noise_metrics';

cast = @(data,format) eval(sprintf('%s(data)', format)); % cast data to format

ip = inputParser;
ip.addParameter('timestamps',[]);
def_disp=1;
ip.addParameter('disp',def_disp,@(x) (x==1 || x==0));
ip.parse(varargin{:});
noreref_dates = ip.Results.timestamps;
disp_ind=ip.Results.disp;


% Chek toolbox
temp = which(DUNGEON_TOOLBOX);
if isempty(temp)
    error('to use processAndReref, %s must be on your path', DUNGEON_TOOLBOX);
end
clear temp;


% sets up input and output directories
subjDir    = fullfile(rootEEGdir,subj); % where subject data is located
noRerefDir = fullfile(subjDir,'eeg.noreref');     % where no-reref data is located
procDir    = fullfile(subjDir,'eeg.processed');   % processed: line noise and detrend
rerefDir   = fullfile(subjDir,'eeg.processedBP'); % where reref data will be saved to
location   = 'local';
% if strfound('FRNU', rootEEGdir) && strfound('volumes', lower(rootEEGdir))
if strfound(rootEEGdir,'FRNU') && strfound(lower(rootEEGdir),'volumes') % needs to be strfound(txt,pattern) otherwise always 0
    location = 'server';
end

%if the new dirs do not exist, make it
if ~exist(procDir,'dir')
    mkdir(procDir)
end
if ~exist(rerefDir,'dir')
    mkdir(rerefDir)
end

% get all session directories (and check noreref nonempty)
if isempty(noreref_dates)
    noreref_dates = lsCell(noRerefDir);
    noreref_dates = noreref_dates(cellfun(@isdir, fullfile(noRerefDir, noreref_dates)));
    noreref_dates = noreref_dates(cellfun(@isExpectedFilename, noreref_dates));
    if isempty(noreref_dates)
        error('noreref directory was empty. You must split raw files first: %s', noRerefDir);
    end
else
    matched = intersect(noreref_dates, lsCell(noRerefDir));
    unmatched = setdiff(noreref_dates, matched);
    if ~isempty(unmatched)
        warning('The following dates passed were not found!!')
        disp(unmatched);
    end
    noreref_dates = matched;
end

reref_dates = lsCell(rerefDir);
reref_dates = reref_dates(cellfun(@isdir, fullfile(noRerefDir, reref_dates)));
reref_dates = reref_dates(cellfun(@isExpectedFilename, reref_dates));

proc_dates = lsCell(procDir);
proc_dates = proc_dates(cellfun(@isdir, fullfile(noRerefDir, proc_dates)));
proc_dates = proc_dates(cellfun(@isExpectedFilename, proc_dates));

all_chans          = getLeads(subj, rootEEGdir, 'leadsType','all', 'chanType','PHYS');
all_subdural_chans = getLeads(subj, rootEEGdir, 'leadsType','all', 'chanType','PHYS', 'hardwareType','subdural');
all_nosubdur_chans = setdiff(all_chans, all_subdural_chans); %- depths and micro-subdural


% Preprocess (single line-noise removal + detrend) and reref every session
samplerate_unique = [];
dataformat_unique = [];
count = 0;
for i = 1:length(noreref_dates)
    timestamp = noreref_dates{i};
    d_proc_sess = fullfile(procDir, timestamp);
    d_reref_sess  = fullfile(rerefDir, timestamp);
    d_noreref_sess = fullfile(noRerefDir, timestamp);
    
    miss_proc_sess  = ~ismember(timestamp, proc_dates) || numel(lsCell(d_proc_sess)) < 2;
    miss_reref_sess = ~ismember(timestamp, reref_dates) || numel(lsCell(d_reref_sess)) < 2;
    
    %- 12/2018, JW added an easy check to see if processing folder exists but was not complete
    tempfile = fullfile(d_proc_sess,'_processingNotDone.txt');
    if ~miss_proc_sess && exist(tempfile,'file'),
        miss_proc_sess = 1; miss_reref_sess = 1;
    end
    
    if miss_proc_sess,
        count = count + 1;
        miss_reref_sess = 1;
        
        % if filestem hasn't already been processed then do it now
        % note that noise_metrics must be called session missing from EITHER folder  << jw 12/2018 says "Huh?"
        
        fprintf('Processing session %s... (%d of %d)\n', timestamp, i, length(noreref_dates));
        
        % In the future, handle stim files in this way:
        %   - bipolar rereference, but no global rereference
        %   - do processing
        % TODO: modify rerefBipolarity to optionally use noreref instead of processed
        % CURRENTLY (06/17): do everything
        is_stim = isStimFile(subj, rootEEGdir, timestamp);
        if is_stim
            % make directories / params.txt copy
            %if ~exist(d_reref_sess,'dir'), mkdir(d_reref_sess); end
            %filename = fullfile(noRerefDir,'params.txt');
            %copyfile(filename, d_reref_sess);
            %continue;
        end
        
        
        % make directories / params.txt copy (will happen separately below for reref folder
        if ~exist(d_proc_sess,'dir'),  mkdir(d_proc_sess); end
        filename = fullfile(noRerefDir,timestamp,'params.txt'); %- jw change 4/2018.. get params from indivial sessions... no longer stored at root because can change between sessions
        copyfile(filename, d_proc_sess);
        
        
        %- jw 12/2018, make it easy to detect when a folder is incomplete because of a processing crash
        fid_temp = fopen(tempfile,'w+');
        fclose(fid_temp);
        
        
        % filter to only those electrodes recorded in this session
        [session_chans] = nkElectrodeFilt(rootEEGdir, subj, timestamp);
        %if isempty(session_chans) % note that nkElectrodeFilt isempty means no filter (all chans present)
        %    session_chans = all_chans;
        %end
        session_subdural_chans = intersect(all_subdural_chans, session_chans);
        session_nosubdur_chans = intersect(all_nosubdur_chans, session_chans); %- depths and micro-subdural
        
        
        % -------------------------------------------------------------------------
        % Processing section
        % -------------------------------------------------------------------------
        % get sampling rate, data format, and gain
        [samplerate,~,dataformat,gain] = GetRateAndFormat(d_proc_sess);
        if ~isempty(session_subdural_chans)
            temp = gete(session_subdural_chans{1},struct('eegfile',fullfile(noRerefDir,timestamp)),0);
            nsamples = length(temp{1});
            clear temp;
            if nsamples < 100 *samplerate % 100 seconds
                fprintf('Error: %s has only %d samples!! It is probably not a legitimate task EEG\n', timestamp, nsamples);
                keyboard;
                continue;
            end
        else
            temp = gete(session_nosubdur_chans{1},struct('eegfile',fullfile(noRerefDir,timestamp)),0);
            nsamples = length(temp{1});
            clear temp;
            if nsamples < 100 *samplerate % 100 seconds
                fprintf('Error: %s has only %d samples!! It is probably not a legitimate task EEG\n', timestamp, nsamples);
                keyboard;
                continue;
            end
        end
        
        
        % Note:  it is imperative that the order of channels returned
        % from eeg_noise_metrics matches the order of session_chans. It may
        % be worth returning the actual names from noise_metrics to verify
        % this...
        fprintf('\n-------------------------------------------------------------------------\n');
        fprintf('Now calling eeg_noise_metrics. This process takes 5-10 minutes (turn on parpool!)\n');
        fprintf('Do not quit early; do not close the figures.\n');
        fprintf('Wait until you see the "ALL CLEAR" message.\nBe patient....\n');
        
        
        if ~isempty(session_subdural_chans)
            [bad_chans_new, bad_chans_old, eeg_raw_proc_subdural, glob_sig_all, glob_sig_clean, mvar_ref_Z_last, mamp_ref_Z_last, line_rel_amp, chan_names ...
                ] = eeg_noise_metrics(subj, timestamp, 0, Inf, Inf, 0, ...
                'outputdir', d_proc_sess, ...
                'z_thresh', Z_THRESH, ...
                'saving', 1, ...
                'loc_detrend', 1, ...
                'rem_sat', 10, ...
                'rmline', 1, ...
                'visual', 1, ...
                'disp', disp_ind, ...
                'location', rootEEGdir);
            % chan_names corresponds to output data. Ensure our list matches!
            assert(size(eeg_raw_proc_subdural, 2) == length(session_subdural_chans), '%d eeg_noise_metrics data points but %d subdural channels', size(eeg_raw_proc_subdural, 2), length(session_subdural_chans));
            chan_ndx = sortBySubstring(session_subdural_chans, chan_names); % sort to chan_name order
            session_subdural_chans = session_subdural_chans(chan_ndx);
            assert(all(strcmpi(chan_names, session_subdural_chans)), 'load_eegfile channel names do not match expected channel names (subdural)');
            bad_chan_names = union(bad_chans_old, bad_chans_new);
            if isempty(bad_chan_names), bad_chan_names = {}; end
            
            % keep track of nans so that they can be replaced with intmax('int16') after casting to int16
            eeg_raw_proc_subdural_nanMask = isnan(eeg_raw_proc_subdural);
            glob_sig_all_nanMask = isnan(glob_sig_all);
            glob_sig_clean_nanMask = isnan(glob_sig_clean);
            
            
            fprintf('ALL CLEAR - eeg_noise_metrics complete for %s!\n\n', timestamp);
            fprintf('-------------------------------------------------------------------------\n');
            
            
            % quick check on nsamples: compare gete to load_eegfile
            assert(size(eeg_raw_proc_subdural,1) == nsamples, 'load_eegfile has different # samples than gete');
            
                        
            % note that eeg_noise_metrics calls load_eegfile which applies
            % gain, so we need to divide by gain in order to follow the
            % convention of raw files not having gain
            glob_sig_all = cast(glob_sig_all ./ gain, 'int16');
            glob_sig_clean = cast(glob_sig_clean ./ gain, 'int16');
            eeg_raw_proc_subdural = cast(eeg_raw_proc_subdural ./ gain, 'int16');
            
            % when casting to int16 the nans become zeros, this replaces those with intmax('int16')
            glob_sig_all(glob_sig_all_nanMask) = intmax('int16');
            glob_sig_clean(glob_sig_clean_nanMask) = intmax('int16');
            eeg_raw_proc_subdural(eeg_raw_proc_subdural_nanMask) = intmax('int16');
            
            clear glob_sig_all_nanMask glob_sig_clean_nanMask eeg_raw_proc_subdural_nanMask
            
            samplerate_unique = unique([samplerate_unique, samplerate]);
            dataformat_unique = unique([dataformat_unique, {dataformat}]);
            
            %assert(length(samplerate_unique) == 1, 'Multiple samplerates across sessions!');
            % M.E: ^ not true anymore (11/18)
            assert(length(dataformat_unique) == 1, 'Multiple data formats across sessions!');
            
            
            % -------------------------------------------------------------------------
            % Write to File section
            % -------------------------------------------------------------------------
            filename = fullfile(d_proc_sess, 'global_avg_all');  %- this was saved in reref, but makes more sense to contain with processed
            fd = fopen(filename,'w');
            fwrite(fd, glob_sig_all, 'int16');
            fclose(fd);
            
            filename = fullfile(d_proc_sess, 'global_avg_good'); %- this was saved in reref, but makes more sense to contain with processed
            fd = fopen(filename,'w');
            fwrite(fd, glob_sig_clean, 'int16');
            fclose(fd);
            
            % save eeg_raw_proc to processed/timestamp
            fprintf('writing processed eeg chanels to %s...', d_proc_sess);
            for j = 1 : size(eeg_raw_proc_subdural,2)
                filename = fullfile(d_proc_sess, session_subdural_chans{j});
                if exist(filename, 'file'), continue; end
                fd = fopen(filename, 'w');
                if fd==-1, fprintf('\n fopen error on %s, trying again>',filename); pause(1); fd = fopen(filename, 'w');  end  %- jw 12/2018 occasionally hitting error on fwrite below. probably OS related max open files
                if fd==-1, fprintf('\n fopen error on %s, keyboard>',filename); keyboard; end
                fwrite(fd, eeg_raw_proc_subdural(:,j), 'int16');
                fclose(fd);
            end
            
            % write variance.csv file
            ndec = 4;
            variance = struct();
            variance.chanName = session_subdural_chans(:);
            variance.is_good  = double(~ismember(session_subdural_chans, bad_chan_names));
            variance.variance = round(mvar_ref_Z_last(:), ndec, 'decimals');
            variance.amplitude= round(mamp_ref_Z_last(:), ndec, 'decimals');
            variance.line_noise = round(line_rel_amp(:), ndec, 'decimals');
            var_table = struct2table(variance);
            filename = fullfile(d_proc_sess, 'variance.csv');
            writetable(var_table,filename);
        else
            
            %- add a note to let user know why variance stuff is missing...
            noSubNote = fullfile(d_proc_sess,'_NoSubdurals_soNoVarOrBadCalc.txt');
            fid_note = fopen(noSubNote,'w+');
            fclose(fid_note);
       
        end
        
        
        % Nonsubdural channels should be processed but not included in
        % global rereferencing scheme and not added to "bad channels"
        warning('off','DTB:GetFilesandChannels:emptyFileList');
        [eeg_raw_proc_nosubdur, ~, ~, chan_names, ~] = load_eegfile(subj,timestamp,0,Inf,...
            'loc_detrend',1,...
            'rmline',1,...
            'rem_sat',10,...
            'visual',0,...
            'ref','raw',...
            'hardwareType','~subdural',...
            'location',rootEEGdir);
        warning('on','DTB:GetFilesandChannels:emptyFileList');
        
        % chan_names corresponds to output data. Ensure our list matches!
        assert(size(eeg_raw_proc_nosubdur, 2) == length(session_nosubdur_chans), '%d eeg_noise_metrics data points but %d nonsubdural channels', size(eeg_raw_proc_nosubdur, 2), length(session_nosubdur_chans));
        chan_ndx = sortBySubstring(session_nosubdur_chans, chan_names); % sort to chan_name order
        session_nosubdur_chans = session_nosubdur_chans(chan_ndx);
        assert(all(strcmpi(chan_names, session_nosubdur_chans)), 'load_eegfile channel names do not match expected channel names (nonsubdural)');
        
        [samplerate,~,dataformat,gain] = GetRateAndFormat(d_proc_sess);
        eeg_raw_proc_nosubdur = cast(eeg_raw_proc_nosubdur ./ gain, 'int16');
        
        %- write the files
        for j = 1 : size(eeg_raw_proc_nosubdur,2)
            filename = fullfile(d_proc_sess, session_nosubdur_chans{j});
            if exist(filename, 'file'), continue; end
            fd = fopen(filename, 'w');
            if fd==-1, fprintf('\n fopen error on %s, trying again>',filename); pause(1); fd = fopen(filename, 'w');  end  %- jw 12/2018 occasionally hitting error on fwrite below. probably OS related max open files
            if fd==-1, fprintf('\n fopen error on %s, keyboard>',filename); keyboard; end
            fwrite(fd, eeg_raw_proc_nosubdur(:,j), 'int16');
            fclose(fd);
        end
        
        %- jw 12/2018, make it easy to detect when a folder is incomplete because of a processing crash
        delete(tempfile);
        if exist(tempfile,'file'),
            fprintf('\n uh oh, temp file still exists but just tried deleting... ');
            keyboard;
        end
        
        fprintf('done!\n');
        
    end
    
    
    
    % save reref averages  - JW moved this outside of the processed condition so rereferencing can be computed without re-doing processing
    if miss_reref_sess || miss_proc_sess, %- but if this session was re-processed, definitely re-ref again
        fprintf('writing rereferec eeg chanels to %s...\n', d_reref_sess);
        
        
        %- make the directory and copy params.txt
        if ~exist(d_reref_sess,'dir'), mkdir(d_reref_sess); end
        filename = fullfile(noRerefDir,timestamp,'params.txt'); %- jw change 4/2018.. get params from indivial sessions... no longer stored at root because can change between sessions
        copyfile(filename, d_reref_sess);
        
        % bipolar rereference (all channels)
        reref_Bipolarity(subjDir, timestamp);
        fprintf('Done!\n');
        
        
        % check
        DBG = 0;
        if DBG
            c=2;
            temp_proc = gete(session_subdural_chans{c},struct('eegfile',d_proc_sess),0);
            temp_raw = gete(session_subdural_chans{c},struct('eegfile',d_noreref_sess),0);
            figure; t1 = round(nsamples/2); t2 = t1+500;
            subplot(4,1,1); plot(glob_sig_all(t1:t2)); title('global all'); ylim([-1000,1000]);
            subplot(4,1,2); plot(glob_sig_clean(t1:t2)); title('global good'); ylim([-1000,1000]);
            subplot(4,1,3); plot(temp_proc{1}(t1:t2)); title([session_subdural_chans{1} ' proc']); ylim([-1000,1000]);
            subplot(4,1,4); plot(temp_raw{1}(t1:t2)); title([session_subdural_chans{1} ' raw']); ylim([-1000,1000]);
        end
    end
end

stemsChecked = length(noreref_dates);
stemsRerefed = count;

end % processAndReref



function yes = isExpectedFilename(filename)
% check if file has date_time format (look for 1 '_' chars) and numeric parts
parts = strsplit(filename, '_');
is_num = @(x) ~isempty(str2double(x)) && ~isnan(str2double(x));
yes = (length(parts) == 2) && is_num(parts{1}) && is_num(parts{2});
end


function yes = isStimFile(subj, rootEEGdir, timestamp)
% There is at this moment,  to Mike's knowledge, no definitive way to check whether a session is STIM
% or not. Best attempt is to look in raw/STIM/ for a timestamp
yes = exist(fullfile(rootEEGdir, subj, 'raw/STIM', timestamp), 'dir');
end