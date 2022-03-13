function [alignedFiles] = alignSpikesWithEcog(subject,ecogDir)
% Looks through behavioral files for subjects,
% in local folder, finds associated spikeInfo.mat files,
% aligns pulses from the two systems, and edits spikeInfo.mat
% and processed.mat such that a copy of spike and LFP data
% is made to produce data that is aligned with the offset provided
% in events.
%
% NOTE: before running this code, put the mountainSort output folders
% into NIHxxx/micro/
%
% code writes out new noreref, processed, and spikeinfo into subject/micro
% folder, as well as description of alignment (alignmentToSpikes.csv),
% which is also the output variable "alignedFiles"
% -- the offset from alignedFiles can be used and applied to the spike timestamps 
% -- and LFP indices produced through the manual spike sorting pipeline, but note
% -- that currently, this does not allow for analysis of spike-ecogPhase
% -- relationships, because there can be a few milliseconds offset that can
% -- only be fixed by applying the stretch factor, rather than only the offset
%
% note: ecog/behavioral dir should be local copy, pre-FRNU upload, so that new
% events won't be overwritten - this code will eventually become
% a part of eegPrepAndAlign, i.e. in the 56_Stage process

% code has 2 branches - either aligns based on DC09 for behavioral tasks,
% or based on DC12 for resting state tasks. for DC09, before processing
% a new subject, populate the cell 'subject_pulseChannels' with a string
% describing the name of the channel we want to align to for that patient
% (name can differ slightly, should open up a spikeInfo.mat and look at
% the variable 'pulses' to see

% there are a few assumptions being made that we should keep in mind
% 1) analog/digital pulses from blackrock are recorded at 30 kHz
% 2) pulses from ecog are at 1 kHz
% 3) noreref/processed LFPs are at 1 kHz


% tweaks needed
% minor - save out new variables that are being written to spikeInfo.mat, processed.mat, and noreref.mat
% major - swap out alignment machinery to incorporate same offset as in double-NSP sync

debug_noWrite = 0;
applyShiftCorrection = 0;

% subject, preferred pulse channel, backup pulse channel, suffix
subject_pulseChannels = {'NIH064' , 'ain1','','_spikeInfo.mat';
    'NIH065', 'ain1','','_spikeInfo.mat';
    'NIH047', 'ainp1','','_spikeInfo.mat';
    'NIH054', 'ain17','','_spikeInfo.mat';
    'NIH057', 'ain17','','_spikeInfo.mat';
    'NIH059', 'ain17','','_spikeInfo.mat';
    'NIH062', 'ain1','ain17','_spikeInfo.mat';
    'NIH050', 'ainp1','','_spikeInfo.mat';
    'NIH066', 'ain1','','_spikeInfo.mat'};


samprate_blackrock = 30000;



for subj = 1:length(subject);
    behavDir = fullfile(ecogDir,'/',subject{subj},'/behavioral/');
    behavFolders =  dir(behavDir);
    behavFolders = behavFolders(find([behavFolders.isdir]==1));
    behavFolders(find(contains({behavFolders.name},'cant'))) = [];
    behavFolders(find(contains({behavFolders.name},'.'))) = [];
    
    
    alignedFiles = cell(1,6);
    alignedFiles{1,1} = 'phys_sess'; alignedFiles{1,2} = 'beh_sess'; alignedFiles{1,3} = 'spikeInfo_loc';
    alignedFiles{1,4} = 't_offset'; alignedFiles{1,5} = 'alignmentSuccess'; alignedFiles{1,6} = 'drift_ms';
    alignedFiles{1,7} = 'drift_correction_applied';
    
    rowIterate = 2;
    
    for exp = 1:length(behavFolders);
        sessionDir = fullfile(behavDir,behavFolders(exp).name,'/');
        sessionFolders = dir(sessionDir);
        sessionFolders = sessionFolders(find(contains({sessionFolders.name},'session')));
        sessionFolders(find(contains({sessionFolders.name},'cant'))) = [];
        
        
        
        
        for sess = 1:length(sessionFolders)
            
            if contains(fullfile(sessionDir,sessionFolders(sess).name),'notReal','IgnoreCase',true) || contains(fullfile(sessionDir,sessionFolders(sess).name),'cantAlign','IgnoreCase',true) || contains(fullfile(sessionDir,sessionFolders(sess).name),'empty','IgnoreCase',true)
                continue
            end
            temp = strsplit(fullfile(sessionDir,sessionFolders(sess).name),'/');
            if strcmp(temp{1,length(temp)}(1),'_')
                continue
            end
            
            if contains(sessionDir,'restingState_jc')
                % special branch, for non-behavioral sessions
                whichDC = 12;
                
                timeStamps = dir(fullfile(sessionDir,sessionFolders(sess).name));
                timeStamps(find(contains({timeStamps.name},'.'))) = [];
                noreref_dir = timeStamps(1).name;
                date_ecog = datenum(noreref_dir,'yymmdd_HHMM');
                
                trigDCDir = fullfile(sessionDir,sessionFolders(sess).name,noreref_dir,'trigDC12.syncBR.txt');
                fileID = fopen(trigDCDir,'r');
                pulses_ecog = textscan(fileID,'%d'); pulses_ecog = double(pulses_ecog{1});
                fclose(trigDCDir);
                
                paramsDir = fullfile(sessionDir,sessionFolders(sess).name,noreref_dir,'params.txt');
                fileID = fopen(paramsDir,'r');
                samprate_ecog = textscan(fileID,'%s'); samprate_ecog = double(samprate_ecog{1,1}{2,1});
                fclose(fileID);
                
                
            else
                % regular branch
                whichDC = 9;
                sess_events = load(fullfile(sessionDir,sessionFolders(sess).name,'/events.mat'));
                sess_events = sess_events.events;
                noreref_dir = sess_events(1).eegfile;
                dashes = strfind(noreref_dir,'/'); noreref_dir = noreref_dir(dashes(end)+1:end);
                
                if ~exist(fullfile(ecogDir,'/',subject{subj},'/eeg.noreref/',noreref_dir))
                    fprintf('\nnoreref directory doesnt contain all files referenced by behavioral folders\n\n')
                    continue
                end
                
                trigDCDir = fullfile(ecogDir,'/',subject{subj},'/eeg.noreref/',noreref_dir,'trigDC09.sync.txt');
                fileID = fopen(trigDCDir,'r');
                pulses_ecog = textscan(fileID,'%d'); pulses_ecog = double(pulses_ecog{1});
                date_ecog = datenum(noreref_dir,'yymmdd_HHMM');
                fclose(fileID);

                
                paramsDir = fullfile(ecogDir,'/',subject{subj},'/eeg.noreref/',noreref_dir,'params.txt');
                fileID = fopen(paramsDir,'r');
                samprate_ecog = textscan(fileID,'%s'); samprate_ecog = double(str2num(samprate_ecog{1,1}{2,1}));
                fclose(fileID);

                
            end
            
            
            % pull the start dates from all the spikeInfos
            
            % this will be rewritten so that it looks in the subject/micro folder to
            % find the matching spikeSession
            
            spikeDir = fullfile(ecogDir,'/',subject{subj},'/micro/');
            spikeSessions = dir(fullfile(ecogDir,'/',subject{subj},'/micro/'));
            spikeSessions(find(contains({spikeSessions.name},'.'))) = [];
            underscores = strfind(spikeSessions(1).name,'_');
            spikeSessionDates = NaN(1,length(spikeSessions));
            for session = 1:length(spikeSessions)
                temp= spikeSessions(session).name(underscores(1)+1:underscores(end)-1);
                spikeSessionDates(1,session) = datenum(temp,'yymmdd_HHMM');
            end
            
            diffs  = (date_ecog-spikeSessionDates);
            % get session closest to start time
            % restrict to positive? negative?
            [~,closest] = find(abs(diffs) == min(abs(diffs)));
            
            if min(abs(diffs)) > 0.04
                fprintf(sprintf('No matching spikeInfo.mat within 1 hour for session %s\n',fullfile(ecogDir,'/',subject{subj},'/eeg.noreref/',noreref_dir)))
                continue
            end
            
            
            if size(closest,2) > 1
                fprintf('two utahs, or duplicate data.. code probably not prepared for this')
                keyboard
            end
            
            indexSubject = find(strcmp(subject{subj},subject_pulseChannels));
            
            %findSpikeInfo = dir(fullfile(spikeDir,spikeSessions(closest(1)).name,'*spikeInfo.mat'));
            
            
            
            
            
            % reorganize folder, if it hanst been organized yet
            % will enter this block of code only if the most recent mountainSort wasnt used (in which
            % case the below rearrangement will have already happened)
            if ~(exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/')))
                
                mkdir(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/'))
                mkdir(fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/'))
                mkdir(fullfile(spikeDir,spikeSessions(closest(1)).name,'sorting/'))
                
                
                
                if exist(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')))
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')),fullfile(spikeDir,spikeSessions(closest(1)).name,'/raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')));
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_sortSummary.csv')),fullfile(spikeDir,spikeSessions(closest(1)).name,'sorting/',sprintf('%s%s',spikeSessions(closest(1)).name,'_sortSummary.csv')));
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeWaveform.mat')),fullfile(spikeDir,spikeSessions(closest(1)).name,'sorting/',sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeWaveform.mat')));
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'sortFigs'),fullfile(spikeDir,spikeSessions(closest(1)).name,'sorting/','sortFigs'));
                end
                
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_noreref.mat')),fullfile(spikeDir,spikeSessions(closest(1)).name,'/raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_noreref.mat')));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,sprintf('%s%s',spikeSessions(closest(1)).name,'_processed.mat')),fullfile(spikeDir,spikeSessions(closest(1)).name,'/raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_processed.mat')));
                
                if (exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'ain_timeseries.mat')))
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'ain_timeseries.mat'),fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/','ain_timeseries.mat'));
                end
                
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'full_ts_raw.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','full_ts_raw.png'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'full_ts_ref.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','full_ts_ref.png'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'global_sigs.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','global_sigs.png'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'mean_spectral.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','mean_spectral.png'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'variance.csv'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','variance.csv'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'vars_amps.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','vars_amps.png'));
                movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'bad_chans.mat'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','bad_chans.mat'));
                
                if (exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw_compare.png')))
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw_compare.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','raw_compare.png'));
                    movefile(fullfile(spikeDir,spikeSessions(closest(1)).name,'reref_compare.png'),fullfile(spikeDir,spikeSessions(closest(1)).name,'cleaning/','reref_compare.png'));
                end
                
                
            end
            
            if exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')))
                suffix = '_spikeInfo.mat';
            else
                suffix = '_noreref.mat';
            end

            spikeInfoFilename = fullfile(spikeDir,spikeSessions(closest(1)).name,'/raw/',sprintf('%s%s',spikeSessions(closest(1)).name,suffix));
            load(spikeInfoFilename);
            spikePulses = pulses;
            % clear spikeInfo
            
            cutString = subject{subj};
            cutIndex = strfind(fullfile(behavDir,behavFolders(exp).name,'/',sessionFolders(sess).name,'/'),cutString);
            cutSpikeInfoIndex = strfind(spikeInfoFilename,cutString);

            
            alignedFiles{rowIterate,1} = noreref_dir;
            temp = fullfile(behavDir,behavFolders(exp).name,'/',sessionFolders(sess).name,'/');
            alignedFiles{rowIterate,2} = temp(cutIndex:end);
            temp = spikeInfoFilename(cutSpikeInfoIndex(1):end);
            temp = strrep(temp,'.mat','_aligned.mat');
            alignedFiles{rowIterate,3} = temp;
            rowIterate = rowIterate + 1;
            
            
            
            % find the exact name of the channel that carries pulses for
            % this patient
            indexPulseType = subject_pulseChannels{indexSubject,2};
            
            if whichDC==9
                try
                    pulses_br = eval(sprintf('spikePulses.%s',indexPulseType));
                catch
                    pulses_br = eval(sprintf('spikePulses.%s',subject_pulseChannels{indexSubject,3}));
                end
            else
                pulses_br = spikePulses.din4;
            end
            
            

            
            
            if whichDC==9;
                
                
                %fullfile(othersDir,sprintf('%s_%s_DC%02d_alignmentFig',subj,ecog_sessName,whichDC));
                pulses_br_cell = cell(1,length(pulses_br));
                for i=1:length(pulses_br); pulses_br_cell{1,i} = pulses_br(i); end
                
                % Use pulseAlign (old way to confirm)
                % Get the minimum time that elapsed between the two systems
                pulses_br = round(pulses_br./(samprate_blackrock/1000));
                br_duration = range(pulses_br)/samprate_ecog/60;
                ecog_duration = range(pulses_ecog)/samprate_ecog/60;
                minPulseDuration = min(br_duration,ecog_duration);
                fprintf('#############################################\n');
                fprintf('Alignment results for DC %d:\n',whichDC);
                % Run pulseAlign
                threshMS = 10;
                mywin = 100;
                foundMatch = 0;
                % keyboard
                while foundMatch == 0 && mywin < 300
                    try
                        [br_match,nk_match] = pulsealign(pulses_br,pulses_ecog,samprate_ecog,threshMS,mywin,0,0,whichDC);
                    catch
                        br_match = NaN;
                        fprintf('no matching start window..\n')
                        mywin = mywin + 20;
                    end
                    
                    if ~isnan(br_match)
                        % 8. Use regression to check for alignment. Get slope and offset for each eeg file
                        bfix = br_match(1);
                        [b,bint,r,rin,stats] = regress(nk_match, [ones(length(br_match),1) br_match-bfix]);
                        b(1) = b(1) - bfix*b(2);
                        
                        % calc max deviation
                        act=[ones(length(br_match),1) br_match]*b;
                        maxdev = max(abs(act - nk_match));
                        meddev = median(abs(act - nk_match));
                        rawDevs = act - nk_match;
                        
                    end
                    
                    
                    if range(br_match) > minPulseDuration*0.9   % Duration of matched pulses isn't too short
                        if abs(1-b(2)) < 0.001   % Regression slope isn't too off
                            if maxdev < 2   % Max deviation isn't too long
                                foundMatch = 1;
                            else
                                fprintf('%s\n','Max deviation is off. Adjusting time window...');
                            end
                        else
                            fprintf('%s\n','Regression slope off. Adjusting time window...');
                        end
                    elseif isempty(br_match)
                        fprintf('%s\n','No matches found. Adjusting time window...');
                    else
                        fprintf('%s\n','Short/incomplete matches. Adjusting time window...');
                    end
                    if foundMatch == 0
                        mywin = mywin + 20;
                    end
                end
                
                if ~isnan(br_match)
                    % 8. Use regression to check for alignment. Get slope and offset for each eeg file
                    bfix = br_match(1);
                    [b,bint,r,rin,stats] = regress(nk_match, [ones(length(br_match),1) br_match-bfix]);
                    b(1) = b(1) - bfix*b(2);
                    
                    % calc max deviation
                    act=[ones(length(br_match),1) br_match]*b;
                    maxdev = max(abs(act - nk_match));
                    meddev = median(abs(act - nk_match));
                    rawDevs = act - nk_match;
                    
                    % save regression stats in string that is dumped at the end of the function
                    fprintf('\tRegression Slope  = %f\n', b(2));
                    fprintf('\tR^2               = %f\n', stats(1));
                    fprintf('\tMedian Deviation  = %f ms\n', meddev);
                    fprintf('\tMax Deviation     = %f ms\n', maxdev);
                    fprintf('\tMatched Pulse range = %.3f minutes\n', range(br_match)/1000/60);
                    fprintf('#############################################\n');
                    
                    % Throw an error if the alignment is off
                    if b(2) < 0.99 || stats(1) < 0.99
                        fprintf(sprintf('Regression slope is way off for %s\n',fullfile(ecogDir,'/',subject{subj},'/eeg.noreref/',noreref_dir)))
                        continue
                        % error('Regression slope is %d',b(2));
                    end
                    
                    alignmentFigSaveDir = fullfile(behavDir,behavFolders(exp).name,'/',sessionFolders(sess).name,'/');
                    
                    pulses_ecog = round(pulses_ecog*1000/samprate_ecog);
                    [timeDiff_xcorr] = alignPulses(pulses_br,pulses_ecog,alignmentFigSaveDir,whichDC,1);
                    
                    timeDiff_regress = br_match(1) - nk_match(1);
                    timeDiff = timeDiff_regress;
                    
                    stretchFactor = [(b(2)*(br_duration*60*1000)) - (br_duration*60*1000)];
                    addOrRemove = '';
                    if b(2) < 1; addOrRemove = 'remove'; 
                    elseif b(2) > 1 addOrRemove = 'add'; 
                    end 
                    
                    b_original = b;
                    alignedFiles{rowIterate-1,6} = round(stretchFactor);
                    debug=1; validCorrection = 0;
                    if (stretchFactor > 5) && (debug==1) && (applyShiftCorrection==1)
                        if debug==1
                            
                            %% create idealized version of timeseries for pulses_br and pulses ecog,
                            % and see if shifted timeseries is a better fit
                            
                            br_analog = zeros(1,max(pulses_br));
                            br_analog(round(pulses_br)) = 1;
                            ecog_analog = zeros(1,max(pulses_ecog));
                            ecog_analog(round(pulses_ecog)) = 1;
                            
                            [xc,lags] = xcorr(br_analog',ecog_analog');
                            maxCorr = max(xc(:));
                            
                            new_br_analog = zeros(1,max(pulses_br) - round(stretchFactor));
                            
                            if b(2) < 1
                                % remove samples
                                remove = round(linspace(1,length(br_analog),stretchFactor));
                                new_br_analog = br_analog;
                                new_br_analog(remove) = [];
                                
                            elseif b(2) > 1
                                
                                % add samples. each sample should be an average of flanking
                                % samples
                                
                                % which samples to change
                                add = round(linspace(1,length(br_analog),round(stretchFactor)));
                                
                                step = 1;
                                iter = 0;
                                for i=1:length(add)-1
                                    % add in existing data
                                    if i==1
                                        new_br_analog(step+iter:iter+add(i+1)) = br_analog(add(i):add(i+1));
                                    else
                                        new_br_analog(step+iter:iter+add(i+1)) = br_analog(add(i):add(i+1)-1);
                                    end                                   
                                    
                                    if i==length(add)-1
                                        % copy the last sample
                                        new_br_analog(add(i+1)+iter) = nanmean([br_analog(add(i+1))]);
                                    else
                                        % take the mean of the flanking samples
                                        new_br_analog(add(i+1)+iter) = nanmean([br_analog(add(i+1)) br_analog(add(i+1)+1)]);
                                    end
                                    % where to start copying next time
                                    step = add(i+1)+1;
                                    iter = iter+1;
                                end
                                
                            end
                            
                            [xc,~] = xcorr(new_br_analog',ecog_analog');
                            maxCorr_shift = max(xc(:));
                            
                            pulses_br = find(new_br_analog==1)';
                            
                            % Use pulseAlign (old way to confirm)
                            % Get the minimum time that elapsed between the two systems
                            br_duration = range(pulses_br)/samprate_ecog/60;
                            ecog_duration = range(pulses_ecog)/samprate_ecog/60;
                            minPulseDuration = min(br_duration,ecog_duration);
                            %                 fprintf('#############################################\n');
                            %                 fprintf('Alignment results for DC %d:\n',whichDC);
                            fprintf('Checking that shift correction was good\n\n')
                            % Run pulseAlign
                            threshMS = 10;
                            mywin = 100;
                            foundMatch = 0;
                            % keyboard
                            while foundMatch == 0 && mywin < 300
                                try
                                    [br_match,nk_match] = pulsealign(pulses_br,pulses_ecog,samprate_ecog,threshMS,mywin,0,0,whichDC);
                                catch
                                    br_match = NaN;
                                    fprintf('no matching start window..\n')
                                    mywin = mywin + 20;
                                end
                                
                                if ~isnan(br_match)
                                    % 8. Use regression to check for alignment. Get slope and offset for each eeg file
                                    bfix = br_match(1);
                                    [b,bint,r,rin,stats] = regress(nk_match, [ones(length(br_match),1) br_match-bfix]);
                                    b(1) = b(1) - bfix*b(2);
                                    
                                    % calc max deviation
                                    act=[ones(length(br_match),1) br_match]*b;
                                    maxdev = max(abs(act - nk_match));
                                    meddev = median(abs(act - nk_match));
                                    rawDevs = act - nk_match;
                                    
                                end
                                
                                
                                if range(br_match) > minPulseDuration*0.9   % Duration of matched pulses isn't too short
                                    if abs(1-b(2)) < 0.001   % Regression slope isn't too off
                                        if maxdev < 2   % Max deviation isn't too long
                                            foundMatch = 1;
                                        else
                                            fprintf('%s\n','Max deviation is off. Adjusting time window...');
                                        end
                                    else
                                        fprintf('%s\n','Regression slope off. Adjusting time window...');
                                    end
                                elseif isempty(br_match)
                                    fprintf('%s\n','No matches found. Adjusting time window...');
                                else
                                    fprintf('%s\n','Short/incomplete matches. Adjusting time window...');
                                end
                                if foundMatch == 0
                                    mywin = mywin + 20;
                                end
                            end
                            
                            if ~isnan(br_match)
                                % 8. Use regression to check for alignment. Get slope and offset for each eeg file
                                bfix = br_match(1);
                                [b_new,bint,r,rin,stats] = regress(nk_match, [ones(length(br_match),1) br_match-bfix]);
                                b_new(1) = b_new(1) - bfix*b_new(2);
                                
                                % calc max deviation
                                act=[ones(length(br_match),1) br_match]*b_new;
                                maxdev = max(abs(act - nk_match));
                                meddev = median(abs(act - nk_match));
                                rawDevs = act - nk_match;
                                
                                if (abs(b_new(2)-1) <= abs(b_original(2)-1)) || (maxCorr_shift < maxCorr)
                                    fprintf('correction will increase deviation of slope from 0, not applying\n')
                                    validCorrection = 0;
                                    keyboard
                                else
                                    fprintf('correction will decrease deviation of slope from 0, applying\n')
                                    validCorrection = 1;
                                end
                                
                            end
                            
                            
                        end
                    end
                    
                    fprintf('finished aligning\n')
                    
                    if (applyShiftCorrection==1) && (stretchFactor>5)
                        alignedFiles{rowIterate-1,7} = 'Y';
                    else
                        alignedFiles{rowIterate-1,7} = 'N';
                    end
                    
                    
                else
                    
                    alignedFiles{rowIterate-1,5} = 'Regression completely fails';
                    microDir = fullfile(ecogDir,'/',subject{subj},'/micro/');
                    alignment_csv_writeOut = fullfile(microDir,'alignmentToSpikes.csv');
                    cell2csv(alignedFiles,alignment_csv_writeOut)
                    continue
                    
                end
                
                
                
                
                
                % Check if old and new metric are consistent. If it's off by more than 50ms, throw an error
                if ((timeDiff_regress - timeDiff_xcorr)/1000) >= 0.05, fprintf('Old and new metrics are inconsistent\n');
                    alignedFiles{rowIterate-1,5} = 'Regression and xcorr are INCONSISTENT (>50 ms)';
                    fprintf('Off by %0.2f seconds\n',((timeDiff_regress - timeDiff_xcorr)/1000))
                    %keyboard;
                else
                    alignedFiles{rowIterate-1,5} = 'Successful - regression and xcorr are CONSISTENT (<50 ms)';
                end
                
            end
            
            
            %time diff is in units of 1 kHz. timestamp is in units of 30 kHz
            if whichDC==12
                
                % if using DC12, let's try the cross-correlation technique
                minuteOffsetViaTimeStamp = (min(abs(diffs)) * 24 * 60);
                
                pulses_br_two = round(pulses_br./30);
                % this '3' in the line below is variable - need to access a
                % params.txt for each session in order to define it
                pulses_ecog_two = round(pulses_ecog./(1)); % assuming din4 is at 30 kHz, and NK file is
                % at 1 kHz
                
                timeseries_ecog = zeros(1,max(pulses_ecog_two));
                timeseries_ecog(pulses_ecog_two) = 1;
                
                timeseries_br = zeros(1,max(pulses_br_two));
                timeseries_br(pulses_br_two) = 1;
                
                
                blackrock_NotDownsampled = pulses_br;
                nk_NotDownsampled = pulses_ecog;
                blackrock_FileTime = spikeSessions(1);
                
                
                fprintf('still under construction\n'); keyboard
                [timeDiff_xcorr,~] = alignPulses_BwSystems(pulses_br_two,pulses_ecog_two,alignmentFigSaveDir,whichDC,'1');
                %[timeDiff,~] = alignPulses_BwSystems(pulses_ecog_two,pulses_ecog_two,alignmentFigSaveDir,whichDC,'1');
                
                timeDiff_regress = timeDiff_regress / 30;
                timeDiff_xcorr = timeDiff_xcorr / 30;
                
                
            end
            
            offsetCalculatedInMs = (timeDiff_regress);
            
            %% record and save changes made to each file in the subject folder
            alignedFiles{rowIterate-1,4} = sprintf('%0.2f seconds',(timeDiff_regress/1000)); %minus 1 because the counter has already iterated
            microDir = fullfile(ecogDir,'/',subject{subj},'/micro/');
            alignment_writeOut = fullfile(microDir,'alignmentToSpikes.mat');
            %save(alignment_writeOut,'alignedFiles')
            
            alignment_csv_writeOut = fullfile(microDir,'alignmentToSpikes.csv');
            cell2csv(alignedFiles,alignment_csv_writeOut)
            
            % check to see if we've already aligned this spike info before
            
            temp_find = spikeInfoFilename(cutSpikeInfoIndex(1):end);
            temp_find = strrep(temp_find,'.mat','_aligned.mat');
            
            
            triedToAlignBefore = find(strcmp(alignedFiles(1:rowIterate-2,3),temp_find));
            resultsFromBefore = alignedFiles(triedToAlignBefore,5);
            alignedThisSpikeInfoAlready  = size(contains(resultsFromBefore,'successful'),1);
            
            %alignedThisSpikeInfoAlready = size(strcmp(alignedFiles(1:rowIterate-2,3),spikeInfoFilename),1);
            
            if debug_noWrite == 0 && alignedThisSpikeInfoAlready == 0
                
                
                
                if exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')))
                    
                    % write new field with adjusted spike times
                    timeStamp_aligned = cell(size(timeStamp,1),1);
                    for unit = 1:size(timeStamp,1)
                        % multiply by 30 b/c correction is applied in 30 kHz
                        % samples
                        timeStamp_aligned{unit,1} = timeStamp{unit,1} - (timeDiff_regress*(samprate_blackrock/1000));
                    end
                    alignedTo = noreref_dir;
                    if whichDC==12; alignmentChan = indexPulseType; else; alignmentChan = indexPulseType; end
                    
                end
                
                
                
                
                % now make timestamp adjustment to the LFPs
                noreref_utah_Filename = strrep(spikeInfoFilename,'spikeInfo','noreref');
                load(noreref_utah_Filename);
                noreref = lfp; clear lfp;
                
                processed_utah_Filename = strrep(noreref_utah_Filename,'noreref','processed');
                load(processed_utah_Filename);
                processed = lfp; clear lfp;
                
                if (applyShiftCorrection==1) && (stretchFactor > 5) && (validCorrection==1);
                    
                    
                    %% shift spike times to align with clock rate of eCog
                    totalTime30kHz = sessDurSec*30000;
                    manipulate = round(linspace(1,totalTime30kHz,round(stretchFactor)));
                    
                    for unit = 1:size(timeStamp,1)
                        tempSpikes = timeStamp{unit,1};
                        hold = tempSpikes;
                        for i=1:length(manipulate)
                            if i==length(manipulate)
                                temp = find(tempSpikes > manipulate(i));
                            else
                                temp = intersect((tempSpikes > manipulate(i)),find(tempSpikes < manipulate(i+1)));
                            end
                            tempSpikes(temp) = tempSpikes(temp) + i;
                            
                        end
                    end
                    timeStamp_aligned{unit,1} = tempSpikes;
                    
                    %% edit noreref/processed to align with clock rate of eCog
                    
                    totalTime1kHz = sessDurSec*1000;
                    manipulate = round(linspace(1,totalTime1kHz,round(stretchFactor)));
                    if strcmp(addOrRemove,'remove')
                        
                        processed(:,manipulate) = [];
                        glob_sig_all(manipulate) = [];
                        glob_sig_good(manipulate) = [];
                        noreref(:,manipulate) = [];
                        
                    elseif strcmp(addOrRemove,'add') % add samples
                        new_processed = zeros(size(processed,1),size(processed,2) + round(stretchFactor));
                        new_noreref = zeros(size(noreref,1),size(noreref,2) + round(stretchFactor));
                        for chan=1:size(processed,1)
                            signal_process = processed(chan,:); iter = 0;
                            signal_noreref = noreref(chan,:);
                            
                            for i=1:length(manipulate)-1 % add in existing data
                                step = manipulate(i)+1;
                                if i==1
                                    signal_process(1+iter:iter+manipulate(i+1)) = processed(chan,manipulate(i):manipulate(i+1));
                                    signal_noreref(1+iter:iter+manipulate(i+1)) = noreref(chan,manipulate(i):manipulate(i+1));
                                else
                                    signal_process(step+iter:iter+manipulate(i+1)) = processed(chan,manipulate(i):manipulate(i+1)-1);
                                    signal_noreref(step+iter:iter+manipulate(i+1)) = noreref(chan,manipulate(i):manipulate(i+1)-1);
                                    
                                end
                                if i==length(add)-1 % copy the last sample
                                    signal_process(manipulate(i+1)+iter) = nanmean([processed(chan,manipulate(i+1))]);
                                    signal_noreref(manipulate(i+1)+iter) = nanmean([noreref(chan,manipulate(i+1))]);
                                else % take the mean of the flanking samples
                                    signal_process(manipulate(i+1)+iter) = nanmean([processed(chan,manipulate(i+1)) processed(chan,manipulate(i+1)+1)]);
                                    signal_noreref(manipulate(i+1)+iter) = nanmean([noreref(chan,manipulate(i+1)) noreref(chan,manipulate(i+1)+1)]);
                                    
                                end % where to start copying next time
                                iter = iter+1;
                            end
                            new_processed(chan,1:size(signal_process,2)) = signal_process;
                            new_noreref(chan,1:size(signal_noreref,2)) = signal_noreref;
                        end
                        processed = new_processed; %clear new_processed
                        noreref = new_noreref; %clear new_noreref
                        
                        %% finally, do the same for global averages
                        
                        
                        
                    end
                    
                    
                    
                    
                end
                
                
                % change in future - a single processed.mat will contain
                % two global_sig_alls and glob_sig_goods (for each array,
                % i.e. 1-64 and 65-128). after that, micros will have their
                % own reference, in which case, look at referencing_info
                
              
                % zero-padding (or removal?)
                if timeDiff_regress < 0 % if blackrock came later, remove data - original method
                    processed_aligned = processed(:,(abs(round(timeDiff_regress))):end); %#ok<*NASGU>
                    glob_sig_all_aligned = glob_sig_all((abs(round(timeDiff_regress))):end);
                    glob_sig_good_aligned = glob_sig_good((abs(round(timeDiff_regress))):end);
                    noreref_aligned = noreref(:,(abs(round(timeDiff_regress))):end);
                elseif timeDiff_regress > 0 % if blackrock came first, add NaNs - original method
                    processed_aligned = cat(2,NaN(size(processed,1),abs(round(timeDiff_regress))),processed); %#ok<*NASGU>
                    glob_sig_all_aligned = cat(1,NaN(abs(round(timeDiff_regress)),1),glob_sig_all); %#ok<*NASGU>
                    glob_sig_good_aligned = cat(1,NaN(abs(round(timeDiff_regress)),1),glob_sig_good); %#ok<*NASGU>
                    noreref_aligned = cat(2,NaN(size(noreref,1),abs(round(timeDiff_regress))),noreref); %#ok<*NASGU>
                elseif timeDiff_regress==0
                    processed_aligned = processed;
                    glob_sig_all_aligned = glob_sig_all;
                    glob_sig_good_aligned = glob_sig_good;
                    noreref_aligned = noreref;
                end
                
                
                %% write data out into the subject directory
                
                microDir = fullfile(ecogDir,'/',subject{subj},'/micro/');
                if ~isfolder(microDir); mkdir(microDir); end
                if ~isfolder(fullfile(microDir,spikeSessions(closest(1)).name)); mkdir(fullfile(microDir,spikeSessions(closest(1)).name)); end
                
                noreref_utah_writeOut = fullfile(microDir,spikeSessions(closest(1)).name,'noreref_aligned.mat');
                processed_utah_writeOut = fullfile(microDir,spikeSessions(closest(1)).name,'processed_aligned.mat');
                
                if exist(fullfile(spikeDir,spikeSessions(closest(1)).name,'raw/',sprintf('%s%s',spikeSessions(closest(1)).name,'_spikeInfo.mat')))
                    spikeInfo_writeOut = fullfile(microDir,spikeSessions(closest(1)).name,'spikeInfo_aligned.mat');
                    save(spikeInfo_writeOut,'applyShiftCorrection','offsetCalculatedInMs','alignedTo','alignmentChan','extractInfoStr','filter_string','metrics','pulses','sessDurSec','sessStr','startTime_mstime','startTime_raw','startTime_str','timeStamp_aligned','uniqueUnitID','-v7.3');
                end
                
                % make edits to save all new variables out 
                lfp_aligned = noreref_aligned; clear noreref_aligned;
                save(noreref_utah_writeOut,'applyShiftCorrection','offsetCalculatedInMs','chan_names','channel_ranges','glob_sig_all_aligned','glob_sig_good_aligned','lfp_aligned','pulses','samplingFreq','-v7.3');
                lfp_aligned = processed_aligned; clear processed_aligned;
                save(processed_utah_writeOut,'applyShiftCorrection','offsetCalculatedInMs','chan_names','channel_ranges','glob_sig_all_aligned','glob_sig_good_aligned','lfp_aligned','pulses','samplingFreq','-v7.3');
                
                
                
            end
            
            
        end
    end
    
    
    
    
    
end


end
