function [alignedFiles] = sam_alignSpikesWithEcog_ManualORAuto(subject,ecogDir,forceAlignAll,varargin)
% Loops through subjects of interest, and through all microphysiology
% data recorded (whether manual or automatic sorts)
% aligns pulses from the two systems (using align_nsps), and outputs aligned versions
% of spike and LFP data such that data is aligned with ecog and thus
% events.eegoffset can be used to index this data
%
% forceAlignAll = true/false; Should be kepy on false, which prevents re-aligning of data
%
% before running this code, put the mountainSort output folders
% into NIHxxx/micro/mountainsort/. if there are files from the manual pipeline,
% put them into NIHxxx/micro/manualsort/
%p
% code writes out new noreref, processed, and spikeinfo into subject/micro
% folder, as well as description of alignment (alignmentToSpikes.csv),
% which is also the output variable "alignedFiles"
%
% ecog/behavioral dir should be local copy, pre-FRNU upload, so that new
% events won't be overwritten
%
% NOTE: if first time processing a subject, populate the subject_pulseChannels
% struct with appropriate info

% created by melkalliny, 1/2019
%
% tweaks needed
% limit code to run only on sessions that havent been aligned already
%
%
%
% 6/23/2020  SJ: changed readtable to readtableSafe & huge overhaul
%               compared to original
% 7/2020     SJ: Added in check to see if micro pulses were a different
%               length than the lfp data. If so, pad the end with zeros. This happens
%               with some of the din_ts from older patients because the ts stops after
%               the last pulse.
% 8/2020     SJ: Added in special cases for zeroNeeded. tested with NIH047 first
% 10/2/2020  SJ: fixing noise channel list to include full name before brackets
% 10/30/2020 SJ: Adding ability to align multiple ecog sessions to 1 micro session -> ie 2nd alignment performed with _SPLIT2 micro session
% 04/01/2021 SJ: A ton of changes, TRC, etc

inp_pars = inputParser;
defaulttaskType = {};
defaultonlyTRC = false;

addParameter(inp_pars,'taskType',defaulttaskType,@iscell);
addParameter(inp_pars,'onlyTRC',defaultonlyTRC,@islogical);

parse(inp_pars,varargin{:})
taskType = inp_pars.Results.taskType;
onlyTRC = inp_pars.Results.onlyTRC;

debug_noWrite = 0;
checkPulses = 1;
%forceAlignAll = true;

up_micro = 'micro';
up_manualsort = 'manualsort';
up_mountain = 'mountainsort'; % IN FUTURE, CHANGE THIS TO 'mountain_sorts'
up_spk = 'spikes';
up_lfp_noreref = 'lfp_noreref';
up_lfp_reref = 'lfp_reref';

% each row: subject, pulse channel string for DC09, pulse channel string
% for DC12, suffix of file that mountainSort outputted (if processing mountainSort
% outputs), type of pulse ('digital' or 'analog' - if analog is available, use
% analog), and empty string
subject_pulseChannels = {...
    'NIH029', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  ''; % 2/24/2020 SJ: added '_ts' to all (incorporated in makePulsesStruct.m on 6/17/19)
    'NIH030', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH034', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH036', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH037', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH039', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH042', {'din1_ts'},    {'din4_ts'},  '_spikeInfo.mat',   {'digital'},  '';
    'NIH046', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH047', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH048', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH049', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH050', {'ainp1_ts','din1_ts'},   {'ainp4_ts','din4_ts'}, '_spikeInfo.mat',   {'analog','digital'},   ''; % Need to use din1_ts for SOME cases (I think only some utahs?)
    'NIH051', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   ''; % SJ added this line because 51 was missing before... not sure why
    'NIH052', {'ainp1_ts'},   {'ainp4_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH053', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH054', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {''},         '';
    'NIH057', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH059', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH060', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {'analog'},   '';
    'NIH062', {'ain17_ts'},   {'ain20_ts'}, '_spikeInfo.mat',   {'analog'},   ''; %- see note below about potential issue with NIH062
    'NIH063', {'ain1_ts'},    {'ain4_ts'},  '_spikeInfo.mat',   {'analog'},   '';
    'NIH064', {'ain1_ts'},    {'ain4_ts'},  '_spikeInfo.mat',   {'analog'},   '';
    'NIH065', {'ain1'},       {''},         '_spikeInfo.mat',   {''},         '';
    'NIH066', {'ain1_ts'},    {'ain4_ts'},  '_spikeInfo.mat',   {'analog'},   '';
    'NIH069', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH071', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH072', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';  %- JW added... just copied from 71.  is it that simple?
    'NIH074', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';  %-
    'NIH076', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH070', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';  %- SJ added... Just copied from 76...
    'NIH079', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH081', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH082', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH086', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   '';
    'NIH088', {'ain1_ts'},    {'ain4_ts'},  '_noreref.mat',     {'analog'},   ''}; 

% Special Cases for alignment, forced Zero-ing applied
              % Subject, ecog_sess, micro_sess, {ecog start, ecog end, micro start, micro end}

zf_subjects =   {'NIH047',           'NIH046',           'NIH052'};
zf_ecog_sess =  {'170311_1312',      '161214_1122',      '170830_1700'};
zf_micro_sess = {'170311_1259',      '161214_1121_beh',  '170830_1625_beh'};
zf_num =        {[0, 0, 1000000, 0], [226000, 0, 0, 0],  [0, 0, 2400000, 0]};

zero_forced = table(zf_subjects', zf_ecog_sess', zf_micro_sess', zf_num', 'VariableNames',{'Subject','ecog_sess','micro_sess','zero'});

% previously, line for NIH062 was: 'NIH062', 'ain1','ain17','_spikeInfo.mat'
% because ain17 was a 'backup' pulse channel, which for some sessions was the only
% valid pulse channel. removing this syntax for now,
% while creating alignment that feeds into align_nsps.m..
% so, expect a bug if this code is run on 062

%manualsort_path = ['micro' filesep 'manualsort' filesep]; %-SJ
%mountainsort_path = ['micro' filesep 'mountain_sorts' filesep]; Not yet! Don't want to do mountain sort yet


if ~iscell(subject), subject = {subject}; end

for sub = 1:numel(subject)
    
    subj = subject{sub};
    
    %- initalize the output structure
    %alignedFiles = {'ecog_sess','beh_sess','spikeInfo_loc','t_offset','alignmentSuccess','eventsEncompassedByMicroFile','lengthOfMicroFileInMin'}; %-SJ replaced from below
    alignedFiles = table('Size',[0 6],'VariableNames',{'beh_sess','ecog_sess','micro_sess','alignmentSuccess','eventsEncompassedByMicroFile','lengthOfMicroFileMin'}, ...
                                      'VariableTypes',{'string',  'string',   'string',    'string',          'double',                      'double'});
    rowIterate = 1;
    
    % Do manualsort and/or mountain_sorts folders exist?
    if      exist(fullfile(ecogDir,subj,up_micro,up_mountain)) &&   exist(fullfile(ecogDir,subj,up_micro,up_manualsort))
        autoOrManual = [1 2];
    elseif  exist(fullfile(ecogDir,subj,up_micro,up_mountain)) &&  ~exist(fullfile(ecogDir,subj,up_micro,up_manualsort))
        autoOrManual = 1;
    elseif ~exist(fullfile(ecogDir,subj,up_micro,up_mountain)) &&   exist(fullfile(ecogDir,subj,up_micro,up_manualsort))
        autoOrManual = 2;
    elseif ~exist(fullfile(ecogDir,subj,up_micro,up_mountain)) && ~exist(fullfile(ecogDir,subj,up_micro,up_manualsort))
        if exist(fullfile(ecogDir,subj,'micro'),'dir')
            fprintf('\n\nfound subj/micro folder, but no mountain or manual subfolder... fix it or skip?\n\n',subj)
            keyboard
            continue
        else
            fprintf('\n\nno micro data found for %s!! skipping!\n\n',subj)
            continue
        end
    end

    %- find all behavioral task folders
    behavDir = fullfile(ecogDir,subj,'behavioral');
    behavFolders = dir(behavDir);
    behavFolders = behavFolders(find([behavFolders.isdir]==1));
    behavFolders(find(contains({behavFolders.name},'.')))  = [];  %- clean up the behavioral folders returned by dir
    behavFolders(find(strncmp({behavFolders.name},'_',1))) = [];
    behavFolders(find(contains({behavFolders.name},'cant'))) = [];
    
    %- Calculate total nubmer of expected sessions (-SJ remove/improve for speed?)
    fullSessCount = 0;
    for exp = 1:length(behavFolders)
        sessionDir     = fullfile(behavDir,behavFolders(exp).name,'/');
        sessionFolders = dir(sessionDir);
        sessionFolders = sessionFolders(find(strncmp({sessionFolders.name},'session',7))); %- rule is the "cant" string preceeds "session" if bad alignment
        fullSessCount  = fullSessCount+length(sessionFolders);
    end
    
    %SJ: rewrite if meant to iterate over both
    if autoOrManual ==1
        microDir = fullfile(ecogDir,subj,up_micro,up_mountain);
        alignment_xlsx_writeOut = fullfile(microDir,'alignmentToSpikes_mountain.xlsx'); %Name for writing out later
    elseif autoOrManual ==2
        microDir = fullfile(ecogDir,subj,up_micro,up_manualsort);
        alignment_xlsx_writeOut = fullfile(microDir,'alignmentToSpikes_manual.xlsx');
    end
    
    csv_alignment_writeOut = strrep(alignment_xlsx_writeOut,'.xlsx','.csv');
    if exist(csv_alignment_writeOut)
        fprintf('%s\n',['ERROR!!! ' csv_alignment_writeOut ' exists!!! There should only be .xlsx files now. Continue to delete.']);
        keyboard
        delete(csv_alignment_writeOut)
    end
    
    if exist(alignment_xlsx_writeOut)
        alignedFiles_old = readtableSafe(alignment_xlsx_writeOut); %SJ from readtable
    end
    
    %- full loop over the behavioral tasks and sessions within
    thisSessCount = 0;
    for exp = 1:length(behavFolders) %SJ CHANGE BACK TO 1: (3 for NIH060)
        sessionDir = fullfile(behavDir,behavFolders(exp).name,'/');
        sessionFolders = dir(sessionDir);
        sessionFolders = sessionFolders(find(strncmp({sessionFolders.name},'session',7))); %- rule is the "cant" string preceeds "session" if bad alignment
        
        %%ADD HERE%%
        %keyboard % Sort sessions in order (use natsort) SJ
        sessionFolders_cell = {sessionFolders.name};
        [~,ndx,dbg] = natsort(sessionFolders_cell);
        sessionFolders = sessionFolders(ndx);
        
        %- loop over viable session folders
        for sess = 1:length(sessionFolders)
            
            thisSessCount = thisSessCount+1;
            fprintf('\n\n >>> ATTEMPTING TO ALIGN: %s/%s/%s  (session %d of %d total sessions) <<<\n', subj,behavFolders(exp).name,sessionFolders(sess).name, thisSessCount,fullSessCount);
            
            pulses_ecog = NaN; pulses_ecog_dc12 = NaN;
            
            % regular branch
            whichDC = 9;
            beh_sess_path = [sessionDir,sessionFolders(sess).name];
            if ~exist(fullfile(beh_sess_path,'events.mat'))
                fprintf('%s\n',['events.mat does not exist in directory: ' beh_sess_path '! Continue only if there is a reason for this.']);
                if contains(beh_sess_path,'attentionFeedback')
                    continue
                end
                keyboard
                continue
            end
            sess_events = load(fullfile(beh_sess_path,'events.mat'));
            beh_sess = char(regexp(beh_sess_path,[subj '.*behavioral.*'],'match'));
            sess_events = sess_events.events;
            if ~isfield(sess_events,'eegfile')
                if ~contains(beh_sess,'CLIN')
                    fprintf('%s\n','Events does not have eegfile as a field, but this is not a CLINICAL stim map session. WHY??');
                    keyboard
                end
                keyboard
                continue
            end
            ecog_sess = sess_events(1).eegfile; %previously noreref_dir
            dashes = strfind(ecog_sess,'/');
            if isempty(dashes), fprintf('\n uh oh... didnt find a backslash in the eegfile: %s',ecog_sess); keyboard; end
            ecog_sess = ecog_sess(dashes(end)+1:end);
            
            if ~exist(fullfile(ecogDir,subj,'eeg.noreref',ecog_sess))
                fprintf('\nnoreref directory doesnt contain all files referenced by behavioral folders\n\n')
                keyboard
                continue
            end
            
            alignedFiles.beh_sess(rowIterate) = beh_sess;
            alignedFiles.ecog_sess(rowIterate) = ecog_sess;

            
            
            %                 for digital pulses
            %                 trigDCDir = fullfile(ecogDir,subject{subj},'/eeg.noreref/',noreref_dir,'trigDC09.sync.txt');
            %                 fileID = fopen(trigDCDir,'r');
            %                 pulses_ecog = textscan(fileID,'%d'); pulses_ecog = double(pulses_ecog{1});
            %                 ecog_start = datenum(noreref_dir,'yymmdd_HHMM');
            %                 fclose(fileID);
            %
            %                 trigDC12Dir = fullfile(ecogDir,subject{subj},'/eeg.noreref/',noreref_dir,'trigDC12.syncBR.txt');
            %                 fileID = fopen(trigDC12Dir,'r');
            %                 pulses_ecog_dc12 = textscan(fileID,'%d'); pulses_ecog_dc12 = double(pulses_ecog_dc12{1});
            %                 fclose(fileID);
            
            % for analog pulses
            trigDCDir = fullfile(ecogDir,subj,'eeg.noreref',ecog_sess,'DC09');
            [fchan,msg] = fopen(trigDCDir,'r+','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', trigDCDir, msg);
            pulses_ecog = fread(fchan,inf,'int16');
            fclose(fchan);
            
            trigDC12Dir = fullfile(ecogDir,subj,'eeg.noreref',ecog_sess,'DC12');
            fileID = fopen(trigDC12Dir,'r'); %SJ - why do we have this here?
            [fchan,msg] = fopen(trigDC12Dir,'r+','l');
            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', trigDC12Dir, msg);
            pulses_ecog_dc12 = fread(fchan,inf,'int16');
            fclose(fchan);
            
            paramsDir = fullfile(ecogDir,subj,'eeg.noreref',ecog_sess,'params.txt');
            fileID = fopen(paramsDir,'r');
            samprate_ecog = textscan(fileID,'%s'); samprate_ecog = double(str2num(samprate_ecog{1,1}{2,1}));
            fclose(fileID);
            
            
            %%ADD HERE%%
            % Get the first and last times of the task (eegoffset 1 and
            % end) SJ
            % Will need to add to ecog_start? is eegoffset a datetime or
            % time relative to beginning of ecog?
            %keyboard
%             eegoffset_start = sess_events(1).eegoffset; %convert samples to day
%             eegoffset_start_day = days(seconds(eegoffset_start/samprate_ecog));
%             eegoffset_end_day = days(seconds(sess_events(end).eegoffset/samprate_ecog));
%             beh_start = ecog_start + eegoffset_start_day;
%             beh_start = ecog_start + days(seconds(sess_events(1).eegoffset/samprate_ecog));
%             beh_end = 1;
            beh_start = epoch2date(sess_events(1).mstime);
            beh_end = epoch2date(sess_events(end).mstime);
            %keyboard %Do we want to do the line below?
            %beh_end = datenum(dateshift(datetime(beh_end,'ConvertFrom','datenum'),'end','minute'));

            %%ADD HERE%%
            %keyboard % We want a # so we can compare > < 
            % Verify this works with CRV too!
            ecog_start = datenum(ecog_sess,'yymmdd_HHMM');
            if numel(pulses_ecog) ~= numel(pulses_ecog_dc12)
                fprintf('%s\n','ERROR!!!! DC09 and DC12 are different lengths!!! Only calculating the duration for DC09!');
                keyboard
            end
            duration_ecog = days(seconds(numel(pulses_ecog)/samprate_ecog)); %in days
            ecog_end = ecog_start + duration_ecog; 
            ecog_end = datenum(dateshift(datetime(ecog_end,'ConvertFrom','datenum'),'end','minute'));% We want to round up to the minute here because we don't know which second of the start minute it started (just assuming 0 sec in)
            %%ADD HERE%%
            % Get the end time of ecog (length of DC09 but take into
            % account sampling rate...?)
            
            
            %% now, cycle through both plexon and mountainSort files
            %   if ecog and micro both use blackrock: then look for
            % double-NSP transform info in raw folder, pull it, and apply it to second
            % file
            
            %   if ecog and micro use different systems, or micro was on second NSP but second
            % NSP contained no ecog: pull all the pulses available and compare the
            % ecog and micro pulses present
            %      - if analog of same kind (dc9 or
            %       12) is present, then feed that into align_nsps and use
            %       transforms info to edit the raw data inside this fx
            %
            %      - if analog from both is not present, then use a version of align_nsps
            %       that is suitable for digital pulses, get transforms, and
            %       edit raw data

            for pipelineType = autoOrManual % not forcing for right now

                
                % pull the start dates from all the spikeInfos looks in the subject/micro folder to
                % find the matching spikeSession    
                microSessions = dir(microDir);
                microSessions(find(contains({microSessions.name},'.'))) = [];
                microSessions(find(contains({microSessions.name},'cant'))) = []; % exclude cant align sessions
                microSessions(find(contains({microSessions.name},'_SPLIT2'))) = []; % SJ: Exclude _SPLIT2 sessions to begin with because if these were already created, we don't want the first ecog session to try to align to both!
                if strcmp(microSessions(1).name(1:3),'NIH')
                    underscores = strfind(microSessions(1).name,'_');
                    dateChar =  underscores(1)+1:underscores(2)+4;
                else
                    dateChar=1:11;
                end
                
                % lets get dates of possible matches, and also grab forced alignment pairs if they exist
                microSessionDates = NaN(1,length(microSessions));
                forced_alignment_pairs = cell(1,length(microSessions));
                for session = 1:length(microSessions)
                    temp= microSessions(session).name(dateChar);
                    microSessionDates(1,session) = datenum(temp,'yymmdd_HHMM');
                    
                    grab_forced = fullfile(microSessions(session).folder,microSessions(session).name,'force_alignment.txt');
                    if exist(grab_forced)
                        %Do we ever get here?? SJ: haven't checked this
                        fid = fopen(grab_forced);
                        tline = fgetl(fid);
                        tlines = cell(0,1);
                        while ischar(tline)
                            tlines{end+1,1} = tline;
                            tline = fgetl(fid);
                        end
                        fclose(fid);
                        forced_alignment_pairs{1,session} = tlines{2,1};
                        clear tlines
                    end
                end
                
                find_forced_pair = find(strcmp(forced_alignment_pairs,ecog_sess));
                if size(find_forced_pair,2) ~= 0
                    keyboard %Do we ever get here?? SJ: haven't checked this
                    closest = find_forced_pair;
                    diffs  = (ecog_start-microSessionDates(closest));
                else
                    % get session closest to start time
                    % SJ: actually, we want to grab any micro sessions that
                    % begin or end within the duration of the ecog file
                    %diffs  = (ecog_start-microSessionDates);
                    %[~,closest] = find(abs(diffs) == min(abs(diffs)));
                    
                    % Let's first grab any micro sessions that start during
                    % the ecog file (don't count one if it starts right at
                    % the end)
                    startWithinEcog = find(microSessionDates >= ecog_start & microSessionDates <= ecog_end);
                    fprintf('%s\n',['Note: There are ' num2str(numel(startWithinEcog)) ' micro sessions that begin within the ecog file.']);
                    %keyboard
                    microIdxInclude = [];
                    % SJ: Changed the following line from >1
                    if numel(startWithinEcog) > 0 % We need to check if we should include more than one
                        if numel(startWithinEcog) > 1
                            %keyboard %We have mostly checked this, but want to see what happens maybe more
                        end
                        if numel(startWithinEcog) > 2
                            keyboard %Do we ever get here?? Just make sure it works properly
                        end
                        
                        %SJ changes with >1 line because we don't want ANY micro file that starts within 5 min of the end... even if only 1 was found.
                        for ww = 1:numel(startWithinEcog)
                            if (ecog_end - microSessionDates(startWithinEcog(ww))) > .0035 % If the micro session starts before 5 min of the eegfile ending
                                microIdxInclude = cat(2, microIdxInclude,startWithinEcog(ww));
                            else
                                if numel(startWithinEcog) == 1
                                    %keyboard %Getting rid of only micro within ecog, make sure we grab another when we reach back an hour!
                                end
                                fprintf('%s\n','Getting rid of 1 micro file because it starts within 5 min of ecog file end.');
                                %keyboard
                            end
                        end
                        
%                     else
%                         %microIdxInclude = microSessionDates(startWithinEcog);
%                         microIdxInclude = startWithinEcog;
                    end
                    %Now, if any of our current micro files don't start
                    %within 5 min of the ecog start time, we need to search
                    %if there are any others that do!
                    if any(abs(microSessionDates(microIdxInclude) - ecog_start) > .0035) || isempty(microIdxInclude)
                        %keyboard % Not tested yet
                        % Time to look 1 hour back (.0417 days)
                        startHrPrior = find(microSessionDates >= (ecog_start-.0417) & microSessionDates < ecog_start);
                        if ~isempty(microIdxInclude) && any(ismember(microIdxInclude, startHrPrior)) %SJ: add isempty check and change to ismember (== fails unless same length)
                            fprintf('%s\n','Error!! Already have this session included! What is up?');
                            keyboard
                        end
                        if numel(startHrPrior) > 0
                            %Need to see if it ends after ecog begins,
                            %otherwise we do not care
                            for ww2 = 1:numel(startHrPrior)
                                % Load in the LFP.mat
                                lfp_check = load(fullfile(microDir,microSessions(startHrPrior(ww2)).name,up_lfp_noreref,'lfp.mat'));
                                %samprate_check = lfp_check.lfpStruct.samplingFreq;
                                microdur_check = days(minutes(lfp_check.lfpStruct.sessDurMin));
                                micro_end_check = microSessionDates(startHrPrior(ww2)) + microdur_check;
                                %Check if the micro file ends at least 5
                                %min into the ecog file
                                if (micro_end_check - ecog_start) > .0035
                                    % You want to include this
                                    %keyboard % Check this
                                    microIdxInclude = cat(2, microIdxInclude,startHrPrior(ww2));
                                end
                            end
                        end
                    end
                    
                    if numel(microIdxInclude) == 0
                        fprintf(sprintf('No matching micro data within the ecog session or 1 hour prior %s\n',fullfile(ecogDir,subj,'eeg.noreref',ecog_sess)))
                        fprintf('\n%s\n','Take a second to verify there are no micro sessions available.');
                        %keyboard %SJ Check
                        %alignedFiles{rowIterate,4} = sprintf('alignment not attempted');
                        %alignedFiles.beh_sess(rowIterate) = beh_sess;
                        %alignedFiles.ecog_sess(rowIterate) = ecog_sess;
                        alignedFiles.micro_sess(rowIterate) = 'no match';
                        alignedFiles.alignmentSuccess(rowIterate) = 'alignment not attempted';
                        alignedFiles.eventsEncompassedByMicroFile(rowIterate) = 0;
                        alignedFiles.lengthOfMicroFileMin(rowIterate) = NaN;
                        %alignedFiles{rowIterate,5} = sprintf('alignment not attempted');
                        %alignedFiles{rowIterate,7} = NaN; % print how long micro file was
                        %alignedFiles{rowIterate,6} = 0; % print how long micro file was
                        %alignedFiles{rowIterate,1} = ecog_sess;
                        %alignedFiles{rowIterate,3} = 'no match';
                        
                        rowIterate = rowIterate+1;
                        continue
                    elseif numel(microIdxInclude) > 1
                        % Need to iterate through loop now
                        %keyboard %Need to check
                        
                    elseif numel(microIdxInclude) == 1
                        %closest = microIdxInclude;
                    else
                        keyboard %What? How is this possible?
                    end                         
                    
                    
                    % Then see if there are any micro sessions starting
                    % within 1 hour before the ecog start and check its
                    % duration -> end time. If the end time occurs after
                    % the ecog_start, then include it (if not right at
                    % beginning)

                    % instead of applying a different threshold for stimMapGui sessions,
                    % lets check to see if the events are encompassed by this
                    % files length
                    % try to find a matching file. if not, go to next session
                    %if min(abs(diffs)) > 0.04 % 4 percent of 1 day is roughly 1 hour

                end
                for cc = 1:numel(microIdxInclude)
                    if cc > 1
                        %keyboard %Time to check this - multiple micro files aligning!
                        alignedFiles.beh_sess(rowIterate) = beh_sess; %SJ: need to add both of these in for consecutive sessions since this usually gets done way above
                        alignedFiles.ecog_sess(rowIterate) = ecog_sess;
                    end
                    closest = microIdxInclude(cc);
                    
                    %micro_sess = microSessions(closest(1)).name;
                    micro_sess = microSessions(closest).name;
                    micro_start = datenum(micro_sess,'yymmdd_HHMM');
                    
                    % check to see if we've already aligned this micro
                    % before - only skip if all beh, ecog, AND micro match
                    %continue_out_of_numPulsesToAlign = false;
                    %triedToAlignBefore = find(strcmp(alignedFiles.micro_sess(1:rowIterate-1),micro_sess));
                    % Maybe here we could do a check for if a micro file is
                    % about to be aligned again to a different ecog file!
                    triedToAlignBefore = find(strcmp(beh_sess,alignedFiles.beh_sess(1:rowIterate-1)) & strcmp(ecog_sess,alignedFiles.ecog_sess(1:rowIterate-1)) & strcmp(micro_sess,alignedFiles.micro_sess(1:rowIterate-1)));
                    triedToAlignBefore_micro_and_ecog = find(strcmp(ecog_sess,alignedFiles.ecog_sess(1:rowIterate-1)) & strcmp(micro_sess,alignedFiles.micro_sess(1:rowIterate-1)));
                    
                    % Micro session is the same, but beh and/or ecog are not:
                    %if ~isempty(find((~strcmp(beh_sess,alignedFiles.beh_sess(1:rowIterate-1)) | ~strcmp(ecog_sess,alignedFiles.ecog_sess(1:rowIterate-1))) & strcmp(micro_sess,alignedFiles.micro_sess(1:rowIterate-1))))
                    if ~isempty(find(~strcmp(ecog_sess,alignedFiles.ecog_sess(1:rowIterate-1)) & strcmp(micro_sess,alignedFiles.micro_sess(1:rowIterate-1))))
                        %STOP!!!!! YOU NEED TO BUILD THIS IN!!!!!!
                        fprintf('\n%s\n',['ATTENTION! 2nd ecog session (' ecog_sess ') aligning to 1 micro session (' micro_sess '). Creating _SPLIT2 micro session and copying over lfp & spike dirs, but not alignment/aligned files, and re-aligning.']);
                        keyboard
                        %!!!!!!!!!!!!!!!!!!!!!!!!
                        micro_sess_prev = micro_sess;
                        micro_sess = [micro_sess '_SPLIT2'];
                        triedToAlignBefore_micro_and_ecog = find(strcmp(ecog_sess,alignedFiles.ecog_sess(1:rowIterate-1)) & strcmp(micro_sess,alignedFiles.micro_sess(1:rowIterate-1)));
                        if isempty(triedToAlignBefore_micro_and_ecog) % Has this ecog session already been aligned to a _SPLIT2 session? If so, continue and go through that pathway instead
                        
                            % Now copy all the non-alignment stuff over
                            if ~exist(fullfile(microDir, micro_sess),'dir') %First time discovering this split
                                mkdir(fullfile(microDir, micro_sess));
                            end
                            fprintf('%s\n','Copying lfp and spike dirs to _SPLIT2');
                            if ~exist(fullfile(microDir,micro_sess_prev,up_lfp_noreref),'dir')
                                keyboard % lfp dir does not exist to be copied to _SPLIT2??
                            end
                            st1 = 100;
                            st2 = 100;
                            if ~exist(fullfile(microDir,micro_sess,up_lfp_noreref),'dir')
                                st1 = copyfile(fullfile(microDir,micro_sess_prev,up_lfp_noreref),fullfile(microDir,micro_sess,up_lfp_noreref));
                            end
                            if ~exist(fullfile(microDir,micro_sess_prev,up_spk),'dir')
                                keyboard % spike dir does not exist to be copied to _SPLIT2??
                            end
                            if ~exist(fullfile(microDir,micro_sess,up_spk),'dir')
                                st2 = copyfile(fullfile(microDir,micro_sess_prev,up_spk),fullfile(microDir,micro_sess,up_spk));
                            end
                            if st1 ~= 1 || st2 ~= 1
                                keyboard %Error copying either lfp (st1) or spike (st2) folder!
                            end
                            %Check SPLIT2 lfp folder and remove lfp_aligned
                            if ~exist(fullfile(microDir,micro_sess,up_lfp_noreref),'dir') 
                                keyboard %Why does this not exist when we just tried to copy it?
                            end
                            lfp_aligned_file = getDirNamesRegexp(fullfile(microDir,micro_sess,up_lfp_noreref),'lfp.*aligned.*.mat');
                            if isempty(lfp_aligned_file)
                                keyboard %How can this be empty if this is the second alignment to this micro sess? Is it named something else? Even if the 1st alignment wasn't successful it still writes it out
                            elseif numel(lfp_aligned_file) > 1
                                keyboard %Why is it picking up more than one aligned lfp .mat?
                            end
                            if st1 ~= 100 %Deleting the aligned lfp only if it comes from the non-split session, not if the split2 session already exists
                                delete(fullfile(microDir,micro_sess,up_lfp_noreref,char(lfp_aligned_file)));
                            end
                            % Now do the same for spike
                            if ~exist(fullfile(microDir,micro_sess,up_spk),'dir') 
                                keyboard %Why does this not exist when we just tried to copy it?
                            end
                            spk_aligned_file = getDirNamesRegexp(fullfile(microDir,micro_sess,up_spk),'.*spikeInfo.*aligned.*.mat');
                            if isempty(spk_aligned_file)
                                % This is okay, just hasn't been sorted yet - just verify that the note is there
                                if numel(getDirNamesRegexp(fullfile(microDir,micro_sess,up_spk),'SPK.*_unsorted.txt')) ~= 1 
                                    keyboard % This is saying that there is no aligned spikeInfo AND it isn't unsorted (where is it since this is the second ecog alignment to the same micro session!)
                                end
                            elseif numel(spk_aligned_file) > 1
                                keyboard %Why is it picking up more than one aligned spikeInfo .mat?
                            elseif st2 ~= 100 %Deleting the aligned spikeInfo only if it comes from the non-split session, not if the split2 session already exists
                                delete(fullfile(microDir,micro_sess,up_spk,char(spk_aligned_file)));
                            end
                            % Okay, so now we have generated our second micro session (_SPLIT2), now we just need to make sure it will work properly
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %SJ: inserting this for TRC sessions misaligned: CVCH
                    %thisIsTRCsess = false;
                    if onlyTRC % Don't bother checking if we're not only doing TRC
                        % is the raw file a blackrock or non-blackrock file?
                        rawEcogDir_TRC = fullfile(ecogDir,subj,'raw');
                        %rawSessions = getDirNamesRegexp(rawEcogDir,'^\d{6}_\d{4}.*');
                        [rawSessions_TRC, rawSessionsisStim_TRC] = getRawSessions(rawEcogDir_TRC);
                        if rawSessionsisStim_TRC(find(strcmp(rawSessions_TRC,ecog_sess))) % STIM_MAP folder
                            rawSession_content_TRC = getDirNamesRegexp(fullfile(rawEcogDir_TRC,'STIM_MAP',ecog_sess),'.*');
                        else
                            rawSession_content_TRC = getDirNamesRegexp(fullfile(rawEcogDir_TRC,ecog_sess),'.*'); % SJ fixed
                        end
                        if any(contains(rawSession_content_TRC,'.TRC')) % SJ fixed
                            [T] = documentCRVsess_micro(subj,ecog_sess,micro_sess,samprate_ecog,'/Users/wittigj/Desktop/TRCsess.xlsx');
                            %thisIsTRCsess = true;
%                         else
%                             %Not a TRC session... continue
%                             continue
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if ~isempty(triedToAlignBefore)
                        keyboard % Check - where does this go after the continue anyway? ** This is is beh session is the same- no idea why/how that would even happen
                        resultsFromBefore = alignedFiles.alignmentSuccess(triedToAlignBefore);
                        % if this pair has already been aligned, dont bother aligning again
                        alignedFiles.micro_sess(rowIterate) = micro_sess;
                        alignedFiles.alignmentSuccess(rowIterate) = alignedFiles.alignmentSuccess(triedToAlignBefore);
                        alignedFiles.eventsEncompassedByMicroFile(rowIterate) = alignedFiles.eventsEncompassedByMicroFile(triedToAlignBefore);
                        alignedFiles.lengthOfMicroFileMin(rowIterate) = alignedFiles.lengthOfMicroFileMin(triedToAlignBefore);
                        
                        rowIterate = rowIterate+1;
                        continue
                    end
                    
                    if ~isempty(triedToAlignBefore_micro_and_ecog)
                        %keyboard % This is the new check I put in place- why the heck does the beh session need to be the same??
                        resultsFromBefore = alignedFiles.alignmentSuccess(triedToAlignBefore_micro_and_ecog(1));%1st to get rid of multiple
                        % if this pair has already been aligned, dont bother aligning again
                        alignedFiles.micro_sess(rowIterate) = micro_sess;
                        alignedFiles.alignmentSuccess(rowIterate) = alignedFiles.alignmentSuccess(triedToAlignBefore_micro_and_ecog(1)); %1st to get rid of multiple
                        % Need to find a better way to calculate this...
%                         perc_encompassed = zeros(1,size(lfpPulses_full,1));
%                         perc_encompassed(numPulsesToAlign) = sum(events_times>=micro_start & events_times<= micro_end)/numel(events_times);
%                         
                        alignedFiles.eventsEncompassedByMicroFile(rowIterate) = NaN; %Do this later
                        %keyboard
                        
                        alignedFiles.lengthOfMicroFileMin(rowIterate) = alignedFiles.lengthOfMicroFileMin(triedToAlignBefore_micro_and_ecog(1));%1st to get rid of multiple

                        rowIterate = rowIterate+1;
                        continue
                    end
    
                    skipLFPalign = false;
                    skipSPKalign = false;

                    if ~forceAlignAll % Just skip this check if you want to force alignment
                        SPKmat_file = char(getDirNamesRegexp(fullfile(microDir, micro_sess, up_spk),'.*spikeInfo.mat'));
                        aligned_SPKmat_file = char(getDirNamesRegexp(fullfile(microDir, micro_sess, up_spk),'.*spikeInfo_aligned.*\.mat'));
                        LFPmat_file = char(getDirNamesRegexp(fullfile(microDir, micro_sess, up_lfp_noreref),'.*lfp.mat'));
                        aligned_LFPmat_file = char(getDirNamesRegexp(fullfile(microDir, micro_sess, up_lfp_noreref),'.*lfp_aligned.*\.mat'));

                        SPKmat = fullfile(microDir, micro_sess, up_spk,SPKmat_file);
                        aligned_SPKmat = fullfile(microDir, micro_sess, up_spk,aligned_SPKmat_file);
                        LFPmat = fullfile(microDir, micro_sess, up_lfp_noreref,LFPmat_file);
                        aligned_LFPmat = fullfile(microDir, micro_sess, up_lfp_noreref,aligned_LFPmat_file);

                        if size(SPKmat_file,1) > 1 || size(aligned_SPKmat_file,1) > 1 || size(LFPmat_file,1) > 1 || size(aligned_LFPmat_file,1) > 1
                            fprintf('%s\n','ERROR!!! There are more than one of SPK/LFP/aligned versions!!!');
                            keyboard
                        end

                        if ~isempty(LFPmat_file) && exist(LFPmat,'file')
                            if ~isempty(aligned_LFPmat_file) && exist(aligned_LFPmat,'file')
                                % lfp.mat exists AND aligned version exists!! Check to see if you need to
                                % replace the aligned version
                                %Check if the lfp.mat has been re-generated since the aligned version was written
                                if strcmp(newerFile(LFPmat,aligned_LFPmat),aligned_LFPmat) % aligned version is newer, don't need to align anything then
                                    %fprintf('%s\n',[aligned_LFPmat ' already exists and is newer than lfp.mat. No need to align.']);
                                    skipLFPalign = true;
                                end
                                % else case means aligned version doesn't exist yet - obvi we want it to
                                % align now - don't stop it!
                            end
                            % else case means there's no lfp?? how? This should catch later on
                        end

                        if ~isempty(SPKmat_file) && exist(SPKmat,'file')
                            if ~isempty(aligned_SPKmat_file) && exist(aligned_SPKmat,'file')
                                % spikeInfo exists AND aligned version exists!! Check to see if you need to
                                % replace the aligned version
                                %Check if the spikeInfo.mat has been re-generated since the aligned version was written
                                if strcmp(newerFile(SPKmat,aligned_SPKmat),aligned_SPKmat) % aligned version is newer, don't need to align anything then
                                    %fprintf('%s\n',[aligned_SPKmat ' already exists and is newer than spikeInfo.mat. No need to align.']);
                                    skipSPKalign = true;
                                end
                                % else case means aligned version doesn't exist yet - obvi we want it to
                                % align now - don't stop it!
                            end
                        else % This means this hasn't been sorted. This means we can skip the alignment!
                            skipSPKalign = true;
                        end
                        if skipSPKalign && skipLFPalign
                            fprintf('%s\n','Aligned spikeInfo and lfp already exist and are newer than un-aligned. No need to align either. Skipping to next beh session.');
                            row_in_old_align = find(strcmp(beh_sess,alignedFiles_old.beh_sess) & strcmp(ecog_sess,alignedFiles_old.ecog_sess) & strcmp(micro_sess,alignedFiles_old.micro_sess));

                            if isempty(row_in_old_align)
                                fprintf('%s\n','ERROR!!! Supposed to be already aligned, but cannot find a match in the alignment csv!!! Continue to re-align.');
                                keyboard
                            else
                                alignedFiles.micro_sess(rowIterate) = micro_sess;
                                alignedFiles.alignmentSuccess(rowIterate) = alignedFiles_old.alignmentSuccess(row_in_old_align);
                                alignedFiles.eventsEncompassedByMicroFile(rowIterate) = alignedFiles_old.eventsEncompassedByMicroFile(row_in_old_align);
                                alignedFiles.lengthOfMicroFileMin(rowIterate) = alignedFiles_old.lengthOfMicroFileMin(row_in_old_align);
                                rowIterate = rowIterate+1;
                                continue
                            end

                        end
                    end

                    % lets make sure the events are encompassed by this files length
                    % SJ- change end to variable
%                     keyboard
%                     if (microSessionDates(closest) - (ecog_start + (sess_events(end).eegoffset / 1000 / 60 / 60 / 24)) > 0) %SJ!! CHange!! Won't work with CRV!!
%                         % last event occurs before the spike session date
%                         fprintf(sprintf('Matching micro data, but events not overlapping with micro file %s\n',fullfile(ecogDir,subj,'eeg.noreref',ecog_sess)))
%                         %alignedFiles{rowIterate,4} = sprintf('alignment not attempted');
%                         %alignedFiles.beh_sess = beh_sess;
%                         %alignedFiles.ecog_sess(rowIterate) = ecog_sess;
%                         alignedFiles.micro_sess(rowIterate) = 'no match';
%                         alignedFiles.alignmentSuccess(rowIterate) = 'alignment not attempted';
%                         alignedFiles.eventsEncompassedByMicroFile(rowIterate) = 0; 
%                         alignedFiles.lengthOfMicroFileMin(rowIterate) = NaN; 
% 
%                         rowIterate = rowIterate+1;
%                         continue
%                     end

                    %alignedFiles{rowIterate,4} = sprintf('see NSP_alignment_summary.txt');
                    alignedFiles.alignmentSuccess(rowIterate) = sprintf('see NSP_alignment_summary.txt'); %SJ: ??
%                     if pipelineType==2
%                         % while we're here, lets look at the duration of the
%                         % spikeSession, and print a couple of things out
%                         % in the future, no need to limit this to pipelineType=2,
%                         % but for now, this jacksheet is only present in the case
%                         % of manually processed micro files
%                         try
%                             jacksheetBR = csv2cell(fullfile(microDir,micro_sess,up_lfp_noreref,'jacksheetBR_withNewChanNamesAndDevNum.csv'));
%                         catch
%                             jacksheetBR = csv2cell(fullfile(microDir,micro_sess,up_lfp_noreref,'jacksheetBR_justSplitChan.csv'));
%                         end
% 
%                         keyboard %SJ: why are we getting duration here?
%                         columnNames = jacksheetBR(1,:);
%                         durationColumn = find(strcmp(columnNames,'DurationMin'));
%                         durationFile = jacksheetBR(2,durationColumn); % just grabbing first channel. even though it could be 1st NSP and not micro data
%                         lengthInMinutes = round(str2num(durationFile{1}));
%                         alignedFiles.lengthOfMicroFileMin(rowIterate) = lengthInMinutes; % print how long micro file was
% 
%                         keyboard %SJ: fix this stuff and the stuff below it too
%                         beginningOfMicro = microSessionDates(closest);
%                         endTimeOfMicro = microSessionDates(closest) + (lengthInMinutes / 60 / 24); %SJ change!!!
% 
%                         lastEventTime = (ecog_start + (sess_events(end).eegoffset / 1000 / 60 / 60 / 24)); %SJ change!!!
%                         firstEventTime = (ecog_start + (sess_events(1).eegoffset / 1000 / 60 / 60 / 24)); %SJ change!!!
%                         if (firstEventTime>beginningOfMicro) && (lastEventTime<endTimeOfMicro)
%                             eventsEncompassed = 1;
%                         else
%                             eventsEncompassed = 0;
%                         end
%                         alignedFiles.eventsEncompassedByMicroFile(rowIterate) = eventsEncompassed; 
%                     end

                    indexSubject = find(strcmp(subj,subject_pulseChannels));

                    % if this is mountainSort block: reorganize folder, if it hanst been organized yet.
                    % this organization is being done at the mountainSort stage now
                    % keeping this block of code only for backwards compatibility

                    % HEY!! WHY ARE WE USING LFP_REREF HERE??? (Doesn't exist yet)
                    if pipelineType==1 && ~(exist(fullfile(microDir,micro_sess,'raw/'))) && ~(exist(fullfile(microDir,micro_sess,up_lfp_reref)))

                        mkdir(fullfile(microDir,micro_sess,'raw/'))
                        mkdir(fullfile(microDir,micro_sess,'cleaning/'))
                        mkdir(fullfile(microDir,micro_sess,'sorting/'))

                        if exist(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_spikeInfo.mat')))
                            movefile(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_spikeInfo.mat')),fullfile(microDir,micro_sess,'/raw/',sprintf('%s%s',micro_sess,'_spikeInfo.mat')));
                            movefile(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_sortSummary.csv')),fullfile(microDir,micro_sess,'sorting/',sprintf('%s%s',micro_sess,'_sortSummary.csv')));
                            movefile(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_spikeWaveform.mat')),fullfile(microDir,micro_sess,'sorting/',sprintf('%s%s',micro_sess,'_spikeWaveform.mat')));
                            movefile(fullfile(microDir,micro_sess,'sortFigs'),fullfile(microDir,micro_sess,'sorting/','sortFigs'));
                        end

                        movefile(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_noreref.mat')),fullfile(microDir,micro_sess,'/raw/',sprintf('%s%s',micro_sess,'_noreref.mat')));
                        movefile(fullfile(microDir,micro_sess,sprintf('%s%s',micro_sess,'_processed.mat')),fullfile(microDir,micro_sess,'/raw/',sprintf('%s%s',micro_sess,'_processed.mat')));

                        if (exist(fullfile(microDir,micro_sess,'ain_timeseries.mat')))
                            movefile(fullfile(microDir,micro_sess,'ain_timeseries.mat'),fullfile(microDir,micro_sess,'raw/','ain_timeseries.mat'));
                        end

                        movefile(fullfile(microDir,micro_sess,'full_ts_raw.png'),fullfile(microDir,micro_sess,'cleaning/','full_ts_raw.png'));
                        movefile(fullfile(microDir,micro_sess,'full_ts_ref.png'),fullfile(microDir,micro_sess,'cleaning/','full_ts_ref.png'));
                        movefile(fullfile(microDir,micro_sess,'global_sigs.png'),fullfile(microDir,micro_sess,'cleaning/','global_sigs.png'));
                        movefile(fullfile(microDir,micro_sess,'nanmean_spectral.png'),fullfile(microDir,micro_sess,'cleaning/','nanmean_spectral.png'));
                        movefile(fullfile(microDir,micro_sess,'variance.csv'),fullfile(microDir,micro_sess,'cleaning/','variance.csv'));
                        movefile(fullfile(microDir,micro_sess,'vars_amps.png'),fullfile(microDir,micro_sess,'cleaning/','vars_amps.png'));
                        movefile(fullfile(microDir,micro_sess,'bad_chans.mat'),fullfile(microDir,micro_sess,'cleaning/','bad_chans.mat'));

                        if (exist(fullfile(microDir,micro_sess,'raw_compare.png')))
                            movefile(fullfile(microDir,micro_sess,'raw_compare.png'),fullfile(microDir,micro_sess,'cleaning/','raw_compare.png'));
                            movefile(fullfile(microDir,micro_sess,'reref_compare.png'),fullfile(microDir,micro_sess,'cleaning/','reref_compare.png'));
                        end

                    end % if strcmp(pipelinePrefix,'mountain') &   [[ clean up mountainsort organization ]]

                    if pipelineType==1
                        %- mountainsort pipeline

                        %- is this spike alignment or LFP aligment?
                        if exist(fullfile(microDir,micro_sess,'raw',sprintf('%s%s',micro_sess,'_spikeInfo.mat')))
                            suffix = '_spikeInfo.mat';
                        else
                            suffix = '_noreref.mat';
                        end
                        microFilename = fullfile(microDir,micro_sess,'raw',sprintf('%s%s',micro_sess,suffix));
                        load(microFilename);
                        lfpPulses_full = pulses;

                    elseif pipelineType==2
                        %- manualsort pipeline
                        microFilename = fullfile(microDir,micro_sess,up_lfp_noreref, 'lfp.mat');
                        if ~exist(microFilename)
                            fprintf('\nfolder exists, but no lfp.mat in %s\n',micro_sess)
                            alignedFiles.alignmentSuccess(rowIterate) = 'no lfp.mat';
                            %cutString = subj;
                            %cutIndex = strfind(fullfile(behavDir,behavFolders(exp).name,sessionFolders(sess).name),cutString);
                            %cutSpikeInfoIndex = strfind(microFilename,cutString);
                            %alignedFiles.ecog_sess(rowIterate) = ecog_sess;
                            %temp = fullfile(behavDir,behavFolders(exp).name,sessionFolders(sess).name,'/');
                            %alignedFiles.beh_sess(rowIterate) = temp(cutIndex:end);
                            %temp = microFilename(cutSpikeInfoIndex(1):end);
                            %temp = strrep(temp,'.mat','_aligned.mat');
                            alignedFiles.micro_sess(rowIterate) = micro_sess;
                            rowIterate = rowIterate + 1;
                            continue
                        end
                        tempLoad  = load(microFilename);
                        lfpStruct = tempLoad.lfpStruct;
                        lfpPulses_full = lfpStruct.pulses;
                    end
                    
                    % Now do the part where you check the duration of the
                    % micro file and see what % of the beh events are
                    % encompassed!
                    % Get rid of the crap earlier
                    %keyboard %Check
%                     duration_micro = days(minutes(lfpStruct.sessDurMin));
%                     micro_end = micro_start + duration_micro;
                    % No, I think it makes more sense to go into the actual
                    % lfp structs to get the lfp data. No idea where that
                    % sessDurMin comes from if there are multiply arrays
                    
                    

                    % lets see if there is a forced alignment file, and load that info if we need it
                    clear forcedPulse forcedFile

                    temp_remove = strfind(microFilename,up_lfp_noreref);
                    temp_filename = microFilename(1:temp_remove-1);
                    temp_filename = fullfile(temp_filename,'force_alignment.txt');
                    fid = fopen(temp_filename);
                    clear tlines
                    if fid ~= -1
                        tline = fgetl(fid);
                        tlines = cell(0,1);
                        while ischar(tline)
                            tlines{end+1,1} = tline;
                            tline = fgetl(fid);
                        end
                        fclose(fid);
                    end
                    if exist('tlines')
                        forcedPulse = tlines{1,1};
                        forcedFile = tlines{2,1};
                    end

                    % lets prepare some variables that will help us keep track of different sets of data that should be aligned separately
                    whichDCused_eachPulseSet = cell(1,size(lfpPulses_full,1));
                    samprate_microLFP_eachPulseSet = cell(1,size(lfpPulses_full,1));
                    signal_noreref_allSets = cell(1,size(lfpPulses_full,1));
                    transform_eachPulseSet = cell(1,size(lfpPulses_full,1));
                    microLFP_eachPulseSet = cell(1,size(lfpPulses_full,1));
                    perc_encompassed = zeros(1,size(lfpPulses_full,1));
                    duration_micro = zeros(1,size(lfpPulses_full,1));
                    samprate_align_eachPulseSet = cell(1,size(lfpPulses_full,1));
                    %micro_zero_eachPulseSet = repelem({[0 0]},1,size(lfpPulses_full,1));
                    %glob_sig_good_allsets  = cell(1,size(lfpPulses_full,1));
                    %glob_sig_all_allsets   = cell(1,size(lfpPulses_full,1));

                    CV_resampled = false;
                    alignmentChanNum = [];

                    whichCellOfData = cell(1,size(lfpPulses_full,1));
                    
                    for numPulsesToAlign = 1:size(lfpPulses_full,1)
                        lfpPulses = lfpPulses_full{numPulsesToAlign,3}; %#ok<*NASGU>

                        if isfield(lfpPulses,'timeseries_freq')
                            samprate_lfpPulses = lfpPulses.timeseries_freq; % SJ: Don't have this in digital
                        else %digital subject?
                            %keyboard
                            samprate_lfpPulses = 1000;
                        end

                        % is the raw file a blackrock or non-blackrock file?
                        rawEcogDir = fullfile(ecogDir,subj,'raw');
                        %rawSessions = getDirNamesRegexp(rawEcogDir,'^\d{6}_\d{4}.*');
                        [rawSessions, rawSessionsisStim] = getRawSessions(rawEcogDir);

                        if (length(find(contains(rawSessions,'ns2'))) > 0) || (length(find(contains(rawSessions,'ns3'))) > 0) %#ok<*ISMT>
                            existTransform = exist(fullfile(rawEcogDir,'NSP_alignment_summary.txt'),'file');
                            if existTransform==0
                                applyTransformOnly=0;
                            else
                                keyboard
                                load(fullfile(rawEcogDir,'NSP_alignment_summary.txt'));
                                % havent hit an example case yet. need to code up the loading of
                                % raw/NSP_alignment_summary.txt, and put it into a cell
                                % variable, 'transforms'
                                applyTransformOnly=1;
                                fprintf('%s\n','Do we get here? The code above needed to be fixed anyway so we probably should not get here. Ask SJ!');
                                keyboard
                                % also, copy NSP_alignment_summary.txt and
                                % NSP_alignment_plot.png into the micro folder
                            end
                        else
                            applyTransformOnly=0;
                        end

                        % raw file is non blackrock, so lets try to get the data into a
                        % format compatible with align_nsps.m
                        %cutString = subj;
                        %cutIndex = strfind(fullfile(behavDir,behavFolders(exp).name,sessionFolders(sess).name),cutString);
                        %cutSpikeInfoIndex = strfind(microFilename,cutString);

                        %alignedFiles.ecog_sess(rowIterate) = ecog_sess;
                        %temp = fullfile(behavDir,behavFolders(exp).name,sessionFolders(sess).name,'/');
                        %alignedFiles.beh_sess(rowIterate) = temp(cutIndex:end);
                        %temp = microFilename(cutSpikeInfoIndex(1):end);
                        %temp = strrep(temp,'.mat','_aligned.mat');
                        alignedFiles.micro_sess(rowIterate) = micro_sess;
                        alignedFiles.alignmentSuccess(rowIterate) = '-';

                        % find the exact name of the channel that carries pulses for this patient
                        if size(indexSubject,1)==0
                            %fprintf('\nneed to define data inputs for this subject in top of file\n')
                            error('need to define data inputs for this subject in top of file')
                        end

                        %indexPulseType = subject_pulseChannels{indexSubject,2}; SJ: don't need
                        pulses_micro_dc09 = NaN;
                        pulses_micro_dc12 = NaN;
                        
                        useAltsig09 = false;
                        useAltsig12 = false;
                        try
                            pulses_micro_dc09 = eval(sprintf('lfpPulses.%s',subject_pulseChannels{indexSubject,2}{1}));
                        catch
                            if ~strcmp(subj,'NIH050')
                                fprintf('%s\n','STOP! There is a session that is not part of NIH050 that does not match the pulses in subject_pulseChannels! Inspect!!!');
                                keyboard 
                            end
                            % SJ: For these particular sessions, the LFP has digital pulses and no analog:
                            %if strcmp(subj,'NIH050') && (strcmp(micro_sess,'170701_1113_beh') || strcmp(micro_sess,'170702_1211_beh') || strcmp(micro_sess,'170630_1320_beh') || strcmp(micro_sess,'170705_1031_beh') || strcmp(micro_sess,'170704_1042_beh'))
                            %old: pulses_micro_dc09 = eval(sprintf('lfpPulses.din1_ts'));
                            pulses_micro_dc09 = eval(sprintf('lfpPulses.%s',subject_pulseChannels{indexSubject,2}{2}));
                            useAltsig09 = true;
                        end
                        removeNegDiffs = find(diff(pulses_micro_dc09) < 0);
                        % remove any series of pulses that occurred inside a
                        % junk segment. if greater than 10, probably not
                        % junk-related and there may have been a clock reset..
                        if removeNegDiffs > 10; fprintf('\nunexpected.. clock reset??\n'); %keyboard;
                        end
                        %pulses_micro_dc09(1:removeNegDiffs(end)+1) = [];
                        
                        try
                            pulses_micro_dc12 = eval(sprintf('lfpPulses.%s',subject_pulseChannels{indexSubject,3}{1}));
                        catch
                            if ~strcmp(subj,'NIH050')
                                fprintf('%s\n','STOP! There is a session that is not part of NIH050 that does not match the pulses in subject_pulseChannels! Inspect!!!');
                                keyboard 
                            end
                            % SJ: For these particular sessions, the LFP has digital pulses and no analog:
                            %if strcmp(subj,'NIH050') && (strcmp(micro_sess,'170701_1113_beh') || strcmp(micro_sess,'170702_1211_beh') || strcmp(micro_sess,'170630_1320_beh') || strcmp(micro_sess,'170705_1031_beh') || strcmp(micro_sess,'170704_1042_beh'))
                            %pulses_micro_dc12 = eval(sprintf('lfpPulses.din4_ts'));
                            pulses_micro_dc12 = eval(sprintf('lfpPulses.%s',subject_pulseChannels{indexSubject,3}{2}));
                            useAltsig12 = true;
                        end
                        removeNegDiffs = find(diff(pulses_micro_dc12) < 0);
                        if removeNegDiffs > 10; fprintf('\nunexpected.. clock reset??\n'); %keyboard;
                        end
                        %pulses_micro_dc12(1:removeNegDiffs(end)+1) = [];

                        if pipelineType==1
                            % prepare the signals
                            noreref_micro_Filename = strrep(microFilename,'spikeInfo','noreref');
                            load(noreref_micro_Filename);
                            keyboard %SJ: yikes... what is going on below?
                            if exist('lfp','var')
                                microLFP = lfp; clear lfp;
                            elseif exist('noreref','var')
                                microLFP = microLFP;
                                fprintf('\n Warning: autosort field names are old and almost out of date!');
                            else
                                fprintf('\n error... field not found'); keyboard;
                            end

                            processed_utah_Filename = strrep(noreref_micro_Filename,'noreref','processed');
                            load(processed_utah_Filename);
                            microLFP_processed = lfp; clear lfp;

                        elseif pipelineType==2
                            % prepare the signals that are to be aligned
                            % by indexing lfpStruct based on which data these pulses
                            % are intended to help in aligning

                            % lets find the filenames corresponding to each of the
                            % lfp data cells
                            temp_filenames_tomatch = cell(1,size(lfpStruct.lfp,2));
                            for j=1:length(temp_filenames_tomatch)
                                allfilenames = lfpStruct.chanIDperNSP{1,j}(:,{'FileName'});
                                allfilenames = unique(allfilenames);
                                if size(allfilenames,1) > 1; fprintf('\njust making sure # unique file names for one NSP == 1\n'); keyboard; end
                                allfilenames = table2cell(allfilenames);
                                temp_filenames_tomatch{1,j} = allfilenames{1,1};
                            end
                            % lets match that filename with these pulses
                            grabThisData = find(strcmp(temp_filenames_tomatch,lfpPulses_full{numPulsesToAlign,2}));
                            if size(grabThisData,1) == 0
                                fprintf('\nweird. pulses point to data from a filename that is not present in the data\n')
                                keyboard
                            end
                            % and grab that data
                            microLFP = lfpStruct.lfp{1,grabThisData};
                            whichCellOfData{1,numPulsesToAlign} = grabThisData;

                            % even if theres no reref, lets pretend
                            reref_micro_Filename = strrep(microFilename,'noreref','reref');
                            if exist(reref_micro_Filename,'file')
                                load(reref_micro_Filename);
                                microLFP_processed = lfp_reref; clear lfp_reref;
                                keyboard % SJ: Haven't gotten here yet. Will we need to go to .lfp? Need to produce another samprate variable to verify that this is 1000?
                                % Also - do we need to get different DC pulses from the processed lfp.mat?
                                % (pulses_micro_DC)
                            else
                                % even if theres no reref, lets pretend there is
                                % for now, to keep compatibility with mountain sort
                                % outputs but if
                                % pipelineType==2, we wont save these out at the end
                                microLFP_processed = microLFP;
                            end

                            % same here.. lets just create NaNs for this variable so that the pipeline
                            % deals with both data in exactly the same way
                            %keyboard % SJ: need to do something with this. If CV, then these assignments must come after the resampling
                            %glob_sig_all = NaN(size(microLFP,2),1);
                            %glob_sig_good = NaN(size(microLFP,2),1);
                            % it only has tLFPdata
                        end
                        
                        % SJ 6/2020: For older patients, if sync pulses timeseries are a different length than the signal
                        % (LFP), then need to pad the end. This is the case because the ts are created from the pulse
                        % timestamps
                        % In the future, may need to keep this within the pipeline #1 b/c pipeline #2 doesn't have these
                        % variables
                        
                        if numel(pulses_micro_dc09) ~= numel(pulses_micro_dc12)
                            fprintf('\n%s\n','ERROR!! pulses_micro_dc09 and pulses_micro_dc12 are not the same length!! Why??');
                            keyboard
                        elseif numel(pulses_micro_dc09) ~= size(microLFP,2)
                            fprintf('\n%s\n','Micro pulse ts are not the same length as lfp data. Padding the end of the pulses with zeros.');
                            if numel(pulses_micro_dc09) > size(microLFP,2)
                                fprintf('\n%s\n','ERROR!!! Micro pulse ts is longer than lfp data!! WHY??');
                                keyboard
                            end
                            % We already know 09 and 12 are the same length, so just calculate with 09
                            pad_amount_pulses = size(microLFP,2) - numel(pulses_micro_dc09);
                            pulses_micro_dc09 = [pulses_micro_dc09 zeros(1,pad_amount_pulses)];
                            pulses_micro_dc12 = [pulses_micro_dc12 zeros(1,pad_amount_pulses)];
                        end
                        
                        samprate_microLFP = lfpStruct.samplingFreq;

                        % if this wasnt a double-nsp ecog patient (in which case the transforms are already
                        % available, go ahead and produce the transforms
                        if applyTransformOnly==0

                            % clear variables, make sure we're not getting
                            % anything left over from previous runs
                            clear signal_noreref_fixed syncPulse syncStartIndex errorStr transforms

                            % prepare a signal1 from the ecog that is the true length
                            % (we're not actually going to use any edits done to signal1,
                            % since that is already fixed through prepAndAlign)
                            split_filename = fullfile(ecogDir,subj,'eeg.noreref',ecog_sess,'DC09'); %-SJ
                            [fchan,msg] = fopen(split_filename,'r+','l');
                            assert(fchan > 0, 'Could not open file %s for writing (error: %s)', split_filename, msg);
                            EEG = fread(fchan,inf,'int16');
                            fclose(fchan);
                            signal1_fake = ones(1,size(EEG,1));
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            % sample rates: [ecog micro]
                            %sampRates = [1000 1000]; %Need to get sample rates from params.txt!!!

                            if rawSessionsisStim(find(strcmp(rawSessions,ecog_sess))) % STIM_MAP folder
                                rawSession_content = getDirNamesRegexp(fullfile(rawEcogDir,'STIM_MAP',ecog_sess),'.*'); 
                            else
                                rawSession_content = getDirNamesRegexp(fullfile(rawEcogDir,ecog_sess),'.*'); % SJ fixed
                            end

                            if any(contains(rawSession_content,'.21E')) % SJ fixed
                                systemType = 'BR'; % need to add new functionality in align_nsps, for NK (pulse height??)
                            elseif any(contains(rawSession_content,'.TRC')) % SJ fixed
                                systemType = 'CV';
                            else
                                if ~any(contains(rawSession_content,'.ns'))
                                    keyboard %SJ check: all folders should have at least one of these endings! SHouldn't get here, take out in future
                                end
                                systemType = 'BR';
                            end
                            % Do a quick check to make sure the sampling rate corresponds to the system type
                            % -SJ
                            if (strcmp(systemType,'BR') && samprate_ecog ~= 1000) || (strcmp(systemType,'CV') && (samprate_ecog ~= 512 && samprate_ecog ~= 1024))
                                fprintf('%s\n',['ERROR!!!!! System type is ' systemType ', but ecog sample rate is ' samprate_ecog '. Ask SJ!']);
                                keyboard
                            end

                            % Do a quick check to make sure the micro pulse and lfp sampling rates are 1000
                            % -SJ
                            if samprate_lfpPulses ~= 1000 || samprate_microLFP ~= 1000
                                fprintf('%s\n%s\n%s\n','ERROR!!! Sample rate for micro LFP and/or pulses are not equal to 1000!',['Sample rate micro lfp Pulses = ' num2str(samprate_lfpPulses)],['Sample rate micro lfp data = ' num2str(samprate_microLFP)]);
                                keyboard
                            end

                            % SJ: now, if CV, need to re-sample ecog pulses ONLY (+signal1_fake) to match  micro (1000 Hz)
                            % Even though we don't have processed yet, let's do it to that as well
                            if strcmp(systemType,'CV')
                                % SJ: It used to be the opposite, where micro was resampled to fit ecog, but commenting that out now
                                fprintf('%s\n','Resampling necessary. Sampling rates of data: ');
                                fprintf('%s\n',['eCOG signal = ' num2str(samprate_ecog)], ...
                                    ['eCOG pulses = ' num2str(samprate_ecog)], ...
                                    ['micro LFP signal = ' num2str(samprate_microLFP)], ...
                                    ['micro LFP pulses = ' num2str(samprate_lfpPulses)], ...
                                    ['micro LFP signal (processed) = ' num2str(samprate_microLFP)], ... %Need to change in future (PROC)
                                    ['micro LFP pulses (processed) = ' num2str(samprate_lfpPulses)]);
                                samprate_align = samprate_microLFP; % Using samprate_align now for everything instead of going back and forth between ecog and micro!

                                if samprate_microLFP ~= samprate_lfpPulses
                                   fprintf('%s\n','ERROR!!!! sample rates of micro LFP and micro LFP pulses are not equal! Need to build this in. Right now, it just uses sample rate of micro LFP for resampling! Ask SJ!');
                                   keyboard                               
                                end

                                keyboard %SJ: Haven't verified this yet!! Can pulses be resampled the same way as signal? Too much smoothing?
                                %resample_from_cell = {microLFPpulses_ecog, microLFP_processed, pulses_micro_dc09, pulses_micro_dc12};
                                resample_from_cell = {signal1_fake, pulses_ecog, pulses_ecog_dc12};
                                resample_to_cell = cell(1,numel(resample_from_cell));

                                fprintf('%s\n',['Preparing to resample ecog pulses to ' num2str(samprate_align) '.']);

                                CV_resampled = true;
                                CV_res_end = ['_' num2str(samprate_ecog)];
                                
                                
                                for rs = 1:numel(resample_from_cell) % Loop over data and pulses to perform the same steps
                                    ecog_mod = resample_from_cell{rs};
                                    if max(ecog_mod(:) >= 32767)
                                        fprintf('%s\n','STOP!! Max ecog_mod is >= intmax16. We have to convert to double when resampling, so how should we treat these values??? (ask SJ)');
                                        keyboard
                                        ecog_mod = ecog_mod.';
                                        ecog_mod = double(ecog_mod);
                                        ecog_mod(ecog_mod == 32767) = NaN;
                                        %micro_mod_resamp = resample(micro_mod,samprate_ecog,samprate_microLFP);
                                        ecog_mod_resamp = resample(ecog_mod,samprate_microLFP,samprate_ecog);
                                        if any(isnan(ecog_mod_resamp(:))) % Putting back the int16max values
                                            ecog_mod_resamp(isnan(ecog_mod_resamp)) = 32767;
                                        end
                                    else
                                        %micro_mod_resamp = resample(double(micro_mod.'),samprate_ecog,samprate_microLFP);
                                        ecog_mod_resamp = resample(double(ecog_mod.'),samprate_microLFP,samprate_ecog);
                                    end
                                    %micro_mod_resamp = resample(double(micro_mod.'),samprate_ecog,samprate_microLFP);
                                    ecog_mod_resamp = int16(ecog_mod_resamp.');
                                    resample_to_cell{rs} = ecog_mod_resamp;
                                    % SJ: faster to do transpose with resample than go row-by-row.
        %                                 for rr = 1:size(microLFP,1)
        %                                     if rr == 1
        %                                         microLFP_resamp = resample(double(microLFP(rr,:)),samprate_ecog,samprate_microLFP);
        %                                     else
        %                                         microLFP_resamp(rr,:) = resample(double(microLFP(rr,:)),samprate_ecog,samprate_microLFP);
        %                                     end
        %                                 end
                                end                                
                                % Now re-assign the variables back
                                signal1_fake = resample_to_cell{1};
                                pulses_ecog = resample_to_cell{2};
                                pulses_ecog_dc12 = resample_to_cell{3};    
                                
%                               OLD WAY (commented out):
                                % SJ: now, if CV, need to re-sample micro to fit ecog
%                                 resample_from_cell = {microLFP, microLFP_processed, pulses_micro_dc09, pulses_micro_dc12};
%                                 resample_to_cell = cell(1,numel(resample_from_cell));
% 
%                                 fprintf('%s\n',['Preparing to resample micro pulses and data to ' num2str(samprate_ecog) '.']);
% 
%                                 CV_resampled = true;
%                                 CV_res_end = ['_' num2str(samprate_ecog)];
%                                 
%                                 
%                                 for rs = 1:numel(resample_from_cell) % Loop over data and pulses to perform the same steps
%                                     micro_mod = resample_from_cell{rs};
%                                     if max(micro_mod(:) >= 32767)
%                                         fprintf('%s\n','STOP!! Max micro_mod is >= intmax16. We have to convert to double when resampling, so how should we treat these values??? (ask SJ)');
%                                         keyboard
%                                         micro_mod = micro_mod.';
%                                         micro_mod = double(micro_mod);
%                                         micro_mod(micro_mod == 32767) = NaN;
%                                         micro_mod_resamp = resample(micro_mod,samprate_ecog,samprate_microLFP);
%                                         if any(isnan(micro_mod_resamp(:))) % Putting back the int16max values
%                                             micro_mod_resamp(isnan(micro_mod_resamp)) = 32767;
%                                         end
%                                     else
%                                         micro_mod_resamp = resample(double(micro_mod.'),samprate_ecog,samprate_microLFP);
%                                     end
%                                     %micro_mod_resamp = resample(double(micro_mod.'),samprate_ecog,samprate_microLFP);
%                                     micro_mod_resamp = int16(micro_mod_resamp.');
%                                     resample_to_cell{rs} = micro_mod_resamp;
%                                     % SJ: faster to do transpose with resample than go row-by-row.
%         %                                 for rr = 1:size(microLFP,1)
%         %                                     if rr == 1
%         %                                         microLFP_resamp = resample(double(microLFP(rr,:)),samprate_ecog,samprate_microLFP);
%         %                                     else
%         %                                         microLFP_resamp(rr,:) = resample(double(microLFP(rr,:)),samprate_ecog,samprate_microLFP);
%         %                                     end
%         %                                 end
%                                 end
% 
%     %                             figure;
%     %                             ct = 1;
%     %                             for sp = 1:numel(resample_to_cell)
%     %                                 if sp == 1
%     %                                     micro_type = 'micro LFP';
%     %                                 elseif sp == 2
%     %                                     micro_type = 'micro LFP Processed';
%     %                                 elseif sp == 3
%     %                                     micro_type = 'micro pulses DC09';
%     %                                 elseif sp == 4
%     %                                     micro_type = 'micro pulses DC12';
%     %                                 else
%     %                                     keyboard %huh?? how do we have more than 4?
%     %                                 end
%     %                                 ct = ct + 1;
%     %                                 subplot(2,4,ct); plot(resample_from_cell{sp}(1,:));
%     %                                 title(['Original ' micro_type ' Chan 1 (sampling rate = ' num2str(samprate_ecog) ')']);
%     % 
%     %                                 ct = ct + 1;
%     %                                 subplot(2,4,ct); plot(resample_to_cell{sp}(1,:),'r'); 
%     %                                 title(['Resampled ' micro_type ' Chan 1 (sampling rate = ' num2str(samprate_microLFP) ')']);
%     %                             end

%                                 % Now re-assign the variables back
%                                 microLFP = resample_to_cell{1};
%                                 microLFP_processed = resample_to_cell{2};
%                                 pulses_micro_dc09 = resample_to_cell{3};
%                                 pulses_micro_dc12 = resample_to_cell{4};

                                fprintf('%s\n','Finished resampling');
                            else
                                CV_res_end = '';
                                samprate_align = samprate_ecog;
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if useAltsig09 && useAltsig12
                                signum = 2;
                            elseif ~useAltsig09 && ~useAltsig12
                                signum = 1;
                            else
                                fprintf('%s\n','ERROR!!!! We are using both analog and digital for either DC09 or DC12 pulses! Not set up for this!');
                                keyboard
                            end
                            alignmentChanNum(numPulsesToAlign) = signum;
                            
                            % prepare the inputs describing pulses
                            if strcmp(subject_pulseChannels{indexSubject,5}{signum},'digital')            
                                sync{1,1} = pulses_ecog_dc12'; sync{1,3} = pulses_ecog';
                                sync{1,2} = double(pulses_micro_dc12); sync{1,4} = double(pulses_micro_dc09);
                            else
                                % if analog, then divide by 4 (gain) to get into proper mV range
                                sync{1,1} = pulses_ecog_dc12'; sync{1,3} = pulses_ecog';
                                sync{1,2} = double(pulses_micro_dc12) / 4; sync{1,4} = double(pulses_micro_dc09) / 4;
                            end

                            % might be a good place to implement a warning about no pulses
                            % found
                            %              if isnan(sync{1,2}) && isnan(sync{1,4})
                            %                 fprintf('\nno pulses pulled from micro data. did you define a chan that is in the pulses struct?\n')
                            %                 keyboard
                            %             end


                            % prepare inputs describing where to output alignment figs
                            sessTitle = ecog_sess;
                            writeFigsHere = fullfile(microDir,micro_sess);
                            
                            %SJ: commented out and moved further up, before
                            %going into numPulsesToAlign

%                             % check to see if we've already aligned this spike info before
%                             %temp_find = microFilename(cutSpikeInfoIndex(1):end);
%                             %temp_find = strrep(temp_find,'.mat','_aligned.mat');
%                             continue_out_of_numPulsesToAlign = false;
%                             temp_find = micro_sess;
%                             triedToAlignBefore = find(strcmp(alignedFiles.micro_sess(1:rowIterate-1),temp_find));
%                             if ~isempty(triedToAlignBefore)
%                                 resultsFromBefore = alignedFiles.alignmentSuccess(triedToAlignBefore);
%                                 %alignedThisFilePairAlready  = size(contains(resultsFromBefore,'alignment'),1); %SJ: I think we want rows?
%                                 % if this pair has already been aligned, dont bother aligning again
%                                 %if (alignedThisFilePairAlready ~= 0) && size(lfpPulses_full,1)==1
%                                 if size(lfpPulses_full,1)==1 %only 1 row
%                                     alignedFiles.alignmentSuccess(rowIterate) = resultsFromBefore{1};
%                                     rowIterate = rowIterate+1;
%                                     continue_out_of_numPulsesToAlign = true;
%                                     continue
%                                 end
%                             end
                            % ok, now lets run alignment
                            differentFileSystems = 1;
                            pulseType = 9;
                            if exist('forcedPulse') && strcmp(forcedPulse,'DC12') %#ok<*EXIST>
                                pulseType = 12;
                            elseif exist('forcedPulse') && strcmp(forcedPulse,'DC09')
                                pulseType = 9;
                            end
                            %  prepare the signals to enter alignment
                            %  in the output, transforms =  {1} = NSP1 add/remove start, end, drift,
                            %    {2} = NSP2 add/remove start, end, drift
                            %    {3} = DRIFT: drift rate, drift added, whichNSP
                            %    {4} = NSP1: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
                            %    {5} = NSP2: originalLength, idxRise, idxFall, deltaRiseFall, finalLengh
                            %    {6} = XCORR correction on indicies: rising diff, falling diff, used correction


                            %% BLOCK OF CODE FOR CORRECTING DATA IN CASES WHERE
                            % ONE FILE SYSTEM WAS STOPPED AT A VERY DIFFERENT
                            % TIME (FOR EXAMPLE, 30 MINUTES OF DATA IN MICRO,
                            % 60 MINUTES OF DATA IN ECOG). STILL IN PROGRESS.
                            % NEEDS TWEAKING IN ORDER TO BE PERFECT AND APPLY
                            % TO MORE CASES

                            % if we know that the micro file ended far before the ecog file, lets
                            % go ahead and NaN pad the later part of the ecog file in hope of making
                            % alignment work out
                            
                            %duration_micro = days(minutes(lfpStruct.sessDurMin));
                            %duration_micro(numPulsesToAlign) = days(seconds(size(microLFP,2)/samprate_ecog)); % Using samprate_ecog bc that's what it was orig resampled to... not anymore!
                            duration_micro(numPulsesToAlign) = days(seconds(size(microLFP,2)/samprate_align)); %SJ changed CVCH
                            micro_end = micro_start + duration_micro(numPulsesToAlign);
                            
                            zfs = find(strcmp(subj,zero_forced.Subject) & strcmp(ecog_sess,zero_forced.ecog_sess) & strcmp(micro_sess,zero_forced.micro_sess));
                            if ~isempty(zfs)
                                % Special case!!
                                %zero_forced_FLAG = true;
                                zeroAmountList = zero_forced.zero{zfs}; % 1) ecog Start, 2) ecog end, 3) micro start, 4) micro end
                            else
                                zeroAmountList = zeros(1,4); % 1) ecog Start, 2) ecog end, 3) micro start, 4) micro end
                            end

                            %Check if there is > 5 min difference between starts and ends of ecog and micro (because in align_nsps, xcorr only looks out 5 min)
                            % zeroNeeded = [padStart, padEnd] (0/1)
                            zeroNeeded = zeros(1,2);
                            if abs(ecog_start - micro_start) > .0035 || (zeroAmountList(1)+zeroAmountList(3)) > 0 %Going to have to zero a part of the START or micro or ecog; or manual zeroing needed for start of either
                                zeroNeeded(1) = 1; 
                            end
                            if abs(ecog_end - micro_end) > .0035 || (zeroAmountList(2)+zeroAmountList(4)) > 0 %Going to have to zero a part of the END or micro or ecog; or manual zeroing needed for end of either
                                zeroNeeded(2) = 1;
                            end
                            %SJ note: in the future, put the zero_forced flag up before zeroNeeded, that way you will go into the different scenarios anyway (done)
                            
%                             if ecog_start - micro_start > .0035 %micro before ecog by 5 min, need to pad ecog
%                                 padNeeded(1) = 1;
%                             elseif ecog_start - micro_start < -.0035 %ecog before micro by 5 min, need to pad micro
%                                 padNeeded(2) = 1;
%                             end
%                             if ecog_end - micro_end < -.0035 %micro ends after ecog by 5 min, need to pad ecog
%                                 padNeeded(3) = 1;
%                             elseif ecog_end - micro_end > .0035 %ecog ends after micro by 5 min, need to pad micro
%                                 padNeeded(4) = 1;
%                             end
                            
                            %SJ added: check if there is a forced zero situation... then apply that INSTEAD of the calculated one if we get there!!

                             
                            if any(zeroNeeded)
                                fprintf('%s\n','Differences > 5 min in start/end times for ecog & micro! About to manually manipulate.');
                                %keyboard
                                for pd = 1:numel(zeroNeeded) %Do 2 loops - one for start and one for end 
                                    if zeroNeeded(pd) == 1 %Need to Zero
                                        if pd == 1 %On the first loop, Zeroing the START
                                            zeroAmount_days = ecog_start - micro_start;
                                            
                                            if zeroAmount_days < 0 || zeroAmountList(1) ~= 0 %zero start of ecog ((SJ- changed to <, somehow was >??? 4/17/20))
                                                % ECOG starts before micro
                                                if zeroAmountList(1) ~= 0
                                                    zeroAmount = zeroAmountList(1);
                                                else
                                                    zeroAmount = ceil((seconds(days(zeroAmount_days))+120)*samprate_align); %convert to samples (while adding 2 min)
                                                end
                                                fprintf('%s\n',['Zeroing the start of the ecog sync pulses (' num2str(abs(zeroAmount)) ' samples).']);
                                                %sync{1,1} = cat(2,NaN(size(sync{1,1},1),zeroAmount),sync{1,1});
                                                %sync{1,3} = cat(2,NaN(size(sync{1,3},1),zeroAmount),sync{1,3});
                                                sync{1,1}(:,1:abs(zeroAmount)) = 0;
                                                sync{1,3}(:,1:abs(zeroAmount)) = 0;
                                                zeroAmountList(1) = zeroAmount;
                                                
                                            elseif zeroAmount_days > 0 || zeroAmountList(3) ~= 0 %Zero start of micro
                                                if zeroAmountList(3) ~= 0
                                                    zeroAmount = zeroAmountList(3);
                                                else
                                                    zeroAmount = ceil((seconds(days(zeroAmount_days))-120)*samprate_align); %convert to samples (while adding 2 min- NEGATIVE)
                                                end
                                                    fprintf('%s\n',['Zeroing the start of the micro sync pulses (' num2str(abs(zeroAmount)) ' samples).']);
                                                %sync{1,2} = cat(2,NaN(size(sync{1,2},1),abs(zeroAmount)),sync{1,2});
                                                %sync{1,4} = cat(2,NaN(size(sync{1,4},1),abs(zeroAmount)),sync{1,4});
                                                sync{1,2}(:,1:abs(zeroAmount)) = 0;
                                                sync{1,4}(:,1:abs(zeroAmount)) = 0;
                                                micro_zero_eachPulseSet{1,numPulsesToAlign}(1) = zeroAmount;
                                                zeroAmountList(3) = zeroAmount;
                                            end
                                            
                                        else %On the second loop, Zeroing the END
                                            
                                            zeroAmount_days = ecog_end - micro_end;
                                            
                                            if zeroAmount_days < 0 || zeroAmountList(4) ~= 0 %zero end of micro ((SJ- changed to <, somehow was >??? 4/17/20))
                                                if zeroAmountList(4) ~= 0
                                                    zeroAmount = zeroAmountList(4);
                                                else
                                                    zeroAmount = ceil((seconds(days(zeroAmount_days))+120)*samprate_align); %convert to samples (while adding 2 min)
                                                end
                                                fprintf('%s\n',['Zeroing the end of the micro sync pulses (' num2str(zeroAmount) ' samples).']);
                                                
                                                sync{1,2}(:,end-abs(zeroAmount)+1:end) = 0;
                                                sync{1,4}(:,end-abs(zeroAmount)+1:end) = 0;
                                                micro_zero_eachPulseSet{1,numPulsesToAlign}(2) = zeroAmount;
                                                zeroAmountList(4) = zeroAmount;
                                                
                                            elseif zeroAmount_days > 0 || zeroAmountList(2) ~= 0 %Zero end of ecog
                                                if zeroAmountList(2) ~= 0
                                                    zeroAmount = zeroAmountList(2);
                                                else
                                                    zeroAmount = ceil((seconds(days(zeroAmount_days))-120)*samprate_align); %convert to samples (while adding 2 min- NEGATIVE)
                                                end
                                                fprintf('%s\n',['Zeroing the end of the ecog sync pulses (' num2str(abs(zeroAmount)) ' samples.)']);
                                                
                                                sync{1,1}(:,end-abs(zeroAmount)+1:end) = 0;
                                                sync{1,3}(:,end-abs(zeroAmount)+1:end) = 0;
                                                zeroAmountList(2) = zeroAmount;
                                            end
                                        end
                                    end
                                end
                            end
                            
%                             %SJ: added special cases list for failed alignments- zero one of them more manually
%                             switch subj
%                                 case 'NIH047'
%                                     if strcmp(ecog_sess,'170311_1312') && strcmp(micro_sess,'170311_1259')
%                                         fprintf('\n%s\n',['Heads up... Manually manipulating special case, ' subj ': ecog: ' ecog_sess ' ;  micro_sess: ' micro_sess]);
%                                         %Beforehand: zeroing start of micro by 11 min (660001 samples), but need to do around 5 min more
%                                         zeroAmountList(3) = zeroAmountList(3) + seconds(minutes(5))*1000;
%                                     end
%                             end

                            
%                             if (size(signal1_fake,2) / size(microLFP,2)) > 1.5
%                                 fprintf('\nbe aware.. code is about to manipulate the pulses before trying to align. may need manual tweaking\n')
%                                 keyboard
% 
%                                 % SJ: Don't do anything actually- how does this perform?
% 
%                                 % ecog is much longer. assuming start at same
%                                 % time, remove big end chunk of ecog
%                                 %SJ: No, instead we will pad the end of the micro!!
%                                 difference = size(signal1_fake,2)-size(microLFP,2);
% 
%                                 % get start and end times of files. if end times
%                                 % are more different than start times (micro was
%                                 % stopped early):
%                                 % if diffs>0, then ecog started later
%                                 proportionOfDifferenceToRemoveFromEnd = 1;
%                                 %sync{1,1}(end-(round(difference*proportionOfDifferenceToRemoveFromEnd)):end) = 0;
%                                 %   sync{1,2} = cat(2,sync{1,2},NaN(size(sync{1,2},1),difference*proportionOfDifferenceToRemoveFromEnd));
%                                 %sync{1,3}(end-(round(difference*proportionOfDifferenceToRemoveFromEnd)):end) = 0;
%                                 %  sync{1,4} = cat(2,sync{1,4},NaN(size(sync{1,4},1),difference*proportionOfDifferenceToRemoveFromEnd));
% 
%                                 %  elseif start times are more different (micro was
%                                 %  started halfway in ecog):
%                                 %                          sync{1,1}(1:(round(difference*proportionOfDifferenceToRemoveFromEnd))) = 0;
%                                 %                          sync{1,3}(1:(round(difference*proportionOfDifferenceToRemoveFromEnd))) = 0;
%                             end




                            % Now check how many events are within that time frame
                            events_times = epoch2date([sess_events.mstime]);
                            perc_encompassed(numPulsesToAlign) = sum(events_times>=micro_start & events_times<= micro_end)/numel(events_times);
                            % Oops, okay the above shit wasn't supposed to happen yet

                            %% now lets align
                            
                            fprintf('\nabout to align micro file %s with noreref %s\n',microSessions(closest).name,ecog_sess)
                            [signal_noreref_fixed, syncPulse, syncStartIndex, errorStr, transforms] = align_nsps(sync, signal1_fake, microLFP,samprate_align,systemType,pulseType,writeFigsHere,sessTitle,differentFileSystems,'keepFirstSignal',1,'isMicro',1,'overwriteLabels',{'ECOG', 'MICRO',lfpPulses_full{numPulsesToAlign,1}},'zeroNeeded',zeroAmountList); %#ok<*ASGLU>
                            
                            if strcmpi(systemType,'CV') %SJ CVCH Has to be done because align_nsps stacks the 2 signals for nsp alignment...
                                signal_noreref_fixed(1,:) = [];
                            end
                            
                            
                            if errorStr == 2
                                whichDCused = 12;
                                sync1 = sync{1,1}; %pulses_ecog_dc12'
                                sync2 = sync{1,2}; %pulses_micro_dc12
                            elseif errorStr == 3
                                whichDCused = 9;
                                sync1 = sync{1,3}; %pulses_ecog'
                                sync2 = sync{1,4}; %pulses_micro_dc09
                            else
                                whichDCused =0;
                            end
                            
%                             % SJ added 8/2020: Create figure similar to NSP alignment showing how the micro, ecog, and beh sessions were aligned, showing:
%                             % - start time and end time of each session
%                             % - any necessary zeroing
%                             % - any necessary trimming/padding
%                             
%                             % Find longest length
%                             m_beh_start   = minutes(days(beh_start));    m_beh_end   = minutes(days(beh_end));
%                             m_ecog_start  = minutes(days(ecog_start));   m_ecog_end  = minutes(days(ecog_end));
%                             m_micro_start = minutes(days(micro_start));  m_micro_end = minutes(days(micro_end));
%                             
%                             max_x = max([m_beh_end - m_beh_start, m_ecog_end - m_ecog_start, m_micro_end - m_micro_start]);
%                             first_x = min([m_beh_start,m_ecog_start,m_micro_start]);
%                             
%                             beh_x = [m_beh_start m_beh_end] - first_x
%                             ecog_x = [m_ecog_start m_ecog_end] - first_x;
%                             micro_x = [m_micro_start m_micro_end] - first_x;
%                             
%                             % (seconds(days(zeroAmount_days))+120)*samprate_ecog);
%                             % minutes(seconds(num_samps/samprate_ecog))
%                             % ecog_edit = [start_orig, start_zero, start_trim, start_align
%                             ecog_micro_x = [ecog_x; micro_x];
%                             em_edit = cell(1,3);
%                             for ii = 1:2
%                                 start_orig = ecog_micro_x(ii,1);
%                                 if zeroAmountList((ii-1)*2+1) == 0
%                                     start_zero = NaN;
%                                 else
%                                     start_zero = zeroAmountList((ii-1)*2+1);
%                                 end
%                                 start_trim = start_orig + minutes(seconds(transforms{ii}(1)/samprate_ecog));
%                                 
%                                 end_orig = ecog_micro_x(ii,2);
%                                 if zeroAmountList((ii-1)*2+2) == 0
%                                     end_zero = NaN;
%                                 else
%                                     end_zero = end_orig-zeroAmountList((ii-1)*2+2);
%                                 end
%                                 end_trim = end_orig - minutes(seconds(transforms{ii}(2)/samprate_ecog));
%                                 em_edit{ii+1} = {[start_orig end_orig], [start_zero end_zero], [start_trim end_trim]};
%                             end
%                             h_line = [4 3 2];
%                             colors_line = {[0 .9 .3], [0 .7 1], [1 .6 .1]};
%                             x_all = {beh_x, ecog_x, micro_x};
%                             
%                             figure
%                             for st = 1:3
%                                 subplot(1,2,1); plot(x_all{st},ones(1,numel(x_all{st}))*h_line(st),'Color',colors_line{st},'LineWidth',15,'Marker','.','MarkerSize',8,'MarkerEdgeColor','k');
%                                 hold on
%                                 xline(x_all{st}(1),'--','Color',colors_line{st});
%                                 xline(x_all{st}(2),'--','Color',colors_line{st});
%                             end
%                             
%                             xlim([-3 max_x+3])
%                             ylim([0 6])   
%                             x_all_aligned = cell(1,3);
%                             
%                             for st2 = 1:3
%                                 if st2 == 1
%                                     subplot(1,2,2); plot(x_all{st2},ones(1,numel(x_all{st2}))*h_line(st2),'Color',colors_line{st2},'LineWidth',15,'Marker','.','MarkerSize',8,'MarkerEdgeColor','k');
%                                     hold on   
%                                     xline(x_all{st2}(1),'--','Color',colors_line{st2});
%                                     xline(x_all{st2}(2),'--','Color',colors_line{st2});
%                                 else
%                                     if ~isnan(em_edit{st2}{2}(1))
%                                         subplot(1,2,2); rectangle('position',[em_edit{st2}{1}(1) h_line(st2)-.05 em_edit{st2}{2}(1)-em_edit{st2}{1}(1) .1],'FaceColor','w','EdgeColor',colors_line{st2});
%                                         %plot([em_edit{st2}{1}(1) em_edit{st2}{2}(1)],ones(1,2)*h_line(st2),'-','Color',colors_line{st2},'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',25,'LineWidth',8); % beginning zero
%                                     end
%                                     if ~isnan(em_edit{st2}{2}(2))
%                                         subplot(1,2,2); rectangle('position',[em_edit{st2}{2}(2) h_line(st2)-.05 em_edit{st2}{2}(2)-em_edit{st2}{1}(2) .1],'FaceColor','w','EdgeColor',colors_line{st2});
%                                         %plot([em_edit{st2}{2}(2) em_edit{st2}{1}(2)],ones(1,2)*h_line(st2),'-','Color',colors_line{st2},'Marker','.','MarkerEdgeColor','k','MarkerEdgeColor','k','MarkerSize',25,'LineWidth',8); %end zero
%                                     end                                    
%                                     % ecog & micro
%                                     subplot(1,2,2); plot(x_all{st2},ones(1,numel(ecog_x))*h_line(st2),'--','Color',colors_line{st2},'Marker','.','MarkerSize',8,'MarkerEdgeColor','k'); % dotted original
%                                     
%                                     % Smaller line with padding (min(trim,orig,zero)) -> max(
%                                     subplot(1,2,2); plot([em_edit{st2}{3}(1) em_edit{st2}{3}(2)],ones(1,2)*h_line(st2),'-','Color',colors_line{st2},'LineWidth',6); % Aligned portion
%                                     % larger line with trim (max(trim,orig)
%                                     subplot(1,2,2); plot([max([em_edit{st2}{3}(1) em_edit{st2}{2}(1) em_edit{st2}{1}(1)]) min([em_edit{st2}{3}(2) em_edit{st2}{2}(2) em_edit{st2}{1}(2)])],ones(1,2)*h_line(st2),'-','Color',colors_line{st2},'LineWidth',15); % Aligned portion
% 
%                                     xline(em_edit{st2}{3}(1),'--','Color',colors_line{st2});
%                                     xline(em_edit{st2}{3}(2),'--','Color',colors_line{st2});
%                                     
%                                     subplot(1,2,2); text(em_edit{st2}{1}(1),h_line(st2)-.1, ['Trim at start = ' num2str(minutes(seconds(transforms{st2-1}(1)/samprate_ecog))) ' min']);
%                                     subplot(1,2,2); text(em_edit{st2}{1}(2),h_line(st2)-.1, ['Trim at end = ' num2str(minutes(seconds(transforms{st2-1}(2)/samprate_ecog))) ' min']);
%                                     subplot(1,2,2); text(em_edit{st2}{1}(1),h_line(st2)-.15, ['Zeroed at start = ' num2str(zeroAmountList((st2-2)*2+1)) ' min']);
%                                     subplot(1,2,2); text(em_edit{st2}{1}(2),h_line(st2)-.15, ['Zeroed at start = ' num2str(zeroAmountList((st2-2)*2+2)) ' min']);
%                                 end
%                                 
%                             end
%                             xlim([-3 max_x+3])
%                             ylim([0 6]) 
%                             
%                             
                            %annotation('textbox',[0 .4 .5 .3],'String',['Trim/Add to start = ' num2str(minutes(seconds(transforms{1}(1)/samprate_ecog))) ' min'],'FitBoxToText','on');
                            
                            % output format:
                            % transforms{1,1} = ecog transforms
                            % [add(-)orRemove(+)fromStart    add(-)orRemove(+)fromEnd   Drift]
                            if checkPulses==1 && errorStr~=1 && errorStr~=4
                                % SJ: this used to call transformSync, which would perform a different set of
                                % operations. Instead, now just call applyAlignCorrection with the sync
                                % pulses
    %                             iMod_lfp = 2;
    %                             driftRate_lfp = 1;
    %                             driftFull_lfp = 1;
    %                             removeDrift_lfp = true;
    %                             keepFirstSignal_lfp = 1;
                                %                   signa1, signal2, iMod, driftRate,            driftFull, removeDrift, keepFirstSignal
                                %applyAlignCorrection(sync1,   sync2,    2,         0, transforms{1,2}(1,3),        true,               1,  

                                % transform original sync using transforms, then run get_triggers on it
    %                             transformedSync = transformSync(sync,transforms);
    %                             if pulseType==12 %SJ: VERIFY that pulsetype matches errorstr output from align_nsps (2 for DC12, 3 for DC09)
    %                                 micro_ms = transformedSync{1,2};
    %                                 ecog_ms = sync{1,1};
    %                             elseif pulseType==9
    %                                 micro_ms = transformedSync{1,1};
    %                                 ecog_ms = sync{1,3};
    %                             end
    %                             triggers_ecog = get_triggers(ecog_ms',1000);
    %                             triggers_micro = get_triggers(micro_ms',1000);
    %                             threshMS = 10; mywin = 200; samplerate = 1000;
    %                             % beh_ms and pulses should be vector of doubles
    %                             [ecog_pulses_use, micro_pulses_use] = pulsealign(triggers_ecog{1,1},triggers_micro{1,1},samplerate,threshMS,mywin,0,0,1);
    %                             ms_field = 'mstime';
    %                             [alignInfo] = logalign_microVsEcog({ecog_pulses_use},{micro_pulses_use},ms_field);
    %                             writeOutAlignInfo = fullfile(microDir,micro_sess,'alignmentStats.txt'); % -SJ
    %                             fid  = fopen(writeOutAlignInfo,'w');
    %                             if fid==-1; fprintf('cant open file..\n'); keyboard; end
    %                             fprintf(fid,'NumPoints=%0.4f\n',alignInfo.reg_numPointFit);
    %                             fprintf(fid,'Intercept=%0.4f\n',alignInfo.reg_intercept);
    %                             fprintf(fid,'Slope=%0.4f\n',alignInfo.reg_slope);
    %                             fprintf(fid,'R^2=%0.4f\n',alignInfo.reg_Rsquare);
    %                             fprintf(fid,'MaxDev=%0.4f\n',alignInfo.reg_maxDev);
    %                             fprintf(fid,'MedianDev=%0.4f\n',alignInfo.reg_medianDev);
    %                             fprintf(fid,'PulsesUsed=%d\n',pulseType);%SJ: VERIFY that pulsetype matches errorstr output from align_nsps (2 for DC12, 3 for DC09)
    %                             fclose(fid);
                            end

                        end

                        % if there is an errorStr indicating not to align,
                        % lets make sure no transforms are executed, and lets
                        % document the csv accordingly

                        if errorStr==1 || errorStr==4
                            alignedFiles.alignmentSuccess(rowIterate) = 'attempted & failed';
                            transforms = []; %apply nothing to the data
                        elseif errorStr==0 || errorStr==2 || errorStr==3
                            alignedFiles.alignmentSuccess(rowIterate) = 'attempted & succeeded';
                        end

                        whichDCused_eachPulseSet{1,numPulsesToAlign} = whichDCused;
                        samprate_microLFP_eachPulseSet{1,numPulsesToAlign} = samprate_microLFP;
                        transform_eachPulseSet{1,numPulsesToAlign} = transforms;
                        microLFP_eachPulseSet{1,numPulsesToAlign} = microLFP;
                        signal_noreref_allSets{1,numPulsesToAlign} = signal_noreref_fixed;
                        samprate_align_eachPulseSet{1,numPulsesToAlign} = samprate_align;
                        %micro_zero_eachPulseSet{1,numPulsesToAlign} = 
                        %glob_sig_good_allsets{1,numPulsesToAlign} = glob_sig_good;
                        %glob_sig_all_allsets{1,numPulsesToAlign} = glob_sig_all;

                        %%%%%%%%%%%%%%%%%%%%%%%%

                        if pipelineType==1
                            keyboard
                        else
                            try %SJ added for special cases for CV
                                lfp_aligned = signal_noreref_fixed{1,2};
                            catch
                                keyboard %SJ: I think this should never happen anymore, right? At least we never got here after fixing NIH050!
                                %if strcmp(subj,'NIH050') && (strcmp(micro_sess,'170701_1113_beh') || strcmp(micro_sess,'170702_1211_beh') || strcmp(micro_sess,'170630_1320_beh') || strcmp(micro_sess,'170705_1031_beh') || strcmp(micro_sess,'170704_1042_beh'))
                                if strcmpi(systemType,'CV') || strcmpi(subj,'NIH050')
                                    if strcmpi(subj,'NIH050') && ~(strcmp(micro_sess,'170701_1113_beh') || strcmp(micro_sess,'170702_1211_beh') || strcmp(micro_sess,'170630_1320_beh') || strcmp(micro_sess,'170705_1031_beh') || strcmp(micro_sess,'170704_1042_beh'))
                                        keyboard
                                    end
                                    lfp_aligned = signal_noreref_fixed;
                                else
                                    keyboard
                                end
                            end
                            lfpStruct.lfp{1,whichCellOfData{1,numPulsesToAlign}} = lfp_aligned;
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%
                    end

                    %SJ: do a few checks to ensure none of these change between pulseSets!
                    if (numel(unique([whichDCused_eachPulseSet{:}])) ~= 1) || (numel(unique([samprate_microLFP_eachPulseSet{:}])) ~= 1) || (numel(unique([samprate_align_eachPulseSet{:}])) ~= 1)
                        fprintf('%s\n','One of the following is different for each of the pulse sets!! This might be a problem: ');
                        fprintf('\t%s\n\t%s\n\t%s\n','whichDCused','samprate_microLFP','samprate_align');
                        keyboard
                    end
                    

                    alignedFiles.eventsEncompassedByMicroFile(rowIterate) = round(max(perc_encompassed),2);
                    alignedFiles.lengthOfMicroFileMin(rowIterate) = round(max(minutes(days(duration_micro))),2);
                    %SJ: Take the max value if there are multiple 
                    
                    % Look thorugh the starting/ending indices of micro
                    % from the transforms and calculate how much of the
                    % beh mstimes are contained within the trimmed micro
                    % file - return the largest % out of the number of
                    % things to align!!!
                    %keyboard
                    %transform_eachPulseSet
                    % If had to zero some of micro, then only count
                    % non-zero portions
                    
                    % Not putting this in now.... SJ

%                     if continue_out_of_numPulsesToAlign
%                         continue
%                     end

                    acn = unique(alignmentChanNum); %SJ added to carry over from numPulsestoAlign, just do a quick check
                    if numel(acn) > 1 || all(acn == 1)
                        % acn > 1 means we are using both analog and digital pulses, as in NIH050 for some sessions. We wil want to keep both in alignmentChan for lfpStruct and spikeInfo!
                        % acn == 1 means we only used analog
                    elseif acn == 2
                        % Only used digital
                    else
                        keyboard % WHAAAAAAT the f is going on here
                    end
                        
                    lfpStruct.alignedTo = ecog_sess;
                    if whichDCused == 9
                        lfpStruct.alignmentChan = subject_pulseChannels{indexSubject,2}(alignmentChanNum);
                    elseif whichDCused == 12
                        lfpStruct.alignmentChan = subject_pulseChannels{indexSubject,3}(alignmentChanNum);
                    else
                        lfpStruct.alignmentChan = 'no alignment';
                    end

                    % while we're outside of pulse-set loop, lets  record and save changes made to each file in the subject folder
                    % also, lets make sure that we havent already aligned this file
                    % old: alignedFiles{rowIterate-1,4} = sprintf('%0.2f seconds',(timeDiff_regress/1000)); %minus 1 because the counter has already iterated

                    %cell2csv(alignedFiles,alignment_xlsx_writeOut)
                    %writetable(alignedFiles,alignment_xlsx_writeOut); %SJ: replaced with after lfp alignment

                    % now, iterate csv
                    rowIterate = rowIterate+1;
                    adjustedSpikes = {};
                    adjustedSpikes_noise = {};
                    


                    for pulseSet=1:size(lfpPulses_full,1)

                        if isempty(transform_eachPulseSet{1,pulseSet}) %SJ : numPulsesToAlign changed to pulseSet
                            % must have aligned this file already
                            continue
                        end

                        whichDCused = whichDCused_eachPulseSet{1,pulseSet};
                        samprate_microLFP = samprate_microLFP_eachPulseSet{1,pulseSet};
                        transforms = transform_eachPulseSet{1,pulseSet}; %SJ : numPulsesToAlign changed to pulseSet
                        microLFP = microLFP_eachPulseSet{1,pulseSet};
                        samprate_align = samprate_align_eachPulseSet{1,pulseSet};
                        %glob_sig_good = glob_sig_good_allsets{1,pulseSet};
                        %glob_sig_all = glob_sig_all_allsets{1,pulseSet};

                        % now, go through each aligned pulse set and apply the
                        % transforms only to chans with nsp suffix described in
                        % lfpPulses_full(#,1)

                        %% now, lets use the transforms to make the necessary corrections, as long as we havent
                        % explicitly defined to not write anything out

                        %remove
                        driftApplyToMicro = transforms{1,2}(1,3);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SPIKES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        manual_spike_folder = 'spikes'; %- was sorted, now spikes

                        spikeInfoFound = 0;
                        if ((pipelineType==1) && exist(fullfile(microDir,micro_sess,'raw',sprintf('%s%s',micro_sess,'_spikeInfo.mat')))) || ((pipelineType==2) && exist(fullfile(microDir,micro_sess,manual_spike_folder)))
                            spikeInfoFound = 1;
                        end
                        if spikeInfoFound

                            % JW asks... how is the data supposed to be organized in here????
                            %   Ohhh. was originally micro/manualsorts/sorted/reref_sortedByDY/xxxspikeInfo.mat
                            %             now we want micro/manualsorts/spikes/xxxspikeInfo.mat
                            %   mountainsort organization is different... might make sense to make them match

                            if pipelineType==2; spikeInfos = dir(fullfile(microDir,micro_sess,manual_spike_folder,'*spikeInfo.mat'));
                            else                fprintf('\n jw probably screwed up the mountainsort path below... check it!'); keyboard; end

                            spikeInfos(find(strncmp({spikeInfos.name},'.',1))) = []; %- if first character is ., then remove it

                            %- heads up... illogical separation between "spikeInfos" and spikeInfoToPull is because JW previously thought multiple sorters might contribute to a sort folder... better to just pick a winner
                            for numSpikeInfo = 1:length(spikeInfos)
                                %spikeInfoToPull = dir(fullfile(microDir,micro_sess,'sorted/',spikeInfos(numSpikeInfo).name,'*_spikeInfo.mat')); %- designed for mountainsort?
                                spikeInfoToPull = dir(fullfile(microDir,micro_sess,manual_spike_folder,spikeInfos(numSpikeInfo).name));
                                if length(spikeInfoToPull)==0
                                    continue
                                end
                                %load(fullfile(microDir,micro_sess,'sorted/',spikeInfos(numSpikeInfo).name,spikeInfoToPull(1).name))
                                tempLoad = load(fullfile(microDir,micro_sess,manual_spike_folder,spikeInfoToPull(1).name)); %- jw hack
                                spikeInfo = tempLoad.spikeInfo;

                                % lets identify the units that these pulses
                                % should be used to transform
                                temp_table = spikeInfo.sessUniqueUnitID(:,{'NSxFileName'});
                                temp_table = table2cell(temp_table);
                                apply_theseChans = find(strcmp(lfpPulses_full{pulseSet,2},temp_table));
                                
                                if ~isfield(spikeInfo,'sessNoiseUnitID')
                                    % Older spikeInfo with no noise units
                                    % built in
                                    spikeInfo.sessNoiseUnitID = {};
                                    if isfield(spikeInfo,'timeStampNoise')
                                        fprintf('%s\n','ERROR!! timeStampNoise exists but sessNoiseUnitID does not??');
                                        keyboard
                                    end
                                    spikeInfo.timeStampNoise = {};
                                    apply_theseChans_noise = [];
                                elseif iscell(spikeInfo.sessNoiseUnitID) && contains(spikeInfo.sessNoiseUnitID,'No isolated noise units')
                                    % No noise units identified. Don't have to do anything
                                    spikeInfo.sessNoiseUnitID = {};
                                    apply_theseChans_noise = [];
                                else
                                    temp_table_noise = spikeInfo.sessNoiseUnitID(:,{'NSxFileName'});
                                    temp_table_noise = table2cell(temp_table_noise);
                                    apply_theseChans_noise = find(strcmp(lfpPulses_full{pulseSet,2},temp_table_noise));
                                end
                                    


                                % if first set of pulses, lets grab the raw data
                                % if > first set, lets grab pulses that
                                % were saved out in previous iterations
                                timeStamp_raw = spikeInfo.timeStamp;
                                
                                %SJ: Now need to do something similar for
                                %noise unit timestamps
                                timeStampNoise = spikeInfo.timeStampNoise;
                                
                                %Check if we even have any noise units...?
                                %(or just begin with empty cell arrays)
                                
                                
    %                             if pulseSet==1
    %                                 timeStamp = timeStamp_raw;
    %                             else
    %                                 timeStamp = adjustedSpikes;
    %                             end

                                timeStamp = timeStamp_raw;
                                % write new field with adjusted spike times.. includes
                                % offset and clock drift rate adjustments..

                                % apply drift correction - if we inserted samples into ecog, lets take samples out
                                % of spike times. if we inserted samples into spikes,
                                % lets insert samples into spikes.. and lets
                                % account that timestamps are in indices at 30 kHz

                                %maxSize = size(signal_noreref_allSets{1,pulseSet}{1,1},2);
                                %correctionInterval = floor((maxSize*SRmatchECOG)/(abs(driftApplyToMicro)*SRmatchECOG));

                                % Shouldn't it be maxSize = size of original micro? (For LFP, it is)

                                SRmatchECOG = samprate_ecog/samprate_microLFP; %1000 is the samprate_micro
                                iMod = 2;
                                keepFirstSignal = 1;

                                % Using 0 drift Rate because it doesn't actually do anything here except put
                                % it in the output transforms (which we don't care about here)
                                % SJ: remove SRmatchECOG!!! CVCH
                                [~, ~, SPKtimestampsFixed, ~, ~] = applyAlignCorrection(signal1_fake, microLFP, iMod, 0, driftApplyToMicro, false, keepFirstSignal, 'SPKtimestamps', timeStamp, 'apply_theseChans',apply_theseChans, 'alignType', 'SPK', 'SPKtrimStart',transforms{1,2}(1,1));
                                
                                if ~isempty(timeStampNoise)
                                    [~, ~, SPKtimestampsFixed_noise, ~, ~] = applyAlignCorrection(signal1_fake, microLFP, iMod, 0, driftApplyToMicro, false, keepFirstSignal, 'SPKtimestamps', timeStampNoise, 'apply_theseChans',apply_theseChans_noise, 'alignType', 'SPK', 'SPKtrimStart',transforms{1,2}(1,1));
                                else
                                    SPKtimestampsFixed_noise = {};
                                end
                                
                                spikeInfo.alignedTo = ecog_sess;
                                if whichDCused == 9
                                    spikeInfo.alignmentChan = subject_pulseChannels{indexSubject,2}(alignmentChanNum); % SJ: Kind of awkward because this writes out twice... just take the alignment chan info (acn) from the lfp alignment, which tells us if 1 or both were used 
                                elseif whichDCused == 12
                                    spikeInfo.alignmentChan = subject_pulseChannels{indexSubject,3}(alignmentChanNum);
                                else
                                    spikeInfo.alignmentChan = 'no alignment';
                                end

                                spikeInfo_writeOut = fullfile(microDir,micro_sess,manual_spike_folder,spikeInfoToPull(1).name);
                                %spikeInfo_writeOut = fullfile(microDir,micro_sess,'sorted/',spikeInfos(numSpikeInfo).name,spikeInfoToPull(1).name); %- jw hack

                                spikeInfo_writeOut = strrep(spikeInfo_writeOut,'.mat',['_aligned' CV_res_end '.mat']);
                                
                                spikeInfo.eegoffsetToMS = 1000/samprate_ecog; %********************* SJ CVCH

                                if CV_resampled
                                    spikeInfo.resampleInfo = ['CV ecog: Ecog pulses resampled from ' num2str(samprate_ecog) ' Hz to ' num2str(samprate_align) ' Hz for alignment.']; %SJ CRCH

                                    spikeInfo_withoutsuffix = strrep(spikeInfo_writeOut,CV_res_end,'');
                                    keyboard
                                    if exist(spikeInfo_withoutsuffix,'file')
                                        fprintf('%s\n','Aligned Spike Info prior to CRV fix exists. Deleting and replacing now.');
                                        keyboard
                                        delete(spikeInfo_withoutsuffix)
                                    end

                                end

                                if pulseSet==1
                                    adjustedSpikes = SPKtimestampsFixed; % Assign to empty adjustedSpikes
                                    adjustedSpikes_noise = SPKtimestampsFixed_noise;
                                else
                                    adjustedSpikes(apply_theseChans) = SPKtimestampsFixed(apply_theseChans); % For subsequent micro names, only apply to those channels
                                    adjustedSpikes_noise(apply_theseChans_noise) = SPKtimestampsFixed_noise(apply_theseChans_noise);
                                end
                                
                                thisNSP = lfpPulses_full{pulseSet,1};
                                % Check if you have any noise units
                                if ~isempty(spikeInfo.sessNoiseUnitID) && ~isempty(spikeInfo.timeStampNoise)

                                    %SJ: Do mask now? Or at end? After pulset
                                    %loop?

                                    micro_num_samp_orig = size(microLFP,2);
                                    micro_num_samp = size(signal_noreref_allSets{1,pulseSet}{1,2},2);


                                    %noiseUnits = adjustedSpikes_noise;
                                    
                                    
                                    % SJ: turns out we need the non-aligned
                                    % timestamps because they just need to
                                    % overlay in with the LFP, NOT macro!!

                                    chanNoise_timestamps = adjustedSpikes_noise(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP)&spikeInfo.sessNoiseUnitID.chanNoise>0)';
                                    %chanNoise_timestamps = spikeInfo.timeStampNoise(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP)&spikeInfo.sessNoiseUnitID.chanNoise>0)';
                                    
                                    thisNSP_arraynum = unique(spikeInfo.sessNoiseUnitID.arrayNum(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP)));
                                    keyboard %Make sure below works for each and every case!!
                                    %SJ fixing noise channel list to include full name before brackets- I don't know why 'thisNSP' would be included in the name- maybe sometimes it is?
                                    %chanNameList = regexp(spikeInfo.sessNoiseUnitID.ChanNameNew,[thisNSP '.*(?=\[)'],'match');
                                    %chanNameList = regexp(spikeInfo.sessNoiseUnitID.ChanNameNew,'.*(?=\[)','match');
                                    chanNameList = regexp(spikeInfo.sessNoiseUnitID.ChanNameNew(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP)),'.*(?=\[)','match');
                                    chanNameList = unique([chanNameList{:}]);

                                    maskArray = cell(1+numel(thisNSP_arraynum),3);
                                    maskArray{1,1} = ['chanNoise[' thisNSP ']'];

                                    if numel(thisNSP_arraynum) == 0
                                        keyboard % Test this
                                        maskArray{2,1} = 'globalNoise';
                                        maskArray{2,2} = [];
                                    elseif numel(thisNSP_arraynum) ~= numel(chanNameList)
                                        % Previously, this was the first check. BUT there are cases where we are on an NSP that has no noise units, but the other does, in which case those units will show up in chanNameList but will not match anything in thisNSP_array
                                        fprintf('%s\n',['ERROR!! For NSP ' thisNSP ', there are ' num2str(numel(thisNSP_arraynum)) ' array numbers and ' num2str(numel(chanNameList)) ' chanNameNews. These numbers should be the same!!']);
                                        keyboard
                                    end

                                    globalNoise_timestamps = cell(1,numel(thisNSP_arraynum));

                                    for aa = 1:numel(thisNSP_arraynum)
                                        %SJ: First find what chan name we are
                                        %working with for each arraynum
                                        thisArray = chanNameList{aa};

                                        maskArray{1+aa} = ['globalNoise[' thisArray ']'];

                                        globalNoise_timestamps{aa} = adjustedSpikes_noise(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP) & (contains(spikeInfo.sessNoiseUnitID.ChanNameNew,thisArray)) & (spikeInfo.sessNoiseUnitID.globalNoise > 0))';
                                        %globalNoise_timestamps{aa} = spikeInfo.timeStampNoise(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP) & (contains(spikeInfo.sessNoiseUnitID.ChanNameNew,thisArray)) & (spikeInfo.sessNoiseUnitID.globalNoise > 0))';

                                        %thisNSP_timestamps = adjustedSpikes_noise(strcmp(spikeInfo.sessNoiseUnitID.NSPsuffix,thisNSP));
                                    end

                                    noiseUnits = cat(2,{chanNoise_timestamps},globalNoise_timestamps);

                                    [timeSeries, mask_stats] = createNoiseMaskSPK(noiseUnits,micro_num_samp);
                                    
%                                     for aa2 = 1:numel(thisNSP_arraynum)+1
%                                         % chanNoise is always the first!!
%                                         maskArray{aa2,2} = [];
%                                         %timeSeries{aa};
%                                     end
                                    % chanNoise is always the first!
                                    maskArray(:,2) = timeSeries';
                                    maskArray(:,3) = mask_stats';
                                    
                                elseif isempty(spikeInfo.sessNoiseUnitID) && isempty(spikeInfo.timeStampNoise)
                                    maskArray = {};
                                else
                                    fprintf('%s\n','ERROR!! sessNoiseUnitID or timeStampNoise is empty but the other is not!! WHY???')
                                    keyboard
                                end
                                    
                                if pulseSet == 1
                                    noiseMask_applytheseChans = {};
%                                 elseif pulseSet == 2
%                                     if ischar(spikeInfo.noiseMask) %means that the first pulseset didnt have any spikes/mask?
%                                         spikeInfo.noiseMask = {};
%                                     end
                                end
                                
                                %SJ 12/21/2020: Changing because spikeInfo is re-loaded each pulseSet loop, so noisemask didn't get saved previously! This will do it
                                noiseMask_applytheseChans{pulseSet,1} = thisNSP;
                                noiseMask_applytheseChans{pulseSet,2} = lfpPulses_full{pulseSet,2};
                                noiseMask_applytheseChans{pulseSet,3}  = maskArray;
                                
%                                 spikeInfo.noiseMask{pulseSet,1} = thisNSP;
%                                 spikeInfo.noiseMask{pulseSet,2} = lfpPulses_full{pulseSet,2};
%                                 spikeInfo.noiseMask{pulseSet,3}  = maskArray;
                                
                                spikeInfo.noiseMask = noiseMask_applytheseChans;
                                
                                spikeInfo.timeStamp = adjustedSpikes; clear SPKtimestampsFixed;
                                spikeInfo.timeStampNoise = adjustedSpikes_noise; clear SPKtimestampsFixed_noise;
                                save(spikeInfo_writeOut,'spikeInfo','-v7.3');
                                %SJ: Shouldn't have to safe spikeInfo until after pulseSet loop since each pulseSet loop loads in the unaligned spikeInfo...

                            end
                        end

                        %keyboard 
                        %continue
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LFP    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    

                        %microLFP = lfpStruct.lfp{1,whichCellOfData{1,pulseSet}}; %SJ: changed because this
                        %may have been resampled if CV. Already assigned from before spikeinfo
                        microLFP_processed = microLFP;
                        maxSize = size(microLFP,2);
                        % insert samples here
                        % manipulateHere = round(linspace(1,maxSize,abs(driftApplyToMicro)));

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % Do all the below crap for processed or mountainsort (leaving commented for now)
                        % Call applyAlignCorrection
                        % Do we have global stuff in this lfp? (only in mountainsort)
    %                     for ll = 1:numel(A)
    % 
    %                     end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    %                     correctionInterval = floor((maxSize)/(abs(driftApplyToMicro)));
    %                     if correctionInterval<1000
    %                         fprintf('\n correction interval is usually hundreds of thousands, but here its %d',correctionInterval);
    %                         keyboard;
    %                     end
    %                     
    %                     %- indicies that will get a second copy (steps of correction interval, starting 1/2 way in correction interval
    %                     manipulateHere = zeros(1,(abs(driftApplyToMicro)));
    %                     for i=1:abs(driftApplyToMicro)
    %                         %-  put it in the middle of the correction interval because that makes sense
    %                         %       and so an extra sample isn't added to the end on even drifts
    %                         manipulateHere(i) = floor(correctionInterval/2)+1 + correctionInterval*(i-1);
    %                     end
    %                     
    %                     %%
    %                     if size(manipulateHere,2) == 0
    %                         keyboard %Does this happen?
    %                         %if there is actually nothing to remove/add
    %                         glob_sig_all_aligned = glob_sig_all;
    %                         glob_sig_good_aligned = glob_sig_good;
    %                         noreref_aligned = microLFP; clear microLFP
    %                         processed_aligned = microLFP_processed; clear microLFP_processed
    %                     else
    %                         
    %                         % if removing data:
    %                         if driftApplyToMicro < 0
    %                             for sample = 1:length(manipulateHere)
    %                                 % lets avg the to-be-removed sample with the next one,
    %                                 % and assign that to the next sample
    %                                 samplesToChangeAcrossChans = nanmean([microLFP(:,manipulateHere(sample)+1),microLFP(:,manipulateHere(sample))],2);
    %                                 microLFP(:,manipulateHere(sample)+1) = samplesToChangeAcrossChans;
    %                                 % do same for reref
    %                                 samplesToChangeAcrossChans = nanmean([microLFP_processed(:,manipulateHere(sample)+1),microLFP_processed(:,manipulateHere(sample))],2);
    %                                 microLFP_processed(:,manipulateHere(sample)+1) = samplesToChangeAcrossChans;
    %                                 % do same for glob_sig_all
    %                                 samplesToChangeAcrossChans = nanmean([glob_sig_all(manipulateHere(sample)+1),glob_sig_all(manipulateHere(sample))],2);
    %                                 glob_sig_all(manipulateHere(sample)+1) = samplesToChangeAcrossChans;
    %                                 % do same for glob_sig_good
    %                                 samplesToChangeAcrossChans = nanmean([glob_sig_good(manipulateHere(sample)+1),glob_sig_good(manipulateHere(sample))],2);
    %                                 glob_sig_good(manipulateHere(sample)+1) = samplesToChangeAcrossChans;
    %                             end
    %                             % now, lets remove all the samples
    %                             microLFP(:,manipulateHere) = [];
    %                             microLFP_processed(:,manipulateHere) = [];
    %                             glob_sig_all(manipulateHere) = [];
    %                             glob_sig_good(manipulateHere) = [];
    %                             
    %                             % and lets prepare this for output
    %                             glob_sig_all_aligned = glob_sig_all; clear glob_sig_all
    %                             glob_sig_good_aligned = glob_sig_good; clear glob_sig_good
    %                             noreref_aligned = microLFP; clear microLFP
    %                             processed_aligned = microLFP_processed; clear microLFP_processed
    %    
    %                         elseif driftApplyToMicro > 0 %adding data
    %                             
    %                             % prepare new series
    %                             noreref_aligned = NaN(size(microLFP,1),size(microLFP,2)+abs(driftApplyToMicro));
    %                             processed_aligned = NaN(size(microLFP_processed,1),size(microLFP_processed,2)+abs(driftApplyToMicro));
    %                             glob_sig_good_aligned = NaN(size(glob_sig_good,1)+abs(driftApplyToMicro),size(glob_sig_good,2));
    %                             glob_sig_all_aligned = NaN(size(glob_sig_all,1)+abs(driftApplyToMicro),size(glob_sig_all,2));
    %                             
    %                             
    %                             for sample = 1:length(manipulateHere)
    %                                 % lets avg the to-be-removed sample with the next one,
    %                                 % and assign that to the next sample
    %                                 samplesToChangeIntoNoreref = nanmean([microLFP(:,manipulateHere(sample)+1),microLFP(:,manipulateHere(sample))],2);
    %                                 % do same for reref
    %                                 samplesToChangeIntoProcessed = nanmean([microLFP_processed(:,manipulateHere(sample)+1),microLFP_processed(:,manipulateHere(sample))],2);
    %                                 % do same for glob_sig_all
    %                                 samplesToChangeIntoGlobSigAll = nanmean([glob_sig_all(manipulateHere(sample)+1),glob_sig_all(manipulateHere(sample))]);
    %                                 % do same for glob_sig_good
    %                                 samplesToChangeIntoGlobSigGood = nanmean([glob_sig_good(manipulateHere(sample)+1),glob_sig_good(manipulateHere(sample))]);
    %                                 
    %                                 % now, we need to carefully stitch
    %                                 % together real data with new insertions..
    %                                 if sample==1
    %                                     % insert the new, stretched data into the first slot (larger than the initiial data by 1 samples)
    %                                     noreref_aligned(:,1:manipulateHere(sample)+sample) = cat(2,microLFP(:,1:manipulateHere(sample)),samplesToChangeIntoNoreref);
    %                                     processed_aligned(:,1:manipulateHere(sample)+sample) = cat(2,microLFP_processed(:,1:manipulateHere(sample)),samplesToChangeIntoProcessed);
    %                                     glob_sig_all_aligned(1:manipulateHere(sample)+sample) = cat(1,glob_sig_all(1:manipulateHere(sample)),samplesToChangeIntoGlobSigAll);
    %                                     glob_sig_good(1:manipulateHere(sample)+sample) = cat(1,glob_sig_good(1:manipulateHere(sample)),samplesToChangeIntoGlobSigGood);
    %                                 elseif sample==length(manipulateHere)
    %                                     % if this is last sample, insert end data first
    %                                     noreref_aligned(:,manipulateHere(sample)+sample:size(noreref_aligned,2)) = microLFP(:,manipulateHere(sample):size(microLFP,2));
    %                                     processed_aligned(:,manipulateHere(sample)+sample:size(noreref_aligned,2)) = microLFP_processed(:,manipulateHere(sample):size(microLFP_processed,2));
    %                                     glob_sig_all_aligned(manipulateHere(sample)+sample:size(glob_sig_all_aligned,1)) = glob_sig_all(manipulateHere(sample):size(glob_sig_all,1));
    %                                     glob_sig_good_aligned(manipulateHere(sample)+sample:size(glob_sig_good_aligned,1)) = glob_sig_good(manipulateHere(sample):size(glob_sig_good,1));
    %                                     % then insert data up to the last manipulation point, along with newly created sample
    %                                     noreref_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,microLFP(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoNoreref);
    %                                     processed_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,microLFP_processed(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoProcessed);
    %                                     glob_sig_all_aligned(manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(1,glob_sig_all(manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoGlobSigAll);
    %                                     glob_sig_good_aligned(manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(1,glob_sig_good(manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoGlobSigGood);
    %                                 else
    %                                     % if this is a middle sample, insert it
    %                                     % along with newly created sample
    %                                     noreref_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,microLFP(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoNoreref);
    %                                     processed_aligned(:,manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(2,microLFP_processed(:,manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoProcessed);
    %                                     glob_sig_all_aligned(manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(1,glob_sig_all(manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoGlobSigAll);
    %                                     glob_sig_good_aligned(manipulateHere(sample-1)+sample:manipulateHere(sample)+sample) = cat(1,glob_sig_good(manipulateHere(sample-1)+1:manipulateHere(sample)),samplesToChangeIntoGlobSigGood);
    %                                 end
    %                             end
    %                         end
    %                     end
    %                     
    %                     % now, lets apply the offset correction (we wont bother
    %                     % trimming the end of the data)
    %                     if transforms{1,2}(1,1)<0
    %                         keyboard %fix the order of concat
    %                         noreref_aligned = cat(2,noreref_aligned,NaN(size(noreref_aligned,1),-transforms{1,2}(1,1)));
    %                         processed_aligned = cat(2,processed_aligned,NaN(size(processed_aligned,1),-transforms{1,2}(1,1)));
    %                         glob_sig_all_aligned = cat(1,glob_sig_all_aligned,NaN(-transforms{1,2}(1,1),1));
    %                         glob_sig_good_aligned = cat(1,glob_sig_good_aligned,NaN(-transforms{1,2}(1,1),1));
    %                     else
    %                         noreref_aligned(:,1:transforms{1,2}(1,1)) = [];
    %                         processed_aligned(:,1:transforms{1,2}(1,1)) = [];
    %                         glob_sig_all_aligned(1:transforms{1,2}(1,1)) = [];
    %                         glob_sig_good_aligned(1:transforms{1,2}(1,1)) = [];
    %                     end
    %                     
    %                     % edit lfp signals based on the transform output
    % %                     glob_sig_all_aligned = glob_sig_all_aligned - transforms{1,2}(1,1);
    % %                     glob_sig_all_aligned(end-transforms{1,2}(1,2):end) = [];
    %                     
    % %                     glob_sig_good_aligned = glob_sig_good_aligned - transforms{1,2}(1,1);
    % %                     glob_sig_good_aligned(end-transforms{1,2}(1,2):end) = [];
    %                     
    %                     % change in future - a single microLFP_processed.mat can contain multiple
    %                     % global_sig_alls and glob_sig_goods (for each array,
    %                     % i.e. 1-64 and 65-128)
    %                     
    %                     %SJ:  what is this even for?
    %                     pad_nan_mask = zeros(1,size(processed_aligned,2)); % SJ: incorrect fix for nan mask in processed
    %                     if transforms{1,2}(1,1) < 0
    %                         pad_nan_mask(1,1:abs(transforms{1,2}(1,1))) = 1;
    %                     end
    %                     
    %                     % need to aggregate this data across pulse sets such
    %                     % that its exactly the same. and need to apply a check
    %                     % that they're all the same length. inevitably, we
    %                     % should run into a case where they arent.. need to
    %                     % figure that out..
    %                     
    %                     % change data from double to int16 for storage size
    %                     noreref_aligned = int16(noreref_aligned);
    %                     processed_aligned = int16(processed_aligned);
    %                     glob_sig_all_aligned = int16(glob_sig_all_aligned);
    %                     glob_sig_good_aligned = int16(glob_sig_good_aligned);
    %                     
    %                     
    % %                     if pipelineType==1
    % %                         keyboard
    % %                     else
    % %                         lfp_aligned = noreref_aligned;
    % %                         lfpStruct.lfp{1,whichCellOfData{1,pulseSet}} = lfp_aligned;
    % %                         lfpStruct.alignedTo = ecog_sess;
    % %                         lfpStruct.alignmentChan = subject_pulseChannels{indexSubject,2};
    % %                         
    % %                     end


                    end
                    
                    %SJ: Try loop now..?
%                     for pulseSet2 = 1:size(lfpPulses_full,1)
%                         
%                         
%                         micro_num_samp = size(microLFP,2);
%                         
%                         
%                         noiseUnits = adjustedSpikes_noise;
%                         
%                         
%                         [timeSeries, mask_stats] = createNoiseMaskSPK(noiseUnits,lengthMicro,varargin);
                        
                        
                        
                        

                    

                    %% write data out into the subject directory
                    if pipelineType==1 %if this is a mountainsort file

                        keyboard

                        if ~isfolder(microDir); mkdir(microDir); end
                        if ~isfolder(fullfile(microDir,micro_sess)); mkdir(fullfile(microDir,micro_sess)); end
                        keyboard % need to add in CV_res_end like below
                        noreref_utah_writeOut = fullfile(microDir,micro_sess,'noreref_aligned.mat');
                        processed_utah_writeOut = fullfile(microDir,micro_sess,'processed_aligned.mat');

                        if spikeInfoFound==1
                            spikeInfo_writeOut = fullfile(microDir,micro_sess,'spikeInfo_aligned.mat');
                            save(spikeInfo_writeOut,'transforms','about_me_spikeInfo','alignedTo','alignmentChan','extractInfoStr','filter_string','metrics','pulses','referencing_info','sessDurSec','sessUniqueUnitID','sessStr','startTime_mstime','timeStamp_aligned','-v7.3');
                        end

                        lfp_aligned = noreref_aligned; clear noreref_aligned;
                        save(noreref_utah_writeOut,'transforms','pad_nan_mask','about_me_noreref','chan_names','channel_ranges','glob_sig_all_aligned','glob_sig_good_aligned','lfp_aligned','pulses','samplingFreq','sessStr','startTime_mstime','-v7.3');
                        lfp_aligned = processed_aligned; clear processed_aligned;
                        save(processed_utah_writeOut,'transforms','pad_nan_mask','about_me_processed','chan_names','channel_ranges','glob_sig_all_aligned','glob_sig_good_aligned','lfp_aligned','pulses','referencing_info','samplingFreq','startTime_mstime','-v7.3');

                    elseif pipelineType==2 %if this is a manually sorted file, write out the LFPs. we've already written out spikeInfos

                        if ~isfolder(microDir); mkdir(microDir); end
                        if ~isfolder(fullfile(microDir,micro_sess)); mkdir(fullfile(microDir,micro_sess)); end

                        noreref_utah_writeOut = fullfile(microDir,micro_sess,up_lfp_noreref,['lfp_aligned' CV_res_end '.mat']);
                        processed_utah_writeOut = fullfile(microDir,micro_sess,up_lfp_reref,['lfp_aligned' CV_res_end '.mat']); %Not ready for this yet -SJ
                        
                        lfpStruct.eegoffsetToMicro = samprate_align/samprate_ecog; %SJ: CVCH
                        
                        if CV_resampled
                            %lfpStruct.resampleInfo = ['CV ecog: Micro timeseries resampled from ' num2str(samprate_microLFP) ' Hz to ' num2str(samprate_ecog) ' Hz.'];
                            lfpStruct.resampleInfo = ['CV ecog: Ecog pulses resampled from ' num2str(samprate_ecog) ' Hz to ' num2str(samprate_align) ' Hz for alignment.']; %SJ CRCH
                            lfpwriteout_withoutsuffix = strrep(noreref_utah_writeOut,CV_res_end,'');
                            keyboard
                            if exist(lfpwriteout_withoutsuffix,'file')
                                fprintf('%s\n','Aligned lfp.mat prior to CRV fix exists. Deleting and replacing now.');
                                keyboard
                                delete(lfpwriteout_withoutsuffix)
                            end
                        end
                        save(noreref_utah_writeOut,'lfpStruct','-v7.3');
                        writetable(alignedFiles,alignment_xlsx_writeOut)
                        %lfp_aligned = processed_aligned; clear processed_aligned;
                        %save(processed_utah_writeOut,'tLFPdata','transforms','pad_nan_mask','lfp_aligned','pulses','-v7.3'); %Not ready for this yet -SJ
                    end
                end
            end
        end
    end % for exp = 1:length(behFolders)

    writetable(alignedFiles,alignment_xlsx_writeOut)
    %cell2csv(alignedFiles,alignment_xlsx_writeOut) %Why are we doing this again?
    fprintf('\n%s done\n',subj)
    
    %- clear the variables if looping over another subject... keep the last set of variables so "alignedFiles" can be returned
    if sub<length(subject)
        clearvars  -except subject subject_pulseChannels ecogDir
    end
    
end %- subject loop
end % end function
