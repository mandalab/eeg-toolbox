
ecogDir = '/Volumes/56PROC/eeg_processing/.update/eeg';

sub_list = getDirNamesRegexp(ecogDir,'NIH.*');
micro_only = true;
rr = 0;
Tmat = {};
for ss = 1:numel(sub_list)
    subj = sub_list{ss};
    rawEcogDir = fullfile(ecogDir,subj,'raw');
    microEcogDir = fullfile(ecogDir,subj,'micro');
    if exist(microEcogDir,'dir')
        micro_align_file = fullfile(ecogDir,subj,'micro','manualsort','alignmentToSpikes_manual.xlsx');
        if exist(micro_align_file,'file')
            alignTable = readtableSafe(micro_align_file);
            [rawSessions, rawSessionsisStim] = getRawSessions(rawEcogDir);
            for sess = 1:numel(rawSessions)
                ecog_sess = rawSessions{sess};
                if rawSessionsisStim(find(strcmp(rawSessions,ecog_sess))) % STIM_MAP folder
                    rawSession_content = getDirNamesRegexp(fullfile(rawEcogDir,'STIM_MAP',ecog_sess),'.*');
                else
                    rawSession_content = getDirNamesRegexp(fullfile(rawEcogDir,ecog_sess),'.*'); % SJ fixed
                end
                if any(contains(rawSession_content,'.TRC')) % SJ fixed
                    if ~strcmp(unique(alignTable.micro_sess(strcmp(alignTable.ecog_sess,ecog_sess))),'no match')
                        % This means there is a micro session that goes with this TRC session
                        rr = rr+1;
                        if rr == 1;
                            Tmat = [{subj},{ecog_sess},unique(alignTable.micro_sess(strcmp(alignTable.ecog_sess,ecog_sess)))];
                        else
                            Tmat = [Tmat ; [{subj},{ecog_sess},unique(alignTable.micro_sess(strcmp(alignTable.ecog_sess,ecog_sess)))]];
                        end
                    else
                    end
                end
            end
        else
            keyboard
            continue
        end
        %rawSessions = getDirNamesRegexp(rawEcogDir,'^\d{6}_\d{4}.*');

    else
        continue
    end

end

T = cell2table(Tmat,'VariableNames',{'subject','ecog_sess','micro_sess'});
