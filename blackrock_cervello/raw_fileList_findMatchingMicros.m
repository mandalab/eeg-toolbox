function [outputArg1] = raw_fileList_findMatchingMicros(subj, ecogDir)

if ~exist(sprintf('%s/%s/behavioral/',ecogDir,subj))
    fprintf('\ninput directory doesnt exist for %s.. aborting\n',subj)
    return
end


temp_dir = dir(sprintf('%s/%s/raw/',ecogDir,subj));
rawFileFound = intersect((find(contains({temp_dir.name},'rawFileList'))),(find(contains({temp_dir.name},'xlsx')))); 
rawFileFound = setdiff(rawFileFound,find(contains({temp_dir.name},'match')));
if isempty(rawFileFound)
    fprintf('\nraw file list doesnt exist for %s.. aborting\n',subj)
    return
end

rawFileList = readtable(fullfile(temp_dir(rawFileFound).folder,temp_dir(rawFileFound).name));

containsBehavior = zeros(size(rawFileList,1),1);

start_end_times = NaN(size(rawFileList,1),2);

% for each of these files, lets go through and assign a start and end time
for file=1:size(rawFileList,1)
timestamp = cellstr(table2cell(rawFileList(file,{'folderName'})));
duration = cell2mat(table2cell(rawFileList(file,{'dur_minutes'})));
micro_start =  datenum(timestamp,'yymmdd_HHMM');
micro_end = micro_start + (round(duration)/60/24); 
start_end_times(file,1) = micro_start;
start_end_times(file,2) = micro_end;
end


% now, lets go through each behavioral events.mat, find whether it overlaps with any micro file, and if so,
% assign that micro file a 1

behavDir = fullfile(ecogDir,subj,'/behavioral/');
behavFolders =  dir(behavDir);
behavFolders = behavFolders(find([behavFolders.isdir]==1));
behavFolders(find(contains({behavFolders.name},'cant'))) = [];
behavFolders(find(contains({behavFolders.name},'.'))) = [];

for exp = 1:length(behavFolders)
    fprintf('%d\n',exp)
    sessionDir = fullfile(behavDir,behavFolders(exp).name,'/');
    sessionFolders = dir(sessionDir);
    sessionFolders = sessionFolders(find(contains({sessionFolders.name},'session'))); %- rule is the "cant" string preceeds
    
    for sess = 1:length(sessionFolders)
        
        
        %if contains(fullfile(sessionDir,sessionFolders(sess).name),'notReal','IgnoreCase',true) || contains(fullfile(sessionDir,sessionFolders(sess).name),'cantAlign','IgnoreCase',true) || contains(fullfile(sessionDir,sessionFolders(sess).name),'empty','IgnoreCase',true) || contains(fullfile(sessionDir,sessionFolders(sess).name),'WONT','IgnoreCase',true)
        if ~strcmp(sessionFolders(sess).name(1:7),'session')
            continue
        end
        temp = strsplit(fullfile(sessionDir,sessionFolders(sess).name),'/');
        if strcmp(temp{1,length(temp)}(1),'_')
            continue
        end
        
        sess_events = load(fullfile(sessionDir,sessionFolders(sess).name,'/events.mat'));
        sess_events = sess_events.events;
        
        length_sess = sess_events(end).mstime - sess_events(1).mstime;
        length_sess = length_sess / 1000 / 60; % in minutes
        
        eegFile = sess_events(1).eegfile; % i think 1st event should always contain accurate pointer..
        eegFileParts = strsplit(eegFile,'/');
        date_start = eegFileParts{1,end};
        
        ecog_start =  datenum(date_start,'yymmdd_HHMM');
        behavior_start = ecog_start + ((double(sess_events(1).eegoffset)/1000/60/60/24));
        behavior_end = ecog_start + ((double(sess_events(end).eegoffset)/1000/60/60/24));
        
        
        % did either the first or last event occur within any micro file
        match_found = 0;
        tic
        while match_found==0
            for file=1:size(start_end_times,1)
                start_file = start_end_times(file,1);
                end_file = start_end_times(file,2);
                if (behavior_start > start_file) && (behavior_start < end_file)
                    match_found = 1;
                    break
                end
                if (behavior_end > start_file) && (behavior_end < end_file)
                    match_found = 1;
                    break
                end
            end
            b = toc;
            if b>10; break; end
        end
        
        % dont overwrite a match with a non-match
        if containsBehavior(file) ~= 1
        containsBehavior(file) = match_found;
        end
        
    end
end



rawFileList = addvars(rawFileList,containsBehavior,'Before','task');
writeInto = fullfile(temp_dir(rawFileFound).folder,temp_dir(rawFileFound).name);
writeInto = strrep(writeInto,'.xlsx','_matchWithMicro.xlsx');
writetable(rawFileList,writeInto);


outputArg1 = [];


end

