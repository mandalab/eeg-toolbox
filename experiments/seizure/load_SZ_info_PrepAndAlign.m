function events = load_SZ_info_PrepAndAlign(subject,rootEEGdir)

% Use main function
events = load_SZ_info(subject);

% Convert to struct
events = table2struct(events);

% go into noreref to grab names
% norerefDir = dir(sprintf('%s%s/eeg.noreref',rootEEGdir,subject));
norerefDir = dir(fullfile(rootEEGdir,subject,'eeg.noreref'));

% Grab the timestamps
iNoreref = 1;
for iTimestamp = 1:length(norerefDir)  
    if ~contains(norerefDir(iTimestamp).name,'.')
        norerefTimes{iNoreref} = norerefDir(iTimestamp).name; 
        iNoreref = iNoreref + 1;
    end
end

%JIC all of these sprintf statements should be replaced with fullfile so that things don't break due to missing '/'
% Grab the directories for each of the raws
for iTimestamp = 1:length(norerefTimes)
    
    % See which folder it's in
    if exist(sprintf('%s%s/raw/%s',rootEEGdir,subject,norerefTimes{iTimestamp}))
        rawStrings{iTimestamp} = sprintf('raw/%s/',norerefTimes{iTimestamp});
    elseif exist(sprintf('%s%s/raw/SEIZURE_BASELINE/%s',rootEEGdir,subject,norerefTimes{iTimestamp}))
        rawStrings{iTimestamp} = sprintf('raw/SEIZURE_BASELINE/%s/',norerefTimes{iTimestamp});
    elseif exist(sprintf('%s%s/raw/STIM_MAP/%s',rootEEGdir,subject,norerefTimes{iTimestamp}))
        rawStrings{iTimestamp} = sprintf('raw/STIM_MAP/%s/',norerefTimes{iTimestamp});
    end
    
end

% Loop through each raw and grab the sample rate in order to replace the
% time fields to be in the ms scale
for iEvent = 1:length(events)
    
    % Loop through each noreref folder to find the one that has this raw
    iTimestamp = 0;
    
    % Set flag to zero
    thisRaw = 0;
    
    % Start while
    while(thisRaw == 0)
        
        % Increment index
        iTimestamp = iTimestamp + 1;
 
        % Grab the directory listing in raw to see if it's the one want
        rawFiles = dir(sprintf('%s%s/%s',rootEEGdir,subject,rawStrings{iTimestamp}));

        % Go through each one to see if the one we want is there
        for iRaw = 1:length(rawFiles)
            
            % Compare the names, set flag if so
            if contains(rawFiles(iRaw).name,events(iEvent).rawFile)
                thisRaw = 1;
            end
            
        end
        
    end
    
    % Grab the params file
    fileID = fopen(sprintf('%s%s/eeg.noreref/%s/params.txt',rootEEGdir,subject,norerefTimes{iTimestamp}),'r');
    params = fgetl(fileID);
    fclose(fileID);
    
    % Split to grab the rate
    rateStr = strsplit(params,' ');
    Fs = str2num(rateStr{2});    

    % Add the new fields
    events(iEvent).eegoffset = events(iEvent).eegoffset_s * Fs;
    events(iEvent).EEG_Change_ms = events(iEvent).EEG_Change_s * Fs;
    events(iEvent).Ictal_Onset_ms = events(iEvent).Ictal_Onset_s * Fs;
    events(iEvent).Clinical_Onset_ms = events(iEvent).Clinical_Onset_s * Fs;
    events(iEvent).Termination_ms = events(iEvent).Termination_s * Fs;
    events(iEvent).eegfile = sprintf('%s%s/eeg.noreref/%s/',rootEEGdir,subject,norerefTimes{iTimestamp});
    events(iEvent).date.Hour = str2double(norerefTimes{iTimestamp}(end-3:end-2));
    events(iEvent).date.Minute = str2double(norerefTimes{iTimestamp}(end-1:end));
    events(iEvent).date = events(iEvent).date + seconds(events(iEvent).EEG_Change_s);
    
end

% Remove the _s fields
events = rmfield(events, 'eegoffset_s');
events = rmfield(events, 'EEG_Change_s');
events = rmfield(events, 'Ictal_Onset_s');
events = rmfield(events, 'Clinical_Onset_s');
events = rmfield(events, 'Termination_s');

end