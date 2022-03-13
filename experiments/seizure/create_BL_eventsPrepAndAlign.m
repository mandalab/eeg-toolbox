function [events, table_fig] = create_BL_eventsPrepAndAlign(subject,sessionDir,rootEEGdir,sessNum)

% Load the csv
str = sprintf('%s/baseline.csv',sessionDir);
baseline = readtable(str,'Delimiter',',');

% Create variables for create_BL_events
stateVec = zeros(height(baseline),1);
timestamps = cell(height(baseline),1);
for iBaseline = 1:height(baseline)
    if strcmp(baseline.State(iBaseline),'awake')
        stateVec(iBaseline) = 1;
    end
    
    timestamps(iBaseline) = baseline.x___TimeStamp(iBaseline);
end

% go into noreref to grab names
norerefDir = dir(sprintf('%s%s/eeg.noreref',rootEEGdir,subject));

% Grab the timestamps
for iTimestamp = 1:length(norerefDir)   
    norerefTimes{iTimestamp} = norerefDir(iTimestamp).name; 
end

% Go through each timestamp to make sure that we have it in noreref
for iTimestamp = 1:length(timestamps)
    
    % Compare
    if ~sum(strcmp(timestamps{iTimestamp},norerefTimes))
        
        % If it wasn't there, throw an error and tell the user to grab it
        fprintf('\n%s was not found in eeg.noreref. Please grab it from FRNU 56 or by using eegRawServerGrabFiles.\n',timestamps{iTimestamp})
        error('MISSING RAW FILE! PLEASE ADD IT AND TRY AGAIN!')
        
    end
end

% Pass variables into create_BL_events
[events, table_fig] = create_BL_events(subject,timestamps,stateVec,rootEEGdir,sessNum);
end