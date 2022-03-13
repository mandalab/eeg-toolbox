function events = extractMultiModal(stimKeyDir,sessLogFile, subject, sessionName)
% Function for extracting behavioral data from multiModal
% Inputs:
% sessLogFile: directory of session.log file
% subject: e.g. 'NIH036'
% sessionName: e.g. 'session_0'

% keyboard

trialIdentifiers = {'IMAGE_ON_SCREEN','TEXT_ON_SCREEN','SOUND_PLAYING'};
% Get the stim modalities for the texts
stimKeyFid = fopen(stimKeyDir,'r');
textStimCat = {};
moreToRead = true;
while moreToRead
    thisLine = fgetl(stimKeyFid);
    if ~ischar(thisLine); moreToRead = false; end
    if moreToRead
        temp = textscan(thisLine,'%s',20);
        temp = temp{1}';
        tempStrArray = cell(1,20-length(temp));
        tempStrArray(:) = {''};
        temp = cat(2,temp,tempStrArray);
        textStimCat = cat(1,textStimCat,temp);
    end
end
fclose(stimKeyFid);

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1
    iKeep=1:find(diff(strNumeric)>1,1,'first');
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep);
end
sessionNum = str2num( sessionName(strNumeric) );
if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...


%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 1;
moreToRead = true;
while moreToRead
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); moreToRead = false; end
    
    if moreToRead
        % Check whether this line is for a trial
        xTOT = textscan(thisLine,'%s',20);
        thisLineStr = {};
        for k = 1:length(xTOT)
            if ~isempty(xTOT{k})
                thisLineStr = cat(2,thisLineStr,xTOT{k});
            end
        end
        for trialMod = trialIdentifiers
            trialMod = char(trialMod);
            i = find(not(cellfun('isempty', strfind(thisLineStr, trialMod))));
            if ~isempty(i)
                mstime = str2double(thisLineStr{1});
                soundDuration = nan;
                
                % Get the stimulus type
                if strcmp(trialMod,'TEXT_ON_SCREEN')
                    modality = 'txt';
                elseif strcmp(trialMod,'IMAGE_ON_SCREEN')
                    modality = 'img';
                elseif strcmp(trialMod,'SOUND_PLAYING')
                    modality = 'aud';
                end
                % Get the actual stimulus
                stimID = thisLineStr{i+1};
                if ~strcmp(trialMod,'TEXT_ON_SCREEN')
                    if strcmp(trialMod,'IMAGE_ON_SCREEN')
                        fileType = strfind(stimID,'.jpg');
                    elseif strcmp(trialMod,'SOUND_PLAYING')
                        fileType = strfind(stimID,'.wav');
                        % Get the sound duration
                        if strcmp(thisLineStr{i+2},'SOUND_DUR') && length(thisLineStr) >= i+3
                            soundDuration = str2double(thisLineStr{i+3});
                        else
                            error('Sound duration missing in session log!');
                        end
                    end
                    slashes = strfind(stimID,'/');
                    stimID = stimID(slashes(end)+1:fileType-1);
                end
                
                % Get stimulus modality
                thisCategory = '';
                if strcmp(trialMod,'TEXT_ON_SCREEN')
                    stimID = lower(stimID);
                    for thisCat = 1:size(textStimCat,1)
                        if ismember(stimID,textStimCat(thisCat,:))
                            thisCategory = textStimCat{thisCat,1};
                        end
                    end
                else
                    % Get rid of number portion
                    thisCategory = stimID;
                    number_i = [];
                    for ii = 1:length(thisCategory)
                        if ~isnan(str2double(thisCategory(ii))) && isreal(str2double(thisCategory(ii)))
                            number_i = cat(2,number_i,ii);
                        end
                    thisCategory(number_i) = [];
                    end
                end
                
                if strcmp(stimID,'audioIcon')
                    break;
                end
                
                clear thisEvent
                thisEvent.experiment        = 'multiModal';
                thisEvent.subject           = subject;
                thisEvent.sessionName       = sessionName;
                thisEvent.sessionNum        = sessionNum;
                thisEvent.mstime            = mstime;
                thisEvent.modality          = modality;
                thisEvent.stimID            = stimID;
                thisEvent.category          = thisCategory;
                thisEvent.soundDuration     = soundDuration;
                
                if (index==1)
                    events        = thisEvent; %- before events defined must convert to structure
                else
                    events(index) = thisEvent;
                end
                
                index = index+1;
            end
        end
    end
end
fclose(fid);  % close session.log


