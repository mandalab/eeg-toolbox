function [events] = paRemap_ExtractEvents_NIH028(sessLogFile, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from paRemap %%%%
%
%   ** NIH028 Version **
%       NIH028 and NIH029 used pssychopy, and the code is throwing errors
%       however, we have old events.mat for them, they just haven't been annotated
%       so the events aren't marked with responseTime or isCorrect
%       this script modifies the *existing* events with these fields
%
%
%%%%%   create an event for every presented word and every response (no words from the training sections should be included)
%
%
%   extraction designed for paRemap v4.2, which is implemented in pyEPL and was used for NIH030 and beyond
%           (earlier versions were implemented in psychoPy) and were used for NIH028 and 029... those session logs will need some tweaking to use with this extraction
%
%   training section NOT saved to session log... for earlier versions this MUST BE DELETED from the session log
%
%
%   NOTE about different versions run with each subject:  
%     NIH28 and NIH029 --- original version of remap in psychopy
%     NIH030 - NIH032  --- separate wav for EVERY word
%     NIH033 and above --- single wav per block
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DBG = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncomment following lines to directly run script
%clear all

rootEEGdir  = '/Users/trottams/eeg/';
%rootEEGdir  = '/Volumes/Shares/FRNU/dataWorking/eeg';
subject     = 'NIH028';   % EEG002  NIH016
sessionName = 'session_0_2015_Feb_26_2010';

sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/paRemap',sessionName);
sessLogFile = fullfileEEG(sessionDir,'session.log');
%eventFile   = fullfileEEG(sessionDir,'events.mat');
%priorEvents = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nOn session: %s', sessionName);

%- For NIH028 and NIH029, copy the session log and add the correct msoffset
sessFolderPath  = sessLogFile(1:strfind(sessLogFile,sessionName)+length(sessionName));
paRemap2Sesslog = [sessFolderPath 'paRemap2_session.log'];
    
if exist(paRemap2Sesslog,'file') % ~exist(sessLogFile,'file')         
    
    dateStrPsycho = sessionName(strfind(sessionName,'2015'):end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert psycho date/time into pyEPL date/time which comes from javaSDF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- imports required to convert "epoch time" saved by pyepl into real date
    import java.lang.System;
    import java.text.SimpleDateFormat;
    import java.util.Date;
    javaSDF = SimpleDateFormat('MM/dd/yyyy HH:mm:ss.SS');  %this java object is used for mstime to date conversion
    
    %- grab the date from the excell file
    dateNumPsycho  = datenum(dateStrPsycho, 'yyyy_mmm_dd_HHMM');
    dateStrPsycho2 = datestr(dateNumPsycho, 'mm/dd/yy HH:MM PM');  %%- attempting to match format of mac info
    
    %- convert matlab datenum into milisecond number used by javaSDF
    % javaSDF.format(Date(0))              --> '12/31/1969 19:00:00.00'
    % javaSDF.format(Date(60000))          --> '12/31/1969 19:01:00.00'     % 60000 is increment of 1 minute in javatime (javatime is in miliseconds, starting at 12/31/1969  1 min = 60000 microsec)
    % javaSDF.format(Date(1424183683378))  --> '02/17/2015 09:34:43.378'    % example mstime from pyepl session.log
    dateNum0java   = datenum(char(cell(javaSDF.format(Date(0)))));     % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    dayJava        = 24 * 60 * 60 * 1000;                                 % number of miliseconds in a day
    dayMatlab      = datenum('01/02/01 1:00')-datenum('01/01/01 1:00');   % number of days in a matlab date num (should be exactly 1)
    daysToAdd      = (dateNumPsycho-dateNum0java)/dayMatlab;
    msStartPyEPL   = round(dayJava*daysToAdd);
    dateNumMSstart = datenum(char(cell(javaSDF.format(Date(msStartPyEPL)))));  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    dateStrMSstart = datestr(dateNumMSstart, 'mm/dd/yy HH:MM PM');  % this magic conversion relies on java 'Date' and 'SimpleDateFormat' imports at the top of the page
    %fprintf('\n >> converting to pyepl time reference: msoffset %d = %s  (should match psychoPy xls date = %s) << ', msStartPyEPL, dateStrMSstart, dateStrPsycho2);  %- uncomment to confirm date conversion is working
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %- copy the psychopy session log with msoffset shifted based on folder name/date
    fidRead  = fopen(paRemap2Sesslog,'r');
    fidWrite = fopen(sessLogFile,'w+'); 
    while true
        thisLine = fgetl(fidRead);
        if ~ischar(thisLine); break; end
        [msTime,pos] = textscan(thisLine,'%f',1);
        fprintf(fidWrite,'%d%s\n',msTime{1}+msStartPyEPL,thisLine(pos:end));
    end
    fclose(fidWrite);    
    fclose(fidRead);  
      
    %- copy the psychopy session log
    paRemap2eeglog = [sessFolderPath 'paRemap2_eeg.eeglog'];
    standardeeglog = [sessFolderPath 'eeg.eeglog'];

    fidRead  = fopen(paRemap2eeglog,'r'); 
    fidWrite = fopen(standardeeglog,'w+'); 
    while true
        thisLine = fgetl(fidRead);
        if ~ischar(thisLine); break; end
        [msTime,pos] = textscan(thisLine,'%f',1);
        fprintf(fidWrite,'%d%s\n',msTime{1}+msStartPyEPL,thisLine(pos:end));
    end
    fclose(fidWrite);    
    fclose(fidRead);  
    fprintf('\n New copies of session.log and eeg.eeglog were created in %s ',  sessFolderPath);
     
    %- Old version, just make the sure copy the first number copy the psychopy session log
    %[SUCCESS,MESSAGE,MESSAGEID] = copyfile(paRemap2Sesslog,sessLogFile);
    
    %[SUCCESS,MESSAGE,MESSAGEID] = copyfile(paRemap2eeglog,standardeeglog);
    
    %fprintf('\n\n Update the new copy of session.log and eeg.eeglog in %s \n --->  REPLACE FIRST TIME POINT WITH % s\n HIT return when done or break with shift F5.',  sessFolderPath, num2str(msStartPyEPL));
    %keyboard
end

fid    = fopen(sessLogFile,'r'); %if loops check to see that sessFile exists
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
name_temp = strsplit(sessionName, '_');
name_temp = name_temp{2};
strNumeric = find( name_temp >= '0' & name_temp <= '9');
if max(diff(strNumeric))>1, iKeep=[1:find(diff(strNumeric)>1,1,'first')]; 
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', ...
        sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep))); strNumeric=strNumeric(iKeep); 
end
sessionNum = str2num( name_temp(strNumeric) );               
if isempty(sessionNum)
    fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  
end %shouldn't need this catch...



%- Read session.log line-by-line and convert to events structure
events      = [];
index       = 0;
probeFound = 0; %- initialize here... 
annFromProbeOn = 0;
type = 0;
isPractice = 0;
while true
    thisLine            = fgetl(fid); % get new line
    if ~ischar(thisLine); break; end  % reached EOF
    
    if DBG, disp(thisLine); end
    
    
    %- Generic text scan to get time, offset, and type
    xTOT                = textscan(thisLine,'%f%d%s');
    msoffset            = xTOT{2}(1);
    eventType           = xTOT{3}{1};
    
    %- default Parameters (details will be filled out/altered based on type)
    experiment          = 'paRemap';
    subject             = subject   ;
    sessionName         = sessionName ;
    sessionNum          = sessionNum  ;  %- store in state var so all events are assigned a sessionNum
    msoffset            = msoffset  ;
%     isCorrect           = nan;
    
    %%- type stores the third column of session.log and loops through row
    %%- by row
    switch eventType
        %%- This case is started by every recording session and ends with
        %%- REC_OFF
        % STORES: - recStart: the absolute time in milliseconds when recording
        % starts
        % - wordType: number from -1, [1,5] which represents the word that
        % was vocalized (or absence of vocalization)
        % - timeAfterRec: the absolute time after recording start, of each word
        % vocalization
        % - vocalizedWord: cell array storing the actual word vocalized
        % - annindex: the number of elements in ann file, for the vocalized
        % words
        case {'REC_ON'} % get name of .ann file and extract the contents
            annindex    = 0;   % initialize going through annotation file index
            recstartTOT = textscan(thisLine,'%f%d%s%s'); % grab the row with REC_ON in session.log file
            
            %- for NIH031 and 32 use annFromProbeOn
            if annFromProbeOn, annFileName = annFromProbeOn;
            else               annFileName = recstartTOT{4}{1};   end
            
            annfilename = fullfile(sessionDir, strcat(annFileName, '.ann')); % get the .ann file path
            recStart    = recstartTOT{1}(1);  % get this block's recording start time
           
            %%
%             if ~strcmp(sessionName, 'session_0b') %%- ONLY NEEDED FOR SESSIONS WITH SEPARATED SESSION (AS A/B)
%                 %%- Error check on opening up annotation .ann file
%                 if (annfid==-1), error('%s.ann not found: \n Exiting.',xTOT{4}(1)); end
%             end
            %%
            
            %disp(['On this file: ' annfilename])
            [annfid, errmsg] = fopen(annfilename, 'r');    % open file to read
            
            % if there does not exist an annotation file yet, or there is
            % an error reading it
            if ~isempty(errmsg)
                % set the response time=0,   isCorrect=-1, vocalizedWord = Nan  -- original version
                % set the response time=nan, isCorrect=-1, vocalizedWord = 'NoAnnFile'  -- JW version
                responseTime  = nan;
                isCorrect     = -1;
                vocalizedWord = 'NoAnnFile'; 
                
                annFileExists = 0;
            else
                annFileExists = 1;
            end
            
            % get lines until reached annotated data
            while true && annFileExists
                try
                    tempLine = fgetl(annfid); % get new line
                catch(e)
                    disp(e)
                    annfilename
                end
                if isempty(tempLine), break; end  % reached EOF
            end
            clear tempTot tempLine
            
            % read rest of annotation
            while true && annFileExists
                line = fgetl(annfid);
                if ~ischar(line); break; end  % reached EOF
                %- Generic text scan to get time, offset, and type
                annTOT                = textscan(line,'%f%d%s');
                
                % add to list these annotation data
                if annindex == 0
                    timeAfterRec          = annTOT{1}(1) + recStart;   %- must add (1) because numbers after the string in the above line cause overflow to first %f
                    wordType              = annTOT{2}(1);
                    vocalizedWord         = {annTOT{3}{1}};
                    annindex = annindex + 1;
                else
                    timeAfterRec          = [timeAfterRec; annTOT{1}(1) + recStart];   %- must add (1) because numbers after the string in the above line cause overflow to first %f
                    wordType              = [wordType; annTOT{2}(1)];
                    vocalizedWord         = [vocalizedWord; annTOT{3}{1}];
                    annindex = annindex + 1;
                end
            end
            
            % store total number of annotations in .ann file
            TOTAL_ANN = annindex;
            
            % variable to help in finding response words/times
            annwordindex = 1; % loop through the annotated words 1:TOTAL_ANN
            %probeFound = 0;  %- JW is commenting this out here... I think this only works for NIH33+... JW's method should work for all though
            
            if annFileExists
                fclose(annfid); % close annotation file
            end
        case {'BLOCK_0', 'BLOCK_1', 'BLOCK_2', 'BLOCK_3', 'BLOCK_4', 'BLOCK_5'}
            blockTOT = textscan(thisLine, '%f%d%s%s%s'); 
            if strcmp(blockTOT{5}(1), 'TEST')
                isPractice = 0;
                blocknumber = eventType;             % store the block number
                miniblocknumber = blockTOT{4}{1};   % store miniblocknumber
                
                %%- modify text to only extract the numbers
                blocknumber = strsplit(blocknumber, '_');
                blocknumber = blocknumber{2};
                miniblocknumber = strsplit(miniblocknumber, '_');
                miniblocknumber = miniblocknumber{2};
                
                newminiblock = 1; % start new counting of pairIndices within a miniblock
            else
                isPractice = 1;
            end
        case {'FIXATION_ON'}
            fixationOnTime = xTOT{1}(1);
        case {'FIXATION_OFF'}
            fixationOffTime = xTOT{1}(1);
        case {'PROBEWORD_ON'}
            if isPractice, continue; end
                
            isProbe    = 1;
            probeTOT   = textscan(thisLine,'%f%d%s%s%s%s%s'); % grab the row of data for probeWordOn
            probeWord  = probeTOT{4}{1};  % store probeWord
            targetWord = probeTOT{6}{1};  % store targetWord
            mstime     = probeTOT{1}(1);  % store absolute time of this row
            type       = probeTOT{3}{1};
            if length(probeTOT{7})>0, annFromProbeOn = probeTOT{7}{1}; else annFromProbeOn=0; end

            % get the pairIndex of this word pair 
            if newminiblock
                newminiblock = 0; % started a new miniblock index count
                pairIndex = 1;
            else
                pairIndex = pairIndex + 1;
            end
            probeFound = 1; % set the flag to make sure probewordon comes first
        case {'MATCHWORD_ON'}
            if isPractice, continue; end
            
            matchOnTime = xTOT{1}(1) - mstime; % get the absolute mstime of matchWord coming on
            
            if annFileExists % make sure there is an annotation file
                %%- When matchword comes on, should have vocalized...
                % determine timerange words can occur from 
                % probewordon -> matchword on (0-timeRange)
                timeRange = matchOnTime;                          % time Range the word can come on
                timeVocalizations = timeAfterRec(:) - mstime;              % convert all vocalizations wrt to the mstime
                validIndices = find(timeVocalizations > 0 & timeVocalizations < timeRange); % find valid indices of word by response times

                %%- found one word that was correct, no other vocalizations
                if length(validIndices) == 1, 
                    if strcmp(targetWord, vocalizedWord{validIndices})
                        isCorrect = 1;

                        % LOG THE EVENT FIELDS and increment index through ann file
                        responseTime = timeVocalizations(validIndices); % responseTime
                        responseWord = vocalizedWord{validIndices};
                    else % incorrect word vocalized
                        isCorrect = 0;
                        responseTime = timeVocalizations(validIndices); % responseTime
                        responseWord = vocalizedWord{validIndices};
                    end
                elseif isempty(validIndices), %%- no word response in this frame period
                    isCorrect = 0;
                    responseTime = 0;
                    responseWord = 'none';
                else %%- more then 1 word vocalization found within time frame
                    isCorrect = 0; % defined since they vocalized more then 1 word
                    if strcmp(targetWord, vocalizedWord{validIndices(1)}) % first try was correct
                        responseTime = timeVocalizations(validIndices(1));
                        responseWord = vocalizedWord{validIndices(1)};
                    elseif strcmp(targetWord, vocalizedWord{validIndices(end)}) % last vocalization was correct
                        responseTime = timeVocalizations(validIndices(end));
                        responseWord = vocalizedWord{validIndices(end)};
                    else % no vocalization was correct, or it was in the middle
                        responseTime = timeVocalizations(validIndices(1)); % get the first response word
                        responseWord = vocalizedWord{validIndices(1)}; % get the first vocalized word
                    end
                end
            else
                % set the response time=0, isCorrect=-1, vocalizedWord = Nan
                responseTime = nan;
                isCorrect    = -1;
                responseWord = 'NoAnnFile';
            end
        case {'PROBEWORD_OFF'} % SAME TIME MATCHWORD TURNS OFF
            if isPractice, continue; end
            probeOffTime = xTOT{1}(1) - mstime; % get the mstime of this line
                        
            if (probeFound == 1),  %- JW's not sure why this is here, but added probeFound=0 to here in an attempt to make the code compatible with <NIH033
                index = index+1; 
                probeFound = 0;  
            else
                fprintf('\n ERROR: probeword_off without probeword_on'); keyboard; 
            end
    end
    
    if isPractice, continue; end
    
    %%%%%% CURRENTLY ONLY APPENDING EVENTS WHEN INDEX ++ed by PROBEWORD ON
    %- asign values to events array
    if index>length(events),
        
        % just making sure all the times are in order
%         if matchOnTime < mstime || probeOffTime < matchOnTime ...
%             || fixationOnTime > probeOffTime ...
%             || fixationOffTime < fixationOnTime
%             disp('error in times')
%         end
        
        
        %- create dummy event structure that is upddated below based on type
        clear thisEvent
        thisEvent.experiment        = experiment  ;
        thisEvent.subject           = subject     ;
        thisEvent.sessionName       = sessionName ;
        thisEvent.sessionNum        = sessionNum  ;   % store in state var so all events are assigned a sessionNum  %% JW updated 2/2015
        thisEvent.type              = type        ;
        thisEvent.msoffset          = msoffset    ;
        thisEvent.mstime            = mstime      ;
        
        thisEvent.pairIndex         = pairIndex;
        thisEvent.probeWord         = probeWord     ;
        thisEvent.targetWord        = targetWord    ;
        thisEvent.isCorrect         = isCorrect     ;
        thisEvent.blocknumber       = blocknumber   ;
        thisEvent.miniblocknumber   = miniblocknumber;
        thisEvent.matchOnTime       = matchOnTime   ;
        thisEvent.probeOffTime      = probeOffTime  ;
        thisEvent.fixationOnTime    = fixationOnTime - mstime;
        thisEvent.fixationOffTime   = fixationOffTime - mstime;
        thisEvent.responseTime      = responseTime;
        thisEvent.vocalizedWord     = responseWord;
        
        if (index==1)
            events        = thisEvent; %- before events defined must convert to structure
        else
            events(index) = thisEvent;
        end
    end
    
end
fclose(fid);  % close session.log