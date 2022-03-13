%function [events] = stimBlast_ExtractEvents(sessLogFile, subject, sessionName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Function for extracting behavioral data from stimulusBlast %%%%
%
%   extraction designed for stimulusBlast
%
%   create an event for every presented word and every response (no words from the training sections should be included)
%
% JHW 3/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncomment following lines to directly run script
clear all
subject     = 'NIH063_preOp';   % EEG002  NIH016
sessionName = 'session_2_v5a';


%- if modifying in FRNU (or local FRNU) use these lines
%rootEEGdir   = '/Volumes/JW24TB/data24TB/eeg_new';                      %office-local
%rootEEGdir   = '/Volumes/jwGDRIVE/data24TB/eeg_new';                    %home -- check first in case work is mounted
%sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/stimulusGuess',sessionName);


%- if modifying in experiments folder use these two lines
rootEEGdir   = '/Users/wittigj/Desktop/experiments/stimulusGuess/trunk/data';                      %office-local
sessionDir  = fullfileEEG(rootEEGdir,subject,sessionName);



%subject     = 'NIH047';   % EEG002  NIH016
%sessionName = 'session_2';
%sessionDir  = fullfileEEG('../data/',subject,sessionName);

sessLogFile = fullfileEEG(sessionDir,'session.log');
eventFile   = fullfileEEG(sessionDir,'events.mat');
priorEvents = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%fprintf('\nOn session: %s ', sessionName)


VERBOSE = 0;

%---------------------------------------------------------------------------------------%
%- Confirm Session File Exists & Create cleaned session number string
fid    = fopen(sessLogFile,'r');
if (fid==-1)
    error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else
    [sessionDir,~,~] = fileparts(sessLogFile); %- used for finding annotation files
    %disp([' The session.log file is located in: '  sessLogFile]);
    sessLogOriginal    = fullfileEEG(sessionDir,'sessionORIGINAL.log');
    sessLogUpdate      = fullfileEEG(sessionDir,'sessionUPDATEimgName.log');
    sessLogStandard    = fullfileEEG(sessionDir,'session.log');
    
    mapStandard =  fullfileEEG(sessionDir,'/../mapOld2NewImages_standardBlur');
    mapNoBlur   =  fullfileEEG(sessionDir,'/../mapOld2NewImages_noBlur');
    
    if ~exist(sessLogOriginal,'file'), copyfile(sessLogFile,sessLogOriginal); end
end

%- Convert session folder name into a number.  Should make sessions like "session_9trim" --> 9
strNumeric = find( sessionName >= '0' & sessionName <= '9');
if max(diff(strNumeric))>1,
    iKeep=[1:find(diff(strNumeric)>1,1,'first')];
    fprintf('\n Possible issue converting session name into a numeric.. %s --> %s; use %s', sessionName, sessionName(strNumeric), sessionName(strNumeric(iKeep)));
    strNumeric=strNumeric(iKeep);
end;
sessionNum = str2num( sessionName(strNumeric) );
if isempty(sessionNum), fprintf('\n ERROR: problem converting session name into a numeric'); keyboard;  end; %shouldn't need this catch...
%---------------------------------------------------------------------------------------%


thisIs90 = nan;
thisIs60 = nan;
mapImg   = {};
lastImgType     = '';
isExemplarBlock = nan;

%- make a new session log with key presses included
FIX_SESS_IMG_OFF = 1;
if FIX_SESS_IMG_OFF == 1,
    
    
    %- now loop over session.log and get strings + mstimes
    if exist(sessLogOriginal,'file'), fclose(fid); fid = fopen(sessLogOriginal,'r'); end %- use the "original" if it exists
    sessLogMS   = [];
    sessLogLine = {};
    iUpdate = 1;
    while 1
        thisLine = fgetl(fid);
        if ~ischar(thisLine); break; end
        
        %- convert string into cell array for each element in string
        xTOT = textscan(thisLine,'%s',20);
        thisLineStr = {};
        for k = 1:length(xTOT),
            if ~isempty(xTOT{k}),
                thisLineStr = cat(2,thisLineStr,xTOT{k});
            end
        end
        
        
        %- Generic text scan to get time, offset, and type
        mstime       = str2num(thisLineStr{1});
        
        
        
        if strcmp(thisLineStr{3},'STIMULUS_LIST'),
            thisIs90 = contains(thisLineStr{4},'90perf');
            thisIs60 = contains(thisLineStr{4},'60perf');
            
            if thisIs90 + thisIs60 ~= 1,
                keyboard;
            else
                sessLogMS(end+1)   = mstime;
                sessLogLine{end+1} = thisLine;
                if thisIs90,
                    thisLine = sprintf('%s\t%s\t%s\t%s',thisLineStr{[1 2]},'BLUR_DIR','ninetyBlurOG->ninetyBlur');
                    mapImg = load(mapNoBlur);
                    mapImg = mapImg.mapStruct.mapList;
                else
                    thisLine = sprintf('%s\t%s\t%s\t%s',thisLineStr{[1 2]},'BLUR_DIR','standardBlur->multiBlur');
                    mapImg = load(mapStandard);
                    mapImg = mapImg.mapStruct.mapList;
                end
            end
        end
        
        
        if strcmp(thisLineStr{3},'IMG_ON_SCREEN'),
            %- find old image, swap in new image name
            
            oldImgPath  = thisLineStr{9};
            if ~contains(oldImgPath,'noise1.jpg'),
                if strcmp(thisLineStr{5},'CLINTON') | strcmp(thisLineStr{5},'KENNEDY') | strcmp(thisLineStr{5},'BUSH') | strcmp(thisLineStr{5},'TRUMP'),
                    isExemplarBlock = 1;
                else
                    isExemplarBlock = 0;
                end
                lastImgType = thisLineStr{5};
                
                oldStimName = thisLineStr{8};
                iOldImg     = find(contains(mapImg(:,1),oldImgPath));
                newImgPath  = mapImg{iOldImg,2};
                [~,newStimName] = fileparts(newImgPath);
                
                %thisLineStr{8} = sprintf('%s->%s',oldStimName,newStimName);
                %thisLineStr{9} = newImgPath;
                thisLineStr{8} = newStimName;
                thisLineStr{9} = sprintf('[OLD/%s]->%s',oldStimName,newImgPath);
                thisLineNew = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', thisLineStr{[1:9]});
                thisLine = thisLineNew;
            end
        end
        
        
        
        if strcmp(thisLineStr{3},'KEYBOARD_PRESS') & ~strcmp(thisLineStr{5},'NONE'),
            keyPress = {'UP' 'RIGHT' 'DOWN' 'LEFT'};
            if isExemplarBlock,
                txtMatch = {'KENNEDY' 'TRUMP' 'BUSH' 'CLINTON'};
            else
                txtMatch = {'PERSON' 'OBJECT' 'ANIMAL' 'PLACE'};
            end
            
            iKey = find(strcmp(keyPress,thisLineStr{5}));
            
            newLineStr = thisLineStr;
            newLineStr{5} = sprintf('%s(%s)',keyPress{iKey},txtMatch{iKey});
            
            if strcmp(txtMatch(iKey),lastImgType),
                newLineStr{6} = 'CORRECT(hit)';
            else
                newLineStr{6} = 'INCORRECT(fa)';
            end
            
            if ~strcmp(thisLineStr{6}(1:7),newLineStr{6}(1:7)),
                %hasErrorInBlock = 1;
                fprintf('\n error in "correct" reporting: (%s) %s -> %s',lastImgType,thisLineStr{6},newLineStr{6});
                %keyboard;
            end
            
            thisLineNew = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', newLineStr{[1:7]});
            thisLine = thisLineNew;
        end
        
        
        sessLogMS(end+1)   = mstime;
        sessLogLine{end+1} = thisLine;
        
        
    end
    fclose(fid);
    
    %- now combine to create "sesion log with keyboard log"
    fprintf('\n writing new file');
    fidNewLog = fopen(sessLogUpdate,'w+');
    
    fprintf(fidNewLog,'%s\n',sessLogLine{:});
    fclose(fidNewLog);
    
    %- filecopy
    copyfile(sessLogUpdate,sessLogStandard)
    
end


