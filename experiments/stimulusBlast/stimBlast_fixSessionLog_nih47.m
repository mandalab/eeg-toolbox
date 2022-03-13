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

%rootEEGdir   = '/Volumes/JW24TB/data24TB/eeg_new';                      %office-local
rootEEGdir   = '/Volumes/jwGDRIVE/data24TB/eeg_new';                    %home -- check first in case work is mounted
subject     = 'NIH047';   % EEG002  NIH016
sessionName = 'session_2';
sessionDir  = fullfileEEG(rootEEGdir,subject,'behavioral/stimulusBlast',sessionName);

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
    sessLogWithVidFile = fullfileEEG(sessionDir,'sessionWithVidLog.log');
    sessLogOriginal    = fullfileEEG(sessionDir,'sessionOriginal.log');
    
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


%- make a new session log with key presses included
FIX_SESS_TIMES = 1;
if FIX_SESS_TIMES == 1 & strcmp(subject,'NIH047'),  
    
    keyLogFile = fullfileEEG(sessionDir,'video.vidlog');
    if exist(keyLogFile),
        
        fidKL = fopen(keyLogFile,'r');
        updateTimeMS = [];
        updateImage  = {''};
        updateLine   = {};
        iSave = 1;
        while 1
            thisLine = fgetl(fidKL);
            if ~ischar(thisLine); break; end
            
            %- convert string into cell array for each element in string
            xTOT = textscan(thisLine,'%s',20);
            thisLineStr = {};
            for k = 1:length(xTOT),
                if ~isempty(xTOT{k}),
                    thisLineStr = cat(2,thisLineStr,xTOT{k});
                end
            end
            
            %-increment before "UPDATE" happens
            if length(thisLineStr)>6,
                
                if strcmp(thisLineStr{7},'UPDATE') & length(updateImage{iSave})>0,
                    if length(updateImage)<2 ||  ~strcmp(updateImage{iSave},updateImage{iSave-1})
                        iSave=iSave+1;
                        updateImage{iSave}='';
                    end
                    
                elseif strcmp(thisLineStr{7},'IMAGE'),
                    %- Generic text scan to get time, offset, and type
                    %- fair way to incorporate error.  Actual event happens sometime between first and second number.  So lets just take 1/2 of second number
                    mstime3 = str2num(thisLineStr{1}) + str2num(thisLineStr{2})/2.0;
                    %- want the lines that preceed "UPDATE", as that says what the top layer of the video is (i.e., the visible image/text)
                    updateLine{iSave}   = thisLine;
                    updateTimeMS(iSave) = mstime3;
                    updateImage{iSave}  = thisLineStr{8};
                end
            end
     
        end
        fclose(fidKL);
        
        
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
            
            if strcmp(thisLineStr{3},'TXT_ON_SCREEN')  | strcmp(thisLineStr{3},'IMG_ON_SCREEN'),
                thisImgFile = thisLineStr{9};
                if strcmp(thisImgFile,updateImage{iUpdate}),
                    mstime = updateTimeMS(iUpdate);
                    thisLine(1:length(thisLineStr{1})) = sprintf('%.0f',mstime);
                    iUpdate=iUpdate+1; 
                else
                    fprintf('\n uh oh');
                end
            end
            
            sessLogMS(end+1)   = mstime;
            sessLogLine{end+1} = thisLine;
                
                
        end
        fclose(fid);
        
        %- now combine to create "sesion log with keyboard log"
        fprintf('\n writing new file');
        fidNewLog = fopen(sessLogWithVidFile,'w+');
        
        fprintf(fidNewLog,'%s\n',sessLogLine{:});
        fclose(fidNewLog);
        
        
       
    end
end

