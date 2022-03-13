function [events] = extractManipMemEventsExp_v02(subject, behDir, sessionName)
% confirm session file exists
sessLogFile = sprintf('%s/session_%d/session.log',behDir,sessionName);

fid    = fopen(sessLogFile,'r');
if (fid==-1); error('Session.log not found: \n %s \n Exiting.',sessLogFile);
else [sessionDir,~,~] = fileparts(sessLogFile); end

% convert session folder name into a number.
%sessionNum = str2num(sessionName(end))
sessionNum = sessionName(end);

% read session.log line-by-line and convert to events structure
events     = [];
index      = 1;

while 1
    thisLine = fgetl(fid);
    if thisLine==-1; break; end
    
    % convert string into cell array for each element in string
    xTOT = textscan(thisLine,'%s',50); thisLineStr = {};
    for k = 1:length(xTOT),
        if ~isempty(xTOT{k});thisLineStr = cat(2,thisLineStr,xTOT{k});end
    end
    
    %- generic text scan to get time, offset, and type
    mstime       = str2num(thisLineStr{1});
    try
        eventType    = thisLineStr{4};
    catch
        eventType    = thisLineStr{3};
    end
    
    makeNewEvent = 0;
    thisEncImg = NaN;
    thisRetImg = NaN;
    subjResp  = NaN;
    corrResp = NaN;
    RT = NaN;
    isAdded = NaN;
    isRepeated = NaN;
    isRemoved = NaN;
    isCorr = NaN;
    xclick = NaN;
    yclick = NaN;
    
    
    switch eventType
        case 'ENCFIX'
            
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            
        case 'ENCIMG'
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            [~,name,ext] = fileparts(thisLineStr{12}); thisEncImg = [name ext];
            [~,name,ext] = fileparts(thisLineStr{14}); thisRetImg = [name ext];
            if strcmp(thisRetImg,thisEncImg); isRepeated = 1; end
            if thisRetImg(1)=='m' & thisEncImg(1)~='m'; isRemoved = 1; end
            if thisRetImg(1)~='m' & thisEncImg(1)=='m'; isAdded = 1; end
            
            
        case 'RETFIX'
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            
        case 'RETIMG'
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            [~,name,ext] = fileparts(thisLineStr{12}); thisEncImg = [name ext];
            [~,name,ext] = fileparts(thisLineStr{14}); thisRetImg = [name ext];
            if strcmp(thisRetImg,thisEncImg); isRepeated = 1; end
            if thisRetImg(1)=='m' & thisEncImg(1)~='m'; isRemoved = 1; end
            if thisRetImg(1)~='m' & thisEncImg(1)=='m'; isAdded = 1; end
            isCorr = str2num(thisLineStr{16});
            RT = str2num( thisLineStr{18});
            
            
        case 'KEYPRE'
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            [~,name,ext] = fileparts(thisLineStr{12}); thisEncImg = [name ext];
            [~,name,ext] = fileparts(thisLineStr{14}); thisRetImg = [name ext];
            if strcmp(thisRetImg,thisEncImg); isRepeated = 1; end
            if thisRetImg(1)=='m' & thisEncImg(1)~='m'; isRemoved = 1; end
            if thisRetImg(1)~='m' & thisEncImg(1)=='m'; isAdded = 1; end
            isCorr = str2num(thisLineStr{16});
            RT = str2num(thisLineStr{18});
        case 'MOUSEON'
            makeNewEvent = 1;
            thisRound = thisLineStr{5};
            thisPos = thisLineStr{7};
            thisTrial = thisLineStr{9};
            [~,name,ext] = fileparts(thisLineStr{12}); thisEncImg = [name ext];
            [~,name,ext] = fileparts(thisLineStr{14}); thisRetImg = [name ext];
            if strcmp(thisRetImg,thisEncImg); isRepeated = 1; end
            if thisRetImg(1)=='m' & thisEncImg(1)~='m'; isRemoved = 1; end
            if thisRetImg(1)~='m' & thisEncImg(1)=='m'; isAdded = 1; end
            isCorr = str2num(thisLineStr{16});
            RT = str2num(thisLineStr{18});
            xclick = str2num(thisLineStr{20});
            yclick = str2num(thisLineStr{22});
            
    end
    
    if makeNewEvent
        clear thisEvent
        thisEvent.experiment = 'manipMemExp';
        thisEvent.subject = subject;
        thisEvent.session = sessionNum;
        thisEvent.type = eventType;
        thisEvent.enc = thisEncImg;
        thisEvent.ret = thisRetImg;
        thisEvent.added = isAdded;
        thisEvent.removed = isRemoved;
        thisEvent.repeated = isRepeated;
        thisEvent.correct = isCorr;
        thisEvent.rt = RT;
        thisEvent.mstime = mstime;
        %thisEvent.str = thisLineStr;
        thisEvent.xclick = xclick;
        thisEvent.yclick = yclick;
        if index==1; events = thisEvent;
        else events(index) = thisEvent; end
        index = index+1;
        
    end
    
end
fclose(fid);



try
    [{events.type}; {events.enc}; {events.ret}; {events.added}; ...
        {events.removed}; {events.repeated}; {events.correct}; {events.rt}; ...
        {events.mstime}];
    ret_events = events(strcmp({events.type},'RETIMG'));
    ret_add = ret_events([ret_events.added]==1);
    c = logical([ret_add.correct]);
    i = logical(~[ret_add.correct]);
    crt = [ret_add(c).rt];
    irt = [ret_add(i).rt];
    
    % fprintf(['ADD CORR: %d \n'],sum(c))
    % fprintf(['ADD INCO: %d \n'],sum(i))
    % fprintf(['ADD PERC: %d \n'],round(sum([ret_add.correct])/length([ret_add])*100))
    % fprintf(['ADD CORR RT: %d +- %d \n'],round(mean(crt)),round(std(crt)))
    % fprintf(['ADD INCO RT: %d +- %d \n'],round(mean(irt)),round(std(irt)))
    % fprintf('\n\n')
    
    
    
    ret_rem = ret_events([ret_events.removed]==1);
    c = logical([ret_rem.correct]);
    i = logical(~[ret_rem.correct]);
    crt = [ret_rem(c).rt];
    irt = [ret_rem(i).rt];
    
    % fprintf(['REM CORR: %d \n'],sum(c))
    % fprintf(['REM INCO: %d \n'],sum(i))
    % fprintf(['REM PERC: %d \n'],round(sum([ret_rem.correct])/length([ret_rem])*100))
    % fprintf(['REM CORR RT: %d +- %d \n'],round(mean(crt)),round(std(crt)))
    % fprintf(['REM INCO RT: %d +- %d \n'],round(mean(irt)),round(std(irt)))
    % fprintf('\n\n')
    
    ret_rep = ret_events([ret_events.repeated]==1);
    c = logical([ret_rep.correct]);
    i = logical(~[ret_rep.correct]);
    crt = [ret_rep(c).rt];
    irt = [ret_rep(i).rt];
    
    % fprintf(['REP CORR: %d \n'],sum(c))
    % fprintf(['REP INCO: %d \n'],sum(i))
    % fprintf(['REP PERC: %d \n'],round(sum([ret_rep.correct])/length([ret_rep])*100))
    % fprintf(['REP CORR RT: %d +- %d \n'],round(mean(crt)),round(std(crt)))
    % fprintf(['REP INCO RT: %d +- %d \n'],round(mean(irt)),round(std(irt)))
catch
end
