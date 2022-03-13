function timeStr = printTiming(ticVal)
thisRawTime = toc(ticVal);
thisTimeMin = floor(thisRawTime/60);
if thisTimeMin > 60, thisTimeHour = floor(thisTimeMin/60); thisTimeMin = mod(thisTimeMin,60); else thisTimeHour = 0; end
thisTimeSec = mod(thisRawTime,60);

thisHrStr = num2str(thisTimeHour); while length(thisHrStr) < 2, thisHrStr = ['0',thisHrStr]; end
thisMinStr = num2str(thisTimeMin); while length(thisMinStr) < 2, thisMinStr = ['0',thisMinStr]; end
thisSecStr = sprintf('%.2f',thisTimeSec);  while length(thisSecStr) < 5, thisSecStr = ['0',thisSecStr]; end

    
timeStr = sprintf('%s:%s:%s',thisHrStr,thisMinStr,thisSecStr);

%keyboard