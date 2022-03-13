function [analogDCfix] = fix_CV_analog(analogDC_mV, chanName)
%
%  fix_CV_analog
%
%     sometimes cervello "analog" traces (DC09, DC10, DC12) do this weird thing with a few very negative samples
%     this screws up our pulse detection algorith, which scales the voltage traces and takes 25th percentile as threshold
%       which would cause the threshold to only catch the random down peaks instead of the real pulses between 0 and 5V
%
%     JW thought this could come from the transform from raw uint16 to an analog value in readCervelloTRC(), but it does not
%       instead, occassioal raw values were in the range of 0 to ~20, whereas ~32000 is equal to 0mV and 64000 is equal to 5000mV
%       so a fix needs to be made to detect the abberant values
%
%     This function finds those abberant values and sets them to the average height of the "good" pulses
%
%     Called by split_CV and optionally called by eegPulseVisualize
%
%  JW 12/2018
%

%- 
SHOW_FIG = 0;  %- debug fig

%- create the output time series
analogDCfix = analogDC_mV;


%- does appear to have pulses that are >1V?
Yrng  = prctile(analogDC_mV,[0 1 99.99]);  %-  Use 99.99 to always catch true pulse peaks (~10 samp/1000 if DC09 on the whole time)
rangeRobust = Yrng(3)-Yrng(2); %- 1st to 99.9th percentile should be >1000mV if pulses are recorded
minVal      = Yrng(1); 
loVal       = Yrng(2);
hiVal       = Yrng(3);


%- if a true "bad" signal, then a very small set of samples should be <-4500, and the rest should be > -500mV.
Yends = prctile(analogDC_mV,[0.01 1 99.99 100],2);  %- check out the percentiles.


%- bad values are usually less than -4500, but in NIH060 they were around -3000
%negThresh = -4500;
negThresh = -1000;


error2see = 0; %- automatically pull up figure if something weird

%- maybe this isn't a pulse recording... if so, the big negative thresh crossing should be paired with smaller negative thres
%- try to make sure this is actually a pulse recording gone bad, and not an analog trace recording on DC09,10,12
%    so if the minimum is less than negative thresh, there should be nothing between -4500 and -500.
if minVal<negThresh & loVal > -500  &  rangeRobust>1000,
    
    numNegSamples = sum(analogDC_mV<negThresh);
    fixedValue    = hiVal;  %- set to one off from max in case max is a fluke.
    analogDCfix(analogDC_mV<negThresh)= fixedValue;
    
    fprintf('\n fix_CV_analog just converted %d samlpes of %s from <%d to %.1f', numNegSamples, chanName, negThresh, fixedValue);

    %- expect fixed value to be above 2V... should that be enforced?
    if fixedValue < 2000,
        fprintf('\n expecting fixedValue to be less than 2000, but it isnt');
        error2see = 1; 
    end
    
    
elseif minVal<negThresh  &  rangeRobust>1000,
    fprintf('\n weird case?  valuves less than -4500mV, robust range>1V, but there are also values between -4500 and -500.  Ask JW to take a look');
    error2see = 1;

end



if SHOW_FIG | error2see,
    figure(987); clf; set(gcf,'color','w','name','Cervello Analog Fix');
    fS = 15;
    subplot(211)
    plot(analogDC_mV);
    set(gca,'fontsize',fS,'tickdir','out','box','off');
    ylabel('analogDC, pre fix');
    title(chanName)
    
    subplot(212)
    plot(analogDCfix);
    set(gca,'fontsize',fS,'tickdir','out','box','off');
    ylabel('analogDC, post fix');
    
    fprintf('\n in fix_CV_analog... ');
    keyboard
end



return


% %%% OTHER WAYS TO DO IT
%
% rngInner = Y(3)-Y(2); %- try to avoid outliers
% rngOuter = Y(4)-Y(1); %- if dont see pulses with 99th percentile, maybe here?
% threshUse = nan;
%
% if rngInner>1000,
%     threshUse = Y(2)+rngInner/2;
% elseif rngOuter>1000,
%     threshUse = Y(1)+rngOuter/2;
% else
%     threshUse = 6000; %- range 5000, so this should result in no pulses;
% end
% %bncPulseCnt = length(find( bncData_mV(iiBNC,2:end)>=threshUse(iiBNC) & bncData_mV(iiBNC,1:end-1)<threshUse(iiBNC)));
%
%
%
%
% %- fill up temp with a count of samples where analogDC is greater than some thresh
% sampAboveThresh = [];  counter = 1;
% for thresh=6000:-100:-500;
%     sampAboveThresh(counter, [1:2]) = [thresh length(find(analogDC(1:end-1) > thresh & analogDC(2:end) <= thresh))]; %- now count transitions
%     %sampAboveThresh(counter, [1:2]) = [thresh length(find(analogDC(:) > thresh))]; %- was counting samples above
%     counter=counter+1;
% end
%
% %- if temp found ANY samples > ANY thresh
% if sum(sampAboveThresh(:,2))>0
%     %- find the biggest jump in super-thresh samples... this should roughly mark the peak of the pulses (i.e., ~5V)
%     diffs      = diff(sampAboveThresh(1:end-1,2)); %- end-1 so the thresh can equal the thresh below max
%     threshJump = sampAboveThresh(find(diffs==max(diffs),1,'first')+1,1);
%     th
%
%     %- in some weird cases ain1 is hovering at 1000mV for 1/4 of the session... that screws up this approach
%     newHiVal =  mean(analogDC(analogDC > threshJump));
%     analogDC(analogDC < 0)          = newHiVal;           %- convert negatives to mean positive.
%     analogDC(analogDC > threshJump) = newHiVal;           %- convert positives to mean positive.
% end
%
%
%
