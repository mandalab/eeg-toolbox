function [sync,signalFs] = getBR_sync(sync1filepath,sync2filepath,signal1filepath)
%
% 12/2018, jw updated so sync returned in mV

% open files, check validity


NSP1 = concatOpenNSx(sync1filepath,0);
NSP2 = concatOpenNSx(sync2filepath,0);

syncFs = NSP1.MetaTags.SamplingFreq;
if syncFs ~= NSP2.MetaTags.SamplingFreq
    fprintf('Sync files at different samprates.. Should never hit this');
    keyboard
    error
end

signalFs = checkFs(signal1filepath);

% prepare output pulses
sync = cell(1,4);

% add in DC 12 pulses, if they're present
temp = strfind({NSP1.ElectrodesInfo.Label},'ain4');
DC12present = sum(cellfun(@(x) sum(x(:)),temp));
if DC12present==1
    DC12_index = find(~(cellfun(@(x) isempty(x(:)),temp)));
    
    NSP1_elec = concatOpenNSx(sync1filepath,0,1,strcat('c:',int2str(DC12_index)));
    NSP2_elec = concatOpenNSx(sync2filepath,0,1,strcat('c:',int2str(DC12_index)));
    
    maxAnalog  = NSP1_elec.ElectrodesInfo(1).MaxAnalogValue;
    maxDigital = NSP1_elec.ElectrodesInfo(1).MaxDigiValue;
    conv2mV    = double(maxAnalog)/double(maxDigital);
        
    sync1mV    = double(downsample(NSP1_elec.Data, syncFs/signalFs))*conv2mV;
    sync2mv    = double(downsample(NSP2_elec.Data, syncFs/signalFs))*conv2mV;
    
    sync{1,1} = sync1mV; sync{1,2} = sync2mv; %- DC12 on 1 and 2
end


% add in DC09 pulses, if they're present
temp = strfind({NSP1.ElectrodesInfo.Label},'ain1');
DC09present = sum(cellfun(@(x) sum(x(:)),temp));
if DC09present==1
    DC09_index = find(~(cellfun(@(x) isempty(x(:)),temp)));
    
    NSP1_elec = concatOpenNSx(sync1filepath,0,1,strcat('c:',int2str(DC09_index)));
    NSP2_elec = concatOpenNSx(sync2filepath,0,1,strcat('c:',int2str(DC09_index)));
    
    maxAnalog  = NSP1_elec.ElectrodesInfo(1).MaxAnalogValue;
    maxDigital = NSP1_elec.ElectrodesInfo(1).MaxDigiValue;
    conv2mV    = double(maxAnalog)/double(maxDigital);
        
    sync1mV    = double(downsample(NSP1_elec.Data, syncFs/signalFs))*conv2mV;
    sync2mv    = double(downsample(NSP2_elec.Data, syncFs/signalFs))*conv2mV;
    
    sync{1,3} = sync1mV; sync{1,4} = sync2mv; %- DC09 on 3 and 4
end


end %- function getBR_sync



%%% HELPER FUNCTION %%%

function [Fs] = checkFs(filename)
ext = filename(length(filename)-3:length(filename));
if strcmp('.ns1',ext)
    Fs = 500;
elseif strcmp('.ns2',ext)
    Fs = 1000;
elseif strcmp('.ns3',ext)
    Fs = 2000;
elseif strcmp('.ns4',ext)
    Fs = 10000;
elseif strcmp('.ns5',ext)
    Fs = 30000;
elseif strcmp('.ns6', ext)
    Fs= 30000;
end
end








%% some draft code in the case of needing to use NEVs for alignment,
% to be used instead of the above code


% NSP1_NEV = openNEV(sync1filepath,'nosave');
% NSP2_NEV = openNEV(sync2filepath,'nosave');
% 
% timeStamp_9_NSP1 = getBlackRockPulsesDC(NSP1_NEV,9);
% timeStamp_9_NSP2 = getBlackRockPulsesDC(NSP2_NEV,9);
% 
% timeStamp_12_NSP1 = getBlackRockPulsesDC(NSP1_NEV,12);
% timeStamp_12_NSP2 = getBlackRockPulsesDC(NSP2_NEV,12);
% 
% 
% pulses_br_dc09analog = zeros(1,size(noreref,2));
% addHere = round(pulses_br_dc09./30);
% for i=1:length(addHere)
%     pulses_br_dc09analog(addHere(i):addHere(i)+100) = 5000;
% end
% pulses_br_dc09 = pulses_br_dc09analog;
% 
% pulses_br_dc12_analog = zeros(1,size(noreref,2));
% addHere = round(pulses_br_dc12./30);
% for i=1:length(addHere)
%     pulses_br_dc12_analog(addHere(i):addHere(i)+50) = 5000;
% end
% pulses_br_dc12 = pulses_br_dc12_analog;
% 
% if size(pulses_br_dc12,2) ~= size(pulses_br_dc09,2)
%     minSize = min(size(pulses_br_dc12,2),size(pulses_br_dc09,2));
%     pulses_br_dc12 = pulses_br_dc12(:,1:minSize);
%     pulses_br_dc09 = pulses_br_dc09(:,1:minSize);
% end
% 
% sync{1,1} = pulses_ecog_dc12'; sync{1,3} = pulses_ecog';
% sync{1,2} = double(pulses_br_dc12); sync{1,4} = double(pulses_br_dc09);
%                         
                        

