function [events] = extractStimEvents(subject, sessName, behDir, eegDir)

%makeNewEventfrmSL
%closest can get to having WN events structures right now

rootDir = '/Volumes/Shares/FRNU/dataWorking/eeg/';
funDir = '/Volumes/FRNU/people/steinhardt/Usefulmatfun';
labfunDir = '/Users/steinhardtcr/Documents/eeg_toolbox_trunk';
addpath(genpath(funDir));
addpath(genpath(labfunDir));

% find eegfiles
%  cd([rootDir subject '/raw/stimMapping/']);
% stimEEGfileFolders=dir;
% files ={stimEEGfileFolders(~strncmp({stimEEGfileFolders.name},'.',1)).name};
dbl_DC = 1;
%typical annotations
annotation_cell = {'Counting','Reading','Naming','Pointing','Other','AD','Seizure'};


cd([behDir '/' sessName])
fName = dir('*stimlog');
stimFile = fName.name;

fid = fopen(stimFile,'r');
header = textscan(fid,'%s%s%s%s%s%s%.1f%s%s',1,'headerLines',1);
fclose(fid);
GUI_ver = header{7};

[offset sites amp1 amp2 freq pulses type annotation]=textread(stimFile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t');

sites=sites(4:end);

sitesbyName{size(sites,1)} = [];

[AllTAGS] = getjckOrder(subject,'/Volumes/FRNU/dataWorking/eeg/');
% get all the site pair names
sitepairs= unique(sites);
for n = 1:length(sitepairs)
sitenums= strsplit(sitepairs{n},'-');
stimNames=strcat([AllTAGS{str2num(sitenums{1})} '-' AllTAGS{str2num(sitenums{2})}]);
[sitesbyName{find(strcmp(sites,sitepairs{n}))}] = deal(stimNames);
end
sitesbyName = sitesbyName(1:size(sites,1));
%% assume all NK uses channels 5 and 6
NKnames=strcat([AllTAGS{(5)} '-' AllTAGS{(6)}]);% pretty sure always use if it's not said

% make events structure
idx= 4;evcnt = 1; tot_pulses =0 ;
while idx < size(offset,1)
events(evcnt).subject=subject;
events(evcnt).exptFolder= sessName;
events(evcnt).amplitude=  str2num(amp1{idx});
events(evcnt).frequency = str2num(freq{idx});
events(evcnt).stimType = type{idx};
events(evcnt).stimLocTag = sitesbyName{evcnt};
events(evcnt).annotation = annotation{idx};
% initialize as no pulses
events(evcnt).DCpulses = 0;
%% getting the dc pulse count

 if strncmp(type{idx},'SM',2)
        
        if GUI_ver > 5          % correction for extra DC pulses due to 1s max zap protection
            events(evcnt).DCpulses = floor(str2num(pulses{idx})/events(evcnt).frequency);
        else
            events(evcnt).DCpulses = 1* dbl_DC;
        end
    elseif strncmp(type{idx},'CC',2)
   
            events(evcnt).monopolar = str2double(type{idx}(3));
        events(evcnt).DCpulses = str2num(pulses{idx});
    elseif strncmp(type{idx},'WN_REP',6)
        %events(idx).stimType = 'WN';      
        events(evcnt).DCpulses = 1;
        if length(type{idx})>6 %- annotation mode
            pieces = strsplit(type{idx},'_');
            isDone = strcmp(pieces{3},'DONE');
            WN_obj.start = ~isDone;
            WN_obj.mod_num = str2double(pieces{3+isDone});
            WN_obj.tot_mod = str2double(pieces{4+isDone});
            WN_obj.type = pieces{5+isDone};
            WN_obj.completed = isDone;
            events(evcnt).DCpulses = 0;
        end
        events(evcnt).isWN_REP = true;
        events(evcnt).WN_REP_obj = WN_obj;
          if str2num(amp1{idx}) < 0
         events(evcnt).monopolar = 0; % 0 is negative 1 is positie
        else
         events(evcnt).monopolar = 1;
        end  
    elseif strncmp(type{idx},'WN',2)
       % events.type = 'WN';
        events(evcnt).DCpulses = 1;
        if str2num(amp1{idx}) < 0
         events(evcnt).monopolar = 0; % 0 is negative 1 is positie
        else
         events(evcnt).monopolar = 1;
        end
    elseif strncmp(type{idx},'Switching',9)
        events(evcnt).stimType = 'SWITCH';
        tmp = strsplit(type{idx},'_');
        elec_loc{1} = tmp{2};
        elec_loc{2} = tmp{3};
        events(evcnt).elec_loc = elec_loc;
    elseif any(strcmp(type{idx},annotation_cell))
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
    elseif strncmp(type{idx},'Telemark_',9)
        events(evcnt).stimType = 'INFO';
        telemark_mode = strncmp(fliplr(type{idx}),fliplr('_ON'),3);
        events(evcnt).tele = telemark_mode;
    else 
      %  bad_ev = [bad_ev typeStr];
      %  bad_ln = [bad_ln i];
        events(evcnt).stimType = 'ANN';
        events(evcnt).annotation = type{idx};
      %  continue
 end

 tot_pulses = tot_pulses+events(evcnt).DCpulses;
    events(evcnt).cumulative_pulses = tot_pulses;
%%

%NEED
%events(idx).eegfile = eegfile;
%isWNREP
%wnRepNum
%cumpulse
%dcpulse
curMatTime = str2num(offset{idx})/1e11;
t= addtodate(curMatTime, 4, 'hour');
tmpdate =datetime(t,'ConvertFrom','datenum');
curTimePos =num2str(round(posixtime(tmpdate)*1000));    

events(evcnt).NKreference = NKnames;
events(evcnt).pulses = pulses{idx};
events(evcnt).offset = offset{idx};
events(evcnt).mstime = curTimePos;
evcnt= evcnt+1;
idx =idx+1;
end
cd([rootDir subject '/behavioral/stimMapping/' sessName ]);
save('eventstmp','events')
%% make the alignment file for later based on the dc and cumulative pulse count:

    pulseTimes= find(([events.DCpulses]>0)); % the events that aren't just information and accounting for multple frequencies of spiking
       fileId =fopen('eeg.eeglog.up','w');  % for now just replacign the eeg.eeglog.up file 
    for n = 1:length(pulseTimes)
    
    fprintf(fileId,'%s\t%d\t%s\n',events(pulseTimes(n)).mstime,0,'CHANNEL_0_UP');
    end
    fclose(fileId);



end
