function processStimLog(logLocation,dbl_DC)
% processStimLog(logLocation,dbl_DC)
%
% dbl_DC    | fixes an issue with some sessions double counting DC pulses
% for some reason


if ~exist('logLocation','var')
    logLocation = uigetdir([],'Find stimlog folder');
end
if ~exist('dbl_DC','var') || dbl_DC == 0
    dbl_DC = 1; %-  multiply num DC pulses for all conditions
end

cd(logLocation)
D = dir;
for i = 1:length(D)
    my(i) = strncmpi(fliplr(D(i).name),fliplr('.stimlog'),7);
end
if sum(my) ~= 1
    fprintf('skipping... folder has %d stimLogs\n',sum(my))
    return
end
mine = D(my).name;

fid = fopen([logLocation '/' mine],'r');

C = textscan(fid,'%d64%s%d%d%d%d%s',-1,'headerLines',2);
fclose(fid);
fid = fopen([logLocation '/' mine],'r');
header = textscan(fid,'%s%s%s%s%s%s%.1f%s%s',1,'headerLines',1);
fclose(fid);
GUI_ver = header{7};

time_ms = C{1}/1000;
annotation_cell = {'Counting','Reading','Naming','Pointing','Other','AD','Seizure'};
bad_ln = []; bad_ev = {};

tot_pulses = 0;

time_start = time_ms(1);
elec_loc = {'',''};
telemark_mode = NaN;
NK_reference = [5 6];

subj = header{2};
sess = header{4};
for i = 1:length(time_ms)
    events = struct();
    events.subject = subj;
    events.session = sess;
    events.time_ms = time_ms(i);
    events.time_elapsed = time_ms(i)-time_start;
    events.location = C{2}{i};
    events.amp1 = C{3}(i);
    events.freq = C{5}(i);
    events.numPulses = C{6}(i);
    typeStr = C{7}{i};
    events.typeStr = typeStr;
    events.elec_loc = elec_loc;
    events.ann = '';
    events.pol = NaN;
    events.tele = telemark_mode;
    events.DCpulses = 0;
    events.GUIver = GUI_ver;
    events.NK_referenece = NK_reference;
    events.isWN_REP = false;
    events.WN_REP_obj = [];
    
    
    if strncmp(typeStr,'SM',2)
        events.type = 'MAPPING';
        if GUI_ver > 5          % correction for extra DC pulses due to 1s max zap protection
            events.DCpulses = floor(events.numPulses/events.freq);
        else
            events.DCpulses = 1* dbl_DC;
        end
    elseif strncmp(typeStr,'CC',2)
        events.type = 'CCEP';
        if length(C{7}{i}) >2
            events.pol = str2double(C{7}{i}(3));
        end
        events.DCpulses = C{6}(i);
    elseif strncmp(typeStr,'WN_REP',6)
        events.type = 'WN_REP';
        WN_obj.start = false;
        WN_obj.completed = false;
        events.DCpulses = 1;
        if length(typeStr)>6 %- annotation mode
            pieces = strsplit(typeStr,'_');
            isDone = strcmp(pieces{3},'DONE');
            WN_obj.start = ~isDone;
            WN_obj.mod_num = str2double(pieces{3+isDone});
            WN_obj.tot_mod = str2double(pieces{4+isDone});
            WN_obj.type = pieces{5+isDone};
            WN_obj.completed = isDone;
            events.DCpulses = 0;
        end
        events.isWN_REP = true;
        events.WN_REP_obj = WN_obj;
        
        
    elseif strncmp(typeStr,'WN',2)
        events.type = 'WN';
        events.DCpulses = 1;
    elseif strncmp(typeStr,'Switching',9)
        events.type = 'SWITCH_ELEC';
        tmp = strsplit(typeStr,'_');
        elec_loc{1} = tmp{2};
        elec_loc{2} = tmp{3};
        events.elec_loc = elec_loc;
    elseif any(strcmp(typeStr,annotation_cell))
        events.type = 'ANN';
        events.ann = typeStr;
    elseif strncmp(typeStr,'Telemark_',9)
        events.type = 'INFO';
        telemark_mode = strncmp(fliplr(typeStr),fliplr('_ON'),3);
        events.tele = telemark_mode;
    else 
        bad_ev = [bad_ev typeStr];
        bad_ln = [bad_ln i];
        continue
    end
    
    %events.DCpulses = events.DCpulses * dbl_DC;
    
    tot_pulses = tot_pulses+events.DCpulses;
    events.cumulative_pulses = tot_pulses;
    event(i) = events;
    
end

if length(bad_ev)>0
    warning('Issue with stimLog, manually check to figure out what!')
    for i = 1:length(bad_ln)
         fprintf('Line %d ... %s\n',bad_ln(i)',bad_ev{i});
    end
    fprintf('Fix, these lines to allow sucessful extraction!\n')
    fprintf('events3.mat NOT created \n');
    error !
    return
end

events3 = event;

        
save('events3','events3');
fprintf('events3.mat successfully created \n');

% maybe good to look at
% line = fgetl(fid); % intentional-- throw out first line
% 
% C = textscan(fid);
% 
% while true
%     this_line = fgetl(fid)