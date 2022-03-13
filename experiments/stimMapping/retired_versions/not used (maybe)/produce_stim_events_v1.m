function produce_stim_events_v1(rootEEGdir,subj)

% produce stim events
%
% [0] stim data from GUI in behavioral & stim data extracted to raw/STIM
% [1] process sessions and identify files [automate]
% [2] make correct events from pulse_HI times
% [4] make final events.mat file
% clear; clc;
task = 'stimMapping'; 
% rootEEGdir = '/Volumes/Shares/FRNU/dataWorking/eeg/';
% subj = 'NIH047_2';

beh_fold = fullfileEEG(rootEEGdir,subj,'behavioral',task);
raw_fold = fullfileEEG(rootEEGdir,subj,'raw','STIM');

eeg_reref = fullfileEEG(rootEEGdir,subj,'eeg.reref');
eeg_noreref = fullfileEEG(rootEEGdir,subj,'eeg.noreref');

%% [1]
% get session folders
sess_tmp = dir(beh_fold);
sess = sess_tmp(strncmp({sess_tmp.name},'session_',8));
n_sess = length(sess);
% process
sess_locs = cell(1,n_sess);
for ind = 1:n_sess
    this_loc = fullfileEEG(beh_fold,sess(ind).name);
    processStimLog(this_loc);
    sess_locs{ind} = this_loc;
end
clear sess_tmp this_loc ind
%% [2]  
% get STIM folders
% raw_tmp = dir(raw_fold);
% for i = 1:length(raw_tmp)
%     is21E(i) = strncmp(fliplr(raw_tmp(i).name),'E12',3);
% end
%stim_file_names = {'161107_1250','161107_1600','161107_1712'}; %-NIH045 %-cheating, need to 
%stim_file_names = {'161214_1529','161215_1539'};
D = rm_hidden(dir(raw_fold)); 
stim_file_names = {D([D.isdir]==1).name};
n_stim_files = length(stim_file_names);
if n_stim_files ~= n_sess
    warning('Incompatable Number of Stim sessions to files saved. May by missing file or need to split session!')
end
% get DC10 times
f_end = 'DC10.updown.txt';
fprintf('Attempting to align stim Sessions \n')
for i = 1:n_stim_files
    fprintf('Session [%d / %d]',i,n_stim_files )
    this_stim_loc = fullfileEEG(eeg_noreref,stim_file_names{i});
    these_files = dir(this_stim_loc);
    my_dc10 = false(1,length(these_files));
    for ii = 1:length(these_files)
        my_dc10(ii) = strncmp(fliplr(these_files(ii).name),fliplr(f_end),length(f_end));
    end
    assert(sum(my_dc10)==1,'NoDC10 File!')
    my_file = these_files(my_dc10).name;
    
    % extract file info and close
    my_file_loc = fullfileEEG(this_stim_loc,my_file);
    fid = fopen(my_file_loc);
    
    f_name = textscan(fid,'%d%s%s',1); f_name = f_name{3};
    f_stem = textscan(fid,'%d%s%s',1); f_stem = f_stem{3};
    pulse_times = textscan(fid,'%d%s%s',-1);
    fclose(fid);
    
    % get PULSE_HI times 
    is_pulse_hi = strcmp(pulse_times{2},'PULSE_HI');
    pulse_hi_times = pulse_times{1}(is_pulse_hi);
    
    % try to align to stimlog
    % assume 1 to 1 correspondence 
    this_sess_fold = sess_locs{i};  %- also where I will save sub events 
    stimLogLoc = fullfileEEG(this_sess_fold,'events3.mat');
    load(stimLogLoc);
    
    % change f_stem to data
    
    f_stem_use = strrep(strrep([eeg_reref filesep f_stem{1}],'/dataWorking/','/data/'),'eeg.noreref','eeg.reref');

    INFO.subject = subj;
    INFO.session = i;
    INFO.eegfile = f_stem_use;% filesep stim_file_names{i}];
    %INFO.eegfile = [f_stem_use{1}];% filesep stim_file_names{i}];
    INFO.eegfolder = f_stem_use;
    

    % align
    if events3(end).cumulative_pulses == length(pulse_hi_times)
        fprintf('Same Number of Pulses... we are done here!! ')
        
    elseif strncmp(subj,'NIH047',6) && i==3 && length(pulse_hi_times)==682
        pulse_hi_times(1:2) = [];
        fprintf('Repair already prepared for this sesion! ')
    else
        
        warning(['NOT Same Number of Pulses... we have some work to do\n\n' ...
        '%d GUI pulses \n%d DC10 pulses'],events3(end).cumulative_pulses,length(pulse_hi_times))
        
        figure(1),clf, hold on
%         plot(diff([events3.time_elapsed]),'b');
%         hold on, plot(diff(pulse_hi_times(3:end)),'r');

        plot(pulse_hi_times([events3([events3.DCpulses]>0).cumulative_pulses]) )
        plot([events3([events3.DCpulses]>0).time_elapsed],'r')
        %plot([events3([events3.DCpulses]>0).time_elapsed],pulse_hi_times([events3([events3.DCpulses]>0).cumulative_pulses]+2) )
        legend DC GUI
        FT(20)
        xlabel Pulse#
        ylabel 'Time (ms)'
        
%         plot([events3.time_elapsed],'b');
%         hold on, plot(pulse_hi_times(1:end),'r');
%         plot(diff([events3.time_elapsed]),diff(pulse_hi_times(3:end)),'k.')
        
        keyboard
        
        %- NIH046 session_1
%         pulse_hi_times(1:end-1);
        %pulse_hi_times = pulse_hi_times(3:end); %-NIH0
    end
    
    
    % make events

    events = create_stim_events_v1(events3,pulse_hi_times,INFO);
    
    % save to folder and prepare to concatenate
    cd(this_sess_fold)
    save('events.mat','events')
    fprintf('Events Saved! \n')
    ev_cell{i} = events;
    clear events


end
fprintf('Combining Events...')

events = [ev_cell{:}];
cd(beh_fold)
save('events.mat','events')
fprintf('Done!\n')
%clc
