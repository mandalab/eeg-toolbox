% produce stim events.m by taking no.reref outputs from eegPrepandAlign.m
%
% [0] stim data from GUI in behavioral & stim data extracted to raw/STIM
% [1] process sessions and identify files [automate]
% [2] make correct events from pulse_HI times
% [4] make final events.mat file
clear; clc;
task = 'stimMapping'; 
%addpath(genpath('/Users/steinhardtcr/Documents'))
addpath(genpath('/Users/steinhardtcr/Documents/OrganizedStimAnalysis/'))

addpath('/Users/steinhardtcr/Dropbox/ZaghloullabMaterial/FunctionsfromTim/processing_data_essential')
rootEEGdir = '/Users/steinhardtcr/Documents/OrganizedStimAnalysis/newEvent3sGood';%%'/Volumes/Shares/FRNU/dataWorking/eeg/';
subj = 'NIH044';

beh_fold = fullfileEEG(rootEEGdir,subj,'behavioral',task);
raw_fold = fullfileEEG(rootEEGdir,subj,'raw');%,'stimMapping');

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
D = (dir(raw_fold)); 
D = D(cellfun(@(x) length(x) == 11,{D.name}));
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
    
    f_stem_use = strrep([eeg_reref filesep f_stem{1}],'eeg.reref','eeg.noreref');

    INFO.subject = subj;
    INFO.session = i;
    INFO.eegfile = f_stem_use;% filesep stim_file_names{i}];
    %INFO.eegfile = [f_stem_use{1}];% filesep stim_file_names{i}];
    INFO.eegfolder = f_stem_use;
    

    % align
    if events3(end).cumulative_pulses == length(pulse_hi_times)
        fprintf('Same Number of Pulses... we are done here!! ')
        figure(2)   
         pulseTimes= find(([events3.DCpulses]>0));
        plot(pulse_hi_times,'g' ); hold on;       
         
        plot([events3( pulseTimes).time_elapsed],'r')
     
    else
        warning(['NOT Same Number of Pulses... we have some work to do\n\n' ...
        '%d GUI pulses \n%d DC10 pulses'],events3(end).cumulative_pulses,length(pulse_hi_times))
        
        figure(2),clf, hold on
         plot([events3( pulseTimes).time_elapsed],'b');
         pht2 = pulse_hi_times(3:end-1);
         hold on, plot((pht2([events3(pulseTimes).cumulative_pulses])),'g');
 
%         pulseTimes= find(([events3.DCpulses]>0));
%         for n = 1: length(events3)
%         [events3(n).time_elapsed] = deal(str2num(events3(n).offset) - str2num(events3(1).offset));
%         end


%          plot([events3(pulseTimes(101:end)).time_elapsed],'b');
         hold on, plot(pulse_hi_times(1:end),'r');
%         plot(diff([events3.time_elapsed]),diff(pulse_hi_times(3:end)),'k.')
        


        keyboard
          pulse_hi_times=pulse_hi_times(3:(end-1));
        legend DC GUI location Northwest
        %FT(20)
        set(gca,'FontSize',20)
        xlabel Pulse#
        ylabel 'Time (ms)'
        
        
       % pulse_hi_times=pulse_hi_times([events3(pulseTimes(2:end)).cumulative_pulses])
        
        %- NIH046 session_1
  %     pulse_hi_times=  pulse_hi_times([events3(pulseTimes).cumulative_pulses]); 
  %      pulse_hi_times= pulse_hi_times(1:end-1);
 %        pulse_hi_times = pulse_hi_times(3:end); %-NIH047 [3/5]
       %pulse_hi_times=pulse_hi_times(1:24)%end-2);
       %nih049 [1/4] (2:end) 3/4 (3:end)
%       pulse_hi_times=pulse_hi_times(3:end);
  
%              events3= events3(pulseTimes(101:(end)));
%              org_cp = events3(1).cumulative_pulses;
%            for n = 1:length(events3)
%            events3(n).cumulative_pulses = events3(n).cumulative_pulses- org_cp+1;
%            end
% % %        pulse_hi_times([events3(pulseTimes(1:(end-100))); %[4/4]
    end
    
    
    % make events
addPathPublicOnly('/Users/steinhardtcr/Dropbox/ZaghloullabMaterial/FunctionsfromTim/processing_data_essential')
    events = create_stim_events_v_crs(events3,pulse_hi_times,INFO);
   
    % save to folder and prepare to concatenate
    cd(this_sess_fold)
   save('events.mat','events')
    fprintf('Events Saved! \n')
    ev_cell{i} = events;
    clear events

% all the session markers disappear because they didn't have dc pulses in
% there but the offsets can be used to reliably find session breaks anyways
end
fprintf('Combining Events...')

for n = 1: length(events)

end
events = [ev_cell{:}];
cd(beh_fold)
save('events.mat','events')
fprintf('Done!\n')
%clc
clearvars -except events