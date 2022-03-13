function [newevents] = align_mangold(events)

root_dir = '/Volumes/Shares/FRNU/dataWorking/eeg/NIH053/';
eye_sess = fullfile(root_dir,'behavioral','eyeTracking');


% load task events
sess_files = unique({events.eegfile});


% load eyetracking sessions
tmp_eye = dir(eye_sess);
eye_sess2 = {tmp_eye(3:5).name}
win_size = 20;

newevents = [];
for i = 1:length(sess_files)
  
  % sync file
  sess_loc = strrep(sess_files{i},'data','dataWorking')
  fid = fopen(fullfile(sess_loc,'trigDC10.syncStim.txt'))
  sync_times = textscan(fid,'%f'); sync_times = sync_times{1};
  sp_time = diff(sync_times);
  
  % eeg log
  fid = fopen(fullfile(eye_sess,eye_sess2{i},'eeg.eeglog'));
  log_times = textscan(fid,'%f %s'); log_times = log_times{1};
  eg_time = diff(log_times);
  
  % align
  mid_ind = round(length(sp_time)/2);
  sp_win = sp_time(mid_ind:mid_ind+win_size);
  
  % regress
  for j = 1:length(eg_time)-win_size
    r(j) = corr(eg_time(j:j+win_size),sp_win);
  end
  
  % check
  [rmax,eg_ind] = max(r);
  [diff(log_times(eg_ind:eg_ind+win_size)) diff(sync_times(mid_ind:mid_ind+win_size))]
  
  % offset
  offset = log_times(eg_ind)-sync_times(mid_ind);
  
  % eye data for each event
 
  % add offset 
  sess_events = events(strcmp({events.eegfile},sess_files(i)));
  for k = 1:length(sess_events)
    sess_events(k).eyeoffset = sess_events(k).eegoffset+offset;  
    sess_events(k).eyefile = fullfile(eye_sess,eye_sess2{i},'eyedata.txt');
  end
  
  % add file 
  newevents = [newevents sess_events];
 
end









