function [events, table_fig, uitable1] = create_BL_events(subj,fileStem,state_ind,rootEEGdir,sessNum)
% ONLY ACCEPTS A SINGLE fileStem
eeg_dir = fullfile(rootEEGdir,subj,'eeg.noreref/'); % HARDCODED
event_dir=fullfile(rootEEGdir,subj,'behavioral','baseline');

if iscell(fileStem)
   fileStem = fileStem{1}; 
end
%% get date and format information
date_vec=zeros(1,6);
date_vec(1,1:5)=str2num(fileStem([1:2;3:4;5:6;8:9;10:11])); % str2num uses eval which may be dangerous...
[~, ~, Fs, Gn ] = GetFilesandChannels(subj,fileStem,'location',rootEEGdir);
date_vec(1)=date_vec(1)+2000;
dt_array=datetime(date_vec);

%% define default values
if state_ind==1
    state_data={'awake'};
elseif state_ind==0
    state_data={'asleep'};
else
    state_data={'unknown'};
end
notes={'none'};

t_sec=0;
eegoffset=0;
eegfile=strcat(fullfile(eeg_dir),fileStem);

% Determine recording system
fidSourceType = fopen(sprintf('%s/sourcetype.txt',eegfile));
sourceTypeString = fscanf(fidSourceType,'%s');
fclose(fidSourceType);

% If statement based on string
if contains(sourceTypeString,'cervello')
    recordingSystem = {'Cervello'};
elseif contains(sourceTypeString,'blackrock')
    recordingSystem = {'Blackrock'};
else
    recordingSystem = {'NK'};
end

% cell format for creating uitable
events_cell=[{subj},cellstr(dt_array),num2cell(sessNum),state_data,notes,t_sec,eegoffset,eegfile,Fs,Gn,recordingSystem]; % cell for now

% events struct
fields={'subjID','date','baselineNum','State', 'notes', 'timeSec','eegoffset','eegfile','Fs','Gn','recordingSystem'};
events=cell2struct(events_cell,fields,2);
events.date=datetime(events.date);

%% populate the uitable
figure('visible','on');
table_fig=gcf;
set(table_fig,'units','normalized','outerposition',[0.175 0.45 0.6475 0.25])

uitable1 = uitable('ColumnName',fields,'RowName',1,'unit','normalized','Position', [0 -0.2 1 1],'FontSize',12);
set(uitable1,'data',events_cell,'fontsize',16);

column_sizes=50*ones(1,size(get(uitable1,'data'),2));
column_sizes([1,3,4,5,7,10])=75;
column_sizes(2)=200;
% column_sizes(5)=250;
column_sizes(8)=750;
column_sizes(11)=125;
% allow for manual editing, may remove this (*)
colmn_edit=false(1,size(get(uitable1,'data'),2));
colmn_edit(5)=true; % so that behavioralProcessing can edit baselineNum

% colmn_edit([4,5])=true;

set(uitable1,'columnwidth',num2cell(column_sizes),'ColumnEditable',colmn_edit);

%% save (NOW THAT THIS IS BEING CALLED FROM 'BehavioralProcessing' NOTHING SHOULD BE SAVED IN THIS FUNCTION)

% % save_button = questdlg('The table has been updated to reflect the user inputs, are you REALLY sure you want to save/overwrite the events structure','save events?','yes','no','no');
% % if btn1_dat && strcmp(save_button,'yes')
% %     if ~exist(event_dir,'dir')
% %         mkdir(event_dir)
% %     end
% %     save(fullfile(event_dir,'events.mat'),'events')
% %     saveas(uitable1,fullfile(event_dir,'events_table.fig'))
% %     hgexport(table_fig, fullfile(event_dir,'events_table'),hgexport('factorystyle'), 'Format', 'png'); % this is a hack so it may break in future releases
% % end
%
% if ~exist(event_dir,'dir')
%     mkdir(event_dir)
% end
% % save(fullfile(event_dir,'events.mat'),'events')
% saveas(table1,fullfile(event_dir,'events_table.fig'))
% hgexport(table_fig, fullfile(event_dir,'events_table'),hgexport('factorystyle'), 'Format', 'png'); % this is an undocumented hack so it may break in future releases

