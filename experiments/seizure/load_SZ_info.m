function [SeizureTable, SeizureCell, events, isSubjFound, sz_filepath]=load_SZ_info(subj, sz_filepath)
%LOAD_SZ_INFO
%Written by: Julio I. Chapeton
%
%This function loads the seizure information from the master spreadsheet in table and cell array format
%Can be run with no inputs to load from the default location, or a path can be provided.
%
%Usage:
%LOAD_SZ_INFO
%LOAD_SZ_INFO(pth)
%
%Inputs:
%    Optional:
%        subj          = string (full path to the seizure table)
%
%Outputs:
%    Required:
%        SeizureTable = table  (table with seizure information)
%        SeizureCell  = cell   (cell array with seizure information, contains an additional row at the top with the field names)
%        events       = cell   (cell with as many elements as subjects, where each element is the events structure for that subject)
%
%
% NOTES:
% 1. 
%% parse inputs
SZ_FILEPATH_DEFAULT = '/Volumes/shares/FRNU/dataWorking/sz/Seizure Table_06-02-17.xlsx';
if nargin < 1, subj = []; end
if nargin < 2, sz_filepath = SZ_FILEPATH_DEFAULT; end

%% Import the data as a cell
[~, ~, SeizureCell] = xlsread(fullfile(sz_filepath),'Sheet1');
SeizureCell = SeizureCell(:,1:11);
SeizureCell(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),SeizureCell)) = {''};

%% Import the data as a table
[~, ~, raw] = xlsread(fullfile(sz_filepath),'Sheet1');
raw = raw(2:end,1:11);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,4,9,10,11]);
raw = raw(:,[3,5,6,7,8]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
SeizureTable = table;

%% Allocate imported array to column variable names
SeizureTable.subjID = cellVectors(:,1);
SeizureTable.rawFile = cellVectors(:,2);
SeizureTable.date = data(:,1);
SeizureTable.date = datetime(SeizureTable.date,'ConvertFrom','excel');
%SeizureTable.Date = datetime(SeizureTable.Date,'ConvertFrom','excel1904');
SeizureTable.State = cellVectors(:,3);
SeizureTable.EEG_Change_s = data(:,2); 
SeizureTable.Ictal_Onset_s = data(:,3);
SeizureTable.Clinical_Onset_s = data(:,4);
SeizureTable.Termination_s = data(:,5);
SeizureTable.SeizureType = cellVectors(:,4);
SeizureTable.OnsetElectrodes = cellVectors(:,5);
SeizureTable.Pulled = cellVectors(:,6);

SeizureTable.eegfile = cell(size(SeizureTable,1),1);
SeizureTable.eegoffset_s = data(:,2); % this times Fs should be eegoffset

%% remove entries with missing files or no subjID
% SeizureTable(cellfun(@isempty,SeizureTable.subjID),:)=[];
SeizureTable(cellfun(@(x) ~strncmp(x,'NIH',3),SeizureTable.subjID),:)=[];
SeizureTable(cellfun(@isempty,SeizureTable.rawFile),:)=[];
SeizureTable=sortrows(SeizureTable,'subjID');
%% Clear temporary variables
clearvars data raw cellVectors R;

% Enforce subjID is length 6
subjID_char = char(SeizureTable.subjID);
subjID_char_numeric_part = subjID_char(:,4:end); % part after 'NIH'
subjID_cell_numeric_part = cellstr(subjID_char_numeric_part);
subjID_int_numeric_part  = cellfun(@str2double, subjID_cell_numeric_part);
subjID_char_numeric_pad  = sprintf('%03d\n', subjID_int_numeric_part);
subjID_cell_numeric_pad  = strsplit(subjID_char_numeric_pad,'\n');
subjID_cell_numeric_pad  = subjID_cell_numeric_pad(1:end-1); % remove empty
SeizureTable.subjID = strcat('NIH', subjID_cell_numeric_pad)';

%% create events structures (need to make this a separate code to get eegdir information)
% [inds names]=findgroups(SeizureTable.subjID);
names=unique(SeizureTable.subjID);

for k=1:length(names)
    events{k}=table2struct(SeizureTable(strcmp(SeizureTable.subjID,names{k}),:));
      
    %%%*** can save event to the right path here  
%     pth='';
%     save(fullfile(path,'events.mat'),events)
end

% Subject filter if passed
isSubjFound = false;
if ~isempty(subj)
    isSubjFound = ismember(subj, names);
    rowMask = strcmpi(SeizureTable.subjID, subj);
    SeizureTable = SeizureTable(rowMask,:);
    eventMask = strcmpi(names, subj);
    events = events(eventMask);
end

