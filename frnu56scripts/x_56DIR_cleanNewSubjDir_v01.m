%
%  Helper Script for FRNU56
%
%   
%   Standard procedure
%
%
%

dbstop if error


%% RUN THIS FUNCTION TO ORGANIZE the SUBJECT data_raw directory

SUBJ = 'NIH090';
%path2data_raw = sprintf('/Volumes/56A/UTAH_A/%s/data_raw', SUBJ); 
path2data_raw = sprintf('/Volumes/56D/UTAH_D/%s/data_raw', SUBJ);     %- office-local (could change this to server)
%path2data_raw = sprintf('/Volumes/72D-1/UTAH_D/%s/data_raw',SUBJ);     %- office-local (could change this to server)

SUCCESS = clean_rawFileDir56(path2data_raw, SUBJ);


%- redo raw file list on FRNU56... sync will copy that to 72 later
MAKE_JACKSHEETS = 1;  %
if MAKE_JACKSHEETS & contains(path2data_raw,'/Volumes/56') & SUCCESS,
    rootFRNU56dir = path2data_raw(1:strfind(path2data_raw,SUBJ)-1); %- path to root UTAH_X where this subject resides
    getBR_rawFileList(SUBJ,'',rootFRNU56dir, 0, 0); %- this call will also update local copies of jacksheet.
    
    
    %- and make a copy of all the potential behavioral/stim sessions in 56E (does this work?)
% The following commented out by SJ:
%     overrideDataPath56 = '/Volumes/56E/public/micro/';
%     processEEGdir = fullfile(overrideDataPath56,SUBJ,'manual_sorts','data_raw'); %- should be pointing to 56PUB/publich/micro or something like that
%     sessList = {};  JUST_CHECK = 0;
%     getBR_grabFileList(SUBJ,processEEGdir,rootFRNU56dir, sessList, JUST_CHECK);
end


%rootEEGdir = '/Volumes/JW24TB/data24TB/eeg_new';
%sessList = {};  JUST_CHECK = 0;
%getBR_grabFileList(subj,rootEEGdir,rootFRNU56dir, sessList, JUST_CHECK);

