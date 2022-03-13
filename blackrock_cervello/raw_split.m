function raw_split(subj, rootEEGdir, EEG_file, BATCH)
%function allTags = nk_split(subj, nk_dir, tagNameOrder, varargin)
% nk_split - Splits an nk .EEG datafile into separate channels into
% the specified directory.
%
% FUNCTION:
%    nk_split(subj,nk_dir,EEG_file,...)
%
% INPUT ARGs:
% subj = 'TJ022'
% rootEEGdir = '/Users/damerasr/damerasr/data/eeg/'
% EEG_file = '/Users/damerasr/damerasr/data/eeg/160614_1008/x.EEG'
%
% Output
%	output is saved out to: output_dir = '/Users/damerasr/damerasr/data/eeg/[SUBJ]/eeg.noreref/' set in line 49
%
%
% Edited from previous version so that manual input of the ouput_dir would not be necessary
%
% 12/2013 ... now uses jacksheetMaster.txt to guide channel number outputs
% 10/2015 ... now can handel "new" EEG-1200 extended file formats
% 08/2016 MST outputs using channel names instead of numbers
% 07/2017 MST Check for shift channels, modify .21E's
% 01/2018 ME converted from nk_split to raw_split that executes different split code based on recording machine
% 12/2018 JW cleaned up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% JW thinks none of this needs to happen here... already happening in each of the split functions
% %%- Load the jacksheetMaster, or create it, or give a warning that it can't be created...
% subjDir             = fullfile(rootEEGdir, subj);
% rawDir              = fullfile(subjDir, 'raw');
% [rawSubDir,~,~]     = fileparts(EEG_file); %- path to the raw file being processed, so FRNU/data/eeg/NIHXXX/raw/121401_1415/
% jackMaster_file_new = fullfile(subjDir, 'docs/jacksheetMaster.csv');
% 
% % create jacksheetMaster if it doesn't exist or hasn't looked at raws
% if ~exist(jackMaster_file_new, 'file')
%     createMasterJack(subj, rootEEGdir);
% end
% 
% jacktable = getJackTable(subj, rootEEGdir);
% jackMaster_names = jacktable.chanName;
% jackMaster_chans = [1:length(jackMaster_names)];
% 
% 
% %- read the raw_info file to identify channels intentionally EXCLUDED from jacksheetmaster... no reason to give a warning if those are found below
% raw_info_file = fullfile(subjDir, 'docs/raw_info.csv');
% raw_info      = readtableSafe(raw_info_file);
% chanDontSplit = raw_info.chanName(find(raw_info.in_jackSheetMaster==0));
% chanDontSplit{end+1} = ''; %- dont report error or try to split channels with no name


% splits raw according to file suffix, ie recording system present
if     contains(EEG_file,'.ns2'), split_BR(subj, rootEEGdir, EEG_file);
    
elseif contains(EEG_file,'.EEG'), split_NK(subj, rootEEGdir, EEG_file);
    
elseif contains(EEG_file,'.TRC'), split_CV(subj, rootEEGdir, EEG_file, BATCH);
    
end


end
% end raw_split function




%#ok<*NOCOL>