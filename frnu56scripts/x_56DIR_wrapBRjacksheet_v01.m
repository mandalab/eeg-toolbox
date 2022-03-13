%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   JW's calls to eeg_toolbox utah processing and output cunctions
%
%   this one is a wrapper to getBR_rawFileList... making them on 56 and copying them to a local EEG folder if desired
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

dbstop if error


FORCE_REMAKE = 1;  %- use this to rename the rawFileList even if one already exists.. good if you dont want to delete the old ones so notes are preserved


rootEEGdir  = '';                        %office-local
%rootEEGdir  = '/Volumes/JW24TB/data24TB/eeg_new';                        %office-local
%rootEEGdir = '/Volumes/JW24TB/dataEEGpristine/update/eeg';


root56dirList = {'/Volumes/56A/UTAH_A/' '/Volumes/56B/UTAH_B/' '/Volumes/56C/UTAH_C/' '/Volumes/56D/UTAH_D/'};
%root56dirList = root56dirList(1);

%%%  FRNU56 DRIVE                        SUBJECTS ON THAT DRIVE             SUBJECTS through First Pass (fixed all errors)
%root56dir = '/Volumes/56A/UTAH_A/';  fullSubjList = [29 30 34 36 37];    subj2process = [];
%root56dir = '/Volumes/56B/UTAH_B/';  fullSubjList = [39 42 46 47 50];   subj2process = [];
%root56dir = '/Volumes/56C/UTAH_C/';  fullSubjList = [38 48 49 51:61];   subj2process = [];
%root56dir = '/Volumes/56E/UTAH_E/';  fullSubjList = [62:69];            subj2process = [69:-1:62];


%-
%subjNum2Process = [1:10 29:69];  %- will only work on these subjects: 1:10 is TRE
subjNum2Process = 76;



for i56=1:length(root56dirList),
%for i56=length(root56dirList):-1:1,
    
    %- find the subject folder, could end in "micro" or something, so use dir to fish it out
    subjNIH = dir(fullfile(root56dirList{i56},'NIH*'));  %-
    subjTRE = dir(fullfile(root56dirList{i56},'TRE*'));  %- no TRE in here now, but they could get added back
    
    subjList = {subjNIH.name subjTRE.name};
    
    
    fprintf('\n Processing %s subjects in %s', length(subjList),root56dirList{i56});
    
    for iSubj=1:length(subjList)
        thisSubj = subjList{iSubj};
        subjNum = str2num(thisSubj(4:6));  %- expecting format NIHXXX
        if ~any(subjNum2Process==subjNum)
            fprintf('\n heads up: subject %s found but will not be processed', thisSubj);
            %keyboard;
            
        else
            fprintf('\n     processing %s',thisSubj);
            
            %rawRootDir = fullfileEEG(root56dirList{i56},thisSubj,'data_raw');
            
            
            if FORCE_REMAKE | ~exist(fullfile(root56dirList{i56},thisSubj,'notes',sprintf('rawFileList_%s.xlsx',thisSubj))),
                
                shortTest = 0;  onlyReturnUtahInfo = 0;
                getBR_rawFileList(thisSubj,rootEEGdir,root56dirList{i56}, shortTest, onlyReturnUtahInfo); %- this call will also update local copies of jacksheet.
             
            else
                fprintf('  --> rawFileList.xlsx already created'); 
            end
            
            %sessList = {};  JUST_CHECK = 0;
            %getBR_grabFileList(subj,rootEEGdir,rootFRNU56dir, sessList, JUST_CHECK);
            
        end
        
        
    end
end

