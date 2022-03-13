function [actualName_ALL, actualCode_ALL]=getChanNames(elecFile,eegFile)
%
% getChanNames.m
% 
% Takes in a pair of raw files and returns the names and codes for
% all recorded channels.
% 
% INPUT ARGs:
% subj = 'TJ022'
%
% Output
%  actualName_ALL --- cell array of strings of channel names
%  actualCode_ALL --- cell array of strings corresponding with .21E row code (JW thinks)
%
%
% 10/2015 - updated to handle NEW EEG-1200 file format
% 07/2017 - MST updated to not return '' entries
% 12/2017 - melkalliny updated to provide capacity for Blackrock files


[parent,filestem] = fileparts(eegFile);
    
if     contains(filestem,'ieeg')     | contains(filestem,'INST'),   
    [actualName_ALL, actualCode_ALL] = getChanNames_BR(elecFile,eegFile);
    
elseif contains(filestem,'cervello') | contains(filestem,'EEG'),
    eegFile = strrep(eegFile,'.EEG','.TRC');
    [actualName_ALL, actualCode_ALL] = getChanNames_CV(eegFile);
    
else % NK
    [actualName_ALL, actualCode_ALL] = getChanNames_NK(elecFile,eegFile);
    
end


end