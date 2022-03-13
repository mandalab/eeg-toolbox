function [actualName_ALL, actualCode_ALL]=getChanNames_BR(elecFileBR,eegFileBR)
%
% getChanNames_BR.m
%
% Takes in a pair of raw blackrock files and returns the names and codes for
% all recorded channels. Can handle multiple raws, if they are labeled
% with suffixes 'ieeg1' and 'ieeg2', rather than only 'ieeg'. Writes them
% out into a ChanNames.txt that is used in lieu of NK's .21E file
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
% 12/2017 - melkalliny adapted getChanNames into this BR version
% 11/2018 - JW tweaked to deal with inconsistent labeling of NSP1 and 2...
%         now a NSP_chanCount.txt is ALWAYS created (even if jsut one NSP), and this code uses that file to figure out if both NSPs have been processed
% 04/2019 - JW converts to jacksheetBR method
% 0225/2020 - CB + SJ: found writing issue where writing jackBR would add another row of blanks when overwriting previous table
%                   - solution: delete jacksheet_localBR before writing out
%                       again.
% 9/2020 - SJ: minor improvements, readtableSafe; Keep ChanNames.txt but delete it later in the pipeline
%VERBOSE = 0 ; %

[thisPath, thisFile, thisSuff] = fileparts(eegFileBR);


%% NEW WAY: use jacksheetBR as the place where channel name changes are stored
%jackBR = dir(fullfile(thisPath,'jacksheetBR_local*.csv'));  %- would like to use complete here, but for some subjects the NSx files were renamed
jackBRfile = fullfile(thisPath,'jacksheetBR_local.csv');  %- would like to use complete here, but for some subjects the NSx files were renamed
if ~exist(jackBRfile,'file')
    jackBR = makeJacksheetBR(thisPath,'',0,1,0);
else
    jackBR = readtableSafe(jackBRfile); %SJ changed from readtable
end

flag_updateJack = 0; %- flag... will get set to 1 below if ChanNameNews are added and jacksheet needs to be resaved to disk


% remove blank spaces (shouldnt have to do this, but OK)
jackBR.ChanName = deblank(jackBR.ChanName);


%- find the channels associated with this NSX file (usualy ns2 for physio)
thisNSX = [thisFile thisSuff];
iNSX2   = [];
iNSX1   = find(strcmp(jackBR.FileName,thisNSX));

if contains(thisNSX,'INST') % SJ 9/20/19, add in functionality for ieeg (used to just be INST)
    INST_ls = {'INST0', 'INST1'};
elseif contains(thisNSX,'ieeg')
    INST_ls = {'ieeg1', 'ieeg2'};
end

if contains(thisNSX,INST_ls{1})
    iNSX2 = find(strcmp(jackBR.FileName,strrep(thisNSX,INST_ls{1},INST_ls{2})));
elseif contains(thisNSX,INST_ls{2})
    iNSX2 = find(strcmp(jackBR.FileName,strrep(thisNSX,INST_ls{2},INST_ls{1})));
end

% if contains(thisNSX,'INST0')
%     iNSX2 = find(strcmp(jackBR.FileName,strrep(thisNSX,'INST0','INST1')));
% elseif contains(thisNSX,'INST1')
%     iNSX2 = find(strcmp(jackBR.FileName,strrep(thisNSX,'INST1','INST0')));
% end
iNSX = union(iNSX1,iNSX2);  %- combine NSP1 and NSP2 indicies into the jacksheet, so all NS2 channels are linked


% Now addd all the ain channels, so they are added even if not in the current NSx filetype (i.e, if physio is NS2 but pulses are NS3)
for ain=1:16
    iChanNSx = find( (strcmp(jackBR.ChanName, sprintf('ain%d',ain)) |  strcmp(jackBR.ChanName, sprintf('ainp%d',ain))) & ~contains(jackBR.FileName,'.nev'));
    if length(iChanNSx)>1
        fprintf('\n shouldnt happen'); keyboard; 
    elseif length(iChanNSx==1)
        newName = sprintf('DC%02d',ain+8); %- push DC names to ChanNameNew
        if ~strcmp(jackBR.ChanNameNew{iChanNSx},newName)
            jackBR.ChanNameNew{iChanNSx} = newName;
            flag_updateJack = 1;
        end
        iNSX = cat(1,iNSX,iChanNSx);
    end
end
    

iNSX = unique(iNSX); %accounts for case in which ain5 is detected both in first auto step and in explicit ain-detection step

if length(iNSX)==0
    fprintf('\n No channels found with expected NSx file... either file isnt present or physio is on a different chan?');
    keyboard;
end


actualName_ALL={}; actualCode_ALL={};
for ii=1:length(iNSX)
    if ~strcmp(jackBR.ChanNameNew{iNSX(ii)},'-')
        actualName_ALL{ii} = jackBR.ChanNameNew{iNSX(ii)}; %- use a new name
    else
        actualName_ALL{ii} = jackBR.ChanName{iNSX(ii)};    %- use the embedded name
    end
    actualCode_ALL{ii} = sprintf('PhysicalChan=%d',jackBR.PhysicalChan(iNSX(ii)));
end



% now, lets write this out into chanNames.txt if it's not already there
if flag_updateJack
    % First, delete it (CB + SJ found writing issue where it would add
    % another row of blanks when overwriting previous table)
    delete(jackBRfile)
    writetable(jackBR,jackBRfile);
end


%% EVENTUALLY CUT chanNames.txt from the blackrock pipeline... for now let it live
USE_CHANNAME_BR = 1; %SJ: For now, still generate it but then delete it later in the pathway
if USE_CHANNAME_BR
    dashSplit = strfind(elecFileBR,'/');
    pathToFolder = elecFileBR(1:dashSplit(end)-1);
    pathToChanCounts = sprintf('%s/NSP_ChanCounts.txt',pathToFolder);
    pathToChanNames = sprintf('%s/chanNames.txt',pathToFolder);
    createTxtFileBR = 0;
    if ~exist(pathToChanNames, 'file') %SJ
        createTxtFileBR = 1;
    end
    
    
    %- need to create or append 2nd NSP to the ChanNames textfile....
    if createTxtFileBR  % first time pulling chan namesç?
        %- first NSP or only NSP
        fid = fopen(pathToChanNames,'w');
        fprintf(fid,'%s\n',actualName_ALL{:});
        fclose(fid);
        fid = fopen(pathToChanCounts,'w+');
        temp = [length(iNSX1)];
        fprintf(fid,'%d\n',temp);
        if length(iNSX2) > 0 %#ok<*ISMT>
            temp = [length(iNSX2)];
            fprintf(fid,'%d\n',temp);
        end
        fclose(fid);
    end
end %if USE_CHANNAME_BR


return;


