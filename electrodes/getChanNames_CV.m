function [actualName_ALL, actualCode_ALL]=getChanNames_CV(eegFileCV)
%
% getChanNames.m
%
% Takes in a pair of raw Cervello files and returns the names and codes for
% all recorded channels. Writes them out into a ChanNames.txt that is used in
% lieu of NK's .21E file.
%
% Note: Can handle multiple raws, but not setup quite yet. Needs a few
% tweaks, mentioned in the code
%
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
% 01/2018 - melkalliny adapted getChanNames_CV into this Cervello version
% 07/08/2020 - SJ added creation of jackSheetCV_local.csv + add in DC in
%           chanNameNew
%         - get rid of ChanNames.txt
% 09/29/2020 - SJ Centralizing place where channel names are extracted from jacksheetBR/CV- now calling new function getChansFromJack
%
% open EEG file
fid = fopen(eegFileCV);
if fid==-1, fprintf('\n ERRROR:  cant open .EEG file %s',eegFileCV); keyboard; end
fclose(fid); %jw added.. open wasn't matched. consider switching to exist

[thisPath, thisFile, thisSuff] = fileparts(eegFileCV);

jackCVfile = fullfile(thisPath,'jacksheetCV_local.csv');
if ~exist(jackCVfile,'file')
    jackCV = makeJacksheetCV(thisPath,0);
else
    jackCV = readtableSafe(jackCVfile);
end

flag_updateJack = 0; %- flag... will get set to 1 below if ChanNameNews are added and jacksheet needs to be resaved to disk

% Now addd all the ain channels, so they are added even if not in the current NSx filetype (i.e, if physio is NS2 but pulses are NS3)
for ain=1:16
    iChan = find((strcmp(jackCV.ChanName, sprintf('ain%d',ain)) |  strcmp(jackCV.ChanName, sprintf('ainp%d',ain))));
    if length(iChan)>1
        fprintf('\n shouldnt happen'); keyboard; 
    elseif numel(iChan) == 1
        newName = sprintf('DC%02d',ain+8); %- push DC names to ChanNameNew
        if ~strcmp(jackCV.ChanNameNew{iChan},newName)
            jackCV.ChanNameNew{iChan} = newName;
            flag_updateJack = 1;
        end
        %iNSX = cat(1,iNSX,iChan);
    end
end

if flag_updateJack
    % First, delete it (CB + SJ found writing issue where it would add
    % another row of blanks when overwriting previous table)
    delete(jackCVfile)
    writetable(jackCV,jackCVfile);
end

% temp = strfind(eegFileCV,'/');
% txtFileCV = eegFileCV(1:temp(end));
% existTxtFileCV = exist(fullfile(txtFileCV,'ChanNames.txt'));

%txtFileCV = fullfile(thisPath,'ChanNames.txt'); %SJ


%SJ new:
% Used to pull in CV create the chanNames.txt straight from that... Since we do that in makeJacksheetCV,
% now we just need to use the jacksheetCV_local :)
% First, iterate over jackCV.ChanName and use ChanNameNew unless it is blank, in which case fill it in
% with ChanName

%exclude = {'din' 'ain17' 'ain18' 'ain19' 'ain20'};

%SJ: This code is replacing the code below because I have centralized it to a new function
%(getChansFromJack)
[actualName_ALL, actualCode_ALL] = getChansFromJack(jackCV);

% actualName_ALL = cell(numel(jackCV.ChanName),1);
% actualCode_ALL = cell(numel(jackCV.ChanName),1);
% 
% for ch = 1:numel(jackCV.ChanName)
%     if ~strcmp(jackCV.ChanNameNew{ch},'-')
%         actualName_ALL{ch} = jackCV.ChanNameNew{ch};
%     else
%         actualName_ALL{ch} = jackCV.ChanName{ch};
%     end
%     actualCode_ALL{ch} = sprintf('PhysicalChan=%d',jackCV.PhysicalChan_guess(ch));
% end
% 
% % Now we need to remove the 'exclude' channels (din*, ain17-20)
% exclude_chans = ~ismember(actualName_ALL,regexpCell(actualName_ALL,'^(din.*|ain1[7-9]|ain20)$'));
% actualName_ALL = actualName_ALL(exclude_chans);
% actualCode_ALL = actualCode_ALL(exclude_chans);





% SJ: The following is commented out because it should be taken care of by the above or by makeJacksheetCV:

% % if existTxtFileCV ~= 2
% if ~exist(txtFileCV,'file') % means its first time pulling chan names
%     chansGrab=[];
%     % pull EEG raw data
%     EEG = readCervelloTRC(eegFileCV,chansGrab);  %usually pulled in uV
%     actualName_ALL = {EEG.electrodes.tagName}';
%     % remove these chans that we won't end up splitting
%     % remember to do so also while actually splitting
%     exclude = {'din' 'ain17' 'ain18' 'ain19' 'ain20'};
%     updatedChans = actualName_ALL;
%     for y= 1:length(exclude)
%         temp = strfind(updatedChans, exclude{y});
%         updatedChans = updatedChans(find(~not(cellfun('isempty', temp))));
%     end
%     % replace with the names that will be found in elementInfo
%     ainToDC = find(not(cellfun('isempty', strfind(lower(updatedChans), 'ain'))));
%     if length(ainToDC) ~=4
%         fprintf('cervello didnt have 4 ain chans. weird. if theres ain1, probably OK')
%         keyboard;
%     end
%     % SJ: This should already be done in makeJacksheetCV now
% %     for y= 1:length(ainToDC)
% %         updatedChansOld = updatedChans{ainToDC(y)};
% %         updatedChansNew = strrep(updatedChansOld,'ain','DC');
% %         updatedChansNew(end:end+1) = num2str(sprintf('%02d',y+8));
% %         updatedChans{ainToDC(y)} = updatedChansNew;
% %     end
%     actualName_ALL = updatedChans;
%     % give em #s
%     actualCode_ALL = {};
%     for name = 1:length(actualName_ALL)
%         actualCode_ALL{1,name} = sprintf('%0.4d',name);
%     end
%     
%     %% this code was used in blackrock to handle multiple NSPs. Will need to be tweaked
%     %     % pull DC names as well
%     %     if isempty(strfind(eegFileBR,'ieeg1'));  %only for single NSP or for 2nd NSP
%     %     NS3FileDir = eegFileBR; NS3FileDir = strrep(NS3FileDir,'ns2','ns3');
%     %     NS3data = openNSx(NS3FileDir); % For analog DC channel
%     %     rawData_NS3 = getNSxData(NS3data);
%     %     allChanNames = {NS3data.ElectrodesInfo.Label}';
%     %     if length(allChanNames) ~=4;
%     %         fprintf('there should be 4 DC channels here. theres not')
%     %         keyboard;
%     %     end
%     %     withoutDcEnd = length(actualCode_ALL);
%     %     for name = 1:length(allChanNames)
%     %         temp = allChanNames{name,:}; temp=temp(find(temp>=' '));
%     %         % replace ain# with DC# that it carries
%     %         temp = sprintf('DC%0.2d',name+8);
%     %         actualName_ALL{1,withoutDcEnd+name} = temp;
%     %         actualCode_ALL{1,withoutDcEnd+name} = sprintf('%0.4d',name+withoutDcEnd);
%     %     end
%     %     end
%     
%     % write out fileNames into .txt
%     % if there is an cervello1 in the name, i.e 2 nsps, write out different
%     % .txts before aggregating into the same one during cervello2
%     
%     %% No such thing as cervello1 or cervello2 files
% %     if ~isempty(strfind(eegFileCV,'cervello1'))
% %         fid = fopen(fullfile(txtFileCV,'ChanNames1.txt'),'w');
% %         fprintf(fid,'%s\n',actualName_ALL{:});
% %         fclose(fid);
% %     elseif ~isempty(strfind(eegFileCV,'cervello2'))
% %         % aggregate chan names from both NSPs into one .txt
% %         fid = fopen(fullfile(txtFileCV,'ChanNames1.txt'),'r');
% %         names = textscan(fid,'%s');
% %         fclose(fid);
% %         actualName_ALL_1 = {};
% %         actualName_ALL_1{1} = names{1}';
% %         actualName_total = [actualName_ALL_1{1} actualName_ALL];
% %         fid = fopen(fullfile(txtFileCV,'ChanNames.txt'),'w');
% %         fprintf(fid,'%s\n',actualName_total{:});
% %         fclose(fid);
% %         delete(fullfile(txtFileCV,'ChanNames1.txt'))
% %         % write out a marker of which chans were in which NSP
% %         fid = fopen(fullfile(txtFileCV,'NSP_ChanCounts.txt'),'w+');
% %         temp = [length(actualName_ALL_1{1}) length(actualName_ALL)-4];
% %         fprintf(fid,'%d\n',temp);
% %         fclose(fid);
% %     else
%         fid = fopen(fullfile(txtFileCV,'ChanNames.txt'),'w');
%         fprintf(fid,'%s\n',actualName_ALL{:});
%         fclose(fid);
% %    end
% else % chan names have been pulled before, so just retrieve them
%     fid = fopen(fullfile(txtFileCV,'ChanNames.txt'),'r');
%     names = textscan(fid,'%s');
%     actualName_ALL = {};
%     actualName_ALL{1} = names{1}';
%     actualCode_ALL = {};
%     for name = 1:length(actualName_ALL{1})
%         actualCode_ALL{1,name} = sprintf('%0.4d',name);
%     end
%     fclose(fid);
%     actualName_ALL = actualName_ALL{1};
% end
% 
% % what are 'good' and 'bad' electrodes, as in getChanNames_NK ??


end
