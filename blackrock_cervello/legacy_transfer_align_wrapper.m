       

% full_subj_list  = {'NIH029', 'NIH030', 'NIH034', 'NIH037', 'NIH039', 'NIH042', 'NIH046', 'NIH047', 'NIH050', 'NIH054', 'NIH057'};

% 39 - 1=posterior, 2=anterior
% 42 - 1=MTG, 9=frontal
% 50 - 1=ATL, 3=MAID, 4=MPID




subj_list  = {'NIH029'};
% about to be done with NIH034
fid = fopen('/Volumes/FRNU_EXT/micro_Ecog_alignment_fullset/progress.txt','w');
for i=1:length(subj_list)
subj = subj_list{1,i};
fprintf('\n\n\n\n\n\n--starting loop for %s--\n\n\n',subj);
fprintf(fid,'\nstarting loop for %s\n',subj);
subjectInfo.dataPath56 = '/Volumes/56E/public/micro/';
sessFolders  = ''; %all sessions in data_processed
grabSortedData(subj,subjectInfo.dataPath56,sessFolders); %- probably should make this an automatic part of extractSpikeInfo...
fprintf(fid,'\ngrabbed sorted data for %s\n',subj);
extractSpikeInfo_v3(subj,subjectInfo.dataPath56,sessFolders);
fprintf(fid,'\nextracted spike infos for %s\n',subj);
raw_fileList_pullMatchingMicros(subj,'/Volumes/FRNU_EXT/micro_ecog_alignment_fullset/','/Volumes/56E/public/micro/')
subjInCell = {subj};
alignSpikesWithEcog_ManualORAuto(subjInCell,'/Volumes/FRNU_EXT/micro_ecog_alignment_fullset/')
fprintf('\n\n\n--ending loop for %s--\n\n\n',subj);
fprintf(fid,'\nending loop for %s\n',subj);
end
fclose(fid);




subj = 'NIH054';
extractSpikeInfo_v3(subj,subjectInfo.dataPath56,sessFolders);

