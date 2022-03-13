function writeSortNotes_readme(path_to_readme)
%
% Writes a sortNotes_readme.txt for the inputted path
% This was previously written straight in extractMicrophysiology_v12b, which now calls this function
% instead.
% Now also called in updateUtahInfoCSV? Or just called in own script?
%
% INPUT:
%   path_to_readme: path to the sortNotes_readme.txt
%       ex) path_to_readme = '/Volumes/56PROC/micro_behavioral/micro_pristine/NIH029/processed_SPK/150607_1542_beh/sortNotes_readme.txt';
%
%
% Created by Samantha Jackson
%
% 1/15/2020 - Created by SJ
%
%



fid = fopen(path_to_readme,'w+');

fprintf(fid,'\n%s\n','Instructions for sorting and filling out sortNotes(000000_0000)_sortedBy??.xlsx');
fprintf(fid,'\n%s\n','STEP ONE: Setting up the sortNotes.xlsx');
fprintf(fid,'\t%s\n','**IMPORTANT NOTE: Do not edit the structure of the excel file, and do not write anything outside of the pre-established columns. ');
fprintf(fid,'\t%s\n','  If you have a general comment about the entire session, please enter it in the SortComment for the first channel.');
fprintf(fid,'\t%s\n','1) Rename sortNotes(000000_0000)_sortedBy??.xlsx to include your initials in place of the ????. Do not use any numbers.');
fprintf(fid,'\t%s\n','2) Open up sortNotes(000000_0000)_sortedByXX.xlsx and look at the columns ?GradeNPlay1? and ?GradeNPlay2?, which should be filled with ?????'); 
fprintf(fid,'\t%s\n','   This is where you will be grading each channel twice for potential sort quality. Delete these question marks now.');
fprintf(fid,'\n%s\n','STEP TWO: Grading in nPlay');
fprintf(fid,'\t%s\n','1) In Parallels Desktop, open nPlay and load in the .ns5 file from the nPlay_visualize folder');
fprintf(fid,'\t%s\n','2) Follow the nPlay set-up instructions on the Wiki');
fprintf(fid,'\t%s\n','3) You will be grading each channel twice: once in the first half of the recording (for GradeNPlay1) and the other in the second half (for GradeNPlay2)');
fprintf(fid,'\t%s\n','4) You should grade the channels upon the following criteria (+/- are okay):');
fprintf(fid,'\t\t%s\n','"A" = definitely one or more units with great SNR that should be easy to sort/isolate');
fprintf(fid,'\t\t%s\n','"B" = definitely one or more units with decent SNR. Definitely a unit, but might not be super easy to isolate');
fprintf(fid,'\t\t%s\n','"C" = very likely one or more units. Best SNR is decent to poor. Likely will get something in plexon, but poor isolation.');
fprintf(fid,'\t\t%s\n','"D" = possibly one or more units, could be worth checking in plexon');
fprintf(fid,'\t\t%s\n','"F" or Blank = definitely NO units, do not bother sorting.');
fprintf(fid,'\t%s\n','5) Use the ?toSort? column to help you indicate which channels you plan to sort.');
fprintf(fid,'\n%s\n','STEP THREE: Sorting in Plexon');
fprintf(fid,'\t%s\n','1) In Parallels Desktop, open Plexon and follow the directions in the Wiki for setting up.');
fprintf(fid,'\t%s\n','2) Load in the binary file for each channel you plan to sort one-by-one from the sort_reref folder and sort. ');
fprintf(fid,'\t%s\n','3) In the ?SortComment? column, include any notes that a future analyzer may find useful, such as the number of units and an indication of the quality');
fprintf(fid,'\t%s\n','   (clear units, noisy, cant sort, etc)');
fprintf(fid,'\t%s\n','Things to note while sorting:');
fprintf(fid,'\t\t%s\n','* The way the grading works is that you should treat each channel?s grade as the highest of the 2 grades. Sort ALL channels with an ?A? or ?B? grade, ');
fprintf(fid,'\t\t%s\n','  no matter what. You should then begin sorting Cs. If the Cs are promising, begin sorting Ds next. Only stop on a letter grade if you have sorted all ');
fprintf(fid,'\t\t%s\n','  the channels with higher grades. For example, you should not be sorting from a D channel if you have not finished sorting all of the Cs first.');
fprintf(fid,'\t\t%s\n','* If ?stimChan? has a number in it and ?dummySort? has ?????');
fprintf(fid,'\t\t\t%s\n','- If it is possible to sort, include a ?dummy? unit in your sort with the stim artifact, then place the number of the unit in place of the ?????');
fprintf(fid,'\t\t%s\n','* If you planned to sort a channel based on nPlay, but in Plexon it is unsortable');
fprintf(fid,'\t\t\t%s\n','- Export a blank .txt file with no units');
fprintf(fid,'\t\t%s\n','* If you come across something you do not understand or the data looks crazy, please ask SJ or JW.');
fprintf(fid,'\n%s\n','STEP FOUR: When you are done sorting');
fprintf(fid,'\t%s\n','1) Append ?_ready? to the end of the session folder name and email SJ.');
fprintf(fid,'\t%s\n','2) When SJ has grabbed the session for processing, it will automatically replace the ?_ready? with ?_(grabbed)?.');
fprintf(fid,'\t%s\n','3) When the spikeInfo is done, you will be able to see it in a newly-created folder called ?reref_sortedByXX? with your initials, which includes:');
fprintf(fid,'\t\t%s\n','* spikeInfo.mat (not aligned)');
fprintf(fid,'\t\t%s\n','* sort_log.log');
fprintf(fid,'\t\t%s\n','* sortFigs');
fprintf(fid,'\t\t%s\n','* sorts_complete/INCOMPLETE.txt: If incomplete, this will indicate which sorts you need to re-sort ');
fprintf(fid,'\t%s\n','4) When alignment is done, you will be able to see the above info as well as the aligned spikeInfo.mat in Shares/FRNU or in 56PUB/readOnly/eeg');

fclose(fid);


end