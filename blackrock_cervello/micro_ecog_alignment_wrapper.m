


subjList = {'NIH029','NIH030','NIH034','NIH036','NIH037','NIH039','NIH042','NIH046','NIH047','NIH048','NIH049','NIH050','NIH052','NIH053','NIH054','NIH057','NIH059','NIH060','NIH062','NIH063','NIH064','NIH066','NIH069','NIH071'};

% align
for j=23:25
    alignSpikesWithEcog_ManualORAuto({subjList{j}},'/Volumes/FRNU_EXT/micro_ecog_alignment_fullset/')
    fprintf('\n\n\n\n\n\n\n\n\n\naligned micro-ecog for %s\n\n\n\n\n\n\n\n',subjList{j})
end


% copy everything into a folder to put on 56
for j=24
    fprintf('\nbeginning to copy %s ......\n',subjList{j})
    if ~exist(sprintf('/Volumes/FRNU_EXT/micro_ecog_finalset/%s/micro',subjList{j}))
        mkdir(sprintf('/Volumes/FRNU_EXT/micro_ecog_finalset/%s/micro',subjList{j}))
    end
    copyfile(sprintf('/Volumes/FRNU_EXT/micro_ecog_alignment_fullset/%s/micro',subjList{j}),sprintf('/Volumes/FRNU_EXT/micro_ecog_finalset/%s/micro',subjList{j}))
    fprintf('done with %s\n',subjList{j})
end