function [T] = documentCRVsess_micro(subj,ecog_sess,micro_sess,samp_rate_TRC,fileOut)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if isempty(fileOut)
    yeswrite = false;
    
    %T = table([],[],[],[],'VariableNames',{'subject','ecog_sess','micro_sess','sampRate'});
else
    yeswrite = true;
    if exist(fileOut,'file')
        Tmat = readtableSafe(fileOut);
        Tmat = table2cell(Tmat);
        
        Tmat = [Tmat ; [{subj},{ecog_sess},{micro_sess},{samp_rate_TRC}]];
    else
        Tmat = [{subj},{ecog_sess},{micro_sess},{samp_rate_TRC}];
    end
    
end

T = cell2table(Tmat,'VariableNames',{'subject','ecog_sess','micro_sess','sampRate'});

if yeswrite
    writetableSafe(T,fileOut);
end

end

