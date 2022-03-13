function subjectInfo = getPalRamSubjInfo(subj,sessNum)
% Get info about patients who did palram with utah arrays
% keyboard
br_fileName = '';
nk_sessName = '';
subjectInfo = struct;

if strcmp(subj,'NIH029') || (isnumeric(subj) && subj == 29)
    adjacentGridElecs = [10,11,18,19];
    utahChanCount = 96;
    palRamSessNums = [1,2,3];
    dataLoc = '56A/UTAH_A';
    suffix = 2;
    if sessNum == 1
        br_fileName = '150608_1133_utah';
        nk_sessName = '150608_1113';
    elseif sessNum == 2
        br_fileName = '150621_124631_utah';
        nk_sessName = '150621_1243';
    elseif sessNum == 3
        br_fileName = '150621_160648_utah2';
        nk_sessName = '150621_1605';
    end
    
    
elseif strcmp(subj,'NIH030') || (isnumeric(subj) && subj == 30)
    adjacentGridElecs = [20,21,28,29];
    utahChanCount = 96;
    palRamSessNums = [0,1];
    dataLoc = '56A/UTAH_A';
    suffix=2;
    if sessNum == 0
        br_fileName = '150628_1758_utah2';
        nk_sessName = '150628_1758';
    elseif sessNum == 1
        br_fileName = '150629_1611_utah2';
        nk_sessName = '150629_1604';
    end
    
    
elseif strcmp(subj,'NIH034') || (isnumeric(subj) && subj == 34)
    adjacentGridElecs = [10,17,19,26];
    utahChanCount = 96;
    palRamSessNums = [0,1,2];
    dataLoc = '56A/UTAH_A';
    suffix=2;
    if sessNum == 0
        br_fileName = '150926_1707_utah2';
        nk_sessName = '150926_1704';
    elseif sessNum == 1
        br_fileName = '150927_1020_utah2';
        nk_sessName = '150927_1021';
    elseif sessNum == 2
        br_fileName = '150929_1043_utah2';
        nk_sessName = '150929_1042';
    end
    
    
elseif strcmp(subj,'NIH037') || (isnumeric(subj) && subj == 37)
    adjacentGridElecs = [3,4,11,12];
    utahChanCount = 96;
    palRamSessNums = [0,1,3,4,5];
    suffix=1;
    %attentionTaskSessNums = 0:21;
    dataLoc = '56A/UTAH_A';
    if sessNum == 0
        br_fileName = {'160126_1011_utah_m_beh.ns6','160126_1011_utah_u_beh.ns6'};
        nk_sessName = '160126_1009';
    elseif sessNum == 1
        br_fileName = {'160126_1535_utah_m_beh.ns6','160126_1535_utah_u_beh.ns6'};
        nk_sessName = '160126_1535';
    elseif sessNum == 3
        br_fileName = '160202_1012_utah_m_beh.ns6';
        nk_sessName = '160202_1008';
    elseif sessNum == 4
        br_fileName = '160205_1103_utah_m_beh.ns6';
        nk_sessName = '160205_1102';
    elseif sessNum == 5
        br_fileName = '160211_1112_utah_m_beh.ns6';
        nk_sessName = '160211_1114';
    end
    
    
elseif strcmp(subj,'NIH039') || (isnumeric(subj) && subj == 39)
    adjacentGridElecs = [];
    utahChanCount = 96;
    palRamSessNums = [0,1,2,3,4];
    dataLoc = '56B/UTAH_B';
    suffix=1;
    if sessNum == 0
        br_fileName = {'160306_1620_utah_m_beh.ns6','160306_1620_utah_u_beh.ns6'};
        nk_sessName = '160306_1619';
    elseif sessNum == 1
        br_fileName = {'160307_1905_utah_m_beh.ns6','160307_1905_utah_u_beh.ns6'};
        nk_sessName = '160307_1904';
    elseif sessNum == 2
        br_fileName = {'160308_1622_utah_m_beh.ns6','160308_1622_utah_u_beh.ns6'};
        nk_sessName = '160308_1619';
    elseif sessNum == 3
        br_fileName = {'160309_1032_utah_m_beh.ns6','160309_1032_utah_u_beh.ns6'};
        nk_sessName = '160309_1032';
    elseif sessNum == 4
        br_fileName = '160312_0928_utah_u_beh.ns6';
        nk_sessName = '160312_0928';
        suffix=2;
    end
    
    
elseif strcmp(subj,'NIH042') || (isnumeric(subj) && subj == 42)
    adjacentGridElecs = [];
    utahChanCount = 96;
    palRamSessNums = [0,1,2];
    dataLoc = '56B/UTAH_B';
    suffix=1;
    nk_sessName = '160615_1116';
    if sessNum == 0
        br_fileName = '160615_1101_utahT_beh.ns6';
        nk_sessName = '160615_1116';
    elseif sessNum == 1
        br_fileName = '160619_1527_utahT_beh.ns6';
        nk_sessName = '160619_1526';
    elseif sessNum == 2
        br_fileName = '160625_1753_utahT_beh.ns6';
        nk_sessName = '160625_1753';
    end

        
else
    error('NOT FOUND');
    
end

subjectInfo.subj = subj;
subjectInfo.allSessNums = palRamSessNums;
subjectInfo.br_fileName = br_fileName;
subjectInfo.nk_sessName = nk_sessName;
subjectInfo.adjacentGridElecs = adjacentGridElecs;
subjectInfo.utahChanCount = utahChanCount;
subjectInfo.dataLoc = dataLoc;
subjectInfo.suffix = suffix;


end

