close all; 
clear all;

homeDir = '/Volumes/shares/FRNU/dataworking/dbs/';


patients(1).subj='DBS050';
patients(1).T_date = '170224';
patients(1).dataFile = {'LT4D1.117', 'RT5D2.254'};
patients(1).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(2).subj='DBS051';
patients(2).T_date = '170510';
patients(2).dataFile = {'LT1D0.327'; 'RT1D1.998'};
patients(2).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(3).subj='DBS053';
patients(3).T_date = '170526';
patients(3).dataFile = {'RT1D0.547'; 'LT2D-1.364'; 'LT2D-0.994'}; % session 1 was stopped and restarted as session 2 at a lower depth
patients(3).tasks = {'MemoryDecision'; 'MemoryDecision'; 'MemoryDecision'};

patients(4).subj='DBS054'; 
patients(4).T_date = '170915';
patients(4).dataFile = {'garbagesession'; 'RT2D-2.067'}; % File names for the Left side are improperly labeled. Could not export left side (session 0). 
patients(4).tasks = {'garbage'; 'MemoryDecision'};

patients(5).subj='DBS055';
patients(5).T_date = '170927';
patients(5).dataFile = {'LT1D-0.226'};
patients(5).tasks = {'MemoryDecision'};

patients(6).subj='DBS056';
patients(6).T_date = '171213';
patients(6).dataFile = {'LT1D0.819'; 'RT1D-0.091'};
patients(6).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(7).subj='DBS057';
patients(7).T_date = '180328';
patients(7).dataFile = {'LT4D0.977'};
patients(7).tasks = {'MemoryDecision'};

patients(8).subj='DBS058';
patients(8).T_date = '180509';
patients(8).dataFile = {'LT2D-1.222'; 'RT1D-0.089'};
patients(8).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(9).subj='DBS059';
patients(9).T_date = '180523';
patients(9).dataFile = {'LT1D2.084';'RT4D1.985'};
patients(9).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(10).subj='DBS061';
patients(10).T_date = '180720';
patients(10).dataFile = {'LT2D1.021'; 'RT1D0.068'};
patients(10).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(11).subj='DBS062';
patients(11).T_date = '181017';
patients(11).dataFile = {'LT3D-2.128'; 'RT2D0.444'};
patients(11).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(12).subj='DBS063'; % GPI recording
patients(12).T_date = '181214';
patients(12).dataFile = {'LT1D-1.281'};
patients(12).tasks = {'MemoryDecision'};

patients(13).subj='DBS064';
patients(13).T_date = '190227';
patients(13).dataFile = {'LT2D-0.243', 'RT2D1.575'};
patients(13).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(14).subj='DBS065';
patients(14).T_date = '190315';
patients(14).dataFile = {'LT3D0.060', 'RT2D2.613'}; 
patients(14).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(15).subj='DBS066';
patients(15).T_date = '190410';
patients(15).dataFile = {'LT2D1.021', 'RT1D2.505'};
patients(15).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(16).subj='DBS067';
patients(16).T_date = '190828';
patients(16).dataFile = {'LT3D-0.831', 'RT2D0.520'}; 
patients(16).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(17).subj='DBS068';
patients(17).T_date = '191002';
patients(17).dataFile = {'RT2D0.905', 'LT3D0.577'};
patients(17).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(18).subj='DBS069';
patients(18).T_date = '191016';
patients(18).dataFile = {'LT8D1.530', 'RT1D2.005'};
patients(18).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(19).subj='DBS070';
patients(19).T_date = '191120';
patients(19).dataFile = {'LT4D-2.261', 'RT2D1.307'};
patients(19).tasks = {'MemoryDecision'; 'MemoryDecision'};

patients(cellfun('isempty',{patients.subj}))=[];
for patient=1:length(patients)
    patients(patient).subj
    DBSPrepAndAlign(homeDir,patients(patient).subj, patients(patient).T_date, patients(patient).dataFile, patients(patient).tasks); % the new method with no_split. DBS041 and beyond
end
