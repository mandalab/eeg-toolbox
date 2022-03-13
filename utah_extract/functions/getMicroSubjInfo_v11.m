function subjectInfo = getUtahSubjInfo(subj, overrideDataPath56)
% Get info about patients who were implanted with utah arrays
% subj: subject number (e.g. 29 or 'NIH029')
%
%  device numbers: %- 1:2 ATL (usually MTG, but also STG),  3:5 depths  7 parietal  9 frontal
%
%
%  optional parameter:
% overrideDataPath56:  if data resides locally, instead of on FRNU56, put the path to the subjects here
%                     e.g. '/Volumes/JW24TB/data24TB/localFRNU56/';
%     alternatively, if you want to work on FRNU72, pass "72" for the override path
%
%
% oldChanName and newChanName should be numLocationsx2 - see NIH050 for
% example. if there is a 128chan utah array in same location, make sure
% both newChanNames contain a conserved string (I.e. AMTG and PMTG), so
% that during creation of spikeInfos, it is recognized that these two
% arrays should go into one spikeInfo
%
%  MicroDeviceNum code: %- 1:2 ATL,  3:5 depths  7 parietal  9 frontal
%
% AJ created 2015-ish
% ME tweaked 2018
% JW tweaked again 2019
%

if nargin<2,
    overrideDataPath56 = '';
end

subjectInfo = struct;
nspSuffix    = {};   %- file suffix is a unique identifier for which ns5/ns6 file came from which NSP.  (e.g., INST0, ieeg, utah, micro)
nspSuffixAlt = {};   %- "B" version for rare cases when patient switched systems and suffix changed (e.g. NIH064 switched from matlab recording to Cervello)
dataPath56   = [];   %- Which "56" drive the raw data resides on.   (e.g., '/Volumes/56A/UTAH_A)
utahChanCnt  = [];   %- number of channels recorded in each device.   96 or 128 for utah arrays, 16 for depths.
oldChanName  = {};   %- name root as found in NSX file  :  (e.g.,   utah1-96)
newChanName  = {};   %- corrected/updated channel string:  (e.g.,   utah1-96 -->  utah_MGG01 -96)
deviceNum    = [];   %- numerical code for utah or depth location: %- 1:2 ATL (usually MTG, but also STG),  3:5 depths  7 parietal  9 frontal
adjacentMacros = []; %- This needs to be updated for most subjects... should be strings of channel names, or grid numbers, but in cell arrays to distinguish multiple arrays
description  = {};

if isnumeric(subj)
    subj = sprintf('NIH%03d',subj);
end


switch subj,
    
    case 'NIH029',
        dataPath56  = '/Volumes/56A/UTAH_A';
        %- NOTE: in the case of the same old-channel name being on two NSPs (e.g., NIH029), the first utah array is the one coming from NSP1 in jacksheetBR
        nspSuffix   = {'ieeg',         'utah'        };  %- NIH029, NSP1 "ieeg" recorded Left AMTG array
        utahChanCnt = [96              96            ];
        oldChanName = {'utah',     1;  'utah',      1};
        newChanName = {'utahAMTG', 1;  'utahSTG',   1};  %- NIH029, NSP1 (ieeg) recorded Left AMTG array, which died on day 2. NSP2 (utah) recorded Superior temporal gyrus array, which became great
        deviceNum   = [1,              2             ];
        adjacentMacros = [10,11,18,19];
        description = {'Cereplex I serial #000101 old style, with ground leads; new and old digital hubs (more anterior and superior). Left Superior MTG',
            'Cereplex I serial #000104 new style, with chassis grounded. Left Anterior MTG'};
        
        
    case 'NIH030',
        dataPath56  = '/Volumes/56A/UTAH_A';
        nspSuffix   = {'utah'       };  %- NIH030, NSP2 "utah" recorded utah array
        utahChanCnt = [96           ];
        oldChanName = {'utah',     1};
        newChanName = {'utahMTG',  1};
        deviceNum   = [1            ];
        adjacentMacros = [20,21,28,29];
        description = {'Cereplex I serial #000102 old style; right MTG (?);'};
        
        
    case 'NIH034',
        dataPath56  = '/Volumes/56A/UTAH_A';
        %- NOTE: in the case of the same old-channel name being on two NSPs (e.g., NIH029), the first utah array is the one coming from NSP1 in jacksheetBR
        nspSuffix   = {'utah1',        'utah2'       };   %- NIH034, utah1 --> ASTG
        utahChanCnt = [96              96            ];
        oldChanName = {'utah',     1;  'utah',      1};
        newChanName = {'utahASTG', 1;  'utahPMTG',  1};   %- NIH034, utah1 --> ASTG (array failed early; first 32 chan never worked),  utah2 --> PMTG
        deviceNum   = [1,              2             ];
        adjacentMacros = [10,17,19,26];
        description = {'Cereplex I serial #14374; posterior (MTG according to KZ), worked well  ? possible mix up on which was which',
            'Cereplex I serial #14377; anterior (STG), failed early'};
        
        
    case  'NIH036'
        dataPath56  = '/Volumes/56A/UTAH_A';
        nspSuffix   = {'utah1'      };  %- NIH036, NSP2 "utah1" recorded utah array
        utahChanCnt = [96           ];
        oldChanName = {'utah',     1};
        newChanName = {'utahMTG',  1};
        deviceNum   = [1            ];
        adjacentMacros = [];
        description = {'Cereplex I serial #14376; MTG placement .  Non-functioning (no single units, possibly all noise)'};
        
        
    case  'NIH037'
        dataPath56  = '/Volumes/56A/UTAH_A';
        nspSuffix   = {'utah_m',       'utah_u'      };  %- NIH037, NSP1=utah_m --> AMTG,  NSP2=utah_u --> PMTG
        utahChanCnt = [96              96            ];
        oldChanName = {'utah',      1; 'utah',      1};
        newChanName = {'utahAMTG',  1; 'utahPMTG',  1};  %- NIH037, utah_m --> AMTG,  utah_u --> PMTG
        deviceNum   = [1,              2             ];
        description = {'Cereplex I serial #14380;  Anterior MTG placement, Marked  worked well  entire time',
            'Cereplex I serial #14375;  Posterior MTG placement; Unmarked  worked well for a few days then failed'};
        
        
    case 'NIH039'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'utah_m',       'utah_u'      };   %- NIH039, utah_m --> PMTG,  utah_u --> AMTG
        utahChanCnt = [96              96            ];
        oldChanName = {'utah',      1; 'utah',      1};
        newChanName = {'utahPMTG',  1; 'utahAMTG',  1};   %- NIH039, utah_m --> PMTG,  utah_u --> AMTG
        deviceNum   = [1,              2             ];
        adjacentMacros = [];
        description = {'Cereplex I serial #14381; Posterior MTG, Marked pigtail.  Worked entire time, but with low unit yield', 
            'Cereplex I serial #14392; Anterior MTG, unmarked pigtail.  Failed after 6 days.  Had a few units before that'};
        
        
    case 'NIH042'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'utahT',        'utahF'       };   %- NIH042, utahT --> MTG,  utahF --> Frontal
        nspSuffixAlt= {'utah_u',       'utah_m',     };  %- NIH042, first few sessions... "marked" pigtail was frontal  utah_u --> MTG,  utah_m --> Frontal. Was it ever INST0?
        utahChanCnt = [96              96            ];
        oldChanName = {'UT', 1 ;       'UF', 1       };
        newChanName = {'utahMTG',  1;  'utahFRNT',  1};
        deviceNum   = [1,              9             ];
        adjacentMacros = [];
        description = {'Cereplex I serial #14398; Frontal, Marked pigtail.  1.0mm array (all previous arrays were 1.5mm)',
            'Cereplex I serial #14397; Anterior MTG, unmarked pigtail. 1.0mm array'};
        
        
    case 'NIH046'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'utah'          };  %- NIH046, NSP2 "utah" recorded utah array
        utahChanCnt = [96              ];
        oldChanName = {'utah',     1   };   %-
        newChanName = {'utahPAR',  1   };
        deviceNum   = [7               ];
        adjacentMacros = [];
        description = {'Cereplex I serial # 14382; right parietal (angular gyrus).  1.0mm array?'};
        
        
    case 'NIH047'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'utah'          };  %- NIH047, NSP2 "utah" recorded utah array
        utahChanCnt = [96              ];
        oldChanName = {'utah',     1   };
        newChanName = {'utahMTG',  1   };
        deviceNum   = [1               ];
        adjacentMacros = [];
        description = {'Cereplex I serial #14388; left MTG.  1.0mm array'};
        
        
    case 'NIH048'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'ieeg',          'ieeg',         'ieeg',         };   %- NIH048, NSP1 'ieeg' captured both depths
        utahChanCnt = [16               11              5               ];
        oldChanName = {'chan',     69;  'chan', 85;     'chan', 97      };  %- first microwire depth; recorded without preAmp.
        newChanName = {'microAD',   1;  'microPD',  1;  'microPD',  12  };  %- PDM was put into clinical jacksheet chan 85-95, 97-101 (skipping 96 for EKG)
        deviceNum   = [3                4               4               ];
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH049'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'micro'       , 'micro'        };  %- NIH049, NSP1 "micro" recorded depths
        utahChanCnt = [16              16             ];
        oldChanName = {'LADM', 1     ; 'chan', 1      };  %- recorded without preAmp; initally named "chan", then LADM. possible unit activity observed during recording (JW)
        newChanName = {'microLAD',  1; 'microLAD',  1 };
        deviceNum   = [3               3              ];  %- same device recorded with two different oldChanNames
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH050'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'utah',          'micro',        'micro',      };   %- NIH050, NSP1 'micro' captured both depths, NSP2 'utah' got the utah,
        utahChanCnt = [96                16              16           ];
        oldChanName = {'utah',     1;   'MAID', 1;      'MPID',      1};
        newChanName = {'utahMTG',  1;   'microAID',  1; 'microPID',  1};
        deviceNum   = [1                 3               4            ];  %- 1:2 ATL,  3:4 depths  7 parietal  9 frontal
        adjacentMacros = [];
        description = {'Cereplex I serial SN 14401;  right MTG, ~3-4cm from pole.  1.0mm array'};
        
        
    case 'NIH051'
        dataPath56  = '/Volumes/56B/UTAH_B';
        nspSuffix   = {'micro'        };   %- NIH051, NSP1 'micro' captured depths, left caudate depth
        utahChanCnt = [16             ];
        oldChanName = {'LLCDM',       1};
        newChanName = {'microLLCD',  1};
        deviceNum   = [5              ];  %- 1:2 ATL,  3:5 depths  7 parietal  9 frontal
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH052'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'micro',          'micro'          };   %- NIH052, NSP1 'micro' captured both depths
        utahChanCnt = [16                 16              ];
        oldChanName = {'RAHDm',      1;  'LAHDm', 1       };
        newChanName = {'microRAHD',  1;  'microLAHD', 1   };
        deviceNum   = [3                  4               ];
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH053'  %-
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1',        'INST1'        };   %- NIH053, NSP2 'INST1' captured both depth micros
        utahChanCnt = [16              16             ];
        oldChanName = {'LADm',     1;  'RADm',1       };
        newChanName = {'microLAD', 1;  'microRAD',1   };
        deviceNum   = [3               4              ];
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH054'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1'      };   %- NIH054, NSP2 'INST1' captured utah
        utahChanCnt = [96           ];
        oldChanName = {'utah',     1};
        newChanName = {'utahMTG',  1};
        deviceNum   = [1            ];
        adjacentMacros = [];
        description = {'Cereplex I serial SN 14402;  right MTG, ~3-4cm from pole.  1.0mm array'};
        
        
    case 'NIH057'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1',     'INST1'      };   %- NIH057, NSP2 'INST1' captured utah
        utahChanCnt = [96           96           ];
        oldChanName = {'utah',     1;'el',    129};   %- at some point ALL channel names lost (macro and micro) and set to "el", with utah starting at el129.
        newChanName = {'utahMTG',  1;'utahMTG',  1};
        deviceNum   = [1            1            ];
        adjacentMacros = [];
        description = {'Cereplex I serial SN 14404;  left MTG, ~3-4cm from pole.  1.0mm array'};
        
        
    case 'NIH059'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1',      'INST1'      };   %- NIH059, NSP2 'INST1' captured utah.   For first 2 sessions utah channels were named "chan", but still on INST1
        utahChanCnt = [96            96           ];
        oldChanName = {'utah',     1;'chan',     1};
        newChanName = {'utahMTG',  1;'utahMTG',  1};
        deviceNum   = [1             1            ];
        adjacentMacros = [];
        description = {'Cereplex I serial SN 14410;  left MTG, ~3-4cm from pole.  1.0mm array'};
        
        
    case 'NIH060'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1',        'INST1',           'INST1',         'INST1',      };   %- NIH060, NSP2 captured all 3 depths
        utahChanCnt = [16              16                  16              16            ];
        oldChanName = {'LADm',      1; 'LAHDm',      1;   'LPHDm',      1; 'LPHD',      1};   %- at some point mistakenly deleted the "m" from LPHD
        newChanName = {'microLAD',  1; 'microLAHD',  1;   'microLPHD',  1; 'microLPHD', 1};
        deviceNum   = [3                4                  5               5             ];
        adjacentMacros = [];
        description = {''};
        
        
    case 'NIH062'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'INST1',         'INST1'       };   %- NIH062, NSP2 'INST1' captured utah.
        utahChanCnt = [64                64           ]; %- beginning of multi-ports (two 64 channel arrays with different referecnes)
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahMTGX',  1;  'utahMTGY',  1}; %- NIH062: not sure which was anterior/posterior
        deviceNum   = [1                 2            ];
        adjacentMacros = [];
        description = {'Cereplex SI serial SN 000002;  left MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan)'};
        
        
    case 'NIH063'
        dataPath56  = '/Volumes/56C/UTAH_C';             %- NIH063s itah only worked for a few session in the begining.  no task i think
        nspSuffix   = {'INST1',         'INST1'       }; %- NIH063, NSP2 'INST1' captured utah.
        nspSuffixAlt= {'utah',          'utah'       };  %- NIH063, also used 'utah'
        utahChanCnt = [64                64           ]; %-
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahMTGX',  1;  'utahMTGY',  1}; %- NIH063: not sure which was anterior/posterior
        deviceNum   = [1                 2            ];
        adjacentMacros = [];
        description = {'Cereplex SI serial SN 000001;  right MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan)'};
        
        
    case 'NIH064'
        dataPath56  = '/Volumes/56C/UTAH_C';
        nspSuffix   = {'utah',         'utah',         };  %- NIH064, NSP2 'utah' captured utah in beginning. Switch to Cervello at the end
        nspSuffixAlt= {'INST0',        'INST0',        };  %- NIH064, NSP2 'utah' captured utah in beginning. Switch to Cervello at the end
        utahChanCnt = [64               64             ];
        oldChanName = {'utah',      1;  'utah',     65 };
        newChanName = {'utahMTGX',  1;  'utahMTGY',  1 };  %- NIH064: not sure which was anterior/posterior
        deviceNum   = [1                 2             ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex I serial SN 14415;  right MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT (2x64 chan)'};
        
        
    case 'NIH066'
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH066, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahPMTG',  1;  'utahAMTG',  1};  %- for NIH066: ELECTRODES 1-64 were anterior (which means chan 66-128), ELECTRDOES 65-128 were posterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN 10658-003;  left MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan)'};
        
        
    case 'NIH067'
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'utah'            };  %- NIH067, implanted old-design single port array.  never worked, only have an entry here for bookkeeping
        utahChanCnt = [96                ];  %-
        oldChanName = {'utah',       1;  };
        newChanName = {'notValid',   1;  };  %- never worked, so label "notValid"
        deviceNum   = [1                 ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {''};
        
        %case 'NIH068'  % no utah, no micros... but mistakenly recorded "utah" data with all zeros for the first day or two.
        
    case 'NIH069',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH069, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];
        oldChanName = {'utah_pmtg', 1;  'utah_amtg', 1};
        newChanName = {'utahAMTG',  1;  'utahPMTG',  1}; %- for NIH069: ELECTRODES 1-64 were posterior (means chan 65-128), ELECTRODES 65-128 were anterior
        deviceNum   = [1                 2            ]; %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN 10658-0007;  right MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan). Chan1-64 is POSTERIOR, 65-128 ANTERIOR'};
        
        
    case 'NIH071',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',        'INST0',           'INST0',         'INST0',        'INST0'        };   %- NIH071, first benke-friend style depths (PMT).  Implanted 3 depths. Only recorded 2 for a bit and missed the "m" in the names for a bit
        utahChanCnt = [10              10                  10              10              10             ];   %-     9 microwires and 1 ref, which is channel 4
        oldChanName = {'RPCDm',     1; 'RACDm',      1;   'ROFDm',      1; 'RACD',      1; 'ROFD',     1  };   %- at some point mistakenly deleted the "m" from LPHD
        newChanName = {'microRPC',  1; 'microRAC',   1;   'microROF',   1; 'microRAC',  1; 'microROF', 1  };
        deviceNum   = [3                4                  5                4               5             ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'PMT Benke Fried elelectrodes with updated Blackrock headstage'};
        
    case 'NIH072',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH072, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahPMTG',  1;  'utahAMTG',  1};  %- for NIH072: string on anterior array, so ELECTRODES 1-64 were anterior (which means chan 65-128), ELECTRDOES 65-128 were posterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN 10658-00011;  right MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan).  String on anterior array'};
              
    case 'NIH074',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH074, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahAMTG',  1;  'utahPMTG',  1};  %- for NIH074: string on posterior array, so ELECTRODES 1-64 were posterior (which means chan 65-128), ELECTRDOES 65-128 were anterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN10658-00005 [note: Ti case label said 0003];  left MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan).  String on posterior array'};

    case 'NIH076',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH076, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahAMTG',  1;  'utahPMTG',  1};  %- for NIH076: *string on posterior array*, so ELECTRODES 1-64 were posterior (which means chan 65-128), ELECTRDOES 65-128 were anterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN10658-00013 [note: Ti case label was blank];  left MTG, ~2-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan).  String on posterior array'};

    case 'NIH070',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH070, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahAMTG',  1;  'utahPMTG',  1};  %- for NIH070: *string on posterior array*, so ELECTRODES 1-64 were posterior (which means chan 65-128), ELECTRDOES 65-128 were anterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [18,           12            ];  %- NIH070, under TG18 for anterior, under TG12 for posterior
        description = {'Cereplex SI serial SN10658-000002 [note: Ti case label was 000005];  right MTG, ~4-6cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan).  String on posterior array'};
        
    case 'NIH079',
        dataPath56  = '/Volumes/56A/UTAH_A';
        nspSuffix   = {'INST0',         'INST0'       };  %- NIH079, NSP1 'INST0' captured utah in the upper channels (129-256)
        utahChanCnt = [64                64           ];  %- started marking the location of array 1 and 2
        oldChanName = {'utah',      1;  'utah',     65};
        newChanName = {'utahAMTG',  1;  'utahPMTG',  1};  %- for NIH070: *string on posterior array*, so ELECTRODES 1-64 were posterior (which means chan 65-128), ELECTRDOES 65-128 were anterior
        deviceNum   = [1                 2            ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [17,           18            ];  %- NIH070, under TG18 for anterior, under TG12 for posterior
        description = {'Cereplex SI serial SN10658-000002 [note: Ti case label was -------];  right MTG, ~3-4cm from pole.  1.0mm array. MULTI-PORT w/STIM (2x64 chan).  String on posterior array'};
        
    case 'NIH081',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0',           'INST0',         'INST0'        };   %- NIH081, second benke-friend style depths (PMT).  Implanted 4 depths. 
        utahChanCnt = [10                10                 10               10            ];   %-     9 microwires and 1 ref, which is channel 4
        oldChanName = {'microRSLP',6;   'microRILP',6;     'microRIMP', 6;   'microRSMP', 6};   %- numbering started at 6 -> 15 instead of 1
        newChanName = {'microSLP', 1;   'microILP', 1;     'microIMP',1;     'microSMP',  1};
        deviceNum   = [3                 4                  5                6             ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'PMT Benke Fried elelectrodes with updated Blackrock headstage'};
        
    case 'NIH082',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0',        };   %- NIH081, second benke-friend style depths (PMT).  Implanted 2 depths. 
        utahChanCnt = [10                10             ];   %-     9 microwires and 1 ref, which is channel 4
        oldChanName = {'micro_PURPLE',6; 'micro_GREEN',6};   %- numbering started at 6 -> 15 instead of 1
        newChanName = {'microPCD',    1; 'microROF',   1};   %- No check after surgery, so not sure which corresponds to which Macro, keep same names
        deviceNum   = [3                 4              ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'PMT Benke Fried elelectrodes with updated Blackrock headstage'};            

    case 'NIH086',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'      };   %- NIH086
        utahChanCnt = [64                64          ];   %-     
        oldChanName = {'utah',        1;'utah',    65};   %- 
        newChanName = {'utahMTGX',    1;'utahMTGY', 1};  %- No string, so not sure which is anterior/posterior
        deviceNum   = [1                 2           ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN10658-000015 [TI case label was the same (000015)]. Left anterior temporal lobe, MTG, ~2-4 cm from pole. 1.0mm array. MULTI-PORT w/STIM (2x64 chan).'};        % (32) Cereplex SI serial SN10658-000015 [TI case label ???].   Left anterior temporal lobe, MTG, ~2-4 cm from pole.    1.0mm array. MULTI-PORT w/STIM (2x64 chan).   
        % NOTE: Pt also had microburst wires BUT they were not recorded from due to noise
        % from Utah summaries: 1 benke-fried PMT implanted in left LPFC.  Revised amp.  Iniitally didn?t plug in
        % SJ: I think it is actually supposed to be LACC (pt does not have LPFC electrodes)

    case 'NIH088',
        dataPath56  = '/Volumes/56D/UTAH_D';
        nspSuffix   = {'INST0',         'INST0'      };   %- NIH081, second benke-friend style depths (PMT).  Implanted 2 depths. 
        utahChanCnt = [64                64          ];   %-     
        oldChanName = {'utah',        1;'utah',    65};   %- 
        newChanName = {'utahMTGX',    1;'utahMTGY', 1};  %- 
        deviceNum   = [1                 2           ];  %- 1 and 2 for ATL,  7 for parietal, 9 for frontal
        adjacentMacros = [];
        description = {'Cereplex SI serial SN8411-014411 [TI case label was 014411]. Right anterior temporal lobe, MTG, ~?-? cm from pole. 1.0mm array. MULTI-PORT (2x64 chan).'};        %  
        % Arrays did not have a string identifying one from the other. Both are in the MTG in the ATL. the one on the right (assuming they did not get twisted when coming out of the digitizer) is more posterior, and should sit right posterior to electrode 17 on the grid. the one on the left sits just anterior to electrode 17 on the grid.
        
    otherwise
        fprintf('\n WARNING: Subject %s does not have microInfo in getMicroSubjInfo --> returning empty variable', subj);
        subjectInfo = [];
        return;
        
end

%- should we force the utahMapping to be done at this point as well?  One issue is not all subjects in this document have utahs
if any(contains(newChanName(:,1),'utah')),
    utahMap = getUtahChannelMap(subj,0);
else
    utahMap = [];
end

%if ~exist('fileSuffix'),
%    fileSuffix = {''};
%end

%- populate output structure
subjectInfo.subj           = subj;
subjectInfo.dataPath56     = dataPath56 ;
subjectInfo.nspSuffix      = nspSuffix;     %- INST0,1
subjectInfo.nspSuffixAlt   = nspSuffixAlt;  %- utah   <-- intent is for this to come into play for subjects that switch NSP filenames
subjectInfo.utahChanCnt    = utahChanCnt;
subjectInfo.oldChanName    = oldChanName;
subjectInfo.newChanName    = newChanName;  %- 1:2 ATL,  3:4 depths  7 parietal;   8 frontal
subjectInfo.deviceNum      = deviceNum;
subjectInfo.adjacentMacros = adjacentMacros;
subjectInfo.description    = description;  %- copied from JW's UTAH update word doc
subjectInfo.utahMap        = utahMap;      %- structure with information about how to map recorded channels to space


%- data quailty check
if length(deviceNum)~=length(unique(deviceNum)),
    %fprintf('\n Heads up... one or more device numbers are repeated. \nOnly ok if both have the same new Chan Str \n(this is a way to deal with changing jacksheet names for one device)');
    
    uniqueDevNums = unique(deviceNum);
    for ii=1:length(uniqueDevNums),
        iDev = find(deviceNum == uniqueDevNums(ii));
        numDev     = length(iDev); %- how many devices with this device num
        numNewName = length( unique( newChanName(iDev,1) ) ); %- now many different newChanName with this device num (should always be 1)
        if length(numNewName)>1,
            fprintf('\n Uh oh... device numbers should probably be unique for each element. This will affect global reref');
            keyboard;
            error('\n fix this');
        end
    end
end


if ~isempty(overrideDataPath56),
    %- JW local hack, shoudn't affect anybody but JW
    %localFRNU56 = '/Volumes/JW24TB/data24TB/localFRNU56/';
    if exist(fullfile(overrideDataPath56,subj),'dir'), %#ok<*NOCOL>
        %fprintf('\n HEADS UP: converting %s to %s in getUtahSubjInfo (working local on JWs machine)\n',subjectInfo.dataPath56,localFRNU56);
        subjectInfo.dataPath56 = overrideDataPath56;
        
    elseif overrideDataPath56==72 | strcmp(overrideDataPath56,'72'),
        fprintf('\n HEADS UP... running processes on FRNU72');
        subjectInfo.dataPath56 = regexprep(dataPath56,'56','72');
    end
end


