function utahMap = getUtahChannelMap(subj, SHOW_PLOT)
%
%   function getUtahChannelMap
%
%   Note: this funtion needs to be maually updated with each new subject so the proper map is pulled for that subject
%
%   inputs:  subject name  e.g. 'NIH029'
%            SHOW_PLOT (1) optional flag to turn off plot
%
%   outputs:  utahMapStruct
% utahMap.subj               = subj;
% utahMap.elecNumSpace1      = elecNumSpace1;      %- "electrode to space mapping"... not very helpful but what blackrock gives us
% utahMap.elecNumSpace1_wire = elecNumSpace1_wire; %- location of gold wire coming from array, required for determining rotation of map onto the brain
% utahMap.elecNumSpace2      = elecNumSpace2;      %- if a "multiport", with two arrays attached to a single digitizer, the second array info goes here
% utahMap.elecNumSpace2_wire = elecNumSpace2_wire;
% utahMap.stimPigtail2Chan   = stimPigtail2Chan;   %- if a Cereplex SI, capable of microstimulation, the mapping between stim pigtail ring and electrode
%
%- Output utahMap struct has these additional fields
% utahMap.chanNumSpace1      = chanNumSpace1;      %- what you want... a mapping between recorded channel and spatial layout
% utahMap.chanNumSpace2      = chanNumSpace2;
% utahMap.chanNumBPlist      -                     %- 2xN list of adjacent channel nummbers for bipolar referencing
% utahMap.stimPigtail2Elec   = stimPigtail2Elec;   %- mapping from stimulation pigtail ring (i.e. the stim channel on cerestim) and recorded electrode
%
%  JW 2/2019
%

if nargin<2, SHOW_PLOT=1; end


switch subj,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'NIH029' 'NIH030' 'NIH034' 'NIH036' 'NIH037' 'NIH039' 'NIH042' 'NIH046' 'NIH047' 'NIH050' 'NIH054' 'NIH057' 'NIH059'}
         %- Cereplex I, serial number 8352-XXXX,8413-XXXX 
         %- single port, 96 channels
        
        %- wiring from stim pigtail ring to electrode number
        stimPigtail2Elec = [nan nan];
        
        %- Electrode Numbering viewing from pad side
        elecNumSpace1 = [
            nan	88	78	68	58	48	38	28	18	nan
            96	87	77	67	57	47	37	27	17	8
            95	86	76	66	56	46	36	26	16	7
            94	85	75	65	55	45	35	25	15	6
            93	84	74	64	54	44	34	24	14	5
            92	83	73	63	53	43	33	23	13	4
            91	82	72	62	52	42	32	22	12	3
            90	81	71	61	51	41	31	21	11	2
            89	80	70	60	50	40	30	20	10	1
            nan	79	69	59	49	39	29	19	9	nan];
        elecNumSpace1_wire = [4 5];
        
        %- second array
        elecNumSpace2 = [];
        elecNumSpace2_wire = [];
        
        %- copy columns "Elec#" from Bank A, B, C, D
        elecNumMap = [
            68	67	36	40	32	31
            66	65	44	48	30	29
            72	71	52	56	28	27
            70	69	60	64	26	25
            76	75	35	39	24	23
            74	73	43	47	22	21
            80	79	51	55	20	19
            78	77	59	63	18	17
            84	83	34	38	16	15
            82	81	42	46	14	13
            88	87	50	54	12	11
            86	85	58	62	10	9
            92	91	33	37	8	7
            90	89	41	45	6	5
            96	95	49	53	4	3
            94	93	57	61	2	1];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'NIH064' 'NIH088'}
        %- Cereplex I, serial number 8411-014XXX
        %- multiport, no stim
        
        %- wiring from stim pigtail ring to electrode number
        stimPigtail2Elec = [nan nan];
        
        %- Electrode Numbering viewing from pad side
        elecNumSpace1 = [
            64	56	48	40	32	24	16	8
            63	55	47	39	31	23	15	7
            62	54	46	38	30	22	14	6
            61	53	45	37	29	21	13	5
            60	52	44	36	28	20	12	4
            59	51	43	35	27	19	11	3
            58	50	42	34	26	18	10	2
            57	49	41	33	25	17	9	1];
        elecNumSpace1_wire = [4 5];
        
        elecNumSpace2 = [
            128	120	112	104	96	88	80	72
            127	119	111	103	95	87	79	71
            126	118	110	102	94	86	78	70
            125	117	109	101	93	85	77	69
            124	116	108	100	92	84	76	68
            123	115	107	99	91	83	75	67
            122	114	106	98	90	82	74	66
            121	113	105	97	89	81	73	65];
        elecNumSpace2_wire = [68 69];
        
        %- copy columns "Elec#" from Bank A, B, C, D
        elecNumMap = [
            100	99	68	72	64	63	32	28
            98	97	76	80	62	61	24	20
            104	103	84	88	60	59	16	12
            102	101	92	96	58	57	8	4
            108	107	67	71	56	55	31	27
            106	105	75	79	54	53	23	19
            112	111	83	87	52	51	15	11
            110	109	91	95	50	49	7	3
            116	115	66	70	48	47	30	26
            114	113	74	78	46	45	22	18
            120	119	82	86	44	43	14	10
            118	117	90	94	42	41	6	2
            124	123	65	69	40	39	29	25
            122	121	73	77	38	37	21	17
            128	127	81	85	36	35	13	9
            126	125	89	93	34	33	5	1];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'NIH062' 'NIH063' 'NIH066' 'NIH069' 'NIH072' 'NIH074' 'NIH076' 'NIH070' 'NIH079' 'NIH086'}
        %- Cereplex SI.  Serial Numbers 10658-XXX; 
        %- multiport with microstimulation
        
        %- Electrode Numbering viewing from pad side
        elecNumSpace1 = [
            64	56	49	42	31	22	15	8
            63	53	48	41	30	21	23	7
            62	55	47	40	33	20	14	6
            61	54	46	39	29	19	13	5
            60	52	36	38	28	25	12	4  %% wire coming out from here
            59	51	45	37	27	18	10	3
            58	35	44	34	26	17	11	2
            57	50	43	32	24	16	9	1];
        elecNumSpace1_wire = [4 5];
        
        elecNumSpace2 = [
            128	120	112	105	97	89	82	73
            127	119	111	104	96	88	78	72
            126	117	110	103	81	87	79	71
            125	116	109	102	94	86	77	70
            124	115	95	101	93	80	76	69 %% wire coming out from here
            123	114	108	100	92	85	75	68
            122	118	107	99	91	84	67	66
            121	113	106	98	90	83	74	65];
        elecNumSpace2_wire = [69 70];
        
        %- copy columns "Elec#" from Bank A, B, C, D
        elecNumMap = [
            84	77	104	97	64	39	13	20
            71	65	118	117	45	51	26	32
            91	85	111	105	57	63	6	12
            78	72	98	119	38	44	19	25
            66	92	120	112	50	56	31	5  %- NOTE, raw xls file for NIH069 says 66;9;120; etc... looks like that 9 should be 92
            86	79	106	99	62	37	11	18
            73	67	125	121	43	49	24	30
            93	87	113	107	55	61	4	10
            80	74	100	126	36	42	17	23
            68	94	122	114	48	54	29	3
            88	81	108	101	60	35	9	16
            75	69	127	123	41	47	22	28
            95	89	115	109	53	59	2	8
            82	76	102	128	34	40	15	21
            70	96	124	116	46	52	27	1
            90	83	110	103	58	33	7	14];
        
        
        %- Define the mapping between stim outputs (chan 1-15) and recording channel :: UNIQUE FOR EACH ARRAY
        stimPigtail2Elec =                                  [1     2     3     4     5     6     7     8     9     10    11    12    13    14    15];
        if strcmp(subj,'NIH062'), stimPigtail2Elec(2,:) =   [26   118    36    11    67    25    23    53    33    78    35    80    81    95   119];  end % SN 000002       (NIH062) [got from map file]
        if strcmp(subj,'NIH063'), stimPigtail2Elec(2,:) =   [23   118    53    67    78    95   119    11    36    80    35    25    26    33    81];  end % SN 000001       (NIH063) [got from map file]
        if strcmp(subj,'NIH066'), stimPigtail2Elec(2,:) =   [78    33    26   118   119    35    11    25    80    36    23    67    95    81    53];  end % SN 10658-000003 (NIH066) [got from map file]
        if strcmp(subj,'NIH069'), stimPigtail2Elec(2,:) =   [67    11    33    78    26   119    35    95   118    80    53    36    81    25    23];  end % SN 10658-000007 (NIH069) [got from map file]
        if strcmp(subj,'NIH072'), stimPigtail2Elec(2,:) =   [81   119    95    11    53    80   118    25    78    33    67    35    36    26    23];  end % SN 10658-000011 (NIH072) [got from map file]
        if strcmp(subj,'NIH074'), stimPigtail2Elec(2,:) =   [26    78    80    33    36    95    53    67    23    35    25    11    81   119   118];  end % SN 10658-00005  (NIH074) [got from map file]
        if strcmp(subj,'NIH076'), stimPigtail2Elec(2,:) =   [118   95    11    23    81    80    53    36    67    78    25    26    33    35   119];  end % SN 10658-00013  (NIH076) [got from map file]
        if strcmp(subj,'NIH070'), stimPigtail2Elec(2,:) =   [78    81    67   119    80    33    23    36    35    26    95    11    25    53   118];  end % SN 10658-000002 (NIH070) [got from map file]
        if strcmp(subj,'NIH079'), stimPigtail2Elec(2,:) =   [78    53    33    35    23    11    25    80    95   119    36    81    26   118    67];  end % SN 10658-000014 (NIH079) [got from map file]
        if strcmp(subj,'NIH086'), stimPigtail2Elec(2,:) =   [26    22	  5	   23     3	   21	  4	    6	 25	    24	  7	   10	 28	   27	  9];  end % SN 10658-000015 (NIH086) [got from map file]
        %if strcmp(subj,'NIH0??'), stimPigtail2Elec(2,:) =   [1     2     3     4     5     6     7     8     9     10    11    12    13    14    15];  end % SN 10658-?????? (NIH0??) [need to get from map file]

        
        if numel(stimPigtail2Elec)~=30, fprintf('\n ERROR: stimPigtail2Elec is unique for each subject and must be defined'); keyboard; error('missing stim pigtail map'); end
        stimPigtail2Elec = stimPigtail2Elec'; %- transpose back to expected form
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        fprintf('ERROR: subject %s NOT FOUND in getUtahChannelMap', subj);
        keyboard;
        
end


%- infostring for packaged readme
infoStr = 'Readme for utahMap data structure (created by getUtahChannelMap):';
infoStr = sprintf('%s\n     NOTE: run getUtahChannelMap with "SHOW_PLOT" flag to see image of channel locations in physical space',infoStr);
infoStr = sprintf('%s\n subj               = subject string (e.g., NIH029)',infoStr);
infoStr = sprintf('%s\n elecNumSpace1      = electrode numbers mapped to physical space (from blackrock xlsx file)',infoStr);
infoStr = sprintf('%s\n elecNumSpace1_wire = electrodes adjacent to the wire bundle in physical space (from blackrock xlxs)',infoStr);
infoStr = sprintf('%s\n elecNumSpace2      = for multi-port (two arrays), elecNumSpace for second array',infoStr);
infoStr = sprintf('%s\n elecNumSpace2_wire = for multi-port (two arrays), elecNumSpace2_wire for second array',infoStr);
infoStr = sprintf('%s\n stimPigtail2Elec   = for Cereplex SI (microstim), mapping between stimulator channels and electrode numbers (from Blackrock xlsx)',infoStr);
infoStr = sprintf('%s\n elecNumMap         = mapping from electrode numbers to output banks on the Cerelex digitizer, which feed recording channels',infoStr);

%-package data
utahMap.readme             = infoStr;
utahMap.subj               = subj;
utahMap.elecNumSpace1      = elecNumSpace1;
utahMap.elecNumSpace1_wire = elecNumSpace1_wire;
utahMap.elecNumSpace2      = elecNumSpace2;
utahMap.elecNumSpace2_wire = elecNumSpace2_wire;
utahMap.stimPigtail2Elec   = stimPigtail2Elec;
utahMap.elecNumMap         = elecNumMap;


%- convert to spatial plot and (optionally) plot the data
utahMap = convertUtahElecs2Chans(utahMap,SHOW_PLOT);

