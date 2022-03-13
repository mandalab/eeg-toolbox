 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   JW's calls to eeg_toolbox alignment functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------------------%
%----------------------  Choose the subject of interest and root EEG directory -----------------------%
%clc; clear all;

subj       = 'NIH078';



%rootEEGdir = '/Volumes/56E/private/.eeg_pristine/update/eeg';
%rootEEGdir = '/Volumes/56PROC/private/.eeg_pristine/update/tempUpdate';

%rootEEGdir = '/Volumes/JW24TB/dataEEGpristine/update/eeg';

%rootEEGdir = '/Volumes/EEG_56STAGE-1/stage/working';
rootEEGdir = '/Volumes/EEG_56STAGE-1/stage/ready';




%% SOME HELPFUL FUNCTIONS THAT CAN BE RUN ON the subj/rootDir... comment in/out depending on need:

%subjectDirCheck(subj,rootEEGdir);

%atlas_simplify(subj, rootEEGdir, 'allow_ties', 0);
%atlas_simplify(subj,rootEEGdir, 'allow_ties', 0, 'fname_in', fullfileEEG(rootEEGdir,subj,'tal/atlas/atlas_bipolar.csv'));
  
%createMasterJack(subj,rootEEGdir);

%oldphrase = '/eeg.reref/';  newphrase = '/eeg.noreref/';
%changeEventsEegFile(subj, rootEEGdir, oldphrase, newphrase);
%changeEventsEegFile(subj, rootEEGdir, '', '');



%-----------------------------------------------------------------------------------------------------%
%----------------------  EEG server directory listing ------------------------------------------------%
%--- NIHON KHODEN PIPELINE FOR GRABBING DATA:  specify the location of files on Inati's server   -----%
%-----------------------------------------------------------------------------------------------------%
nkArchive = '';

switch subj,
    %case 'NIH018', prefixMachine = 'DA8661S';	strDateStart = '10/09/13';	strDateEnd = '10/10/13';
    case 'NIH005', prefixMachine = 'MA888';     strDateStart = '02/21/12';	strDateEnd = '04/26/12';  pulseChan = {'EKG2' 'EKG1'};
    case 'NIH009', prefixMachine = 'XA866';     strDateStart = '09/16/12';	strDateEnd = '09/19/12';  pulseChan = {'EKG2' 'EKG1'};
    case 'NIH016', prefixMachine = 'XA866';     strDateStart = '06/15/13';	strDateEnd = '06/30/13';  pulseChan = {'EKG2' 'EKG1'};
    case 'NIH017', prefixMachine = 'XA866';     strDateStart = '07/15/13';	strDateEnd = '07/24/13';  pulseChan = {'EKG2' 'EKG1'};   pulseChan = {'TP5'  'TP6'};
    case 'NIH018', prefixMachine = 'DA866';     strDateStart = '09/25/13';	strDateEnd = '10/10/13';  pulseChan = {'EKG2' 'EKG1'};
    case 'NIH019', prefixMachine = 'BA693';     strDateStart = '10/10/13';	strDateEnd = '10/16/13';  pulseChan = {'EKG2' 'EKG1'};
    case 'NIH020', prefixMachine = 'DA866';     strDateStart = '10/25/13';	strDateEnd = '11/01/13';  pulseChan = {'DC09'};
    case 'NIH021', prefixMachine = 'EA877';     strDateStart = '11/06/13';	strDateEnd = '11/16/13';  pulseChan = {'DC09'};
        
    case 'NIH022', prefixMachine = 'DA866';     strDateStart = '01/20/14';	strDateEnd = '02/03/14';  pulseChan = {'DC09'};
    case 'NIH023', prefixMachine = 'EA877';     strDateStart = '05/19/14';	strDateEnd = '05/29/14';  pulseChan = {'DC09'};
    case 'NIH024', prefixMachine = 'EA877';     strDateStart = '07/07/14';	strDateEnd = '07/19/14';  nkArchive='Archive2';  pulseChan = {'DC09'};  behavSess={'0Z4','0ZJ','0ZW','101','10G','10X'};  stimSess={'11D','11H'};
    case 'NIH025', prefixMachine = 'DA866';     strDateStart = '09/04/14';	strDateEnd = '09/15/14';  pulseChan = {'DC09'};  behavSess={'GH','GI','GY','H2','H6'};  stimSess={''};
    case 'NIH026', prefixMachine = 'DA8662';    strDateStart = '09/23/14';	strDateEnd = '10/07/14';  pulseChan = {'DC09'};  behavSess={'JD','JS','JW','KD','KO','KT','KY','L9','LG','LI','LM','M3','MP'};  stimSess={''};
    case 'NIH027', prefixMachine = 'EA8771';    strDateStart = '10/21/14';	strDateEnd = '11/03/14';  pulseChan = {'DC09'};  behavSess={'7S','85','87','89','8R','9P','C3','CC','AP','AR','CE'};  stimSess={'AT'};
        
    case 'NIH028', prefixMachine = 'DA866';     strDateStart = '02/12/15';	strDateEnd = '02/28/15';  pulseChan = {'DC09'};  behavSess={'2WD','2WK','2WY','2X3','2XF','2XK','2XZ','2Y7','2YC','2YQ','30H','30N','311','32V'};  stimSess={'2YW','32I','32N','32R'};
    case 'NIH029', prefixMachine = 'CA211';     strDateStart = '06/02/15';	strDateEnd = '06/25/15';  pulseChan = {'DC09'};  behavSess={'1GI','1GO','1H1','1H6','1HJ','1HO','1HZ','1I4','1IL','1IN','1IZ','1J3','1K5','1KB','1KN','1KU','1L7','1LF','1LZ','1MC','1MS','1MX','1NN','1NS','1O5','1O7','10J','10L','1ON','1OQ','1OU','1P8'};  stimSess={'1IX','1J1'};
    case 'NIH030', prefixMachine = 'CA211';     strDateStart = '06/26/15';	strDateEnd = '07/07/15';  pulseChan = {'DC09'};  behavSess={'1RE','1SG','1UF','1UK','1US','1VN','1WF','1XE','1XL','1YK'};  stimSess={'1VB','1VD'};
    case 'NIH031', prefixMachine = 'CA211';     strDateStart = '07/11/15';	strDateEnd = '07/20/15';  pulseChan = {'DC09'};  behavSess={'20T','21E','21J','222','22A','22E','23B','23F','23H','23S','23X','24P'};  stimSess={'22G','22K','22O','233','237'};
    case 'NIH032', prefixMachine = 'CA211';     strDateStart = '07/20/15';	strDateEnd = '08/31/15';  pulseChan = {'DC09'};  behavSess={'26A','26G','26M','26X','273','276','282','298','29I'};  stimSess={'28N','294'};
    case 'NIH033', prefixMachine = 'HA866';     strDateStart = '09/01/15';  strDateEnd = '09/20/15';  pulseChan = {'DC09'};  behavSess={''}; stimSess={''};
    case 'NIH034', prefixMachine = 'HA866';     strDateStart = '09/23/15';  strDateEnd = '10/01/15';  pulseChan = {'DC09'};  behavSess={'0AE','0BE','0BI','0BX','0C3'}; stimSess={''};
    case 'NIH035', prefixMachine = 'HA866';     strDateStart = '11/12/15';  strDateEnd = '11/21/15';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH036', prefixMachine = 'CA211';     strDateStart = '12/09/15';  strDateEnd = '12/23/15';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
        
    case 'NIH037', prefixMachine = 'HA866';     strDateStart = '01/21/16';  strDateEnd = '02/17/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH038', prefixMachine = 'CA211';     strDateStart = '02/15/16';  strDateEnd = '02/27/16';  pulseChan = {'DC09'};  behavSess={'315'}; stimSess={};
    case 'NIH039', prefixMachine = 'HA866';     strDateStart = '03/03/16';  strDateEnd = '03/15/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH040', prefixMachine = 'HA866';     strDateStart = '04/11/16';  strDateEnd = '04/26/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH041', prefixMachine = 'CA211';     strDateStart = '05/17/16';  strDateEnd = '05/26/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={'3A7','3DV','3Dz','3E1'};
    case 'NIH042', prefixMachine = 'CA211';     strDateStart = '06/09/16';  strDateEnd = '06/30/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={'3G5','3IL','3IN','3L0','3L2','3L3','3LF','3LG','3LI','3MY','3NH','3NQ','3NS'};
    case 'NIH043', prefixMachine = 'HA866';     strDateStart = '07/07/16';  strDateEnd = '07/30/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'TRE003', prefixMachine = 'CA211';     strDateStart = '07/10/16';  strDateEnd = '08/12/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH044', prefixMachine = 'CA211';     strDateStart = '10/06/16';  strDateEnd = '10/21/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH045', prefixMachine = 'CA211';     strDateStart = '10/30/16';  strDateEnd = '11/11/16';  nkArchive='Archive4';  pulseChan = {'DC09'};  behavSess={'47U','489','48F'}; stimSess={'49X','49Z','4A3','4A5','4A7'};
    case 'NIH046', prefixMachine = 'HA866';     strDateStart = '12/05/16';  strDateEnd = '12/17/16';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
        
    case 'NIH047', prefixMachine = 'HA866';     strDateStart = '03/06/17';  strDateEnd = '03/19/17';  nkArchive='Archive4';  pulseChan = {'DC09'};  behavSess={'1YG','1YL','1YZ','1ZB','1ZN','1ZV','20D','20I','20U','20Y','21O','21U'}; stimSess={'212','216','21I','21S','21W'};
    case 'NIH048', prefixMachine = 'HA866';     strDateStart = '05/18/17';  strDateEnd = '06/01/17';  pulseChan = {'DC09'};  behavSess={'24Y','255','25H','25L','25N','25Z','26I','272','27F','27H','285','28k'}; stimSess={'26Q','276','27C'};
    case 'NIH049', prefixMachine = 'HA866';     strDateStart = '05/31/17';  strDateEnd = '06/19/17';  pulseChan = {'DC09'};  behavSess={'29I','29U','29Z','2A8','2AD','2AM','2AR','2B1','2B6','2CB','2CI','2DC','2DV'}; stimSess={'2DM','2DR','2DX','2DZ'};
    case 'NIH050', prefixMachine = 'HA866';     strDateStart = '06/20/17';  strDateEnd = '07/15/17';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH051', prefixMachine = 'CA211';     strDateStart = '07/04/17';  strDateEnd = '07/22/17';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
        
    case 'NIH058', prefixMachine = 'CA211';     strDateStart = 'XX/XX/18';  strDateEnd = 'XX/XX/18';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
        
    case 'NIH064', prefixMachine = 'HA866';     strDateStart = '09/20/18';  strDateEnd = '10/20/18';  pulseChan = {'DC09'};  behavSess={'3H0','3HH','3I4','3LN','3LS'}; stimSess={'3FY','3G0'};
    case 'NIH064', prefixMachine = 'HA866';     strDateStart = '09/20/18';  strDateEnd = '10/20/18';  pulseChan = {'DC09'};  behavSess={'3H0','3HH','3I4','3LN','3LS'}; stimSess={'3FY','3G0'};
    case 'NIH064', prefixMachine = 'HA866';     strDateStart = '09/20/18';  strDateEnd = '10/11/18';  pulseChan = {'DC09'};  behavSess={'3H0','3HH','3I4','3LN','3LS'}; stimSess={'3FY','3G0'};
    case 'NIH065', prefixMachine = 'CA211';     strDateStart = '10/01/18';  strDateEnd = '10/13/18';  pulseChan = {'DC09'};  behavSess={}; stimSess={};
    case 'NIH067', prefixMachine = 'CA211';     strDateStart = '10/29/18';  strDateEnd = '11/15/18';  nkArchive='Archive7';  pulseChan = {'DC09'};  behavSess={}; stimSess={};  
    otherwise,     prefixMachine = '???';       strDateStart = '???';       strDateEnd = '???';
end
%eegRawServerDir(subj, rootEEGdir, prefixMachine, strDateStart, strDateEnd, nkArchive);  %- try none, archive3, archive4

%%eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, rawFileList, typeListOptional, nktdir, stimFlag, szFlag )
%eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, behavSess, {'21E','EEG'});              %rawPath = fullfileEEG(rootEEGdir,subj,'raw/_grabbed');  eegPulseVisualize(rawPath, {'DC09','DC10'});   % STIM {'21E','EEG','LOG'}
%eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, stimSess,  {'21E','EEG','LOG'},'',1,0);
%rawPath = fullfileEEG(rootEEGdir,subj,'raw/_grabbed');  eegPulseVisualize(rawPath, {'DC09','DC10'});   % STIM {'21E','EEG','LOG'}

%eegRawServerGrabFiles(subj, rootEEGdir, prefixMachine, {'5AY'}, {'21E','EEG','LOG'}); rawPath = fullfileEEG(rootEEGdir,subj,'raw/_grabbed');  eegPulseVisualize(rawPath, {'DC09','DC10'});   % STIM {'21E','EEG','LOG'}

%%%%%%%%
%eegRawServerGrabFiles('NIH064', rootEEGdir, 'HA866', {'3LN','3LS'}, {'21E','EEG'});  % STIM {'21E','EEG','LOG'}
%eegRawServerGrabFiles('NIH027', rootEEGdir, 'EA8771', {'9P'},{'VOR'}); %- video output for missing audio session  EA87719P




%-----------------------------------------------------------------------------------------------------%
%----------------------  EEG server directory listing ------------------------------------------------%
%--- BLACKROCK PIPELINE FOR GRABBING DATA:  data is saved on Zaghloul's server                   -----%
%-----------------------------------------------------------------------------------------------------%

%%%%%% Blackrock version of file listing
%rootFRNU56dir = '/Volumes/56C/UTAH_C/';   shortTest=0;
%getBR_rawFileList(subj,rootEEGdir,rootFRNU56dir, shortTest);

%sessList = {};  JUST_CHECK = 0;
%getBR_grabFileList(subj,rootEEGdir,rootFRNU56dir, sessList, JUST_CHECK);





%-----------------------------------------------------------------------------------------------------%
%----------------------  pulse visualize -------------------------------------------------------------%

%%<<<>>> USE THIS TWO LINES FOR THE REAL DEAL <<<<>>>>>
%rawPath = fullfileEEG(rootEEGdir,subj,'raw/');              eegPulseVisualize(rawPath, {'EKG1','EKG2'});
%rawPath = fullfileEEG(rootEEGdir,subj,'raw');               eegPulseVisualize(rawPath, {'DC09','DC10'});
%
%rawPath = fullfileEEG(rootEEGdir,subj,'raw/STIM_MAP');       eegPulseVisualize(rawPath, {'DC09','DC10'});



%- Look for open files
if 0,
    fIDs = fopen('all');
    if length(fIDs)>0,
        fprintf('\n uh oh: %s left %d files open:\n', subj, length(fIDs));
        for iF=1:length(fIDs), filename = fopen(fIDs(iF)); fprintf(' %d) %s\n',iF,filename); end;
    end
end




%-----------------------------------------------------------------------------------------------------%
%----------------------  Extract Behavioral Events for the Attention Task ----------------------------%
taskListFound = dir(fullfileEEG(rootEEGdir,subj,'behavioral')); taskListCln = {}; for iT=1:length(taskListFound), if ~strcmp(taskListFound(iT).name(1),'.') & taskListFound(iT).isdir, taskListCln{end+1} = taskListFound(iT).name; end; end
%fprintf('\n Task folders found:'); disp(taskListCln);

taskList = {'stimMapGUI' 'manipMemImp' 'palRam' 'palRamStim' 'PAL' 'stimulusBlast'  'stimulusGuess' 'attentionTask' 'attentionZap' 'stimMapAnn'  'stimulusBlast' 'semanticSpan' 'zapList'  'paRemap' 'languageTask' 'moveTask' 'pa3' 'paRepeat'  'playPass'  'auditoryLexicalDecision', 'auditoryVowels', 'goAntiGo', 'SequenceMem', 'theDoors'};

%taskList = taskListCln; for iTask = [1:length(taskListCln)],
for iTask = [],
    %behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral_preOp');   % extract preOp tasks
    %behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral_postOp');  % extract postOp tasks
    
    fclose('all');
    behavioralProcessing(subj, rootEEGdir, taskList{iTask}, 'behavioral');         % extract iEEG tasks
    fIDs = fopen('all');  if length(fIDs)>0, fprintf('\n uh oh: %s left %d files open', taskList{iTask}, length(fIDs)); keyboard; end
end%





%-----------------------------------------------------------------------------------------------------%
%---------------------- Prep for and Align; Create alignmentSummary.txt  -----------------------------%
%-----------------------------------------------------------------------------------------------------%



%-STANDARD VERSION
%eegPrepAndAlign(subj, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',0, 'batch',0, 'skipProcessReref', 0);%

%eegPrepAndAlign(subj, rootEEGdir, 'batch',0);%
%eegPrepAndAlign(subj, rootEEGdir,  'freshJackandSplit',1, 'freshExtractAlign',1, 'batch',1, 'skipProcessReref', 1);%
%eegPrepAndAlign(subj, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',1, 'batch',0, 'skipProcessReref', 1);%
%eegPrepAndAlign(subj, rootEEGdir,  'freshJackandSplit',1, 'freshExtractAlign',1, 'batch',1, 'skipProcessReref', 0);%
%eegPrepAndAlign(subj, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',0, 'batch',1, 'skipProcessReref', 1);%


%subjectDirCheck(subj,rootEEGdir);


%scanRoot  = '/Volumes/JW24TB/data24TB/eeg_new';
%checklistScan(scanRoot,'writeScanRes',1,'verifyChecklists',1,'debug',0)
% rootPath - path where NIHXXX subject folders are present
% writeScanRes - binary flag for writing out master warning list
% verifyChecklists - binary flag that will trigger a re-run of subjectDirCheck on each subject
% debug - binary flag, all written files are output to debug directory instead of rootPath
%





%-----------------------------------------------------------------------------------------------------%
%----------------------  REMAP ALL EVENTS EEGFILE FIELDS, ALL SUBJECTS
% e0ecute the code within the if statement to remap all local copies from server to local
if 1,
    %- problems, needs loving:  [5 7 9 21 35]   haven't checked 53+?
    thisSubj=subj;
    
    % JW tried to do a fresh build on 64:76.   crash on 64, 65
    %rebuildList  = [66 67 68 69   71 72 73 74 75];  rebuildFail = [64 65 68];  rerefProcessFail = [66 55 61];
    %rebuildList2 = [54 55 56 57 58 59 60 61 62 63]; rebuildPotentialissue = [56];  rebuildError = [62]; %note NIH056 used NEV for alignment... best to check this
    
    rebuildList = [64 65 68];
    %- rebuildList2 covers remaining patients from first cervello (54) to 64
    for subjNum = rebuildList %  %
        
        thisSubj = sprintf('NIH%03d', subjNum);
        %thisSubj = sprintf('TRE%03d', subjNum);
        %thisSubj = sprintf('BEH%03d', subjNum);
        
        fprintf('\n\n PROCESSING %s:',thisSubj);
        
        
        
        %taskList = {'stimMapGUI' 'palRam' 'palRamStim' 'PAL'};
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{iTask}, 'behavioral');  % attentionTask languageTask moveTask pa3 paRepeat playPass
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{1}, 'behavioral');  % attentionTask languageTask moveTask pa3 paRepeat playPass
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{2}, 'behavioral');  % attentionTask languageTask moveTask pa3 paRepeat playPass
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{3}, 'behavioral');  % attentionTask languageTask moveTask pa3 paRepeat playPass
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{iTask}, 'behavioral_preOp');  % attentionTask languageTask moveTask pa3 paRepeat playPass
        %behavioralProcessing(thisSubj, rootEEGdir, taskList{iTask}, 'behavioral_postOp');  % attentionTask languageTask moveTask pa3 paRepeat playPass

        
        
        %createMasterJack(thisSubj,rootEEGdir);
        eegPrepAndAlign(thisSubj, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',0, 'batch',1, 'skipProcessReref', 0);%  %- just process
        %eegPrepAndAlign(thisSubj, rootEEGdir,  'freshJackandSplit',0, 'freshExtractAlign',1, 'batch',1, 'skipProcessReref', 1);%  %- just process
        %eegPrepAndAlign(thisSubj, rootEEGdir,  'freshJackandSplit',1, 'freshExtractAlign',1, 'batch',0, 'skipProcessReref', 1);%  %- just process
        
        
        %atlas_simplify(thisSubj, rootEEGdir, 'allow_ties', 0, 'force_save', 1);
        
        
        %- check to see if atlas simple was made, if not, make it now
%         atlasDir = fullfileEEG(rootEEGdir,thisSubj,'tal/atlas'); %- if tal/atlas directory doesn't exist, don't bother checking for simples!
%         if exist(atlasDir,'dir'),
%             atlasSimpleMP = fullfileEEG(atlasDir,'atlas_monopolar_simple.csv');
%             if ~exist(atlasSimpleMP,'file'),
%                 fprintf('\n\n Missing atlas_monopolar_simple.csv ... making it now:\n');
%                 atlas_simplify(thisSubj, rootEEGdir, 'allow_ties', 0);
%             end
%             atlasSimpleBP = fullfileEEG(atlasDir,'atlas_bipolar_simple.csv');
%             if ~exist(atlasSimpleBP,'file'),
%                 fprintf('\n\n Missing atlas_bipolar_simple.csv ... making it now:\n');
%                 atlas_simplify(thisSubj,rootEEGdir,'allow_ties', 0, 'fname_in', fullfileEEG(atlasDir,'atlas_bipolar.csv'));
%             end
%         end
        
        
        %checkTal(subj, root, [ isBatch=0, rerun=0 ]) outputs metrics and summary figure files... usually run when tal folder created
        %checkTal(this%Subj, rootEEGdir, 1 , 1); close all;
        
        %verifyElementInfo(thisSubj, rootEEGdir);
        %rawPath = fullfileEEG(rootEEGdir,thisSubj,'raw/');      eegPulseVisualize(rawPath, {'DC09','DC10'});
        


        drawnow();
    end
end




