function myout = metaproc_elmt2jack(varargin)
% Desc: 
%   Given an element_info.csv file, this function expands the rows so that
%   each row has 1 electrode, outputing jacksheetMaster.csv.
%
%   Note that this function does NOT cross-reference the channels with
%   those found in EEG raws. That is done in
%   readElementInfo->eegPrepAndAlign->nk_split->createMasterJack. Thus the
%   file output by this function is not the final form of
%   jacksheetMaster.csv
%
% Inputs (pairs):
%   element_csv - path to file input [Default: ./docs/element_info.csv]
%   jack_csv    - path to file for output [Default: ./jacksheetMaster.csv]

% Outputs:
%   creates jacksheetMaster.csv

% REVISION HISTORY
%   ??/16 Cocjin - Created
%   07/16    MST - Modified to work with prepAndAlign
%

%%
if ~isempty(varargin) && istable(varargin{1})
    myout = elmt2jack_table(varargin{1});
    
    % sort by element_info order
    tagNames = varargin{1}.tagName;
    [~,presortOrder] = sort(cellfun(@(x) sscanf(x,'%*[^0-9]%d'),myout.chanName));
    chanNames = myout.chanName(presortOrder);
    sortOrder = sortBySubstring(chanNames, tagNames);
    myout = myout(sortOrder,:);
    %myout = myout(ismember(upper(myout.chanType), 'PHYS'),:); % get PHYS chans
    
else
    myout = elmt2jack_dirstuff(varargin{:});
end

return

function [fullstruct] = elmt2jack_dirstuff(varargin)
%% extract inputs, also considering defaults
default_cell = {'element_csv','./element_info.csv','parameter';
                'jack_csv','./jack_info.csv','parameter';
                'run',false,'parameter'};
inputs = util_extract_inputs(default_cell, varargin);
%% read in the input tables/csv
if inputs.run
    verifyElementInfo('','','pathToElmtInfo',inputs.element_csv);
    element_table = readtable(inputs.element_csv);
    jack_table = metaproc_elmt2jack(element_table);
    writetable(jack_table, inputs.jack_csv);
end
%%
outputs.jack_csv = inputs.jack_csv;
outputs.log = '';
%%
fullstruct.inputs = inputs;
fullstruct.outputs = outputs;
%%
return
function jackInfo = elmt2jack_table(elmtTable)
%% format class information
classTypeList = {};
for myField = elmtTable.Properties.VariableNames
    if strcmp(myField{1}(1:2),'is')        
        % convert the is[] fields to numerics
        elmtTable.(myField{1}) = cellfun(@(x) (str2num(x)), elmtTable.(myField{1}), 'UniformOutput',false);
        
        % extract out class information and convert binary
        rawClassTable = elmtTable(:,{'tagName',myField{1}});
        rawClassTable.Properties.VariableNames = {'tagName', 'absElmtIdx'};
        classInfo.(myField{1}) = pd_extract_class_table(rawClassTable, myField{1});
        classTypeList = [classTypeList; myField{1}];
    end
end

%% generate augmented element info table
% find which is nk recorded
isNKrecord = elmtTable.nkRecOrder > 0;
isSyncRef = strcmpi(elmtTable.chanType, 'SYNC') | strcmpi(elmtTable.chanType, 'CLIN');

% split into nk and notnk recorded tables
nkTable = elmtTable(isNKrecord & ~isSyncRef,:);
refTable = elmtTable(~isNKrecord & isSyncRef, :);notnkTable = elmtTable(~isNKrecord & ~isSyncRef,:);


% sort and expand
nkTable = nkTable(sort(nkTable.nkRecOrder),:);
nkInfo = pd_expandNKtable(nkTable);

%% create jackInfo from nkTableCut and class information
refTable{:,'isCut'} = repmat({[]},height(refTable),1);
refInfo = pd_expandNKtable(refTable);

notnkTable{:,'isCut'} = repmat({[]},height(notnkTable),1);
notnkInfo = pd_expandNKtable(notnkTable);

jackInfo = [nkInfo; refInfo; notnkInfo];

% if class table populated at all, then try to merge, if not, then set to 0
for iClass = 1:length(classTypeList)% merge
    myClassType = classTypeList{iClass};
    myClassTable = classInfo.(myClassType);
    myIsClass = myClassTable.Properties.VariableNames{3};    
    if ~isempty(myClassTable) 
        jackInfo = outerjoin(jackInfo, myClassTable, 'MergeKeys',1, 'Keys',{'tagName','absElmtIdx'});
        isCol = jackInfo{:,myIsClass} == 1;
        jackInfo{~isCol,myIsClass} = 0;
    else
        jackInfo{:,myIsClass} = 0;
    end
end
%%
tagName = jackInfo.tagName;
absElmtIdx = arrayfun(@num2str, jackInfo.absElmtIdx, 'UniformOutput',false);

% DC channels should always be 2 digits
twoDigitIdxs = find(strcmp('DC', jackInfo.tagName));
if ~isempty(twoDigitIdxs)
    for i = twoDigitIdxs
        absElmtIdx{i} = sprintf('%02s',absElmtIdx{i});
    end
end

jackInfo{:,'chanName'} = strcat(tagName, absElmtIdx);
jackInfo.tagName = [];
jackInfo.absElmtIdx = [];

%%
return

function expTable = pd_expandNKtable(nkTable)
% initialize variables 
subjId = {};
tagName = {};
absElmtIdx = [];
recordIdx = [];
whichHemi = {};
locType = {};
hardwareType = {};
chanType = {};

% counter for recordIdx
iRecord = 1;

% expand out to electrodes
for iRow = 1:height(nkTable)  % expand each element
    % how many elements?
    if nkTable{iRow,'bigDim'} > 0
        nObj = nkTable{iRow,'bigDim'}*nkTable{iRow,'smallDim'};
    else 
        nObj = 1;
    end
    
    % build each row
    for iObj = 1:nObj % add rows...
        % build most columns
        subjId = [subjId; nkTable{iRow, 'subjId'}];
        tagName = [tagName; nkTable{iRow,'tagName'}];
        absElmtIdx = [absElmtIdx; iObj];
        whichHemi = [whichHemi; nkTable{iRow,'whichHemi'}];
        locType = [locType; nkTable{iRow,'locType'}];
        hardwareType = [hardwareType; nkTable{iRow, 'hardwareType'}];
        chanType = [chanType; nkTable{iRow,'chanType'}];
        
        % build recordIdx column
        if nkTable{iRow,'nkRecOrder'} > 0 
            if ~ismember(iObj,nkTable{iRow,'isCut'}{:})
                recordIdx = [recordIdx; iRecord];
                iRecord = iRecord + 1;
            else
                recordIdx = [recordIdx; -2];
            end
        else
            recordIdx = [recordIdx; nkTable{iRow,'nkRecOrder'}];
        end
    end
end
expTable = table(subjId, tagName, absElmtIdx, recordIdx, whichHemi, locType, hardwareType, chanType);
return

function outputTable = pd_extract_class_table(jackInfo, classType)
% expand out to table
tagName = {};
absElmtIdx = [];

iItr = 1;
for iRow = 1:height(jackInfo)
    myTagName = jackInfo.tagName{iRow};
    idxList = jackInfo.absElmtIdx{iRow};
    for iIdx = 1:length(idxList)
        tagName = [tagName; myTagName];
        absElmtIdx = [absElmtIdx; idxList(iIdx)];
        iItr = iItr + 1;
    end
end

% add classType column
isColumn = ones(size(tagName));

% make table out of this shit
outputTable = table(tagName, absElmtIdx, isColumn);
outputTable.Properties.VariableNames = {'tagName', 'absElmtIdx', classType};
return
