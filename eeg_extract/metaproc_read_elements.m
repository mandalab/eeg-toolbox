
%%
function elmtTable = metaproc_read_elements(elPath)
elmtTable = readtable(elPath);
classTypeList = {};
for myField = elmtTable.Properties.VariableNames
    if strcmp(myField{1}(1:2),'is')        
        % convert the is[] fields to numerics
        elmtTable.(myField{1}) = cellfun(@(x) (str2num(x)), elmtTable.(myField{1}), 'UniformOutput',false);
        
        % extract out class information
        rawClassTable = elmtTable(:,{'tagName',myField{1}});
        rawClassTable.Properties.VariableNames = {'tagName', 'absElmtIdx'};
        classInfo.(myField{1}) = pd_extract_class_table(rawClassTable, myField{1});
        classTypeList = [classTypeList; myField{1}];
    end
end
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