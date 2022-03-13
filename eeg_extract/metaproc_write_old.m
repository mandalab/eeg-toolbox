function [fullstruct] = metaproc_write_old(varargin)
% METAPROC_WRITE_OLD makes electrodes.m, tagNames, good_leads, bad_leads, leads txt's


%%
default_cell = {...
    'elements_csv','./element_info.csv','parameter';
    'jack_csv','','parameter';
    'tal_dir','','parameter';
    'docs_dir','./docs/','parameter';
    'run', true','parameter';
    'makeOldSource', false, 'parameter'};

inputs = util_extract_inputs(default_cell, varargin);
%%
outputs.elec_m = fullfile(inputs.docs_dir, 'electrodes.m');
if inputs.makeOldSource
    outputs.tag_txt = fullfile(inputs.docs_dir, 'tagNames.txt');
    outputs.good_leads_txt = fullfile(inputs.tal_dir, 'good_leads.txt');
    outputs.bad_leads_txt = fullfile(inputs.tal_dir,'bad_leads.txt');
    outputs.leads_txt = fullfile(inputs.tal_dir, 'leads.txt');
end
%%
if inputs.run
    elements = metaproc_read_elements(inputs.elements_csv);
    nk_elements = metaproc_extract_nk_elements(elements);
    [r, gridLayout, missingElecs] = metaproc_extract_elec_em(nk_elements);
    metaproc_write_elec_m(outputs.elec_m,r,gridLayout,missingElecs, nk_elements.tagName, nk_elements.fullTagName, height(nk_elements));
    if inputs.makeOldSource
        jack_info = readtable(inputs.jack_csv);
        metaproc_write_leads(outputs,jack_info);
        metaproc_write_tags(outputs.tag_txt, nk_elements)
    end
end
%%
fullstruct.inputs = inputs;
fullstruct.outputs = outputs;
%%
return

function nkTable = metaproc_extract_nk_elements(elmtTable)
%% generate augmented element info table
% find which is nk recorded
isNKrecord = elmtTable.nkRecOrder > 0;
isRef = elmtTable.nkRecOrder == 0;

% split into nk and notnk recorded tables
nkTable = elmtTable(isNKrecord,:);
notnkTable = elmtTable(~isNKrecord & ~ isRef,:);
refTable = elmtTable(isRef, :);

% sort and expand
nkTable = nkTable(sort(nkTable.nkRecOrder),:);
return

function metaproc_write_tags(tagPath, nkTable)
% tagNames 
tagNames = [nkTable{:,'tagName'}; 'Z';'DC'];
tagTable = table(tagNames);

% ensure that directory exists
[tagDir fname ext] = fileparts(tagPath);
if ~isdir(tagDir)
    mkdir(tagDir);
end

writetable(tagTable, tagPath, 'Delimiter', ' ','WriteVariableNames',false);
return
function [r, gridLayout, missingElecs] = metaproc_extract_elec_em(nkTableMissing)
missingElecs = nkTableMissing{:,'isCut'};

% generate gridLayout
gridLayout = nkTableMissing{:,{'bigDim','smallDim'}};

% generate r
numElecsRaw = gridLayout(:,1).*gridLayout(:,2);
numMissing = arrayfun(@(x) numel(x{:}), missingElecs);
numElecs = numElecsRaw-numMissing;
for iRow = 1:length(numElecs)
    if iRow == 1
        r(iRow,1) = 1;
    else
        r(iRow,1) = r(iRow-1,2) + 1;
    end
    r(iRow,2) = numElecs(iRow)-1+r(iRow,1);
end

return
function metaproc_write_elec_m(elecPath, r, gridLayout, missingElecs, tagName, fullTagName, numElecTags)
[elecDir fname ext] = fileparts(elecPath);
if ~isdir(elecDir)
    mkdir(elecDir);
end

% r
fid = fopen(elecPath,'w');
fprintf(fid,'%s\n','r = [');
for iTag = 1:numElecTags
    if iTag == numElecTags
        fprintf(fid,'%d, %d \t',r(iTag,:));
    else
        fprintf(fid,'%d, %d;\t',r(iTag,:));
    end
    fprintf(fid,'%s\t',['%' tagName{iTag}]);
    fprintf(fid,'%s',['(' fullTagName{iTag} ')']);
    fprintf(fid,'\n');
end
fprintf(fid,'%s\n','];');

% gridLayout
fprintf(fid,'\n');
fprintf(fid,'%s\n','gridLayout = [');
for iTag = 1:numElecTags
    if iTag == numElecTags
        fprintf(fid,'%d, %d ',gridLayout(iTag,:));
    else
        fprintf(fid,'%d, %d; ',gridLayout(iTag,:));
    end
    fprintf(fid,'\n');    
end
fprintf(fid,'%s\n','];');

% missingElecs
fprintf(fid,'\n');
fprintf(fid,'%s\n','missingElecs = {');
for iTag = 1:numElecTags
    if iTag == numElecTags
        fprintf(fid,'[%s] \n',num2str(missingElecs{iTag}));
    else
        fprintf(fid,'[%s];\n',num2str(missingElecs{iTag}));
    end
end
fprintf(fid,'%s\n','};');
fclose(fid);
return
function metaproc_write_leads(oldPathSubj, leadsInfo)
% leadsdir
[leadsPath fname ext] = fileparts(oldPathSubj.leads_txt);
if ~isdir(leadsPath)
    mkdir(leadsPath);
end


% sort by recordIdx
[~, sortRec]= sort(leadsInfo.recordIdx);
leadsInfo = leadsInfo(sortRec,:);

% filter into leads, good_leads, and bad_leads
leadsFilt = leadsInfo(leadsInfo.recordIdx > 0, :);
isIctal = leadsFilt.isIctal == 1;

leadsErr = leadsFilt{:,{'recordIdx'}};
goodErr = leadsFilt{~isIctal,{'recordIdx'}};
badErr = leadsFilt{isIctal,{'recordIdx'}};

% write out to file
dlmwrite(oldPathSubj.leads_txt,leadsErr)
dlmwrite(oldPathSubj.good_leads_txt,goodErr)
dlmwrite(oldPathSubj.bad_leads_txt,badErr)
return