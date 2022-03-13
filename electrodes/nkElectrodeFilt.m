function [chan_names,leads_bp] = nkElectrodeFilt(rootEEGdir,subj,filestem, varargin)
%function [chan_nums,chan_names,leads_bp] = nkElectrodeFilt(rootEEGdir,subj,filestem, varargin)
% Function [filestem,elecLabels,lBp,glBp,blBp] = nkElectrodeFilt(eegDir,subj,filestem)
%
%   Description: Given a filestem, this function looks at the relevant
%   jacksheet and identifies those electrodes that were actually recorded
%   from for it. The function then returns those electrodes that were
%   actually recorded from as well as the different bipolar pairs for that
%   filestem.
%
%   **NOTE 1** if all channels are recorded in all sessions, this
%       function returns empty - so this condition should be checked!
%
%   **Note 2** chan_names are based on the jacksheet's second column and
%       chan_nums are based on that name's chan2num translation based on
%       jacksheetMaster (NOT the first column!)
%
%   INPUT:
%   	egDir=directory where data is held
%                  (ex. '/Users/damerasr/Sri/data/eeg/')
%         --subj=subject name (ex. 'NIH001')
%         --filestem = '121211_0933'
%
%   Key-Value Pair Input
%         --useOldSource, [true/false] = whether to use old leads.txt files
%               (old) or new csv files (true, default)
%   OUTPUT:
%       **NOTE**: empty output means all chans recorded in this session!

%       chan_names: list of recording electrodes (or empty)
%       leads_bp:  names of leads_bp filtered by recording electrodes (or empty)
%
%  REVISION HISTORY:
%   07/16 MST - Works without needing good/bad/leads txt's
%   07/17 MST - remove useOldSource option
%   08/17 JHW - remove "chan_nums" output... not using numerical indexing of channels anymore
%   04/18 JHW - changed location of leads_bp from tal to docs
%   12/18 JHW - return full chan list and BP chan list in case where no channels are missing

subjDir     = fullfile(rootEEGdir,subj); % where subject data is located
docsDir     = fullfile(subjDir,'docs'); % leads bp now stored in docs folder...
noRerefDir  = fullfile(subjDir,'eeg.noreref'); % where no-reref data is located
phys        = getLeads(subj, rootEEGdir, 'chanType', 'PHYS');


allChanExist = 0;
% check whether any channel is missing in any raw
filename = fullfile(docsDir, 'raw_info.csv');
if exist(filename, 'file')
    allChanExist = 1;
    rawInfo = readtableSafe(filename);
    %raws = setdiff(rawInfo.Properties.VariableNames, {'chanName','in_element_info','raw_info'});
    raws = setdiff(rawInfo.Properties.VariableNames, {'chanName','in_jackSheetMaster','raw_info'}); %- JW fixed 12/2018. column name was wrong
    for i = 1 : length(raws)
        if sum(rawInfo.(raws{i})) ~= height(rawInfo) % ie all 1's in column
            allChanExist = 0;
            break;
        end
    end
end



chan_names  = []; %#ok<*NASGU>
leads_bp    = [];

if allChanExist
    fprintf('ALL CHANNELES ARE RECORDED IN ALL SESSIONS (nkElecrodeFilt full chan list)\n')
    chan_names = rawInfo.chanName;  %- 12/2018, JW uncommented the following two lines
    leads_bp   = getLeads(subj, rootEEGdir, 'isBP', true);
else
    
    % reads in jacksheet
    jackFile = fullfile(noRerefDir, filestem, 'jacksheet.txt');
    fd = fopen(jackFile);
    jackSheet = textscan(fd,'%d%s');
    fclose(fd);
    chan_names = jackSheet{2};
    
    % only use PHYS electrodes for re-referencing
    ref_mask = ismember(upper(chan_names), upper(phys));
    chan_names = chan_names(ref_mask);
    
    % reads in all bipolar pairs
    bpFile = fullfile(docsDir,'leads_bp.txt');
    if exist(bpFile,'file')
        fd = fopen(bpFile);
        lines = textscan(fd, '%s');
        lines = lines{1};
        fclose(fd);
        leads_bp = cellfun(@strsplit, lines, repmat({'-'}, length(lines), 1), 'uniformoutput',0);
        
        col = @(x,i) cellfun(@(c) c{i}, x, 'uniformOutput',0);
        leads_bp = [col(leads_bp, 1), col(leads_bp, 2)];
        
        
        % outputs relevant leads bipolar
        mask = ismember(leads_bp(:,1),chan_names) & ismember(leads_bp(:,2),chan_names);
        leads_bp = leads_bp(mask,:);
        
    else
        error('%s does not exist!', bpFile);
    end
    
end
