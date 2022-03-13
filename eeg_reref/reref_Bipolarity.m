function varargout = reref_Bipolarity(subjDir, eegFilestem, varargin)
% function r = reref_Bipolarity(subjDir,eegFilestem)
%
% Creates bipolar rereerenced files for the given filestem in processedBP
%
% ** Note ** Bipolar files are overwritten if they already exist
%
% INPUT:
%         --subjDir      = subject directory. For example, '/Users/dongj3/Jian/data/eeg/NIH003/'
%         --eegFilestem  = 'NIH001_121211_0933'
%
%
% REVISION HISTORY
%  07/16 MST - Added nested functions to clean up code
%              Adapted for new csv source
%  05/17 MST - Adapted to pull from eeg.processed/ instead of eeg.noreref
%  04/18 JHW - changed reref folder name and location of leads_bp (from tal to docs)
%  05/18 MST - Changed responsibilities: no longer creates/checks leads_bp

SOURCE_DIR = 'eeg.processed'; % eeg.noreref; *MST 05/17
OUTPUT_DIR = 'eeg.processedBP'; % JW change
%SOURCE_CHECK = false; % parse leads_bp again to double-check getLeads function

[rootEEGdir,subj] = fileparts(subjDir);
docsDir     = fullfile(subjDir,'docs');  %- JW switched from tal to docs 4/2018
bipolarFile = fullfile(docsDir,'leads_bp.txt');
assert(exist(bipolarFile, 'file') > 0, 'File not found: %s', bipolarFile)


%% JW cut 12/2018. No good reason to do this here.
%% It will get done within nkElectrodesFilt, which additionally confirms that only recorded channels are returned
% leads = getLeads(subj, rootEEGdir, ...
%     'leadsType',    'all', ...
%     'chanType', 'PHYS', ...
%     'hardwareType', 'subdural');
%
% bipolarPairs = getLeads(subj, rootEEGdir, 'isBP', true);
%
% subdural_match = contains(bipolarPairs,leads);
% keep = find(sum(subdural_match,2)==2);
% bipolarPairs = bipolarPairs(keep,:);



% get leadsBP from session list of bp elecs (actually recorded elecs)
[~,bipolarPairs] = nkElectrodeFilt(rootEEGdir, subj, eegFilestem); %- JW 12/2018, now this function always returns channels/bp_pairs
%if ~isempty(bipolarPairs2)
%    bipolarPairs = bipolarPairs2;
%end

makeBipolarFiles();
% end reref_Bipolarity code

%%%%%%%%%%%%%%%%%%%%%%% nested function definitions %%%%%%%%%%%%%%%%%%%%%%%


%Actually creates the bipolar files...
    function makeBipolarFiles()
        
        if ~exist(fullfile(subjDir, OUTPUT_DIR))
            mkdir(fullfile(subjDir, OUTPUT_DIR))
        end
        
        for k = 1:size(bipolarPairs,1)
            pair = bipolarPairs(k,:);
            
            if ~exist(fullfile(subjDir, OUTPUT_DIR, eegFilestem))
                mkdir(fullfile(subjDir, OUTPUT_DIR, eegFilestem))
            end
            
            chanfile = sprintf('%s/%s-%s', fullfile(subjDir, OUTPUT_DIR, eegFilestem),bipolarPairs{k,1},bipolarPairs{k,2});
            if isnumeric(pair(1)) && isnumeric(pair(2))
                file1 = sprintf('%s/%03i', eegFilestem, pair(1));
                file2 = sprintf('%s/%03i', eegFilestem, pair(2));
            else
                file1 = sprintf('%s/%s', eegFilestem, pair{1});
                file2 = sprintf('%s/%s', eegFilestem, pair{2});
            end
            
            file1Dir = fullfile(subjDir,SOURCE_DIR,file1);
            file2Dir = fullfile(subjDir,SOURCE_DIR,file2);
            
            dataFormat = 'int16';
            file1Handle = fopen(file1Dir, 'r','l');
            if file1Handle == -1
                fprintf('  ERROR: reref_Bipolarity can''t open %s\n',file1Dir);
                continue; % skip pair
            end
            data1 =  fread(file1Handle, inf, dataFormat);
            fclose(file1Handle);
            
            file2Handle = fopen(file2Dir, 'r','l');
            if file2Handle == -1
                fprintf('  ERROR: reref_Bipolarity cant open %s\n', file2Dir);
                continue; % skip pair
            end
            data2 =  fread(file2Handle, inf, dataFormat);
            fclose(file2Handle);
            
            if length(data1) ~= length(data2)
                fprintf('SEVERE ERROR: %s has size %d while %s has size %d! Cannot create %s.\n',...
                    file1, length(data1), file2, length(data2), chanfile);
                fprintf('  One possibility is that one split channel has a number/chanel name mismatch. Check your jacksheetMaster\n');
                keyboard;
                return;
            end
            
            
            disp(['Making file:' chanfile]);
            chanHandle = fopen(chanfile,'w','l');
            assert(chanHandle > 0, 'Cannot open file %s for writing (fopen code: %d)', chanfile, chanHandle)
            fwrite(chanHandle,data1-data2,dataFormat);
            fclose(chanHandle);
        end
    end % makeBipolarFiles

end % function reref_Bipolarity
