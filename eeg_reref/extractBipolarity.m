function bipolarPairChans = extractBipolarity(subj, rootEEGdir, write_leads_bp)
    % function [bipolarPairs,bipolarPairNums] = extractBipolarity(subjDir)
    %
    % Finds the appropriate pairs of bipolar electrodes 
    % 
    % INPUT:
    %   subj        - subject e.g. NIH003
    %   rootEEGdir  - parent of subject directory (For example, '/Users/dongj3/Jian/data/eeg/)
    %  
    % Optional Input:
    %   rewrite_leads_bp (default 0) - If 1 write/overwrite docs/leads_bp.txt
    %
    % OUTPUT:
    %   bipolarPairChans     = Array of n by 2 elements denoting the electrode numbers, where n is the number of pairs
    %   bipolarPairNums      = Corresponding channel numbers  %%-- JW removed this output 8/2017... no longer supporting numerical index of channels
    %
    % NOT OUTPUT?: bipolarDistances = Array of n elements denoting the distance, where n is the number of pairs
    %
    % JHW cleaned up a little 11/2013
    % MST get hardware type in order to exclude microwires 07/17
    % MST rewrote function based on findNeighbors
    % JHW switched leads_bp from tal to docs folder
    % SJ added 'subduralHD' to hwType_bp list (2/5/2020)
    % SJ removed 'subduralHD' from hwType_bp list (2/11/2020)

    if nargin < 3 || isempty(write_leads_bp), write_leads_bp = 0; end
    
    hwType_bp   = {'subdural','depth','micro-subdural'};
    subjDir     = fullfile(rootEEGdir, subj);  %- jw switched leads_bp.txt to docs folder (tal stuff is done independently)
    fname_leads = fullfile(subjDir, 'docs/leads_bp.txt');

    info = getElementInfo(subj, rootEEGdir);
    tagnames_bp = getTagNames(subj, rootEEGdir, 'hardwareType',hwType_bp, 'elementInfo',info);

    bipolarPairChans = [];
    for i = 1 : length(tagnames_bp)
        tag = tagnames_bp{i};
        pairs_tag = findNeighbors(subj, rootEEGdir, tag, 1, 'elementInfo', info);
        pairs_chan= util_combine_strnum(tag, pairs_tag);
        bipolarPairChans = [bipolarPairChans; pairs_chan];
    end
    
    
    if write_leads_bp
        %if ~exist(tal, 'dir'), mkdir(tal); end
        % write
        fd = fopen(fname_leads, 'w');
        for i = 1 : length(bipolarPairChans)
            fprintf(fd, '%s-%s\n', bipolarPairChans{i,1}, bipolarPairChans{i,2});
        end
        fclose(fd);
    end
    
end
