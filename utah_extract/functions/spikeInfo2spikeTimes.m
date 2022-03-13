function [chans, unitIDs, spikeTimes] = spikeInfo2SpikeTimes(fname)
    % spikeInfoToSpikeTimes examines spikeInfo.mat, processes it, and extracts each good unit's spikeTimes
    % 
    % The purpose of this function is to adapt to the structure of SpikeInfo.mat, if its structure changes, and to perhaps
    % pre-process or filter it.
    %
    % For example, spikeInfo as output by the mountain sort algorithm may "find" 20 units on a channel, but we may believe this
    % is erroneous. Any examination of metrics and filtering out of these units would happen in this function.
    %
    % It was written to be fed to spike-kernel extraction for real-time processing. The output is conventient for loading a 
    % and finding a waveform from the raw ns6 file for each unit in the n units output from this function.
    %
    %
    % INPUT
    %   fname   - filepath to spikeInfo
    %
    % OUTPUT
    %   chans   - Length n array containing numeric indices into utah's channels (e.g. between 1-96). Not necessarily unique.
    %   unitIDs - Length n cell string array containing times of i'th unit's ID. Unique.
    %   spikeTimes - Length n cell array containing of timestamps of i'th unit's array of spike times (in samples)
    %
    % NOTES
    %   - spikeTimes are trough-aligned (output of mountain sort), ie if spikeTime=100, a spike trough occurs at the 100'th sample
    %       but the waveform (50 samples) might begin earlier
    %
    %   - It is also possible this function should be written to examine some of the other fiels, like the summary CSVs
    %
    % REVISION HISTORY
    %   08/18 MST - Created
    %
    % See also: extractSpikeInfo, trainSpikeKernel
    
    
    DBG = 0;
    if nargin < 1
        fname = '/Volumes/56D/UTAH_RAFI/spikesForJohn/NIH057/NIH057_180226_1516_utah/NIH057_180226_1516_utah_spikeInfo.mat';
        DBG = 1;    
    end
    
    %% Load
    if exist('spikeInfo', 'var')
        DBG = 1;
        disp(spikeInfo);
    else
        spikeInfo = load(fname);
        %spikeInfo = data.spikeInfo;
        %clear data;
    end
    
    %% New format (08/31/2018)
    unitIDs = spikeInfo.uniqueUnitID;
    spikeTimes = spikeInfo.timeStamp;
    chans = {};
    for i = 1:numel(unitIDs)
        parts = strsplit(unitIDs{i}, '_');
        chans(i) = parts(1);
    end
    
    
    return

    
end


