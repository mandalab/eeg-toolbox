function pulses = makePulsesStruct(varargin)
%
% Function created by Chris Zawara to nicely package all the pulse info for subsequent alignment steps
%
% created 2018 by CZ
%         2019/6  CZ added digital time series
% tweaked and renamed by JW
%

%output struct
pulses = struct;

VERBOSE = 0; %- optional flag for outputs that help with debugging


p = inputParser;

p.addParameter('ns3_fpath', '', @ischar);
p.addParameter('nev_fpath', '', @ischar);
p.addParameter('postProc', []);

parse(p, varargin{:});

%record filenames passed in
ns3_fpath = p.Results.ns3_fpath;
nev_fpath = p.Results.nev_fpath;
postProc  = p.Results.postProc;

ns3_fpath_splits = strsplit(ns3_fpath, '/');
nev_fpath_splits = strsplit(nev_fpath, '/');

pulses.ain_filename = ns3_fpath_splits{end};
pulses.din_filename = nev_fpath_splits{end};

%first read ns3

if ~isequal(ns3_fpath, '') && exist(ns3_fpath, 'file')
    
    micro_freq = 30000; % sampling rate for micros
    output_timeseries_freq = 1000;
    
    pulses.timeseries_freq = output_timeseries_freq;
    pulses.uptimes_freq = micro_freq;
    
    rawData_NS3 = concatOpenNSx(ns3_fpath);
    pulses.pulse_nsx_postProc = rawData_NS3.postProc;
    
    ns3_freq = rawData_NS3.MetaTags.SamplingFreq;
    
    %create vector of indices for downsampling analog pulse timeseries
    downsampled_step_factor = ns3_freq/output_timeseries_freq;
    
    
    
    electrodeLabels = { rawData_NS3.ElectrodesInfo.Label };
    
    %remove null spaces from electrodeLabels
    for iElec = 1:length(electrodeLabels)
        
        currentElec = electrodeLabels{iElec};
        endIdx = 1;
        
        while (endIdx + 1) <= length(currentElec) && double(currentElec(endIdx + 1)) ~= 0
            
            endIdx = endIdx + 1;
        end
        
        
        electrodeLabels{iElec} = electrodeLabels{iElec}(1:endIdx);
    end
    
    
    if VERBOSE,
        fprintf('electrode labels\n');
        disp(electrodeLabels);
    end
    
    
    % look for ain channels
    ains = electrodeLabels( find(cellfun( @(x) contains(x, 'ain'), electrodeLabels)) );
    
    if VERBOSE,
        fprintf('ain electrodes\n');
        disp(ains);
    end
    
    %loop looking for ain
    for iAin = 1:length(ains)
        
        current_ain = ains{iAin};
        channelIdx = find(cellfun( @(x) isequal(x, current_ain), electrodeLabels));
        
        if length(channelIdx) > 1
            error('more than one channel in %s with name %s. Is this possible?\n', ns3_fpath, current_ain);
        end
        
        if VERBOSE,
            fprintf('channelIdx %d\n', channelIdx);
            size(rawData_NS3.Data(channelIdx,:)')
        end
        
        % call get_triggers for current channel, set pulse times to physio sample rate and store in pulses struct
        thisAinPulses = get_triggers(double(rawData_NS3.Data(channelIdx,:)'),ns3_freq);
        eval( sprintf('pulses.%s_uptimes = thisAinPulses{1}*micro_freq/ns3_freq;', current_ain ) );
        
        % downsample current channel and store in the pulses struct
        downsampled_channel = rawData_NS3.Data(channelIdx,1:downsampled_step_factor:size(rawData_NS3.Data(channelIdx,:),2));
        eval( sprintf('pulses.%s_ts = downsampled_channel;', current_ain ) );
        
        
    end
    
    
end


pulses.din1_uptimes = [];
pulses.din2_uptimes = [];
pulses.din3_uptimes = [];
pulses.din4_uptimes = [];

pulses.din1_ts = [];
pulses.din2_ts = [];
pulses.din3_ts = [];
pulses.din4_ts = [];

pulses.din1_IPIViolationFlag = [];
pulses.din2_IPIViolationFlag = [];
pulses.din3_IPIViolationFlag = []; 
pulses.din4_IPIViolationFlag = [];
pulses.what_is_IPIViolationFlag = 'tldr: when correctSplitNEV is called on din uptimes AND a clock reset is detected an estimated timestamp is calculated marking the end of a segment. If the time between that timestamp and the last preceeding pulse is greater than the channel-specifc IPI, this flag is set to 1. read correctSplitNEV';

if ~isequal(nev_fpath, '') && exist(nev_fpath, 'file')

    NEVdata = openNEV( nev_fpath , 'nosave' , 'nomat');

    [pulses.din1_uptimes, pulses.din1_ts, pulses.din1_IPIViolationFlag] = getBlackRockPulsesDC(NEVdata, 9, postProc);
    [pulses.din2_uptimes, pulses.din2_ts, pulses.din2_IPIViolationFlag] = getBlackRockPulsesDC(NEVdata, 10, postProc);
    [pulses.din3_uptimes, pulses.din3_ts, pulses.din3_IPIViolationFlag] = getBlackRockPulsesDC(NEVdata, 11, postProc);
    [pulses.din4_uptimes, pulses.din4_ts, pulses.din4_IPIViolationFlag] = getBlackRockPulsesDC(NEVdata, 12, postProc);

    %- zero-pad the DC pulses so they are all the same length, and matching the ain_ts if that is longer/shorter
    %   alternative would be to force them to be same size as ain as long as ain exists and is not empty... (i.e., possible truncation)
    if isfield(pulses,'ain1_ts'), ainDuration = length(pulses.ain1_ts); else ainDuration = 0; end
    longestBin = max([ainDuration length(pulses.din1_ts) length(pulses.din2_ts) length(pulses.din3_ts) length(pulses.din4_ts)]);
    for iDIN=[1:4],
        thisLength = length(pulses.(sprintf('din%d_ts',iDIN)));
        pulses.(sprintf('din%d_ts',iDIN))(thisLength+1:longestBin) = 0;  %- zeropad so all are same length
    end
end 

end