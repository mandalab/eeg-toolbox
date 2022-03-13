function events = create_stim_events_v1(events3,eegoffset,INFO)
% function create_stim_events_v1(events3,eegoffset,INFO)
% INFO.subject = 'NIH045';
% INFO.session = 0;
% INFO.f_stem  = '..'
%
% eegfile: '/Volumes/Shares/FRNU/data/eeg/NIH042/eeg.reref/NIH042_160615_1116'
% eegoffset: 277093
%
% function to generate events for subject built with GUI 

%- initialize
fields_needed = {
    'experiment'        %- always stimMapping
    'subject'           %- Subj name eg. NIH0XX
    'sessionName'       %- '0a'
    'session'           %- 0
    
    'stimType'          % start w/ SM, CC, or WN
    'amplitude'         % btw 100-8000
    'polarity'          % 0 or 1
    'frequency'         % as defined by stimulator
    'stimLocTag'        % eg 'UNKNOWN' or 'ALT5-ALT6'
    'stimLocNum'        % eg [ 10, 11] or [ NaN, NaN]
    
    'annotation'         %- more for old data, could be useful
    
    'isStimEvent'       %- 0 or 1
    
    'numPulses'         % for older CCEP, maybe don't need
    
    'isWN_REP'          % 0 or 1
    'WN_REP_obj'        % [] or wn eobj
    
    'GUIVer'            %- 0 for pre-GUI, 1- if GUI but unknown
    'isSingleWire'      %- 1 if preGUI or telemark
    'NKreference'       %- NaN, [5,6], or [other,other]
    
    'cumulative_pulses' %- a key element for alignment, best not to eliminate
    'num_dc_pulses'
    
    'eegfolder' 
    'eegfile'           %- essential
    'eegoffset'         %- essential
    };


b_event = struct();
for f = 1:length(fields_needed)
    b_event.(fields_needed{f}) = [];
end

ex_ev = events3(2);
%- get constants
C.experiment    = 'stimMapping';
C.subject       = INFO.subject;
C.session       = INFO.session;
C.eegfile       = INFO.eegfile;
C.eegfolder     = INFO.eegfolder;
C.sessionName   = ex_ev.session{1};
C.GUIVer        = ex_ev.GUIver;

% create events
b_event.experiment  = C.experiment;
b_event.subject     = C.subject;
b_event.session     = C.session;
b_event.eegfile     = C.eegfile;
b_event.eegfolder   = C.eegfolder;
b_event.sessionName = C.sessionName;
b_event.GUIVer      = C.GUIVer;



%- create events, asign to eegoffset times based on cumulative pulses!
% cum_pulses =    [events3.cumulative_pulses];
% dc_pulses =     [events3.DCpulses];

for i = 2:length(events3)
    this_ev3 = events3(i);
    ev = b_event;   %- initialize
    
    if this_ev3.numPulses >= 1 %-only do if pulses actually administered
        pulse_st_num = this_ev3.cumulative_pulses - this_ev3.DCpulses + 1;
        ev.eegoffset = eegoffset(pulse_st_num); %-make double? [NO!]
    else
        ev.eegoffset = NaN; %- otherwise last wn_rep event will surpass # of DC10 pulses (eegoffset)
    end
        
    ev.num_dc_pulses = this_ev3.DCpulses;  
    ev.cumulative_pulses = this_ev3.cumulative_pulses;


    ev.annotation = this_ev3.ann;
    ev.isSingleWire = this_ev3.tele;
    ev.NKreference = this_ev3.NK_referenece;
    ev.isStimEvent = true; %- given how we are building events, always true
    
    ev.stimType = this_ev3.typeStr;
    ev.amplitude = this_ev3.amp1;
    ev.polarity = this_ev3.pol;
    ev.frequency = this_ev3.freq;
    ev.stimLocTag = this_ev3.elec_loc;  %- blank
    ev.stimLocNum = split_tag(this_ev3.location);
    ev.numPulses = this_ev3.numPulses;
    ev.isWN_REP = this_ev3.isWN_REP;
    ev.WN_REP_obj = this_ev3.WN_REP_obj;
    
    if ~exist('events','var')
        events = ev;
    else
        events = [events ev];
    end
    
end

function tag = split_tag(in_tag)

if ischar(in_tag)
    ees = strsplit(in_tag,{'-','_'});
    tag = [0 0];
    for i = 1:length(ees)
        tag(i) = str2double(ees{i});
    end
    
    return
else
    tag = in_tag;
end

