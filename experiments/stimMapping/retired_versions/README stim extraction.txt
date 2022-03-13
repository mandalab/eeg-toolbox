Notes on stim extraction codes:



%%%%%—— NIH020 to 35… possibly 36-40-ish as well  —%%%%%

I) createStimMapSessFromAnn_v01   written by JW 7/2017.  
This creates behaviora/stimMapping/session_0_(XXXXX_XXXX) folders with annotations and updown files use to convert old-school stim mapping sessions into events.   
Valid for subjects <NIH040 or so.  Will have problems with NIH020 to NIH024, as they didn’t use Cocjin’s annotation shorthand (seems OK in practice tho)

II) createStimEvents_v3a.m,  modified by TS 1/10/2016, possibly valid for all v1b subjects. JW cleaned up 8/2017 to fit better into prepAndAlign



%%%%%——  NIH036-ish to 40-ish?   —%%%%%
createStimEvents_v4.m, modified by TS 12/2/2016,  eliminates annotation file processing… 
perhaps intended to be used with Cerestim file to figure out which channel was stimulated?  not obvious if that is the case



%%%%%——  NIH044 on (currently at NIH052)   —%%%%%
extractStimEvents_v2.m  created by TS, modified by CS and JW
This extracts events from Tim’s stim mapping gui’s outputs.  This is how it should be done!

Thia says here process is:
ProcessCurrentSubj
produce_stim_events
create_stim_events_v1



%- jw just found these …. what do they do?  what subj are they for?  looks like maybe thia wrote them?
processStimLog.m ???
produce_stim_events_v1.m ???




————————————————————————————————————————————————————————————————————————————————————————

RETIRED FOLDER


eegStimPrep_v3a.m     written by JW… just a hack of an older version of prepAndAlign that specifically splits raw files in the raw/STIM folder. This should not be used anymore.    If you want to process old stim mapping data (before cerestim and Tim’s skimping GUI), then write a new function (createStimMapSessFromAnn_v01) based on pulseViusalize that doesn’t actually split the stim files, but instead just gets the pulses and annotation files and puts them in behavioral/stimMapping

createStimEvents_v1b.m    written by JW in 2015 or so, valid for subj 24-35?
used to create events from Nihon Khoden manual annotation combined with recorded DC10/11/12.  Requires eegStimPrep_v3a or something like that to open the NK files and split out just the DC’s and annotations into a text file that is manually edited so it: 
(1) conforms with John Cocjin’s annotation style (e.g., “i” means current level)… that started with subj NIH024 
(2) so annotations are correctly assigned to stimulation events (e.g., if annotation of channel change occurs after stim starts, manually correct so it occurs before stim starts so event has correct channel assignment) 

