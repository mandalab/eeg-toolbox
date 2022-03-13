Creating RAM_FR Events:

==========================================================================

**UPDATE 07/08/16: 
Modified so that CreateAllEvents can be called by behavioralProcessing.m; 
**functionality is no longer exactly as described below**

==========================================================================

Step 1:
Run the matlab function RAM_FR_CreateAllEvents. This function calls RAM_FR_CreateTASKEvents and RAM_FR_CreateMATHEvents, which create the recall task and math events structs, respectively. RAM_FR_CreateAllEvents operates on one session at a time and saves the events.mat and MATH_events.mat files in the given session directory.

RAM_FR_CreateAllEvents takes in three required inputs:
--subject: (string) code for subject, e.g. 'TJ083'
--expDir: (string) directory where [session_X] directories live, e.g. '/data/eeg/TJ083/behavioral/RAM_FR/nonstim_data/TJ083/'
--session: (numeric) the session number for which to create events, e.g. 0

Example usage:
RAM_FR_CreateAllEvents('TJ083','/data/eeg/TJ083/behavioral/RAM_FR/nonstim_data/TJ083/',0);


Step 2:
Align the events with the eeg data (using AlignTool method).


Step 3:
Run the matlab function RAM_FR_AddEventsToDatabase. This function takes the events for a specific session and adds them to the 'full' events structure that is contained in /data/events/RAM_FR/.

RAM_FR_AddEventsToDatabase takes four required inputs:
--subject: (string) code for subject, e.g. 'TJ083'
--expDir: (string) directory where the events.mat file is for the session you want to add, e.g. '/data/eeg/TJ083/behavioral/RAM_FR/nonstim_data/TJ083a/session_0'
--session: (numeric) the session number that you want to assign to the added events in the full events structure in the database. NOTE: THIS IS NOT NECESSARILY THE SAME AS THE NUMBER IN ~/session_X
--RAMFRver: (string) the version of RAM_FR to which the events belong. ex: 'RAM_FR1' (passive record task) or 'RAM_FR2' (stimulation task)


