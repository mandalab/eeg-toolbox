function test_FR( events )
%TEST_FR(events) uses an events structure of RAM_FR to ensure that it meets
%   requirements
outfile = 'TEST_FR_OUTPUT.txt';
% % 
if exist(outfile,'file')
    system(['rm ', outfile]);
end
diary(outfile)
disp([' FR VERSION ' events(2).expVersion]);
%% 12 word lists
WORD_events = events(strcmp({events.type},'WORD'));
[subjEvents, ~] = splitEventsBy(WORD_events, 'subject');
listLength = [];
for s = 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    sessionEvents = splitEventsBy(this_subjEvents, 'session');
    for se = 1:length(sessionEvents)
        this_sessEvents = sessionEvents{se};
        listEvents = splitEventsBy(this_sessEvents, 'list');
        for l = 1:length(listEvents)
            this_listEvents = listEvents{l};
            listLength = [listLength length(this_listEvents)];
        end
    end
end
fprintf(['****** LIST LENGTH *******' stats(listLength)]);
if all(listLength==12)
    fprintf('TEST: LISTLENGTH=12; PASS\n')
else
    fprintf('TEST: LISTLENGTH=12; FAIL\n')
end
%% 25 lists per session 

WORD_events = events(strcmp({events.type},'WORD'));
[subjEvents, ~] = splitEventsBy(WORD_events, 'subject');
nLists = [];
for s = 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    un_sessions = unique([this_subjEvents.session]);
    for se = 1:length(sessionEvents)
        this_sess_events = sessionEvents{se};
        nLists = [nLists length(unique([this_sess_events.list]))];
    end
end
fprintf(['****** 25 lists per session *******' stats(nLists)]);

if all(nLists==25)
    fprintf('TEST: LISTS/SESSION=25; PASS\n')
else
    fprintf('TEST: LISTS/SESSSION=25; FAIL\n')
end
%% 2 stim sessions (REMOVED DUE TO VARIABLE # STIM SESSIONS)


% WORD_events = events(strcmp({events.type},'WORD'));
% [subjEvents, ~] = splitEventsBy(WORD_events, 'subject');
% nSess = [];
% for s = 1:length(subjEvents)
%     this_subjEvents = subjEvents{s};
%     sessionEvents = splitEventsBy(this_subjEvents, 'session');
%     this_stimSess = 0;
%     for se = 1:length(sessionEvents)
%         if any([sessionEvents{se}.isStim])
%             this_stimSess = this_stimSess + 1;
%         end
%     end
%     nSess = [nSess this_stimSess];
% end
% fprintf(['****** FOUR STIM SESSIONS *******' stats(nSess)]);

%% UNIQUE LISTS BETWEEN SESSIONS 1 and 2

WORD_events = events(strcmp({events.type},'WORD'));
[subjEvents, ~] = splitEventsBy(WORD_events, 'subject');
isUnique = true;
for s = 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    sessionEvents = splitEventsBy(this_subjEvents, 'session');
    overlap = 0;
    for i = 1:length(sessionEvents)
        this_sessionEvents_i = sessionEvents{i};
        listEvents_i = splitEventsBy(this_sessionEvents_i, 'list');
        for list_i = 1:length(listEvents_i)
            this_listEvents_i = listEvents_i{list_i};
            for j = 1:length(sessionEvents)
                if i==j
                    continue
                end
                this_sessionEvents_j = sessionEvents{j};
                listEvents_j = splitEventsBy(this_sessionEvents_j, 'list');
                for list_j = 1:length(listEvents_j)
                    this_listEvents_j = listEvents_j{list_j};
                    if all([this_listEvents_i.itemno]==[this_listEvents_j.itemno])
                        isUnique = false;
                        break
                    end
                end
                if ~isUnique
                    break
                end
            end
            if ~isUnique
                break
            end
        end
        if ~isUnique
            break
        end
    end  
    if ~isUnique
        break
    end
end
if isUnique
    isUnique = 'YES';
else
    isUnique = 'NO';
end
fprintf('****** ARE LISTS UNIQUE BETWEEN SESSIONS? *******\n\t%s\n',isUnique);

if strcmp(isUnique,'YES')
    fprintf('TEST: UNIQUE LISTS BETWEEN SESSIONS; PASS\n')
else
    fprintf('TEST: UNIQUE LISTS BETWEEN SESSIONS; FAIL\n')
end
%% STIM WORDS SESS 1 not STIM WORDS SESS 2

WORD_events = events(strcmp({events.type},'WORD'));
[subjEvents] = splitEventsBy(WORD_events, 'subject');
stimRepeat = false;
foundSess12 = false;
for s = 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    sessionEvents = splitEventsBy(this_subjEvents,'session');
    sessions = cellfun(@(x)unique([x.session]),sessionEvents);
    sessIsStim = cellfun(@(x)[x.isStim]==1,sessionEvents,'uniformOutput',false);
    sessWords = cellfun(@(x)[x.itemno],sessionEvents,'uniformOutput',false);
    sessStimWords = cellfun(@(x,y)x(y),sessWords,sessIsStim,'uniformOutput',false);

    try
     sess1StimWords = sessStimWords{sessions==0};
    sess2StimWords = sessStimWords{sessions==1};
    
    stimRepeat = stimRepeat || any(ismember(sess1StimWords, sess2StimWords));
    foundSess12 = true;
    catch
        continue
    end
end
if stimRepeat
    stimRepeat = 'YES';
else
    stimRepeat = 'NO';
end
fprintf('******* DO STIM WORDS REPEAT BETWEEN SESSIONS 1 and 2? *******\n\t%s\n',stimRepeat);
if ~foundSess12
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 1 AND 2; SKIPPED (no sessions 1&2)\n')
elseif strcmp(stimRepeat,'NO')
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 1 AND 2; PASS\n')
else
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 1 AND 2; FAIL\n')
end
%% STIM WORDS SESS 3 not STIM WORDS SESS 4

WORD_events = events(strcmp({events.type},'WORD'));
[subjEvents] = splitEventsBy(WORD_events, 'subject');
stimRepeat = false;
foundSess34 = false;
for s = 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    sessionEvents = splitEventsBy(this_subjEvents,'session');
    sessions = cellfun(@(x)unique([x.session]),sessionEvents);  
    sessIsStim = cellfun(@(x)[x.isStim]==1,sessionEvents,'uniformOutput',false);
    sessWords = cellfun(@(x)[x.itemno],sessionEvents,'uniformOutput',false);
    sessStimWords = cellfun(@(x,y)x(y),sessWords,sessIsStim,'uniformOutput',false);
    try
    sess1StimWords = sessStimWords{sessions==2};
    sess2StimWords = sessStimWords{sessions==3};
    stimRepeat = stimRepeat || any(ismember(sess1StimWords, sess2StimWords));
    foundSess34 = true;
    catch
        continue
    end
end
if stimRepeat
    stimRepeat = 'YES';
else
    stimRepeat = 'NO';
end
fprintf('******* DO STIM WORDS REPEAT BETWEEN SESSIONS 3 and 4? *******\n\t%s\n',stimRepeat);
if ~foundSess34
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 3 AND 4; SKIPPED (no sessions 3&4)\n')
elseif strcmp(stimRepeat,'NO')
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 3 AND 4; PASS\n')
else
    fprintf('TEST: NO STIM REPEATS BETWEEN SESSIONS 3 AND 4; FAIL\n')
end
%% DOES STIM ALTERNATE?
fprintf('******* DO STIM WORDS ALTERNATE? ********\n')
WORD_events = events(strcmp({events.type},'WORD'));
subjEvents = splitEventsBy(WORD_events, 'subject');
inarow = [];
for s= 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    sessionEvents = splitEventsBy(this_subjEvents,'session');
    for se = 1:length(sessionEvents);
        this_sessionEvents = sessionEvents{se};
        listEvents = splitEventsBy(this_sessionEvents,'list');
        for l = 1:length(listEvents);
            list_inarow = [];
            prev_isStim = 3;
            this_listEvents = listEvents{l};
            for i=1:length(this_listEvents);
                this_item = this_listEvents([this_listEvents.serialpos]==i);
                fprintf('%d',this_item.isStim);
                if prev_isStim == this_item.isStim
                    list_inarow = list_inarow+1;
                else
                    inarow=[inarow, list_inarow];
                    prev_isStim = this_item.isStim;
                    list_inarow = 1;
                end
            end
            inarow=[inarow, list_inarow];
            fprintf('\n');
        end
        fprintf('------------\n')
    end
    fprintf('------------\n')
end

if all(inarow==2 | inarow==12)
    fprintf('TEST: STIM WORDS ALTERNATE; PASS\n');
else
    fprintf('TEST: STIM WORDS ALTERNATE; FAIL\n');
end
%% ON/OFF FIRST LISTS PER SESSION

WORD_events = events(strcmp({events.type},'WORD'));
subjEvents = splitEventsBy(WORD_events, 'subject');
all_onFirst = [];
all_offFirst = [];
for s= 1:length(subjEvents)
    this_subjEvents = subjEvents{s};
    fprintf('\t%s:\n',this_subjEvents(1).subject);
    sessionEvents = splitEventsBy(this_subjEvents,'session');
    for se = 1:length(sessionEvents);
        onFirst = 0;
        offFirst = 0;       
        this_sessionEvents = sessionEvents{se};
        listEvents = splitEventsBy(this_sessionEvents,'list');
        for l = 1:length(listEvents);
            this_listEvents = listEvents{l};
            this_item = this_listEvents([this_listEvents.serialpos]==1);
            if this_item(1).isStim
                onFirst = onFirst + 1;
            else
                offFirst = offFirst + 1;
            end
        end
        all_onFirst(end+1) = onFirst;
        all_offFirst(end+1) = offFirst;
        fprintf('\t%d: ON:%d \tOFF:%d\n',se,onFirst, offFirst)
    end
end

if all(all_onFirst==10) && all(all_offFirst==15)
    fprintf('TEST: 10 ON/OFF FIRST LISTS PER SESSION; PASS\n');
else
    fprintf('TEST: 10 ON/OFF FIRST LISTS PER SESSION; FAIL\n');
end
%% FOR EACH STIM CONDITION, 3 ON LISTS IN ONE SESSION, 2 IN THE OTHER
% 
% WORD_events = events(strcmp({events.type},'WORD'));
% subjEvents = splitEventsBy(WORD_events, 'subject');
% TRIAL_events = events(strcmp({events.type},'TRIAL'));
% trialNumbers = [TRIAL_events.list];
% sessNumbers = [TRIAL_events.session];
% stimTypes = arrayfun(@(x)[x.stimLoc, x.stimAmp], TRIAL_events,'uniformoutput',false);
% unique_stimTypes = unique(stimTypes);
% unique_sessions = unique([events.session]);
% unique_subjects = unique({events.subject});
% onFirsts = zeros(length(unique_stimTypes),length(unique_sessions), length(unique_subjects));
% for s= 1:length(subjEvents)
%     this_subjEvents = subjEvents{s};
%     sessionEvents = splitEventsBy(this_subjEvents,'session');
%     for se = 1:length(sessionEvents); 
%         this_sessionEvents = sessionEvents{se};
%         listEvents = splitEventsBy(this_sessionEvents,'list');
%         for l = 1:length(listEvents);
%             this_listEvents = listEvents{l};
%             this_listNum = this_listEvents(1).list;
%             this_sessNum = this_listEvents(1).session;
%             this_stimType = stimTypes(trialNumbers==this_listNum & sessNumbers==this_sessNum );
%             this_item = this_listEvents([this_listEvents.serialpos]==1);
%             if this_item(1).isStim
%                 onFirsts(strcmp(unique_stimTypes,this_stimType),se, s) = ...
%                     onFirsts(strcmp(unique_stimTypes,this_stimType),se,s) + 1;
%             end
%         end
%     end
% end
% 
% fprintf('****** ONFIRSTS DISTRIBUTED EVENLY ******\n');
% disp(unique_stimTypes)
% disp(onFirsts)
% 
% if all(all(onFirsts(2:end,:)==5)) || ...
%     all(xor(onFirsts(2:end,2:2:end)==3,onFirsts(2:end,2:2:end)==2))
%     
%     fprintf('TEST: FOR GIVEN STIM CONDITION, 3 IN 1 SESSION, 2 IN OTHER: PASS\n');
% else
%     fprintf('TEST: FOR GIVEN STIM CONDITION, 3 IN 1 SESSION, 2 IN OTHER: FAIL\n');
% end
%% STIM TRIGGERED 200ms prior to onset 
STIM_events_mask =  strcmp({events.type},'STIM_ON');
delays = [];
for i=find(STIM_events_mask)
    delays(end+1) = events(i+1).mstime - events(i).mstime;
end

fprintf(['******* STIM BEFORE WORD ONSET (ms) *******',stats(delays)]);
if all(delays>150 & delays<250)
    fprintf('TEST: STIM TRIGGERED 200MS PRIOR TO WORD ONSET: PASS\n')
else
    fprintf('TEST: STIM TRIGGERED 200MS PRIOR TO WORD ONSET: FAIL\n');
end
%% 1600+750 to 1000 (2350-2600) between word presentations
WORD_events_mask = strcmp({events.type},'WORD');
delays =[];
for i=find(WORD_events_mask)
    if strcmp(events(i+1).type,'WORD')
        delays(end+1) = events(i+1).mstime - events(i).mstime;
    end
end

fprintf(['******* TIME BTWN WORDS (ms) *******',stats(delays)]);
if all(delays>2250 & delays<2700)
    fprintf('TEST: 2350-2600ms BETWEEN WORDS: PASS\n')
else
    fprintf('TEST: 2350-2600ms BETWEEN WORDS: FAIL\n');
end
%% 20(+) second distractor
DISTRACT_events = events(strcmp({events.type},'DISTRACT_START') | strcmp({events.type},'DISTRACT_END'));
delays = [];
for i=1:length(DISTRACT_events)
    if strcmp({DISTRACT_events(i).type},'DISTRACT_START') && ...
            strcmp({DISTRACT_events(i+1).type},'DISTRACT_END')
        delays(end+1) = DISTRACT_events(i+1).mstime - DISTRACT_events(i).mstime;
    end
end
fprintf(['******* DISTRACTOR TIME (ms) *******',stats(delays)]);
if all(delays>20000)
    fprintf('TEST: 20000+ms DISTRACTOR: PASS\n')
else
    fprintf('TEST: 20000+ms DISTRACTOR: FAIL\n');
end
%% 30s retrieval
REC_events = events(strcmp({events.type},'REC_START') | strcmp({events.type},'REC_END'));
delays = [];
for i=1:length(REC_events)
    if strcmp({REC_events(i).type},'REC_START') && ...
            strcmp({REC_events(i+1).type},'REC_END')
        delays(end+1) = REC_events(i+1).mstime - REC_events(i).mstime;
    end
end
fprintf(['******* RECALL TIME (ms) *******',stats(delays)]);
if all(delays>29900 & delays<30100)
    fprintf('TEST: 30000ms RETRIEVAL: PASS\n')
else
    fprintf('TEST: 30000ms RETRIEVAL: FAIL\n');
end   
diary off
function strstats = stats(list)
strstats = sprintf('\n\tMEAN:%.1f\n\tMIN:%d\n\tMAX:%d\n',...
    nanmean(list(~isnan(list))),...
    min(list(~isnan(list))),...
    max(list(~isnan(list))));
    
    
function [ events_split, unique_values ] = splitEventsBy( events, field )
%SPLITBY Summary of this function goes here
%   Detailed explanation goes here
if ischar(events(1).(field))
    unique_values = unique({events.(field)});
else
    unique_values = unique([events.(field)]);
end
events_split = cell(size(unique_values));
for i=1:length(events)
    if iscell(unique_values)
        index = strcmp(unique_values,events(i).(field));
         events_split{index} = [events_split{index} events(i)];
    else
        index = unique_values==events(i).(field);
         events_split{index} = [events_split{index} events(i)];
    end
    
   

end

