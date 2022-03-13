function [memUsedGB, memFreeGB, memTotalGB, memComp] = getMemoryUsage()
%-  This function poles the system to figure out how much free memory is available
%-   really only important so we can figure out if/when "extractMicrophys" is going to crash
%
%- 2/2019  JW adapted from "getMemUseage" and "monitor_memory_whos()" from matlab central
%
% https://superuser.com/questions/197059/mac-os-x-sysctl-get-total-and-free-memory-size
% https://apple.stackexchange.com/questions/81581/why-does-free-active-inactive-speculative-wired-not-equal-total-ram


%- total system memory (hardware RAM)
[~, out]=system('sysctl hw.memsize | awk ''{print $2}''');
mem=sscanf(out,'%f');  %- hardware memory
memTotalGB = mem / 2^30;


%- JW trying to figure out what to pull from vm_stat to figure out memory stuff
%   page size is 4096 bytes.  use sscanf to pull the number then scale to GB
[~,out]=system('vm_stat | grep "Pages free"');
mem=sscanf(out,'Pages free: %f.');
memFreeGB=mem*4096/2^30;  %- DOES NOT MATCH actiity monitory (but close to total-memUsed)


%%- additional outputs from VM, don't seem to help nail down usage, so dont bother parsing for now
[~,out]=system('vm_stat | grep "Pages speculative"');
mem=sscanf(out,'Pages speculative: %f.');
memSpec=mem*4096/2^30;


%- App mem = Pages Active + Pages Speculative?
[~,out]=system('vm_stat | grep "Pages active"');
mem=sscanf(out,'Pages active: %f.');
memActive=mem*4096/2^30;

[~,out]=system('vm_stat | grep "Pages inactive"');
mem=sscanf(out,'Pages inactive: %f.');
memInactive=mem*4096/2^30;


[~,out]=system('vm_stat | grep "Pages wired down"');
mem=sscanf(out,'Pages wired down: %f.');
memWired=mem*4096/2^30;  %- matches Wire3dMemory


[~,out]=system('vm_stat | grep "Pages occupied by compressor"');
mem=sscanf(out,'Pages occupied by compressor: %f.');
memComp=mem*4096/2^30;  %- matches Compressed


% %- try different combinations to find the thing that matches activity monitor (still haven't found it)
% memUsed0 = memActive+memWired+memComp %- more than this, but doesn't include spec
% memUsed1 = memTotalGB - memFreeGB  
% memUsed2 = memTotalGB - memFreeGB - memInactive 
% memUsed3 = memTotalGB - memFreeGB - memInactive - memWired - memComp


memUsedGB  = memActive+memWired+memComp;



%- total memory: free + active + specualtive + cahsed + wired + compressed
%- memory used = app + wired + compressed


%% OUTPUT of system('vm_stat')
% 'Mach Virtual Memory Statistics: (page size of 4096 bytes)
% Pages free:                             8295673.
% Pages active:                           1097711.
% Pages inactive:                         3532101.
% Pages speculative:                         2582.
% Pages throttled:                              0.
% Pages wired down:                       1270730.
% Pages purgeable:                             14.
% "Translation faults":                 273977359.
% Pages copy-on-write:                     429192.
% Pages zero filled:                    156947391.
% Pages reactivated:                     36470737.
% Pages purged:                             20055.
% File-backed pages:                      3479077.
% Anonymous pages:                        1153317.
% Pages stored in compressor:             4482651.
% Pages occupied by compressor:           2572478.
% Decompressions:                        98295570.
% Compressions:                         154843230.
% Pageins:                               89252576.
% Pageouts:                                 13041.
% Swapins:                               10920893.
% Swapouts:                              12842498.

