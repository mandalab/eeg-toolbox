function [chanArrayMask,arrayList,firstIdxArray] = getArrayFromSortNoteTable(sortNoteTable)
% getArrayFromSortNoteTable: Get the separate arrays and channel info from a sortNoteTable
%
% Outputs:
%       chanArrayMask       =   Array corresponding to each channel with a number indicating which array that channel
%                               belongs to. The number of the array corresponds to the array name in arrayList
%                               that is at that number index.
%       arrayList           =   List of unique array names, where the position of each array name in arrayList
%                               corresponds to the number in chanArrayMask
%       firstIdxArray       =   Index of the first channel in each array, in which each index in firstIdxArray corresponds 
%                               with the indices in arrayList 
%
% Called by extract_SpiekInfo_v3b (and eventually when sortNoteTable gets created)
%
%
% Created by Samantha N. Jackson 4/2/2020
%
%
% 4/2/2020 SJ: Created function
% 4/1/2021 SJ: Fixed if statement from || to &&

sortNotes = sortNoteTable;

chanNames = sortNotes.ChanNameNew;

% Look for certain conditions:
% 1) micro: any length (usually 16, sometimes 10?)
% 2) utah: 96 chans
% 3) utah: 64 chans (128 total)

tot_num_arrays = 0;
tot_array_names = {};
tot_chan_mask = zeros(numel(chanNames),1);
tot_chan_first_idx = [];

% First do Micro
%micro_mask = contains(chanNames,'micro');
micro_chans = chanNames(contains(chanNames,'micro'));

if ~isempty(micro_chans)
    % Find out how many micro channel arrays there are...?
    % First, strip off micro and numbers
    microArrays = unique(cellstr(string(regexp(micro_chans,'(?<=micro)\D*(?=\d*)','match'))));
    for ii = 1:numel(microArrays)
        tot_num_arrays = tot_num_arrays + 1;
        thisArrayStr = microArrays{ii};
        idx_thisArrayStr = contains(chanNames,'micro') & contains(chanNames,thisArrayStr);
        % Check to see if the number of channels makes sense
        if sum(idx_thisArrayStr) ~= 10 && sum(idx_thisArrayStr) ~= 16
            fprintf('%s\n',['Usually a micro array will have 16 or 10 channels. This one (' thisArrayStr ') has ' num2str(sum(idx_thisArrayStr)) '. Check it out before continuing! Ask SJ!']);
            keyboard
        end

        if sum(tot_chan_mask(idx_thisArrayStr)) ~= 0
            fprintf('%s\n','ERROR!!! This is to check to make sure we are not assigning channels twice!! What is going on?? Ask SJ!');
            keyboard
        end
        tot_chan_mask(idx_thisArrayStr) = tot_num_arrays;
        tot_array_names{tot_num_arrays} = ['micro' thisArrayStr];
        tot_chan_first_idx(tot_num_arrays) = find(idx_thisArrayStr,1);
    end
end

% Now do Utah
utah_chans = chanNames(contains(chanNames,'utah'));
if ~isempty(utah_chans)
    utahArrays = unique(cellstr(string(regexp(utah_chans,'(?<=utah)\D*(?=\[\w*\])','match'))));
    
    if numel(utahArrays) == 1
        % Only 1 array... That means this is most likely a 96 channel array. Do a quick check to make sure the numbers make sense
        switch numel(utah_chans)
            case 96
                % Good to go - this makes the most sense
            case 64
                fprintf('%s\n','There is only one 64 chan Utah array? Double check please, may be okay. Ask SJ!');
                keyboard
            case 128
                fprintf('%s\n','ERROR!!!! Only picking up one Utah array even though there are 128 channels! There should be two arrays!! Ask SJ!');
                keyboard
            otherwise
                fprintf('%s\n',['ERROR!!!! Picking up a utah array, but the number of channels (' num2str(numel(utahArrays)) ') makes no sense! Investigate! (Ask SJ)']);
                keyboard
        end
    elseif numel(utahArrays) == 2
        % Damn... more than 1 array, gotta split it... Can be 2 64s or 2 96s
        switch numel(utah_chans)
            case 192
                % Okay, so two 96 arrays
            case 128
                % Okat, so two 64 arrays
            otherwise
                fprintf('%s\n','ERROR!!! We have 2 Utah arrays but the number of total channels does not indicate two 96 chan or two 64 chan arrays! Ask SJ!');
                keyboard
        end
    else
        fprintf('%s\n','ERROR!!! We have more than 2 Utah arrays?? What is going on? May be okay, but SJ will have to write this into the code!');
        keyboard
    end
    
    for jj = 1:numel(utahArrays)
        tot_num_arrays = tot_num_arrays + 1;
        thisArrayStr = utahArrays{jj};
        idx_thisArrayStr = contains(chanNames,'utah') & contains(chanNames,thisArrayStr);
        
        if sum(tot_chan_mask(idx_thisArrayStr)) ~= 0
            fprintf('%s\n','ERROR!!! This is to check to make sure we are not assigning channels twice!! What is going on?? Ask SJ!');
            keyboard
        end        
        
        tot_chan_mask(idx_thisArrayStr) = tot_num_arrays;
        tot_array_names{tot_num_arrays} = ['utah' thisArrayStr];
        tot_chan_first_idx(tot_num_arrays) = find(idx_thisArrayStr,1);
    end
    
end

% Check to make sure the micro mask doesn't have any zeros left!
if any(tot_chan_mask == 0)
    fprintf('%s\n','ERROR!!! Not all channels accounted for by array mask!');
    keyboard
end

chanArrayMask = tot_chan_mask;
arrayList = tot_array_names;
firstIdxArray = tot_chan_first_idx;


end

