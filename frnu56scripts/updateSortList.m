
PUBpath = '/Volumes/56PUB/readWrite/micro_forSorting';

sub_list = getDirNamesRegexp(PUBpath,'NIH.*');

for ii = 32:numel(sub_list)
    subpath = fullfile(PUBpath,sub_list{ii});
    sess_list = getDirNamesRegexp(subpath,'\d{6}_\d{4}.*');
    omit_list = getDirNamesRegexp(subpath,'\d{6}_\d{4}.*(ready|grabbed).*');
    sess_dir = sess_list(~ismember(sess_list,omit_list));
    for jj = 1:numel(sess_dir)
        sesspath = fullfile(subpath,sess_dir{jj});
        if numel(getDirNamesRegexp(sesspath,'sortNotes.*sortedBy.*')) ~= 1
            keyboard
            continue
        end
        sorttablepath = fullfile(sesspath,char(getDirNamesRegexp(sesspath,'sortNotes.*sortedBy.*')));
        sortTable = readtableSafe(sorttablepath);
        if any(strcmp('dummySort',sortTable.Properties.VariableNames))
            keyboard
        end
        if any(strcmp('chanNoise',sortTable.Properties.VariableNames)) && any(strcmp('globalNoise',sortTable.Properties.VariableNames))
            if iscell(sortTable.globalNoise)
            else
                if all(isnan(sortTable.globalNoise)) && all(isnan([sortTable.chanNoise]))
                    % Nothing to do- both all nan
                    fprintf('%s\n',sesspath);
                    continue
                end
            end
            if any(contains(sortTable.globalNoise,'???'))
                if ~all(strcmp(sortTable.globalNoise,'') | strcmp(sortTable.globalNoise,'???')) || ~all(isnan(sortTable.chanNoise)) %~all(strcmp(sortTable.chanNoise,'') | strcmp(sortTable.chanNoise,'???'))
                    keyboard % Contains a number
                elseif all(contains(sortTable.globalNoise(1:3),'???'))
                    keyboard
                end
                chanNoiseNEW = sortTable.globalNoise;
                sortTable.globalNoise = sortTable.chanNoise;
                sortTable.chanNoise = chanNoiseNEW;
            end
        else
            keyboard
        end

%         if any(strcmp('dummySort',sortTable.Properties.VariableNames))
%             if all(~strcmp('chanNoise',sortTable.Properties.VariableNames)) && all(~strcmp('globalNoise',sortTable.Properties.VariableNames))
%                 sortTable.Properties.VariableNames{'dummySort'} = 'globalNoise';
%                 chanNoise = NaN(numel(sortTable.globalNoise),1);
%                 sortTable = addvars(sortTable,chanNoise,'Before','globalNoise');
%             else
%                 keyboard
%             end
%         else
%             keyboard
%         end
        writetable(sortTable,sorttablepath);
        fprintf('%s\n',sesspath);
    end
    fprintf('%s\n',[subpath ' complete']);
end














