function newTable = addBlankTableRow(origTable)
% ADDBLANKTABLEROW adds a blank row to a table
%   
% Inputs:
%   origTable - table to add
%
% Outputs:
%   newTable

% REVISION HISTORY:
%  10/16 MST - Created

    assert(istable(origTable));
    assert(height(origTable) > 0);
    
    templateRow = origTable(1, :);    
    newRow = struct();
    
    for i = 1 : length(origTable.Properties.VariableNames)
        col = origTable.Properties.VariableNames{i};
        val = templateRow{1, i};
        if iscell(val)
            newRow.(col) = {''};
        elseif isnumeric(val)
            newRow.(col) = NaN;
        elseif isdatetime(val)
            newRow.(col) = NaT;
        else
            fprintf('Value type: %s\n', class(val));
            error('addBlankTableRow Unrecognized column type');
        end
         
        %newTable{1, i} = val;
    end
    
    newTable = [origTable; struct2table(newRow, 'AsArray', true)];
    

end