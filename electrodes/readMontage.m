function [all_index, all_jacks, all_label] = readMontage(subj, rootEEGdir_or_montagePath)
% [all_index, all_jacks, all_label] = readMontage(subj, rootEEGdir_or_montagePath)

    all_index = [];
    all_jacks = [];
    all_label = [];

    if exist(rootEEGdir_or_montagePath, 'file') && ~isempty(strfind(rootEEGdir_or_montagePath, 'xls')),
        % passed path to montage
        filename = rootEEGdir_or_montagePath;
    else
        % try to find a montage
        docs = fullfile(rootEEGdir_or_montagePath, subj, 'docs');
        assert(exist(docs, 'dir') > 0, 'No docs folder: %s\n', docs);
        filename = fullfile(docs, 'clinical_montage.xlsx');
        assert(exist(filename) > 0, 'Montage not found: %s', filename);
%         docs_contents = lsCell(docs);
%         ndx = menuText('Choose the montage file', docs_contents, 'multiselect',0);
%         if isempty(ndx), return; end
%         filename = fullfile(docs, docs_contents{ndx});
    end
    if contains(filename, 'xlsb')
        error('You need to save the xlsb file as another filetype');
    end
    disp(filename)
    
    [STATUS,SHEETS] = xlsfinfo(filename); %- get the list of sheets
    
    sheets_remain = 1;
    i_sheet = 1;
    while sheets_remain
        try
            [num,txt,raw] = xlsread(filename, SHEETS{i_sheet});
            if isempty(raw) 
                i_sheet=i_sheet+1; 
                continue;
            end
        catch e
            sheets_remain = ~strcmp(e.identifier, 'MATLAB:xlsread:WorksheetNotFound');
            if sheets_remain, rethrow(e); % rethrow other errors
            else, continue; end 
        end
        
        disp(raw);
        is_montage_sheet = inputYN('Is this a montage sheet?');
        
        if is_montage_sheet
            index = [];
            jacks = [];
            label = [];
            
            % returned from xlsread are num, txt, and raw variables
            % num may be smaller than txtx and raw because left and top-most
            % non-numeric rows/cols are trimmed. Try to fix this by finding
            % how many rows got trimmed off of "num"
            numeric_mask = cellfun(@(x) isscalar(x) && isnumeric(x) && ~isnan(x), raw);
            [dx,dy] = find(numeric_mask,1);
            dx = dx - 1;
            dy = dy - 1;
            
            assert(~isempty(num), 'No numeric labels');
            nums = sort(num(~isnan(num)));
            
            sz = size(num,1);
            
            for ndx_num = min(nums):max(nums)
                [i,j]  = ind2sub(sz, find(num==ndx_num, 1));
                jack   = raw{i+dx,j+1+dy}; % bank+ndx e.g. B12
                lbl  = raw{i+dx,j+2+dy};
                
                fprintf('%d\t%s\t%s\n', num(i,j), jack, lbl);
                index = [index num(i,j)];
                jacks = [jacks {jack}];
                label = [label {lbl}];
                
            end
            
            % put in cells if there's more than one sheet
            if ~isempty(all_index)
                index = {index};
                jacks = {jacks};
                label = {label};
            end
            
            all_index = [all_index index];
            all_jacks = [all_jacks jacks];
            all_label = [all_label label];
        end
        
        i_sheet = i_sheet + 1;
    end

end