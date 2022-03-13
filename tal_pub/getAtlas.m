function atlas = getAtlas(subj, rootEEGdir, chanName, varargin)
% atlas = GETATLAS(subj, rootEEGdir, chanName, ...)
% 
% Optional Input (key-value pairs)
%
%
    ip = inputParser;
    ip.addParameter('chanNames',[]);
    parse(ip, varargin{:});
    chanNames = ip.Results.chanNames;
    
    d_atlas = fullfile(rootEEGdir, subj, 'tal/zloc/atlas');
    fname_atlas = 

    filename = fullfile(rootEEGdir, subj, 'tal/monopolar.mat');
    if ~exist(filename,'file')
        % for Mike
        dirs = evalin('base','dirs');
        d = {dirs.finalRun, dirs.dataWorking, dirs.data};
        for i = 1:3
            rootEEGdir = d{i};
            filename = fullfile(rootEEGdir, subj, 'tal/monopolar.mat');
            if exist(filename,'file'), break; end
        end
        
    end
    assert(exist(filename,'file') > 0, 'File not found: %s\n', filename);
    data = load(filename);
    monopolar = data.monopolar;
    
    % order correctly (according to leads order)
    order = sortBySubstring({monopolar.chanName}, getLeads(subj, rootEEGdir));
    monopolar = monopolar(order);
    
    % filter
    if ~isempty(chanNames)
        monopolar = monopolar(ismember(upper({monopolar.chanName}), upper(chanNames)));
    end
end
