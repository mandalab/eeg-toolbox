function outDist = klGetDist(obs1,obs2,varargin)

% Set defaults
type = 'euc';
catVars = zeros(1,length(obs1));

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'cat'}
            catVars = varargin{varStrInd(iv)+1};
    end
end

    switch type,
        case 'corr',
            nanCrit = ~isnan(obs1) & ~isnan(obs2);
            if sum(nanCrit) == 0,
                outDist = nan;
            else
                outDist = 1-corr(obs1(nanCrit)',obs2(nanCrit)');
            end
        case 'euc'
            outDist = sqrt(nansum((obs1-obs2).^2));
        case 'sqeuc',
            outDist = (nansum((obs1-obs2).^2));
        case 'prod'
            allSim = nan(1,size(obs1,2));
            if sum(catVars) > 0,
                allSim(catVars) = abs(double(obs1(catVars)) - double(obs2(catVars))) < .00001;
            end
            allSim(~catVars) = 1-abs(obs1(~catVars) - obs2(~catVars));
%                     thisDist = 1-nanprod(1-allDist);
            outDist = 1-(nanprod(allSim)^(1/length(allSim)));
        case 'exprod',
            allSim = nan(1,size(obs1,2));
            if sum(catVars) > 0,
                allSim(catVars) = obs1(catVars) ~= obs2(catVars);
            end
            allSim(~catVars) = abs(obs1(~catVars) - obs2(~catVars));
            outDist = 1-nanprod(exp(-allSim));
    end
end