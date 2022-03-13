function [p0] = pAdjust(p)
%------------------
%%- Conversion of R function to get list of adjusted p-values
%%- Pass in a vector of p values that represent multiple comparisons
%
% just focus on "BH" technique first

if size(p,1)>1 & size(p,2)==1, p=p'; end

n = length(p);
p0 = p ;

if (n <= 1)
    return
end

if sum(p)==n | sum(isnan(p))==n, return; end; %- if all 1's or nan's passed in, then pass them back out


%%% bonferroni
%p0 = pmin(1, n * p);

%%% BH
i = n:-1:1 ;
[y,o] =sort(p,'descend');
[y,ro]=sort(o);

%baby 'cummin' function: min value up to this point in the array
for ii=1:n,
    pTemp(ii) = min(1, min(n./i(1:ii) .* p(o(1:ii))));
end
p0(1:n) = pTemp(ro);

return