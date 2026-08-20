function d = cohend(x, y)
% COHEND  Cohen's d for a paired/one-sample comparison:
%   mean difference divided by the SD of the difference.
if nargin < 2 || isempty(y)
    d = mean(x,'omitnan') / std(x,'omitnan');       % one-sample (vs 0)
else
    d = mean(x - y,'omitnan') / std(x - y,'omitnan');% paired
end
end
