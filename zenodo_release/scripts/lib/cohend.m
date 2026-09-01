function d = cohend(x, y)
% cohen's d: mean difference over sd of the difference (paired or one-sample)
if nargin < 2 || isempty(y)
    d = mean(x,'omitnan') / std(x,'omitnan');        % one-sample (vs 0)
else
    d = mean(x - y,'omitnan') / std(x - y,'omitnan'); % paired
end
end
