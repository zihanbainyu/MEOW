function r = rm_oneway(label, xc, xi, xn)
% one-way rm-anova, 3-level condition; integer df, gg p, partial eta^2
X = [xc(:) xi(:) xn(:)];
X = X(all(~isnan(X),2),:);
n = size(X,1);
t = array2table(X, 'VariableNames', {'c1','c2','c3'});
w = table(categorical({'compared';'isolated';'novel'}), 'VariableNames', {'Cond'});
a = ranova(fitrm(t, 'c1-c3~1', 'WithinDesign', w));
F = a.F(1); df1 = a.DF(1); df2 = a.DF(2);
etap = a.SumSq(1) / (a.SumSq(1) + a.SumSq(2));
e = gg_eps(X);
p = 1 - fcdf(F, df1*e, df2*e);
fprintf('%-22s F(%d,%d) = %.2f, %s (GG, eps = %.3f), eta_p^2 = %.3f  [n = %d]\n', ...
        label, df1, df2, F, pstr(p), e, etap, n);
r = struct('F',F,'df1',df1,'df2',df2,'p',p,'eps',e,'etap',etap,'n',n);
end
