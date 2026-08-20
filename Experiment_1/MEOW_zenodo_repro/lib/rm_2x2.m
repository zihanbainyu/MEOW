function r = rm_2x2(label, cc, ci, ic, ii)
% RM_2x2  2x2 repeated-measures ANOVA (condition x correctness).
%   Every effect has 1 numerator df, so Greenhouse-Geisser is a no-op
%   (epsilon = 1): integer df, ordinary p, partial eta^2.
%   Column order: cc = comp/correct, ci = comp/incorrect,
%                 ic = iso/correct,  ii = iso/incorrect.
M = [cc(:) ci(:) ic(:) ii(:)];
M = M(all(~isnan(M),2),:);
n = size(M,1);
t  = array2table(M, 'VariableNames', {'cc','ci','ic','ii'});
w  = table({'comp';'comp';'iso';'iso'}, {'cor';'inc';'cor';'inc'}, ...
           'VariableNames', {'Cond','Acc'});
rm = fitrm(t, 'cc-ii~1', 'WithinDesign', w);
a  = ranova(rm, 'WithinModel', 'Cond*Acc');
rn = string(a.Properties.RowNames);
term = {'(Intercept):Cond','(Intercept):Acc','(Intercept):Cond:Acc'};
lab  = {'condition','accuracy','interaction'};
fprintf('%s  [n = %d]\n', label, n);
r = struct();
for i = 1:3
    k = find(rn == term{i}, 1);
    F = a.F(k); df1 = a.DF(k); df2 = a.DF(k+1);
    etap = a.SumSq(k) / (a.SumSq(k) + a.SumSq(k+1));
    p = 1 - fcdf(F, df1, df2);                  % df1 = 1 -> no GG
    fprintf('   %-12s F(%d,%d) = %.2f, %s, eta_p^2 = %.3f\n', lab{i}, df1, df2, F, pstr(p), etap);
    r.(lab{i}) = struct('F',F,'df1',df1,'df2',df2,'p',p,'etap',etap);
end
end
