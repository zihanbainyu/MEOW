function r = rm_3x2(label, cc, ci, ic, ii, nc, ni)
% 3x2 rm-anova (condition[3] x accuracy[2]); gg p for the 2-df effects, partial eta^2
% columns: comp/cor comp/inc iso/cor iso/inc nov/cor nov/inc
M = [cc(:) ci(:) ic(:) ii(:) nc(:) ni(:)];
M = M(all(~isnan(M),2),:);
n = size(M,1);
cc=M(:,1); ci=M(:,2); ic=M(:,3); ii=M(:,4); nc=M(:,5); ni=M(:,6);
t = array2table(M, 'VariableNames', {'cc','ci','ic','ii','nc','ni'});
w = table({'comp';'comp';'iso';'iso';'nov';'nov'}, ...
          {'cor';'inc';'cor';'inc';'cor';'inc'}, 'VariableNames', {'Cond','Acc'});
a = ranova(fitrm(t, 'cc-ni~1', 'WithinDesign', w), 'WithinModel', 'Cond*Acc');
rn = string(a.Properties.RowNames);
cond_eps = [(cc+ci)/2, (ic+ii)/2, (nc+ni)/2];   % condition means, for cond epsilon
int_eps  = [cc-ci, ic-ii, nc-ni];               % correct-incorrect, for interaction epsilon
fprintf('%s  [n = %d]\n', label, n);
r = struct();
r.condition   = term(a, rn, '(Intercept):Cond',     'condition',   cond_eps);
r.accuracy    = term(a, rn, '(Intercept):Acc',      'accuracy',    []);
r.interaction = term(a, rn, '(Intercept):Cond:Acc', 'interaction', int_eps);
end

function s = term(a, rn, name, label, epsData)
k = find(rn == name, 1);
F = a.F(k); df1 = a.DF(k); df2 = a.DF(k+1);
etap = a.SumSq(k) / (a.SumSq(k) + a.SumSq(k+1));
if df1 > 1 && ~isempty(epsData)
    e = gg_eps(epsData);
    p = 1 - fcdf(F, df1*e, df2*e);
    fprintf('   %-12s F(%d,%d) = %.2f, %s (GG, eps = %.3f), eta_p^2 = %.3f\n', label, df1, df2, F, pstr(p), e, etap);
    s = struct('F',F,'df1',df1,'df2',df2,'p',p,'eps',e,'etap',etap);
else
    p = 1 - fcdf(F, df1, df2);
    fprintf('   %-12s F(%d,%d) = %.2f, %s, eta_p^2 = %.3f\n', label, df1, df2, F, pstr(p), etap);
    s = struct('F',F,'df1',df1,'df2',df2,'p',p,'eps',1,'etap',etap);
end
end
