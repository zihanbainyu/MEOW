function r = rm_3x2(label, cc, ci, ic, ii, nc, ni)
% RM_3x2  3x2 repeated-measures ANOVA (condition[3] x correctness[2]).
%   condition & interaction have 2 numerator df -> Greenhouse-Geisser p + eps.
%   accuracy has 1 df -> ordinary p. Effect size = partial eta^2.
%   Column order: comp/cor comp/inc  iso/cor iso/inc  nov/cor nov/inc.
M = [cc(:) ci(:) ic(:) ii(:) nc(:) ni(:)];
M = M(all(~isnan(M),2),:);
n = size(M,1);
cc=M(:,1); ci=M(:,2); ic=M(:,3); ii=M(:,4); nc=M(:,5); ni=M(:,6);
t  = array2table(M, 'VariableNames', {'cc','ci','ic','ii','nc','ni'});
w  = table({'comp';'comp';'iso';'iso';'nov';'nov'}, ...
           {'cor';'inc';'cor';'inc';'cor';'inc'}, 'VariableNames', {'Cond','Acc'});
rm = fitrm(t, 'cc-ni~1', 'WithinDesign', w);
a  = ranova(rm, 'WithinModel', 'Cond*Acc');
rn = string(a.Properties.RowNames);
condData = [(cc+ci)/2, (ic+ii)/2, (nc+ni)/2];   % condition means -> epsilon for cond
intData  = [cc-ci, ic-ii, nc-ni];               % correct-minus-incorrect -> epsilon for interaction
fprintf('%s  [n = %d]\n', label, n);
r = struct();
r.condition   = report_term(a, rn, '(Intercept):Cond',     'condition',   condData);
r.accuracy    = report_term(a, rn, '(Intercept):Acc',      'accuracy',    []);
r.interaction = report_term(a, rn, '(Intercept):Cond:Acc', 'interaction', intData);
end

function s = report_term(a, rn, term, label, epsData)
k = find(rn == term, 1);
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
