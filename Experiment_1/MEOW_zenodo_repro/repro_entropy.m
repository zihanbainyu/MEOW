function repro_entropy()
% REPRO_ENTROPY  Spatial gaze entropy (Supplementary Information).
% 2x2 RM-ANOVA: item position (A, B) x encoding condition (compared, isolated).
% Source data: results/spatial_entropy_results.mat
%   (entropy_[ab]_[comp|iso]_subj, common_subjs_[ab])
% All effects have 1 numerator df, so no Greenhouse-Geisser is needed;
% effect size = partial eta^2. Follow-up t-tests use BH-FDR.

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
S = load(fullfile(DATA,'spatial_entropy_results.mat'));

[~, ia, ib] = intersect(S.common_subjs_a, S.common_subjs_b);
ea_c = S.entropy_a_comp_subj(ia); ea_i = S.entropy_a_iso_subj(ia);
eb_c = S.entropy_b_comp_subj(ib); eb_i = S.entropy_b_iso_subj(ib);
n = numel(ia);

fprintf('\n================ SPATIAL ENTROPY (Supplementary) ================\n');
fprintf('n = %d common subjects\n', n);

t = table(ea_c, ea_i, eb_c, eb_i, 'VariableNames', {'A_comp','A_iso','B_comp','B_iso'});
within = table(categorical({'A';'A';'B';'B'}), categorical({'comp';'iso';'comp';'iso'}), ...
    'VariableNames', {'Item','Condition'});
rm  = fitrm(t, 'A_comp-B_iso ~ 1', 'WithinDesign', within);
tbl = ranova(rm, 'WithinModel', 'Item*Condition');
etap = @(k) tbl.SumSq(k) / (tbl.SumSq(k) + tbl.SumSq(k+1));
fprintf('  Item:           F(%d,%d) = %.2f, %s, eta_p^2 = %.3f\n', tbl.DF(3), tbl.DF(4), tbl.F(3), pstr(tbl.pValue(3)), etap(3));
fprintf('  Condition:      F(%d,%d) = %.2f, %s, eta_p^2 = %.3f\n', tbl.DF(5), tbl.DF(6), tbl.F(5), pstr(tbl.pValue(5)), etap(5));
fprintf('  Item*Condition: F(%d,%d) = %.2f, %s, eta_p^2 = %.3f\n', tbl.DF(7), tbl.DF(8), tbl.F(7), pstr(tbl.pValue(7)), etap(7));

% follow-up paired t-tests (BH-FDR over 4)
[~,p1,~,s1]=ttest(ea_c,ea_i); [~,p2,~,s2]=ttest(eb_c,eb_i);
[~,p3,~,s3]=ttest(ea_c,eb_c); [~,p4,~,s4]=ttest(ea_i,eb_i);
q=bh_fdr([p1 p2 p3 p4]);
fprintf('  follow-ups (BH-FDR):\n');
fprintf('    A comp vs iso: t(%d) = %.2f, d = %.2f, %s, q = %.3f\n', s1.df, s1.tstat, cohend(ea_c,ea_i), pstr(p1), q(1));
fprintf('    B comp vs iso: t(%d) = %.2f, d = %.2f, %s, q = %.3f\n', s2.df, s2.tstat, cohend(eb_c,eb_i), pstr(p2), q(2));
fprintf('    comp A vs B:   t(%d) = %.2f, d = %.2f, %s, q = %.3f\n', s3.df, s3.tstat, cohend(ea_c,eb_c), pstr(p3), q(3));
fprintf('    iso  A vs B:   t(%d) = %.2f, d = %.2f, %s, q = %.3f\n', s4.df, s4.tstat, cohend(ea_i,eb_i), pstr(p4), q(4));

% Bayes factors (optional)
try
    subj=(1:n)'; y=[ea_c;ea_i;eb_c;eb_i]; sj=repmat(subj,4,1);
    itm=categorical([repmat({'A'},2*n,1); repmat({'B'},2*n,1)]);
    cnd=categorical([repmat({'comp'},n,1); repmat({'iso'},n,1); repmat({'comp'},n,1); repmat({'iso'},n,1)]);
    B=table(categorical(sj),itm,cnd,y,'VariableNames',{'subj','Item','Condition','y'});
    bf_full=bf.anova(B,'y ~ Item*Condition','treatAsRandom',{'subj'},'verbose',false);
    bf_add =bf.anova(B,'y ~ Item+Condition','treatAsRandom',{'subj'},'verbose',false);
    fprintf('  Bayes: full BF10 = %.2f ; interaction (full/additive) = %.2f\n', bf_full, bf_full/bf_add);
catch
    fprintf('  [BF skipped: bayesFactor toolbox not on path]\n');
end
end
