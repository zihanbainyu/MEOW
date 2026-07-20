clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%
%% setup (mirrors S_behavioral_analysis.m, BH-FDR replaces Holm-Bonferroni)
%%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631];
addpath(genpath('/Users/bai/Documents/GitHub/MEOW/toolbox/bayesFactor-master'));
base_dir = '..';
res_dir = fullfile(base_dir, 'results');
cohend = @(x,y) mean(x-y,'omitnan') / std(x-y,'omitnan');
load(fullfile(res_dir, 'all_subjs_stats.mat'), 'all_subjs');

get_v = @(f1, f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs);
b1_acc_sam = get_v('one','acc_same')'; b1_acc_sim = get_v('one','acc_sim')'; b1_acc_new = get_v('one','acc_new')';
b1_rt_sam  = get_v('one','rt_same')';  b1_rt_sim  = get_v('one','rt_sim')';
b2_ldi_c   = get_v('two','ldi_comp')'; b2_ldi_i   = get_v('two','ldi_iso')'; b2_ldi_n   = get_v('two','ldi_nov')';
b2_dp_c    = get_v('two','dprime_comp')'; b2_dp_i = get_v('two','dprime_iso')'; b2_dp_n   = get_v('two','dprime_nov')';
b2_rt_l_c  = get_v('two','rt_AB_comp')';  b2_rt_l_i = get_v('two','rt_AB_iso')';  b2_rt_l_n = get_v('two','rt_AB_nov')';
b2_rt_t_c  = get_v('two','rt_AA_comp')';  b2_rt_t_i = get_v('two','rt_AA_iso')';  b2_rt_t_n = get_v('two','rt_AA_nov')';
rec_d_c    = get_v('rec','d_comp'); rec_d_i = get_v('rec','d_iso');

stat_results = struct();
within = table(categorical({'compared';'isolated';'novel'}), 'VariableNames', {'Condition'});

%% 1-back (encoding verification)  — Benjamini-Hochberg FDR @ q=0.05
fprintf('\n--- Encoding (1-back) — paired t-tests, BH-FDR q=0.05 ---\n');
[~,p1,~,s1] = ttest(b1_acc_sam, b1_acc_sim);
[~,p2,~,s2] = ttest(b1_acc_sim, b1_acc_new);
[~,p3,~,s3] = ttest(b1_rt_sam, b1_rt_sim);
[~,p4,~,s4] = ttest(b1_acc_sam, b1_acc_new);
q = bh_fdr([p1 p2 p3 p4]);
fprintf('Acc same v similar: t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s1.df, s1.tstat, cohend(b1_acc_sam,b1_acc_sim), p1, q(1));
fprintf('Acc similar v new:  t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s2.df, s2.tstat, cohend(b1_acc_sim,b1_acc_new), p2, q(2));
fprintf('RT  same v similar: t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s3.df, s3.tstat, cohend(b1_rt_sam,b1_rt_sim),  p3, q(3));
fprintf('Acc same v new:     t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s4.df, s4.tstat, cohend(b1_acc_sam,b1_acc_new), p4, q(4));

%% 2-back RM-ANOVAs (post-hocs use BH-FDR within each family of 3)
fprintf('\n--- Retrieval (2-back) — RM-ANOVA (GG-corrected), post-hocs BH-FDR ---\n');
run_rm_anova_bh('LDI',       b2_ldi_c, b2_ldi_i, b2_ldi_n, within);
run_rm_anova_bh('d-prime',   b2_dp_c,  b2_dp_i,  b2_dp_n,  within);
run_rm_anova_bh('RT lure',   b2_rt_l_c,b2_rt_l_i,b2_rt_l_n,within);
run_rm_anova_bh('RT target', b2_rt_t_c,b2_rt_t_i,b2_rt_t_n,within);

%% recognition
fprintf('\n--- Post-task recognition (paired/one-sample t-tests, BH-FDR across 3) ---\n');
[~,pp1,~,ss1] = ttest(rec_d_c);
[~,pp2,~,ss2] = ttest(rec_d_i);
[~,pp3,~,ss3] = ttest(rec_d_c, rec_d_i);
qq = bh_fdr([pp1 pp2 pp3]);
d1 = mean(rec_d_c,'omitnan')/std(rec_d_c,'omitnan');
d2 = mean(rec_d_i,'omitnan')/std(rec_d_i,'omitnan');
d3 = cohend(rec_d_c', rec_d_i');
fprintf('compared > 0:        t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', ss1.df, ss1.tstat, d1, pp1, qq(1));
fprintf('isolated > 0:        t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', ss2.df, ss2.tstat, d2, pp2, qq(2));
fprintf('compared v isolated: t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', ss3.df, ss3.tstat, d3, pp3, qq(3));

%% WM-EM correlations (one-tailed, as in original)
fprintf('\n--- WM-EM correlations (one-tailed) ---\n');
wm_full = (b2_dp_c + b2_dp_i + b2_dp_n)/3;
em_full = (rec_d_c' + rec_d_i')/2;
em_benefit_full = (b2_ldi_c + b2_ldi_i)/2 - b2_ldi_n;
v1 = ~isnan(wm_full) & ~isnan(em_benefit_full);
[r,p] = corr(wm_full(v1), em_benefit_full(v1), 'Tail','right');
fprintf('WM-d'' x EM-benefit: r=%.3f, p=%.4f, n=%d\n', r, p, sum(v1));
v2 = ~isnan(em_benefit_full) & ~isnan(em_full);
[r,p] = corr(em_benefit_full(v2), em_full(v2), 'Tail','right');
fprintf('EM-benefit x EM-d'':  r=%.3f, p=%.4f, n=%d\n', r, p, sum(v2));
v3 = ~isnan(wm_full) & ~isnan(em_full);
[r,p] = corr(wm_full(v3), em_full(v3), 'Tail','right');
fprintf('WM-d'' x EM-d'':      r=%.3f, p=%.4f, n=%d\n', r, p, sum(v3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = bh_fdr(p)
    % Benjamini-Hochberg FDR adjusted p-values (q-values)
    p = p(:)'; n = numel(p);
    [ps, idx] = sort(p);
    q_sorted = ps .* n ./ (1:n);
    for i = n-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    q_sorted = min(q_sorted, 1);
    q = zeros(size(p));
    q(idx) = q_sorted;
end

function run_rm_anova_bh(lbl, d1, d2, d3, within)
    t = table(d1, d2, d3, 'VariableNames', {'compared','isolated','novel'});
    rm = fitrm(t, 'compared-novel~1', 'WithinDesign', within);
    tbl = ranova(rm); eps_ = epsilon(rm);
    F = tbl.F(1); df1 = tbl.DF(1); df2_gg = tbl.DF(2)*eps_.GreenhouseGeisser;
    p_gg = tbl.pValueGG(1);
    eta2_p = tbl.SumSq(1) / (tbl.SumSq(1) + tbl.SumSq(2));
    fprintf('\n%s: F(%d, %.1f) = %.2f, p = %.4f, eta2_p = %.3f, eps = %.3f\n', ...
        lbl, df1, df2_gg, F, p_gg, eta2_p, eps_.GreenhouseGeisser);

    n = length(d1); subj = (1:n)';
    bf_tbl = table(subj, d1, d2, d3, 'VariableNames', {'subj','c1','c2','c3'});
    bf_long = stack(bf_tbl, {'c1','c2','c3'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj); bf_long.cond = categorical(bf_long.cond);
    bf10 = bf.anova(bf_long, 'y~cond', 'treatAsRandom', {'subj'}, 'verbose', false);
    if bf10 < 1, fprintf('  BF01 = %.2f\n', 1/bf10); else, fprintf('  BF10 = %.2f\n', bf10); end

    cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    [~,p1,~,s1] = ttest(d1,d2);
    [~,p2,~,s2] = ttest(d2,d3);
    [~,p3,~,s3] = ttest(d1,d3); %#ok<ASGLU>
    q = bh_fdr([p1 p2 p3]); df = s1.df;
    fprintf('  post-hoc (BH-FDR q-values):\n');
    fprintf('    comp v iso:  t(%d) = %.2f, d = %.2f, p = %.4f, q = %.4f\n', df, s1.tstat, cd(d1,d2), p1, q(1));
    fprintf('    iso  v nov:  t(%d) = %.2f, d = %.2f, p = %.4f, q = %.4f\n', df, s2.tstat, cd(d2,d3), p2, q(2));
    fprintf('    comp v nov:  t(%d) = %.2f, d = %.2f, p = %.4f, q = %.4f\n', df, s3.tstat, cd(d1,d3), p3, q(3));
end
