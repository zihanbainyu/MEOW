clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single-subject behavioral analysis
%
% the group script (S_behavioral_analysis_fixed.m) treats each subject as one
% observation and tests across subjects (paired t, rm-anova, correlation).
% with n = 1 none of that is defined, so here the unit of observation is the
% TRIAL:
%   - proportions come with wilson 95% CIs and explicit trial counts
%   - contrasts are two-proportion z / chi-square over disjoint trial sets
%   - RT contrasts are rank-sum / kruskal-wallis (within-subject RTs are skewed)
%   - d' and LDI get parametric-bootstrap CIs
% counts are printed next to every estimate so thin cells are visible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%
subj_id  = 601;
base_dir = '..';
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');
if ~exist(res_dir,'dir'), mkdir(res_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
min_rt   = 0.150;
n_boot   = 5000;
rng(1);   % reproducible bootstraps

c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255;
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;

subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', subj_id));
load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', subj_id)), 'final_data_output');
r1 = final_data_output.results_1_back_all;
r2 = final_data_output.results_2_back_all;
has_rec = isfield(final_data_output, 'results_recognition');
if has_rec, rec = final_data_output.results_recognition; end

diary_file = fullfile(res_dir, sprintf('sub%03d_behavioral_report.txt', subj_id));
if exist(diary_file,'file'), delete(diary_file); end
diary(diary_file);
fprintf('================================================================\n');
fprintf(' SINGLE-SUBJECT BEHAVIORAL REPORT  --  sub%03d   (%s)\n', subj_id, datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('================================================================\n');

%% recode responses
r1.resp_key = cellstr(r1.resp_key); r1.resp_key(strcmp(r1.resp_key,'NA')) = {'none'};
r1.correct  = strcmp(cellstr(r1.corr_resp), r1.resp_key);
r2.resp_key = cellstr(r2.resp_key); r2.resp_key(strcmp(r2.resp_key,'NA')) = {'none'};
r2.correct  = strcmp(cellstr(r2.corr_resp), r2.resp_key);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. data quality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n=== 0. DATA QUALITY ==============================================\n');
blocks = unique([r1.block; r2.block])';
qc = struct();
for tsk = 1:2
    if tsk==1, T = r1; nm = '1-back'; else, T = r2; nm = '2-back'; end
    miss = strcmp(T.resp_key,'none') & ~strcmp(cellstr(T.corr_resp),'none');   % should have responded, did not
    fprintf('\n%s: %d trials, blocks [%s]\n', nm, height(T), num2str(unique(T.block)'));
    fprintf('  %-8s %7s %10s %12s %10s %10s\n','block','n','no-resp','fix broken','gate t/o','RT<min');
    for b = blocks
        m = T.block==b; if ~any(m), continue; end
        fprintf('  %-8d %7d %9.1f%% %11.1f%% %9.1f%% %9.1f%%\n', b, sum(m), ...
            100*mean(miss(m)), 100*mean(T.fix_broken(m)), 100*mean(T.gate_timeout(m)), ...
            100*mean(T.rt(m) < min_rt & ~isnan(T.rt(m))));
    end
    fprintf('  %-8s %7d %9.1f%% %11.1f%% %9.1f%% %9.1f%%\n','ALL', height(T), ...
        100*mean(miss), 100*mean(T.fix_broken), 100*mean(T.gate_timeout), ...
        100*mean(T.rt < min_rt & ~isnan(T.rt)));
    qc.(sprintf('task%d',tsk)) = struct('n',height(T),'miss',mean(miss), ...
        'fix_broken',mean(T.fix_broken),'gate_timeout',mean(T.gate_timeout));
end
fprintf('\nRT (>%.0f ms): 1-back median %.3f s [IQR %.3f-%.3f] | 2-back median %.3f s [IQR %.3f-%.3f]\n', ...
    1000*min_rt, median(r1.rt(r1.rt>min_rt),'omitnan'), quantile(r1.rt(r1.rt>min_rt),0.25), quantile(r1.rt(r1.rt>min_rt),0.75), ...
    median(r2.rt(r2.rt>min_rt),'omitnan'), quantile(r2.rt(r2.rt>min_rt),0.25), quantile(r2.rt(r2.rt>min_rt),0.75));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. encoding (1-back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1     = r1.rt > min_rt;
i1_sim = r1.condition=="compared" & r1.identity=="B";                       % expect 'k' (similar)
i1_sam = r1.condition=="repeat" & strcmp(r1.corr_resp,'j');                 % expect 'j' (same)
i1_new = (r1.condition=="compared" & r1.identity=="A") | r1.condition=="isolated" | ...
         (r1.condition=="repeat" & strcmp(r1.corr_resp,'none'));            % expect 'none'

fprintf('\n\n=== 1. ENCODING (1-back) =========================================\n');
fprintf('\n-- accuracy (wilson 95%% CI) --\n');
one = struct();
[one.acc_same, one.n_same] = prop_report('same    (repeat 2nd) ', r1.correct, i1_sam);
[one.acc_sim,  one.n_sim ] = prop_report('similar (compared B) ', r1.correct, i1_sim);
[one.acc_new,  one.n_new ] = prop_report('new     (all others) ', r1.correct, i1_new);

fprintf('\n-- accuracy contrasts (two-proportion z; disjoint trial sets) --\n');
prop_test('same vs similar', r1.correct(i1_sam), r1.correct(i1_sim));
prop_test('similar vs new ', r1.correct(i1_sim), r1.correct(i1_new));
prop_test('same vs new    ', r1.correct(i1_sam), r1.correct(i1_new));

fprintf('\n-- RT on correct trials --\n');
rt_sam = r1.rt(i1_sam & r1.correct & v1);
rt_sim = r1.rt(i1_sim & r1.correct & v1);
rt_report('same   ', rt_sam); rt_report('similar', rt_sim);
rank_test('same vs similar', rt_sam, rt_sim);
one.rt_same = median(rt_sam,'omitnan'); one.rt_sim = median(rt_sim,'omitnan');

fprintf('\n-- confusions (row = presented, col = response; proportion [n]) --\n');
[m1, c1] = conf_matrix(r1.resp_key, {i1_sam, i1_sim, i1_new}, {'j','k','none'});
print_conf(m1, c1, {'same','similar','new'}, {'->same','->similar','->new'});

fprintf('\n-- accuracy by block --\n');
acc1_blk = block_acc(r1, blocks, {i1_sam, i1_sim, i1_new}, {'same','similar','new'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. retrieval (2-back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an A-N trial is the 2nd item after an 'A-N' goal item. blocks are concatenated
% here, so the +2 pointer is additionally required to stay inside its own block.
real2 = ~contains(r2.goal, "JUNK");
v2    = r2.rt > min_rt;
pan   = false(height(r2),1);
for i = 1:height(r2)-2
    if strcmp(r2.goal(i),'A-N') && r2.block(i+2)==r2.block(i), pan(i+2) = true; end
end
aa = real2 & strcmp(r2.goal,'A-A') & strcmp(r2.corr_resp,'j');
ab = real2 & strcmp(r2.goal,'A-B') & strcmp(r2.corr_resp,'k');
an = real2 & pan & strcmp(r2.corr_resp,'none');
cnd    = {real2 & strcmp(r2.condition,'compared'), real2 & strcmp(r2.condition,'isolated'), real2 & strcmp(r2.condition,'novel')};
cnd_nm = {'compared','isolated','novel'};
goal_m = {aa, ab, an}; goal_nm = {'AA (same / target)','AB (similar / lure)','AN (new / foil)'};

fprintf('\n\n=== 2. RETRIEVAL (2-back) ========================================\n');
fprintf('\n-- accuracy by goal x condition (wilson 95%% CI) --\n');
acc2 = nan(3,3); ci2 = nan(3,3,2); n2 = zeros(3,3);
for g = 1:3
    fprintf('\n  %s\n', goal_nm{g});
    for c = 1:3
        m = goal_m{g} & cnd{c};
        [acc2(g,c), n2(g,c), lo, hi] = prop_report(sprintf('    %-9s        ', cnd_nm{c}), r2.correct, m);
        ci2(g,c,:) = [lo hi];
    end
    fprintf('    -> across conditions: %s\n', chi2_str(r2.correct, {goal_m{g}&cnd{1}, goal_m{g}&cnd{2}, goal_m{g}&cnd{3}}));
end

fprintf('\n-- LDI  ( p(k|AB) - p(k|AN), bootstrap 95%% CI ) --\n');
ldi = nan(1,3); ldi_ci = nan(3,2);
for c = 1:3
    nh = sum(ab & cnd{c}); nf = sum(an & cnd{c});
    h  = mean(strcmp(r2.resp_key(ab & cnd{c}),'k'));
    f  = mean(strcmp(r2.resp_key(an & cnd{c}),'k'));
    [ldi(c), ldi_ci(c,:)] = boot_index(@(a,b) a-b, h, nh, f, nf, n_boot);
    fprintf('  %-9s LDI = %+.3f  [%+.3f, %+.3f]   (lure-k %.3f of %d, foil-k %.3f of %d)\n', ...
        cnd_nm{c}, ldi(c), ldi_ci(c,1), ldi_ci(c,2), h, nh, f, nf);
end
fprintf('  pairwise differences (bootstrap):\n');
boot_pair_diff(r2, ab, an, 'k', @(a,b) a-b, cnd, cnd_nm, n_boot);

fprintf('\n-- d''  ( z(p(j|AA)) - z(p(j|AN)), bootstrap 95%% CI ) --\n');
dp = nan(1,3); dp_ci = nan(3,2);
for c = 1:3
    nh = sum(aa & cnd{c}); nf = sum(an & cnd{c});
    h  = mean(strcmp(r2.resp_key(aa & cnd{c}),'j'));
    f  = mean(strcmp(r2.resp_key(an & cnd{c}),'j'));
    [dp(c), dp_ci(c,:)] = boot_index(@(a,b) calc_d(a,b,nh,nf), h, nh, f, nf, n_boot);
    fprintf('  %-9s d'' = %.3f  [%.3f, %.3f]   (hit %.3f of %d, FA %.3f of %d)\n', ...
        cnd_nm{c}, dp(c), dp_ci(c,1), dp_ci(c,2), h, nh, f, nf);
end
fprintf('  pairwise differences (bootstrap):\n');
boot_pair_diff_d(r2, aa, an, cnd, cnd_nm, n_boot);

fprintf('\n-- RT on correct trials --\n');
rt_l = cell(1,3); rt_t = cell(1,3);
for c = 1:3
    rt_l{c} = r2.rt(ab & cnd{c} & r2.correct & v2);
    rt_t{c} = r2.rt(aa & cnd{c} & r2.correct & v2);
end
fprintf('  lure (AB):\n');   for c=1:3, rt_report(sprintf('  %-9s', cnd_nm{c}), rt_l{c}); end
fprintf('    -> %s\n', kw_str(rt_l));
fprintf('  target (AA):\n'); for c=1:3, rt_report(sprintf('  %-9s', cnd_nm{c}), rt_t{c}); end
fprintf('    -> %s\n', kw_str(rt_t));

fprintf('\n-- confusions per condition (row = presented, col = response) --\n');
m2 = cell(1,3); cc2 = cell(1,3);
for c = 1:3
    fprintf('\n  %s:\n', cnd_nm{c});
    [m2{c}, cc2{c}] = conf_matrix(r2.resp_key, {aa&cnd{c}, ab&cnd{c}, an&cnd{c}}, {'j','k','none'});
    print_conf(m2{c}, cc2{c}, {'same','similar','new'}, {'->same','->similar','->new'});
end

fprintf('\n-- accuracy by block --\n');
acc2_blk = block_acc(r2, blocks, {aa, ab, an}, {'AA same','AB similar','AN new'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. post-task recognition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rec_cond = {}; rec_d = []; rec_dci = []; rec_hit = []; rec_n = []; far = NaN;
fprintf('\n\n=== 3. POST-TASK RECOGNITION =====================================\n');
if has_rec
    rec.resp_key = cellstr(rec.resp_key); rec.resp_key(strcmp(rec.resp_key,'NA')) = {'none'};
    rec.correct  = strcmp(cellstr(rec.corr_resp), rec.resp_key);
    is_old = rec.trial_type=="old"; is_new = ~is_old;
    n_new  = sum(is_new); n_fa = sum(strcmp(rec.resp_key(is_new),'j') & rec.rt(is_new)>min_rt);
    far    = n_fa/n_new;
    fprintf('\n  foils: %d trials, false-alarm rate %.3f (%d "old" responses)\n', n_new, far, n_fa);
    fprintf('  no response: %.1f%% of %d recognition trials\n', 100*mean(strcmp(rec.resp_key,'none')), height(rec));
    fprintf('\n  %-10s %6s %22s %26s\n','condition','n','hit rate [95% CI]','d'' [bootstrap 95% CI]');
    oc = unique(cellstr(rec.condition(is_old)))';
    for k = 1:numel(oc)
        m  = is_old & strcmp(cellstr(rec.condition), oc{k});
        nh = sum(m); h = sum(rec.correct(m) & rec.rt(m)>min_rt)/nh;
        [lo, hi]  = wilson(h, nh);
        [d, dci]  = boot_index(@(a,b) calc_d(a,b,nh,n_new), h, nh, far, n_new, n_boot);
        fprintf('  %-10s %6d %10.3f [%.3f, %.3f] %12.3f [%.3f, %.3f]\n', oc{k}, nh, h, lo, hi, d, dci(1), dci(2));
        rec_cond{end+1} = oc{k}; rec_hit(end+1) = h; rec_n(end+1) = nh; %#ok<SAGROW>
        rec_d(end+1) = d; rec_dci(end+1,:) = dci; %#ok<SAGROW>
    end
    fprintf('\n  hit-rate contrasts:\n');
    for a = 1:numel(oc)-1
        for b = a+1:numel(oc)
            ma = is_old & strcmp(cellstr(rec.condition), oc{a});
            mb = is_old & strcmp(cellstr(rec.condition), oc{b});
            prop_test(sprintf('  %-9s vs %-9s', oc{a}, oc{b}), rec.correct(ma), rec.correct(mb));
        end
    end
else
    fprintf('\n  no recognition data for this subject.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indiv = struct('subj_id',subj_id,'qc',qc,'one',one, ...
    'acc2',acc2,'acc2_ci',ci2,'acc2_n',n2, ...
    'ldi',ldi,'ldi_ci',ldi_ci,'dprime',dp,'dprime_ci',dp_ci, ...
    'conf_1back',m1,'conf_1back_n',c1,'conf_2back',{m2},'conf_2back_n',{cc2}, ...
    'acc1_blk',acc1_blk,'acc2_blk',acc2_blk, ...
    'rec_cond',{rec_cond},'rec_hit',rec_hit,'rec_n',rec_n,'rec_far',far, ...
    'rec_d',rec_d,'rec_d_ci',rec_dci);
save(fullfile(res_dir, sprintf('sub%03d_indiv_stats.mat', subj_id)), 'indiv');
fprintf('\n\nsaved: %s\n', fullfile(res_dir, sprintf('sub%03d_indiv_stats.mat', subj_id)));
fprintf('saved: %s\n', diary_file);
diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 1: data quality / time on task
f = figure('color','w','Position',[60 60 1250 800],'Name','QC');
subplot(2,2,1);
qmat = [arrayfun(@(b) mean(r1.fix_broken(r1.block==b)), blocks); ...
        arrayfun(@(b) mean(r1.gate_timeout(r1.block==b)), blocks); ...
        arrayfun(@(b) mean(r2.fix_broken(r2.block==b)), blocks); ...
        arrayfun(@(b) mean(r2.gate_timeout(r2.block==b)), blocks)]';
bar(blocks, 100*qmat); ylabel('% of trials'); xlabel('block');
title('fixation / gate failures','FontSize',13); set(gca,'XTick',blocks);
legend({'1b fix broken','1b gate t/o','2b fix broken','2b gate t/o'},'Location','best','Box','off'); box off;

subplot(2,2,2); hold on;
histogram(r1.rt(r1.rt>min_rt), 30, 'FaceColor', c_same, 'FaceAlpha', .55, 'EdgeColor','none','Normalization','pdf');
histogram(r2.rt(r2.rt>min_rt), 30, 'FaceColor', c_comp, 'FaceAlpha', .55, 'EdgeColor','none','Normalization','pdf');
xline(median(r1.rt(r1.rt>min_rt),'omitnan'),'--','Color',c_same,'LineWidth',2);
xline(median(r2.rt(r2.rt>min_rt),'omitnan'),'--','Color',c_comp,'LineWidth',2);
xlabel('RT (s)'); ylabel('density'); title('RT distribution','FontSize',13);
legend({'1-back','2-back'},'Box','off'); box off; hold off;

subplot(2,2,3); hold on;
plot(blocks, 100*acc1_blk', '-o','LineWidth',2,'MarkerFaceColor','w');
yline(100/3,'k:','chance'); ylim([0 105]); xlabel('block'); ylabel('accuracy (%)');
title('1-back accuracy over blocks','FontSize',13);
legend({'same','similar','new'},'Location','southwest','Box','off');
set(gca,'XTick',blocks); box off; hold off;

subplot(2,2,4); hold on;
plot(blocks, 100*acc2_blk', '-o','LineWidth',2,'MarkerFaceColor','w');
yline(100/3,'k:','chance'); ylim([0 105]); xlabel('block'); ylabel('accuracy (%)');
title('2-back accuracy over blocks','FontSize',13);
legend({'AA same','AB similar','AN new'},'Location','southwest','Box','off');
set(gca,'XTick',blocks); box off; hold off;
save_fig(f, fig_dir, sprintf('sub%03d_qc', subj_id));

%% fig 2: 1-back
f = figure('color','w','Position',[60 60 1350 420],'Name','1-back');
subplot(1,3,1);
bar_ci([one.acc_same one.acc_sim one.acc_new], [one.n_same one.n_sim one.n_new], ...
    {c_same,c_sim,c_new}, {'same','similar','new'}, 'accuracy', '1-back accuracy');
ylim([0 1.1]); yline(1/3,'k:','chance','LineWidth',1.2);
subplot(1,3,2);
raincloud_cell({rt_sam, rt_sim}, {c_same,c_sim}, {'same','similar'}, 'RT (s)', '1-back RT (correct)');
subplot(1,3,3);
draw_matrix(m1, c1, {c_same,c_sim,c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
title('1-back confusions','FontSize',13);
save_fig(f, fig_dir, sprintf('sub%03d_1back', subj_id));

%% fig 3: 2-back
f = figure('color','w','Position',[40 40 1350 900],'Name','2-back');
subplot(2,2,1);
grouped_bar_ci(acc2, ci2, n2, {c_comp,c_iso,c_nov}, {'AA same','AB similar','AN new'}, cnd_nm, 'accuracy', '2-back accuracy');
yline(1/3,'k:','chance','LineWidth',1.2);
subplot(2,2,2);
bar_ci_explicit(ldi, ldi_ci, {c_comp,c_iso,c_nov}, cnd_nm, 'LDI', 'lure discrimination');
yline(0,'k-');
subplot(2,2,3);
bar_ci_explicit(dp, dp_ci, {c_comp,c_iso,c_nov}, cnd_nm, 'd''', 'target detection');
yline(0,'k-');
subplot(2,2,4);
raincloud_cell([rt_l, rt_t], {c_comp,c_iso,c_nov,c_comp,c_iso,c_nov}, ...
    {'AB comp','AB iso','AB nov','AA comp','AA iso','AA nov'}, 'RT (s)', '2-back RT (correct)');
xtickangle(30);
save_fig(f, fig_dir, sprintf('sub%03d_2back', subj_id));

%% fig 4: confusion matrices
f = figure('color','w','Position',[40 40 1400 400],'Name','confusions');
for c = 1:3
    subplot(1,3,c);
    draw_matrix(m2{c}, cc2{c}, {c_same,c_sim,c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
    title(sprintf('2-back: %s', cnd_nm{c}), 'FontSize', 13);
end
save_fig(f, fig_dir, sprintf('sub%03d_confusion', subj_id));

%% fig 5: recognition
if has_rec && ~isempty(rec_d)
    f = figure('color','w','Position',[100 100 950 420],'Name','recognition');
    cmap = containers.Map({'compared','isolated','novel'}, {c_comp, c_iso, c_nov});
    cols = cell(1,numel(rec_cond));
    for k = 1:numel(rec_cond)
        if isKey(cmap, rec_cond{k}), cols{k} = cmap(rec_cond{k}); else, cols{k} = c_new; end
    end
    subplot(1,2,1);
    bar_ci(rec_hit, rec_n, cols, rec_cond, 'hit rate', 'recognition hit rate');
    ylim([0 1.1]); yline(far,'r--',sprintf('FA = %.2f',far),'LineWidth',1.5);
    subplot(1,2,2);
    bar_ci_explicit(rec_d, rec_dci, cols, rec_cond, 'd''', 'recognition d''');
    yline(0,'r--','chance','LineWidth',1.5);
    save_fig(f, fig_dir, sprintf('sub%03d_recognition', subj_id));
end

fprintf('\nfigures written to %s\n', fig_dir);

%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%
function [lo, hi] = wilson(p, n)
% wilson score interval: stays inside [0 1] at ceiling/floor and small n,
% unlike the normal approximation
    if n == 0 || isnan(p), lo = NaN; hi = NaN; return; end
    z = 1.959963985; d = 1 + z^2/n; ctr = (p + z^2/(2*n))/d;
    hw = z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d;
    lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end

function [p, n, lo, hi] = prop_report(lbl, correct, mask)
    n = sum(mask);
    if n == 0, p = NaN; lo = NaN; hi = NaN; fprintf('  %s  (no trials)\n', lbl); return; end
    p = mean(correct(mask)); [lo, hi] = wilson(p, n);
    fprintf('  %s  %.3f  [%.3f, %.3f]   %d/%d\n', lbl, p, lo, hi, sum(correct(mask)), n);
end

function prop_test(lbl, a, b)
% two-proportion z-test; the two trial sets are disjoint so independence holds
    n1 = numel(a); n2 = numel(b);
    if n1==0 || n2==0, fprintf('  %s: insufficient trials\n', lbl); return; end
    p1 = mean(a); p2 = mean(b); pp = (sum(a)+sum(b))/(n1+n2);
    se = sqrt(pp*(1-pp)*(1/n1+1/n2));
    if se == 0, fprintf('  %s: identical proportions (%.3f)\n', lbl, p1); return; end
    z = (p1-p2)/se; p = 2*(1-normcdf(abs(z)));
    h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2));   % cohen's h
    fprintf('  %s: %.3f vs %.3f, diff = %+.3f, z = %.2f, h = %.2f, p = %.4f %s\n', ...
        lbl, p1, p2, p1-p2, z, h, p, stars(p));
end

function rt_report(lbl, x)
    x = x(~isnan(x));
    if isempty(x), fprintf('  %s  (no trials)\n', lbl); return; end
    fprintf('  %s  median %.3f s [IQR %.3f-%.3f], mean %.3f, SD %.3f, n = %d\n', ...
        lbl, median(x), quantile(x,0.25), quantile(x,0.75), mean(x), std(x), numel(x));
end

function rank_test(lbl, a, b)
    a = a(~isnan(a)); b = b(~isnan(b));
    if numel(a)<3 || numel(b)<3, fprintf('  %s: too few trials\n', lbl); return; end
    [p, ~, st] = ranksum(a, b);
    if isfield(st,'zval'), z = st.zval; else, z = NaN; end
    r = abs(z)/sqrt(numel(a)+numel(b));   % z-based effect size
    fprintf('  %s: rank-sum z = %.2f, r = %.2f, p = %.4f %s\n', lbl, z, r, p, stars(p));
end

function s = chi2_str(correct, masks)
% correct/incorrect x group contingency test
    obs = zeros(numel(masks),2);
    for k = 1:numel(masks)
        obs(k,:) = [sum(correct(masks{k})), sum(masks{k})-sum(correct(masks{k}))];
    end
    obs(sum(obs,2)==0,:) = [];
    if size(obs,1) < 2 || any(sum(obs,1)==0), s = 'chi2 not computable'; return; end
    exp_ = sum(obs,2)*sum(obs,1)/sum(obs(:));
    x2 = sum(sum((obs-exp_).^2 ./ exp_));
    df = (size(obs,1)-1)*(size(obs,2)-1); p = 1-chi2cdf(x2, df);
    v  = sqrt(x2/(sum(obs(:))*min(size(obs)-1)));   % cramer's V
    s  = sprintf('chi2(%d) = %.2f, V = %.2f, p = %.4f %s', df, x2, v, p, stars(p));
    if any(exp_(:) < 5), s = [s ' (expected count < 5, treat with care)']; end
end

function s = kw_str(cells)
    g = []; x = [];
    for k = 1:numel(cells)
        v = cells{k}(~isnan(cells{k})); x = [x; v(:)]; g = [g; k*ones(numel(v),1)]; %#ok<AGROW>
    end
    if numel(unique(g)) < 2 || numel(x) < 6, s = 'kruskal-wallis not computable'; return; end
    [p, tbl] = kruskalwallis(x, g, 'off');
    s = sprintf('kruskal-wallis H(%d) = %.2f, p = %.4f %s', tbl{2,3}, tbl{2,5}, p, stars(p));
end

function d = calc_d(h, f, nh, nf)
% 1/(2N) clamp keeps d' finite at ceiling / floor
    d = norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));
end

function [est, ci, bs] = boot_index(fn, h, nh, f, nf, n_boot)
% parametric bootstrap: resample hit and FA counts from their binomials.
% with one subject the sampling variability that matters is over trials.
    if nh==0 || nf==0 || isnan(h) || isnan(f), est = NaN; ci = [NaN NaN]; bs = nan(n_boot,1); return; end
    est = fn(h, f);
    bh  = binornd(nh, h, n_boot, 1)/nh;
    bf  = binornd(nf, f, n_boot, 1)/nf;
    bs  = arrayfun(fn, bh, bf);
    ci  = quantile(bs, [0.025 0.975]);
end

function boot_pair_diff(r2, hit_mask, fa_mask, key, fn, cnd, cnd_nm, n_boot)
    bs = cell(1,3); est = nan(1,3);
    for c = 1:3
        nh = sum(hit_mask & cnd{c}); nf = sum(fa_mask & cnd{c});
        h  = mean(strcmp(r2.resp_key(hit_mask & cnd{c}), key));
        f  = mean(strcmp(r2.resp_key(fa_mask  & cnd{c}), key));
        [est(c), ~, bs{c}] = boot_index(fn, h, nh, f, nf, n_boot);
    end
    report_pairs(est, bs, cnd_nm);
end

function boot_pair_diff_d(r2, aa, an, cnd, cnd_nm, n_boot)
    bs = cell(1,3); est = nan(1,3);
    for c = 1:3
        nh = sum(aa & cnd{c}); nf = sum(an & cnd{c});
        h  = mean(strcmp(r2.resp_key(aa & cnd{c}),'j'));
        f  = mean(strcmp(r2.resp_key(an & cnd{c}),'j'));
        [est(c), ~, bs{c}] = boot_index(@(x,y) calc_d(x,y,nh,nf), h, nh, f, nf, n_boot);
    end
    report_pairs(est, bs, cnd_nm);
end

function report_pairs(est, bs, cnd_nm)
    pr = [1 2; 2 3; 1 3];
    for k = 1:size(pr,1)
        a = pr(k,1); b = pr(k,2);
        if isnan(est(a)) || isnan(est(b)), continue; end
        d  = bs{a} - bs{b}; ci = quantile(d, [0.025 0.975]);
        p  = 2*min(mean(d<=0), mean(d>=0));
        fprintf('    %-9s - %-9s = %+.3f  [%+.3f, %+.3f], p_boot = %.4f %s\n', ...
            cnd_nm{a}, cnd_nm{b}, est(a)-est(b), ci(1), ci(2), p, stars(p));
    end
end

function [mat, cnt] = conf_matrix(resp, masks, keys)
    nk = numel(keys); mat = nan(numel(masks), nk); cnt = zeros(numel(masks), nk);
    for r = 1:numel(masks)
        n = sum(masks{r}); if n==0, continue; end
        for c = 1:nk
            cnt(r,c) = sum(strcmp(resp(masks{r}), keys{c}));
            mat(r,c) = cnt(r,c)/n;
        end
    end
end

function print_conf(mat, cnt, rlbl, clbl)
    fprintf('  %-12s', ''); fprintf('%16s', clbl{:}); fprintf('%8s\n','n');
    for r = 1:size(mat,1)
        fprintf('  %-12s', rlbl{r});
        for c = 1:size(mat,2), fprintf('%10.3f [%3d]', mat(r,c), cnt(r,c)); end
        fprintf('%8d\n', sum(cnt(r,:)));
    end
end

function acc = block_acc(T, blocks, masks, lbls)
    acc = nan(numel(masks), numel(blocks));
    fprintf('  %-12s', 'block'); fprintf('%12d', blocks); fprintf('\n');
    for r = 1:numel(masks)
        for b = 1:numel(blocks)
            m = masks{r} & T.block==blocks(b);
            if any(m), acc(r,b) = mean(T.correct(m)); end
        end
        fprintf('  %-12s', lbls{r}); fprintf('%12.3f', acc(r,:)); fprintf('\n');
    end
end

function s = stars(p)
    s = repmat('*', 1, (p<0.05)+(p<0.01)+(p<0.001));
end

function bar_ci(vals, ns, cols, xlbls, ylbl, ttl)
    lo = nan(size(vals)); hi = nan(size(vals));
    for i = 1:numel(vals), [lo(i), hi(i)] = wilson(vals(i), ns(i)); end
    bar_ci_explicit(vals, [lo(:) hi(:)], cols, xlbls, ylbl, ttl);
    for i = 1:numel(vals)
        text(i, 0.03, sprintf('n=%d', ns(i)), 'HorizontalAlignment','center','FontSize',9,'Color',[.3 .3 .3]);
    end
end

function bar_ci_explicit(vals, ci, cols, xlbls, ylbl, ttl)
    hold on;
    for i = 1:numel(vals)
        bar(i, vals(i), 0.65, 'FaceColor', cols{i}, 'EdgeColor', [.3 .3 .3], 'LineWidth', 1);
        if all(~isnan(ci(i,:)))
            plot([i i], ci(i,:), 'k-', 'LineWidth', 1.5);
            plot([i-0.1 i+0.1], [ci(i,1) ci(i,1)], 'k-', 'LineWidth', 1.5);
            plot([i-0.1 i+0.1], [ci(i,2) ci(i,2)], 'k-', 'LineWidth', 1.5);
        end
        top = max([vals(i), ci(i,2)]); if isnan(top), top = vals(i); end
        text(i, top, sprintf('%.2f', vals(i)), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',10);
    end
    set(gca,'XTick',1:numel(vals),'XTickLabel',xlbls,'FontSize',12);
    xlim([0.4 numel(vals)+0.6]); ylabel(ylbl,'FontSize',13);
    if ~isempty(ttl), title(ttl,'FontSize',13); end
    box off; hold off;
end

function raincloud_cell(cells, cols, xlbls, ylbl, ttl)
% single-subject variant of the group raincloud: each cloud is built from that
% cell's TRIALS (unequal n), so there are no across-subject connecting lines
    hold on; n_grps = numel(cells); all_v = [];
    for i = 1:n_grps
        d = cells{i}(:); d = d(~isnan(d)); all_v = [all_v; d]; %#ok<AGROW>
        if isempty(d), continue; end
        if numel(d) > 3
            [fd, xi] = ksdensity(d); fd = fd/max(fd)*0.35;
            patch([i+fd, i*ones(1,numel(fd))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',0.45);
        end
        jit = -0.12 - rand(numel(d),1)*0.22;
        scatter(i+jit, d, 12, cols{i}, 'filled', 'MarkerFaceAlpha', 0.45);
        q = quantile(d, [0.25 0.5 0.75]);
        rectangle('Position', [i-0.05, q(1), 0.10, q(3)-q(1)], 'FaceColor', cols{i}, 'EdgeColor','k','LineWidth',1.1);
        plot([i-0.05 i+0.05], [q(2) q(2)], 'k-', 'LineWidth', 2);
        text(i+0.40, q(2), sprintf('n=%d', numel(d)), 'FontSize', 8, 'Color', [.3 .3 .3]);
    end
    set(gca,'XTick',1:n_grps,'XTickLabel',xlbls,'FontSize',12);
    xlim([0.4 n_grps+0.7]); ylabel(ylbl,'FontSize',13);
    if ~isempty(ttl), title(ttl,'FontSize',13); end
    if ~isempty(all_v)
        rg = range(all_v); if rg==0, rg = 1; end
        ylim([min(all_v)-0.1*rg, max(all_v)+0.1*rg]);
    end
    box off; hold off;
end

function grouped_bar_ci(vals, ci, ns, cols, glbls, clbls, ylbl, ttl)
% vals / ns: [n_goal x n_cond]; ci: [n_goal x n_cond x 2]
    hold on; [ng, nc] = size(vals); w = 0.8/nc; hs = gobjects(1,nc);
    for g = 1:ng
        for c = 1:nc
            x = g - 0.4 + w*(c-0.5);
            hb = bar(x, vals(g,c), w*0.9, 'FaceColor', cols{c}, 'EdgeColor',[.3 .3 .3], 'LineWidth', 0.8);
            if g == 1, hs(c) = hb; end
            lohi = squeeze(ci(g,c,:));
            if all(~isnan(lohi)), plot([x x], lohi, 'k-', 'LineWidth', 1.2); end
            text(x, 0.03, sprintf('%d', ns(g,c)), 'HorizontalAlignment','center','FontSize',8,'Color',[.3 .3 .3]);
        end
    end
    set(gca,'XTick',1:ng,'XTickLabel',glbls,'FontSize',12);
    xlim([0.4 ng+0.6]); ylim([0 1.1]); ylabel(ylbl,'FontSize',13);
    if ~isempty(ttl), title(ttl,'FontSize',13); end
    legend(hs, clbls, 'Location','southoutside','Orientation','horizontal', ...
        'Box','off','AutoUpdate','off');   % AutoUpdate off: the chance line must stay out
    box off; hold off;
end

function draw_matrix(mat, cnt, cols, ylbl, xlbl)
    hold on;
    for r = 1:3
        for c = 1:3
            v = mat(r,c); if isnan(v), v = 0; end
            t_c = (v > 0.5)*[1 1 1] + (v <= 0.5)*[0.2 0.2 0.2];
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], cols{r}, 'EdgeColor','none','FaceAlpha',v);
            text(c, r-0.08, sprintf('%.2f', v), 'HorizontalAlignment','center','Color',t_c,'FontWeight','bold','FontSize',10);
            text(c, r+0.18, sprintf('[%d]', cnt(r,c)), 'HorizontalAlignment','center','Color',t_c,'FontSize',8);
            if r == 1, text(c, 0.4, xlbl{c}, 'HorizontalAlignment','center','Color',cols{c},'FontWeight','bold','FontSize',10); end
        end
        text(0.4, r, ylbl{r}, 'HorizontalAlignment','right','Color',cols{r},'FontWeight','bold','FontSize',10);
    end
    axis ij equal; xlim([0.2 3.8]); ylim([0.3 3.5]);
    set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none'); hold off;
end

function save_fig(f, fig_dir, name)
    set(f, 'Renderer', 'painters', 'PaperPositionMode', 'auto');
    pos = get(f, 'Position');                       % match page to figure, else pdf clips
    set(f, 'PaperUnits', 'points', 'PaperSize', pos(3:4));
    print(f, fullfile(fig_dir, [name '.pdf']), '-dpdf', '-vector');
    print(f, fullfile(fig_dir, [name '.png']), '-dpng', '-r150');
end
