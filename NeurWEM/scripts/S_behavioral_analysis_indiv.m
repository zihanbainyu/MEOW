clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%
% NeurWEM (fMRI) single-subject behavioural report. Same analysis as the
% NeurWEM_plt version, with the ISOLATED condition added to the 2-back
% condition set (compared / isolated / novel) and the response keys mapped to
% this task ('1' = same/old, '2' = similar, 'none'/no-press = new/different).
% Section 3 analyses the post-task MST (old / lure / new).
subj_id  = 101;
base_dir = '..';
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');
if ~exist(res_dir,'dir'), mkdir(res_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
min_rt   = 0.150;
n_boot   = 5000;
rng(1);   % reproducible bootstraps

c_comp = [180 174 211]/255; c_iso = [245 205 150]/255; c_nov = [183 210 205]/255;
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255;  c_new = [219 219 219]/255;

% 2-back / MST condition set (order fixes colours + reporting)
cnd_nm   = {'compared','isolated','novel'};
cnd_cols = {c_comp,    c_iso,     c_nov};
nC       = numel(cnd_nm);

subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', subj_id));
load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', subj_id)), 'final_data_output');
r1 = final_data_output.results_1_back_all;
r2 = final_data_output.results_2_back_all;
has_mst = isfield(final_data_output, 'results_mst') && ~isempty(final_data_output.results_mst);
if has_mst, rm = final_data_output.results_mst; end

diary_file = fullfile(res_dir, sprintf('sub%03d_behavioral_report.txt', subj_id));
if exist(diary_file,'file'), delete(diary_file); end
diary(diary_file);
fprintf('================================================================\n');
fprintf(' SINGLE-SUBJECT BEHAVIORAL REPORT  --  sub%03d   (%s)\n', subj_id, datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(' NeurWEM fMRI  (conditions: %s)\n', strjoin(cnd_nm, ' / '));
fprintf(' keys: 1 = same/old,  2 = similar,  none = new/different\n');
fprintf('================================================================\n');

%% recode responses ('NA' -> 'none'; correct = required == pressed)
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
    fprintf('  %-8s %7s %10s %10s\n','block','n','no-resp','RT<min');
    for b = blocks
        m = T.block==b; if ~any(m), continue; end
        fprintf('  %-8d %7d %9.1f%% %9.1f%%\n', b, sum(m), ...
            100*mean(miss(m)), 100*mean(T.rt(m) < min_rt & ~isnan(T.rt(m))));
    end
    fprintf('  %-8s %7d %9.1f%% %9.1f%%\n','ALL', height(T), ...
        100*mean(miss), 100*mean(T.rt < min_rt & ~isnan(T.rt)));
    qc.(sprintf('task%d',tsk)) = struct('n',height(T),'miss',mean(miss));
end
fprintf('\nRT (>%.0f ms): 1-back median %.3f s [IQR %.3f-%.3f] | 2-back median %.3f s [IQR %.3f-%.3f]\n', ...
    1000*min_rt, median(r1.rt(r1.rt>min_rt),'omitnan'), quantile(r1.rt(r1.rt>min_rt),0.25), quantile(r1.rt(r1.rt>min_rt),0.75), ...
    median(r2.rt(r2.rt>min_rt),'omitnan'), quantile(r2.rt(r2.rt>min_rt),0.25), quantile(r2.rt(r2.rt>min_rt),0.75));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. encoding (1-back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% response types keyed off the required response. isolated items (both A and
% B) require no press, so they fall under 'new' together with compared-A and
% the first repeat item.
v1     = r1.rt > min_rt;
i1_sam = strcmp(cellstr(r1.corr_resp),'1');    % same    (repeat 2nd item)
i1_sim = strcmp(cellstr(r1.corr_resp),'2');    % similar (compared B)
i1_new = strcmp(cellstr(r1.corr_resp),'none'); % new     (compared A / repeat 1st / isolated)
i1_iso = strcmp(cellstr(r1.condition),'isolated'); % isolated subset (all 'new')

fprintf('\n\n=== 1. ENCODING (1-back) =========================================\n');
fprintf('\n-- accuracy (wilson 95%% CI) --\n');
one = struct();
[one.acc_same, one.n_same] = prop_report('same    (repeat 2nd) ', r1.correct, i1_sam);
[one.acc_sim,  one.n_sim ] = prop_report('similar (compared B) ', r1.correct, i1_sim);
[one.acc_new,  one.n_new ] = prop_report('new     (all others) ', r1.correct, i1_new);
[one.acc_iso,  one.n_iso ] = prop_report('  of which isolated  ', r1.correct, i1_iso);

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
[m1, c1] = conf_matrix(r1.resp_key, {i1_sam, i1_sim, i1_new}, {'1','2','none'});
print_conf(m1, c1, {'same','similar','new'}, {'->same','->similar','->new'});

fprintf('\n-- accuracy by block --\n');
acc1_blk = block_acc(r1, blocks, {i1_sam, i1_sim, i1_new}, {'same','similar','new'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. retrieval (2-back)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an A-N trial is the 2nd item after an 'A-N' goal item. blocks are
% concatenated here, so the +2 pointer must stay inside its own block.
real2 = ~contains(r2.goal, "JUNK");
v2    = r2.rt > min_rt;
pan   = false(height(r2),1);
for i = 1:height(r2)-2
    if strcmp(r2.goal(i),'A-N') && r2.block(i+2)==r2.block(i), pan(i+2) = true; end
end
aa = real2 & strcmp(r2.goal,'A-A') & strcmp(cellstr(r2.corr_resp),'1');
ab = real2 & strcmp(r2.goal,'A-B') & strcmp(cellstr(r2.corr_resp),'2');
an = real2 & pan & strcmp(cellstr(r2.corr_resp),'none');
cnd = cellfun(@(c) real2 & strcmp(r2.condition, c), cnd_nm, 'UniformOutput', false);
goal_m = {aa, ab, an}; goal_nm = {'AA (same / target)','AB (similar / lure)','AN (new / foil)'};

fprintf('\n\n=== 2. RETRIEVAL (2-back) ========================================\n');
fprintf('\n-- accuracy by goal x condition (wilson 95%% CI) --\n');
acc2 = nan(3,nC); ci2 = nan(3,nC,2); n2 = zeros(3,nC);
for g = 1:3
    fprintf('\n  %s\n', goal_nm{g});
    for c = 1:nC
        m = goal_m{g} & cnd{c};
        [acc2(g,c), n2(g,c), lo, hi] = prop_report(sprintf('    %-9s        ', cnd_nm{c}), r2.correct, m);
        ci2(g,c,:) = [lo hi];
    end
    fprintf('    -> across conditions: %s\n', chi2_str(r2.correct, cellfun(@(cm) goal_m{g}&cm, cnd, 'UniformOutput', false)));
end

fprintf('\n-- LDI  ( p(2|AB) - p(2|AN), bootstrap 95%% CI ) --\n');
ldi = nan(1,nC); ldi_ci = nan(nC,2);
for c = 1:nC
    nh = sum(ab & cnd{c}); nf = sum(an & cnd{c});
    h  = mean(strcmp(r2.resp_key(ab & cnd{c}),'2'));
    f  = mean(strcmp(r2.resp_key(an & cnd{c}),'2'));
    [ldi(c), ldi_ci(c,:)] = boot_index(@(a,b) a-b, h, nh, f, nf, n_boot);
    fprintf('  %-9s LDI = %+.3f  [%+.3f, %+.3f]   (lure-2 %.3f of %d, foil-2 %.3f of %d)\n', ...
        cnd_nm{c}, ldi(c), ldi_ci(c,1), ldi_ci(c,2), h, nh, f, nf);
end
fprintf('  pairwise differences (bootstrap):\n');
boot_pair_diff(r2, ab, an, '2', @(a,b) a-b, cnd, cnd_nm, n_boot);

fprintf('\n-- d''  ( z(p(1|AA)) - z(p(1|AN)), bootstrap 95%% CI ) --\n');
dp = nan(1,nC); dp_ci = nan(nC,2);
for c = 1:nC
    nh = sum(aa & cnd{c}); nf = sum(an & cnd{c});
    h  = mean(strcmp(r2.resp_key(aa & cnd{c}),'1'));
    f  = mean(strcmp(r2.resp_key(an & cnd{c}),'1'));
    [dp(c), dp_ci(c,:)] = boot_index(@(a,b) calc_d(a,b,nh,nf), h, nh, f, nf, n_boot);
    fprintf('  %-9s d'' = %.3f  [%.3f, %.3f]   (hit %.3f of %d, FA %.3f of %d)\n', ...
        cnd_nm{c}, dp(c), dp_ci(c,1), dp_ci(c,2), h, nh, f, nf);
end
fprintf('  pairwise differences (bootstrap):\n');
boot_pair_diff_d(r2, aa, an, cnd, cnd_nm, n_boot);

fprintf('\n-- RT on correct trials --\n');
rt_l = cell(1,nC); rt_t = cell(1,nC);
for c = 1:nC
    rt_l{c} = r2.rt(ab & cnd{c} & r2.correct & v2);
    rt_t{c} = r2.rt(aa & cnd{c} & r2.correct & v2);
end
fprintf('  lure (AB):\n');   for c=1:nC, rt_report(sprintf('  %-9s', cnd_nm{c}), rt_l{c}); end
fprintf('    -> %s\n', kw_str(rt_l));
fprintf('  target (AA):\n'); for c=1:nC, rt_report(sprintf('  %-9s', cnd_nm{c}), rt_t{c}); end
fprintf('    -> %s\n', kw_str(rt_t));

fprintf('\n-- confusions per condition (row = presented, col = response) --\n');
m2 = cell(1,nC); cc2 = cell(1,nC);
for c = 1:nC
    fprintf('\n  %s:\n', cnd_nm{c});
    [m2{c}, cc2{c}] = conf_matrix(r2.resp_key, {aa&cnd{c}, ab&cnd{c}, an&cnd{c}}, {'1','2','none'});
    print_conf(m2{c}, cc2{c}, {'same','similar','new'}, {'->same','->similar','->new'});
end

fprintf('\n-- accuracy by block --\n');
acc2_blk = block_acc(r2, blocks, {aa, ab, an}, {'AA same','AB similar','AN new'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. post-task MST (old / lure / new)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old  = seen item (correct = same, '1');  lure = unseen pairmate (correct =
% similar, '2');  new = fresh foil (correct = withhold).
%   LDI = p(2|lure) - p(2|new)   (lure discrimination / pattern separation)
%   REC = p(1|old)  - p(1|new)   (recognition)
mst = struct();
fprintf('\n\n=== 3. POST-TASK MST =============================================\n');
if has_mst
    rm.resp_key = cellstr(rm.resp_key); rm.resp_key(strcmp(rm.resp_key,'NA')) = {'none'};
    rm.correct  = strcmp(cellstr(rm.corr_resp), rm.resp_key);
    tt  = string(rm.trial_type);
    old = tt=="old"; lure = tt=="lure"; new = tt=="new";
    nO  = sum(old); nL = sum(lure); nN = sum(new);
    vM  = rm.rt > min_rt;
    fprintf('\n  %d trials (old/lure/new = %d/%d/%d), no-response %.1f%%\n', ...
        height(rm), nO, nL, nN, 100*mean(strcmp(rm.resp_key,'none')));

    hit = pkey(rm,old,'1');  lsim = pkey(rm,lure,'2'); cr = pkey(rm,new,'none'); fa = pkey(rm,new,'1');
    fprintf('\n-- rates (wilson 95%% CI) --\n');
    prop_rate('hit   p(1|old)     ', hit,  nO);
    prop_rate('lure  p(2|lure)    ', lsim, nL);
    prop_rate('CR    p(none|new)  ', cr,   nN);
    prop_rate('FA    p(1|new)     ', fa,   nN);

    fprintf('\n-- indices (bootstrap 95%% CI) --\n');
    [mst.ldi, ldci] = boot_index(@(a,b) a-b, lsim, nL, pkey(rm,new,'2'), nN, n_boot);
    [mst.rec, rcci] = boot_index(@(a,b) a-b, hit,  nO, fa,               nN, n_boot);
    [mst.dpr, dpci] = boot_index(@(a,b) calc_d(a,b,nO,nN), hit, nO, fa, nN, n_boot);
    fprintf('  LDI  = %+.3f  [%+.3f, %+.3f]\n', mst.ldi, ldci(1), ldci(2));
    fprintf('  REC  = %+.3f  [%+.3f, %+.3f]\n', mst.rec, rcci(1), rcci(2));
    fprintf('  d''   = %.3f  [%.3f, %.3f]\n',    mst.dpr, dpci(1), dpci(2));

    fprintf('\n-- RT on correct trials --\n');
    rt_report('old ', rm.rt(old  & rm.correct & vM));
    rt_report('lure', rm.rt(lure & rm.correct & vM));

    fprintf('\n-- confusion (row = presented, col = response) --\n');
    [mM, cM] = conf_matrix(rm.resp_key, {old, lure, new}, {'1','2','none'});
    print_conf(mM, cM, {'old','lure','new'}, {'->old','->similar','->new'});
    mst.rate = [hit lsim cr]; mst.fa = fa; mst.conf = mM; mst.conf_n = cM;
    mst.ldi_ci = ldci; mst.rec_ci = rcci; mst.dpr_ci = dpci; mst.n = [nO nL nN];
else
    fprintf('\n  no MST data for this subject.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indiv = struct('subj_id',subj_id,'conditions',{cnd_nm},'qc',qc,'one',one, ...
    'acc2',acc2,'acc2_ci',ci2,'acc2_n',n2, ...
    'ldi',ldi,'ldi_ci',ldi_ci,'dprime',dp,'dprime_ci',dp_ci, ...
    'conf_1back',m1,'conf_1back_n',c1,'conf_2back',{m2},'conf_2back_n',{cc2}, ...
    'acc1_blk',acc1_blk,'acc2_blk',acc2_blk,'mst',mst);
save(fullfile(res_dir, sprintf('sub%03d_indiv_stats.mat', subj_id)), 'indiv');
fprintf('\n\nsaved: %s\n', fullfile(res_dir, sprintf('sub%03d_indiv_stats.mat', subj_id)));
fprintf('saved: %s\n', diary_file);
diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 1: data quality
f = figure('color','w','Position',[60 60 1250 800],'Name','QC');
subplot(2,2,1); hold on;
miss1 = arrayfun(@(b) mean(strcmp(r1.resp_key(r1.block==b),'none') & ~strcmp(cellstr(r1.corr_resp(r1.block==b)),'none')), blocks);
miss2 = arrayfun(@(b) mean(strcmp(r2.resp_key(r2.block==b),'none') & ~strcmp(cellstr(r2.corr_resp(r2.block==b)),'none')), blocks);
bar(blocks, 100*[miss1; miss2]'); ylabel('% of trials'); xlabel('block');
title('missed responses (should have responded)','FontSize',13); set(gca,'XTick',blocks);
legend({'1-back','2-back'},'Location','best','Box','off'); box off; hold off;

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
raincloud_cell({rt_sam, rt_sim}, {c_same,c_sim}, {'same','similar'}, 'RT (s)', '1-back RT (correct)', {[1 2]});
subplot(1,3,3);
draw_matrix(m1, c1, {c_same,c_sim,c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
title('1-back confusions','FontSize',13);
save_fig(f, fig_dir, sprintf('sub%03d_1back', subj_id));

%% fig 3: 2-back
f = figure('color','w','Position',[40 40 1350 900],'Name','2-back');
subplot(2,2,1);
grouped_bar_ci(acc2, ci2, n2, cnd_cols, {'AA same','AB similar','AN new'}, cnd_nm, 'accuracy', '2-back accuracy');
yline(1/3,'k:','chance','LineWidth',1.2);
subplot(2,2,2);
bar_ci_explicit(ldi, ldi_ci, cnd_cols, cnd_nm, 'LDI', 'lure discrimination');
yline(0,'k-');
subplot(2,2,3);
bar_ci_explicit(dp, dp_ci, cnd_cols, cnd_nm, 'd''', 'target detection');
yline(0,'k-');
subplot(2,2,4);
rc_lbl = [cellfun(@(c) ['AB ' c], cnd_nm, 'UniformOutput', false), ...
          cellfun(@(c) ['AA ' c], cnd_nm, 'UniformOutput', false)];
% bracket within goal type (AB across conditions, then AA across conditions)
rc_pairs = {[1 2],[1 3],[2 3],[4 5],[4 6],[5 6]};
raincloud_cell([rt_l, rt_t], [cnd_cols, cnd_cols], rc_lbl, 'RT (s)', '2-back RT (correct)', rc_pairs);
xtickangle(30);
save_fig(f, fig_dir, sprintf('sub%03d_2back', subj_id));

%% fig 4: 2-back confusion matrices
f = figure('color','w','Position',[40 40 500*nC 400],'Name','confusions');
for c = 1:nC
    subplot(1,nC,c);
    draw_matrix(m2{c}, cc2{c}, {c_same,c_sim,c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
    title(sprintf('2-back: %s', cnd_nm{c}), 'FontSize', 13);
end
save_fig(f, fig_dir, sprintf('sub%03d_confusion', subj_id));

%% fig 5: MST
if has_mst
    f = figure('color','w','Position',[100 100 1300 420],'Name','MST');
    subplot(1,3,1);
    bar_ci(mst.rate, mst.n, {c_same,c_sim,c_new}, {'hit (old)','lure->sim','CR (new)'}, 'rate', 'MST rates');
    ylim([0 1.1]); yline(mst.fa,'r--',sprintf('FA = %.2f',mst.fa),'LineWidth',1.5);
    subplot(1,3,2);
    bar_ci_explicit([mst.ldi mst.rec mst.dpr], [mst.ldi_ci; mst.rec_ci; mst.dpr_ci], ...
        {c_sim,c_same,[.5 .5 .5]}, {'LDI','REC','d'''}, 'index', 'MST indices');
    yline(0,'k-');
    subplot(1,3,3);
    draw_matrix(mst.conf, mst.conf_n, {c_same,c_sim,c_new}, {'Exp Old','Exp Lure','Exp New'}, {'Resp Old','Resp Sim','Resp New'});
    title('MST confusions','FontSize',13);
    save_fig(f, fig_dir, sprintf('sub%03d_mst', subj_id));
end

fprintf('\nfigures written to %s\n', fig_dir);

%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%
function [lo, hi] = wilson(p, n)
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

function prop_rate(lbl, p, n)
    if n==0 || isnan(p), fprintf('  %s  (no trials)\n', lbl); return; end
    [lo, hi] = wilson(p, n);
    fprintf('  %s  %.3f  [%.3f, %.3f]   n=%d\n', lbl, p, lo, hi, n);
end

function p = pkey(T, mask, key)
    if sum(mask)==0, p = NaN; else, p = mean(strcmp(T.resp_key(mask), key)); end
end

function prop_test(lbl, a, b)
    n1 = numel(a); n2 = numel(b);
    if n1==0 || n2==0, fprintf('  %s: insufficient trials\n', lbl); return; end
    p1 = mean(a); p2 = mean(b); pp = (sum(a)+sum(b))/(n1+n2);
    se = sqrt(pp*(1-pp)*(1/n1+1/n2));
    if se == 0, fprintf('  %s: identical proportions (%.3f)\n', lbl, p1); return; end
    z = (p1-p2)/se; p = 2*(1-normcdf(abs(z)));
    h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2));
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
    r = abs(z)/sqrt(numel(a)+numel(b));
    fprintf('  %s: rank-sum z = %.2f, r = %.2f, p = %.4f %s\n', lbl, z, r, p, stars(p));
end

function s = chi2_str(correct, masks)
    obs = zeros(numel(masks),2);
    for k = 1:numel(masks)
        obs(k,:) = [sum(correct(masks{k})), sum(masks{k})-sum(correct(masks{k}))];
    end
    obs(sum(obs,2)==0,:) = [];
    if size(obs,1) < 2 || any(sum(obs,1)==0), s = 'chi2 not computable'; return; end
    exp_ = sum(obs,2)*sum(obs,1)/sum(obs(:));
    x2 = sum(sum((obs-exp_).^2 ./ exp_));
    df = (size(obs,1)-1)*(size(obs,2)-1); p = 1-chi2cdf(x2, df);
    v  = sqrt(x2/(sum(obs(:))*min(size(obs)-1)));
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
    d = norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));
end

function [est, ci, bs] = boot_index(fn, h, nh, f, nf, n_boot)
    if nh==0 || nf==0 || isnan(h) || isnan(f), est = NaN; ci = [NaN NaN]; bs = nan(n_boot,1); return; end
    est = fn(h, f);
    bh  = binornd(nh, h, n_boot, 1)/nh;
    bf  = binornd(nf, f, n_boot, 1)/nf;
    bs  = arrayfun(fn, bh, bf);
    ci  = quantile(bs, [0.025 0.975]);
end

function boot_pair_diff(r2, hit_mask, fa_mask, key, fn, cnd, cnd_nm, n_boot)
    nc = numel(cnd); bs = cell(1,nc); est = nan(1,nc);
    for c = 1:nc
        nh = sum(hit_mask & cnd{c}); nf = sum(fa_mask & cnd{c});
        h  = mean(strcmp(r2.resp_key(hit_mask & cnd{c}), key));
        f  = mean(strcmp(r2.resp_key(fa_mask  & cnd{c}), key));
        [est(c), ~, bs{c}] = boot_index(fn, h, nh, f, nf, n_boot);
    end
    report_pairs(est, bs, cnd_nm);
end

function boot_pair_diff_d(r2, aa, an, cnd, cnd_nm, n_boot)
    nc = numel(cnd); bs = cell(1,nc); est = nan(1,nc);
    for c = 1:nc
        nh = sum(aa & cnd{c}); nf = sum(an & cnd{c});
        h  = mean(strcmp(r2.resp_key(aa & cnd{c}),'1'));
        f  = mean(strcmp(r2.resp_key(an & cnd{c}),'1'));
        [est(c), ~, bs{c}] = boot_index(@(x,y) calc_d(x,y,nh,nf), h, nh, f, nf, n_boot);
    end
    report_pairs(est, bs, cnd_nm);
end

function report_pairs(est, bs, cnd_nm)
    pr = nchoosek(1:numel(cnd_nm), 2);
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

function raincloud_cell(cells, cols, xlbls, ylbl, ttl, pairs)
    if nargin < 6 || isempty(pairs)
        pairs = arrayfun(@(k) [k k+1], 1:numel(cells)-1, 'UniformOutput', false);
    end
    hold on; n = numel(cells); all_v = [];
    bw = 0.30; jw = 0.11;
    for i = 1:n
        d = cells{i}(:); d = d(~isnan(d)); all_v = [all_v; d]; %#ok<AGROW>
        if isempty(d), continue; end
        jit = (rand(numel(d),1)-0.5)*2*jw;
        scatter(i+jit, d, 24, cols{i}, 'filled', 'MarkerFaceAlpha',0.55, ...
            'MarkerEdgeColor',[.25 .25 .25], 'LineWidth',0.25);
        q = quantile(d,[0.25 0.5 0.75]); iqrv = q(3)-q(1);
        lo_w = max(min(d), q(1)-1.5*iqrv); hi_w = min(max(d), q(3)+1.5*iqrv);
        plot([i i],[lo_w q(1)],'k-','LineWidth',0.9);
        plot([i i],[q(3) hi_w],'k-','LineWidth',0.9);
        plot(i+[-.06 .06],[lo_w lo_w],'k-','LineWidth',0.9);
        plot(i+[-.06 .06],[hi_w hi_w],'k-','LineWidth',0.9);
        rectangle('Position',[i-bw/2, q(1), bw, max(iqrv,eps)], ...
            'EdgeColor','k','LineWidth',1.1,'FaceColor','none');
        plot(i+[-bw/2 bw/2],[q(2) q(2)],'k-','LineWidth',1.8);
        text(i+bw/2+0.03, q(2), sprintf('n=%d',numel(d)),'HorizontalAlignment','left', ...
            'VerticalAlignment','middle','FontSize',7.5,'Color',[.5 .5 .5]);
    end
    if ~isempty(all_v)
        yr = range(all_v); if yr==0, yr = 1; end
        base = max(all_v) + 0.12*yr; step = 0.11*yr; k = 0;
        for pp = 1:numel(pairs)
            ij = pairs{pp};
            a = cells{ij(1)}(:); a = a(~isnan(a));
            b = cells{ij(2)}(:); b = b(~isnan(b));
            if numel(a)<3 || numel(b)<3, continue; end
            pv = ranksum(a,b); s = stars(pv); if isempty(s), s = 'ns'; end
            y = base + k*step; k = k + 1;
            plot([ij(1) ij(1) ij(2) ij(2)], [y-0.02*yr y y y-0.02*yr], 'k-','LineWidth',0.9);
            text(mean(ij), y, s, 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize',11);
        end
        ylim([min(all_v)-0.06*yr, base + max(k,1)*step + 0.05*yr]);
    end
    set(gca,'XTick',1:n,'XTickLabel',xlbls,'FontSize',12);
    xlim([0.4 n+0.6]); ylabel(ylbl,'FontSize',13);
    if ~isempty(ttl), title(ttl,'FontSize',13); end
    box off; hold off;
end

function grouped_bar_ci(vals, ci, ns, cols, glbls, clbls, ylbl, ttl)
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
        'Box','off','AutoUpdate','off');
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
    pos = get(f, 'Position');
    set(f, 'PaperUnits', 'points', 'PaperSize', pos(3:4));
    print(f, fullfile(fig_dir, [name '.png']), '-dpng', '-r150');
end
