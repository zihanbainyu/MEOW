clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [618, 619, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631];
base_dir = '..';
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');
min_rt   = 0.150;

% colors
c_comp = [160 174 211]/255; c_iso = [176 230 255]/255; c_nov = [163 210 205]/255;
c_sim  = [255 191 205]/255; c_same = [97 125 164]/255; c_new = [219 219 219]/255;

%%%%%%%%%%%%%%%%%%%%%%%
%% load & compute subject-level stats
%%%%%%%%%%%%%%%%%%%%%%%
all_subjs = struct();
trl_seq   = {};  % trial-level table for sequential dependency LME only

for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('processing subject %d...\n', curr_id);
    subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', curr_id));
    load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', curr_id)), 'final_data_output');
    r1 = final_data_output.results_1_back_all;
    r2 = final_data_output.results_2_back_all;
    stats = [];

    %% 1-back
    r1.resp_key = cellstr(r1.resp_key); r1.resp_key(strcmp(r1.resp_key,'NA')) = {'none'};
    r1.correct  = strcmp(cellstr(r1.corr_resp), r1.resp_key);
    idx_sim = r1.condition=="compared" & r1.identity=="B";
    idx_sam = r1.condition=="repeat"   & strcmp(r1.corr_resp,'j');
    idx_cr  = (r1.condition=="compared" & r1.identity=="A") | r1.condition=="isolated" | ...
              (r1.condition=="repeat"   & strcmp(r1.corr_resp,'none'));
    stats.one.acc_sim  = mean(r1.correct(idx_sim));
    stats.one.acc_same = mean(r1.correct(idx_sam));
    stats.one.acc_new  = mean(r1.correct(idx_cr));
    v_rt = r1.rt > min_rt;
    stats.one.rt_sim  = median(r1.rt(idx_sim & r1.correct & v_rt),'omitnan');
    stats.one.rt_same = median(r1.rt(idx_sam & r1.correct & v_rt),'omitnan');
    n_sim = sum(idx_sim); n_sam = sum(idx_sam); n_cr = sum(idx_cr);
    stats.one.err_sim_as_same = sum(strcmp(r1.resp_key(idx_sim),'j'))/n_sim;
    stats.one.err_sim_as_new  = sum(strcmp(r1.resp_key(idx_sim),'none'))/n_sim;
    stats.one.err_same_as_sim = sum(strcmp(r1.resp_key(idx_sam),'k'))/n_sam;
    stats.one.err_same_as_new = sum(strcmp(r1.resp_key(idx_sam),'none'))/n_sam;
    stats.one.err_new_as_same = sum(strcmp(r1.resp_key(idx_cr),'j'))/n_cr;
    stats.one.err_new_as_sim  = sum(strcmp(r1.resp_key(idx_cr),'k'))/n_cr;

    %% 2-back
    r2.resp_key = cellstr(r2.resp_key); r2.resp_key(strcmp(r2.resp_key,'NA')) = {'none'};
    r2.correct  = strcmp(cellstr(r2.corr_resp), r2.resp_key);
    pan = zeros(height(r2),1);
    for i = 1:height(r2)-2
        if strcmp(r2.goal(i),'A-N'), pan(i+2) = 1; end
    end
    real     = ~contains(r2.goal,"JUNK"); v_rt = r2.rt > min_rt;
    ab_idx   = real & strcmp(r2.goal,'A-B');
    aa_idx   = real & strcmp(r2.goal,'A-A');
    an_idx   = logical(pan);
    comp_idx = real & strcmp(r2.condition,'compared');
    iso_idx  = real & strcmp(r2.condition,'isolated');
    nov_idx  = real & strcmp(r2.condition,'novel');
    j_idx    = real & strcmp(r2.corr_resp,'j');
    k_idx    = real & strcmp(r2.corr_resp,'k');
    n_idx    = real & strcmp(r2.corr_resp,'none');
    calc_d   = @(h,f,nh,nf) norminv(max(1/(2*nh),min(1-1/(2*nh),h))) - ...
                              norminv(max(1/(2*nf),min(1-1/(2*nf),f)));

    for cond = {'comp','iso','nov'}
        c = cond{1};
        switch c
            case 'comp', cidx = comp_idx;
            case 'iso',  cidx = iso_idx;
            case 'nov',  cidx = nov_idx;
        end
        n_AA = sum(aa_idx & cidx & j_idx);
        n_AB = sum(ab_idx & cidx & k_idx);
        n_AN = sum(an_idx & cidx & n_idx);
        stats.two.(sprintf('acc_AA_%s',c)) = mean(r2.correct(aa_idx & cidx & j_idx));
        stats.two.(sprintf('acc_AB_%s',c)) = mean(r2.correct(ab_idx & cidx & k_idx));
        stats.two.(sprintf('acc_AN_%s',c)) = mean(r2.correct(an_idx & cidx & n_idx));
        stats.two.(sprintf('rt_AA_%s',c))  = median(r2.rt(aa_idx & cidx & j_idx & r2.correct & v_rt),'omitnan');
        stats.two.(sprintf('rt_AB_%s',c))  = median(r2.rt(ab_idx & cidx & k_idx & r2.correct & v_rt),'omitnan');
        stats.two.(sprintf('err_AA_%s_as_k',c)) = sum(strcmp(r2.resp_key(aa_idx & cidx & j_idx),'k'))/n_AA;
        stats.two.(sprintf('err_AA_%s_as_n',c)) = sum(strcmp(r2.resp_key(aa_idx & cidx & j_idx),'none'))/n_AA;
        stats.two.(sprintf('err_AB_%s_as_j',c)) = sum(strcmp(r2.resp_key(ab_idx & cidx & k_idx),'j'))/n_AB;
        stats.two.(sprintf('err_AB_%s_as_n',c)) = sum(strcmp(r2.resp_key(ab_idx & cidx & k_idx),'none'))/n_AB;
        stats.two.(sprintf('err_AN_%s_as_j',c)) = sum(strcmp(r2.resp_key(an_idx & cidx & n_idx),'j'))/n_AN;
        stats.two.(sprintf('err_AN_%s_as_k',c)) = sum(strcmp(r2.resp_key(an_idx & cidx & n_idx),'k'))/n_AN;
        stats.two.(sprintf('dprime_%s',c))  = calc_d(...
            stats.two.(sprintf('acc_AA_%s',c)), ...
            stats.two.(sprintf('err_AN_%s_as_j',c)), n_AA, n_AN);
    end

    stats.two.ldi_comp = stats.two.acc_AB_comp - stats.two.err_AN_comp_as_k;
    stats.two.ldi_iso  = stats.two.acc_AB_iso  - stats.two.err_AN_iso_as_k;
    stats.two.ldi_nov  = stats.two.acc_AB_nov  - stats.two.err_AN_nov_as_k;

    %% sequential dependency trial table (A-B trials only)
    for cond = {'compared','isolated','novel'}
        c = cond{1};
        switch c
            case 'compared', cidx = comp_idx;
            case 'isolated', cidx = iso_idx;
            case 'novel',    cidx = nov_idx;
        end
        ab_cond_rows = find(ab_idx & cidx & k_idx);
        for ti = 2:length(ab_cond_rows)
            curr_t = ab_cond_rows(ti);
            prev_t = ab_cond_rows(ti-1);
            trl_seq{end+1} = table(curr_id, {c}, ...
                double(r2.correct(curr_t)), double(r2.correct(prev_t)), ...
                'VariableNames', {'subj','condition','correct','prev_correct'});
        end
    end

    %% recognition
    if isfield(final_data_output,'results_recognition')
        rec = final_data_output.results_recognition;
        rec.correct = strcmp(cellstr(rec.corr_resp), cellstr(rec.resp_key));
        old_r = rec(rec.trial_type=="old",:);
        new_r = rec(rec.trial_type~="old",:);
        n_new = height(new_r);
        n_fa  = sum(strcmp(cellstr(new_r.resp_key),'j') & new_r.rt>min_rt);
        far   = max(1/(2*n_new), min(1-1/(2*n_new), n_fa/n_new));
        for cond = {'compared','isolated'}
            c = cond{1}; tc = old_r(old_r.condition==c,:); nc = height(tc);
            hc = sum(tc.correct & tc.rt>min_rt)/nc;
            stats.rec.(sprintf('d_%s',c(1:3))) = calc_d(hc, far, nc, n_new);
        end
    else
        fprintf('  recognition missing for subject %d.\n', curr_id);
        stats.rec.d_com = NaN; stats.rec.d_iso = NaN;
    end

    all_subjs(s).id    = curr_id;
    all_subjs(s).stats = stats;
end

% finalize sequential dependency table
trl_seq           = vertcat(trl_seq{:});
trl_seq.subj      = categorical(trl_seq.subj);
trl_seq.condition = categorical(trl_seq.condition, {'novel','isolated','compared'});

save(fullfile(res_dir,'all_subjs_stats.mat'),'all_subjs');

%%%%%%%%%%%%%%%%%%%%%%%
%% extract group-level vectors
%%%%%%%%%%%%%%%%%%%%%%%
get_v = @(f1,f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs)';

b1_acc_sam = get_v('one','acc_same'); b1_acc_sim = get_v('one','acc_sim'); b1_acc_new = get_v('one','acc_new');
b1_rt_sam  = get_v('one','rt_same'); b1_rt_sim  = get_v('one','rt_sim');

b2_ldi_c  = get_v('two','ldi_comp');    b2_ldi_i  = get_v('two','ldi_iso');    b2_ldi_n  = get_v('two','ldi_nov');
b2_dp_c   = get_v('two','dprime_comp'); b2_dp_i   = get_v('two','dprime_iso'); b2_dp_n   = get_v('two','dprime_nov');
b2_rt_l_c = get_v('two','rt_AB_comp');  b2_rt_l_i = get_v('two','rt_AB_iso');  b2_rt_l_n = get_v('two','rt_AB_nov');
b2_rt_t_c = get_v('two','rt_AA_comp');  b2_rt_t_i = get_v('two','rt_AA_iso');  b2_rt_t_n = get_v('two','rt_AA_nov');
b2_hr_c   = get_v('two','acc_AB_comp'); b2_hr_i   = get_v('two','acc_AB_iso'); b2_hr_n   = get_v('two','acc_AB_nov');
b2_far_c  = get_v('two','err_AN_comp_as_k'); b2_far_i = get_v('two','err_AN_iso_as_k'); b2_far_n = get_v('two','err_AN_nov_as_k');

rec_d_c = get_v('rec','d_com'); rec_d_i = get_v('rec','d_iso');
d_tot   = (rec_d_c + rec_d_i)/2;
ldi_adv = b2_ldi_c - b2_ldi_i;

%%%%%%%%%%%%%%%%%%%%%%%
%% statistical tests
%%%%%%%%%%%%%%%%%%%%%%%
within = table(categorical({'compared';'isolated';'novel'}), 'VariableNames',{'Condition'});
stat_results = struct();

%% 1-back
fprintf('\n=== 1-BACK STATS ===\n');
[~,p,~,s] = ttest(b1_rt_sam,  b1_rt_sim);  stat_results.b1_rt     = report_t(s,p,'RT: same vs sim');
[~,p,~,s] = ttest(b1_acc_sam, b1_acc_sim); stat_results.b1_acc_sv = report_t(s,p,'Acc: same vs sim');
[~,p,~,s] = ttest(b1_acc_sim, b1_acc_new); stat_results.b1_acc_sn = report_t(s,p,'Acc: sim vs new');

%% 2-back: primary outcomes
fprintf('\n=== 2-BACK: LURE DISCRIMINATION (LDI) ===\n');
[stat_results.ldi, stat_results.ldi_ph] = run_rm_anova('LDI', b2_ldi_c, b2_ldi_i, b2_ldi_n, within);

fprintf('\n=== 2-BACK: RT LURE TRIALS ===\n');
[stat_results.rt_lure, stat_results.rt_lure_ph] = run_rm_anova('RT lure', b2_rt_l_c, b2_rt_l_i, b2_rt_l_n, within);

%% 2-back: selectivity controls (expect null + BF01)
fprintf('\n=== 2-BACK: TARGET DETECTION d-prime — selectivity control ===\n');
[stat_results.dprime, stat_results.dprime_ph] = run_rm_anova('d-prime AA', b2_dp_c, b2_dp_i, b2_dp_n, within);

fprintf('\n=== 2-BACK: RT TARGET TRIALS — selectivity control ===\n');
[stat_results.rt_targ, stat_results.rt_targ_ph] = run_rm_anova('RT target', b2_rt_t_c, b2_rt_t_i, b2_rt_t_n, within);

%% 2-back: response bias control (expect null + BF01)
fprintf('\n=== 2-BACK: FALSE ALARM RATE — response bias control ===\n');
fprintf('  Mean FAR: comp=%.4f, iso=%.4f, nov=%.4f\n', ...
    mean(b2_far_c,'omitnan'), mean(b2_far_i,'omitnan'), mean(b2_far_n,'omitnan'));
[stat_results.far, stat_results.far_ph] = run_rm_anova('FAR', b2_far_c, b2_far_i, b2_far_n, within);

%% recognition
fprintf('\n=== RECOGNITION ===\n');
[~,p,~,s] = ttest(d_tot);            stat_results.rec_overall = report_t(s,p,'Recognition d'' vs 0');
[~,p,~,s] = ttest(rec_d_c, rec_d_i); stat_results.rec_ci      = report_t(s,p,'Recognition: compared vs isolated');

%% individual differences
fprintf('\n=== INDIVIDUAL DIFFERENCES ===\n');
[r_enc, p_enc] = corr(b1_acc_sim, ldi_adv, 'rows','complete');
n_valid = sum(~isnan(b1_acc_sim) & ~isnan(ldi_adv));
fprintf('  1-back sim acc vs LDI advantage (comp-iso): r=%.3f, p=%.4f, n=%d\n', r_enc, p_enc, n_valid);
stat_results.indiv.r = r_enc; stat_results.indiv.p = p_enc;

%% sequential dependency — trial-level LME
% rm-ANOVA on subject means cannot answer this question because the
% question is inherently trial-level: does trial n-1 outcome predict
% trial n outcome, and does this differ between conditions?
% The interaction term (prev_correct × condition) is the critical test.
fprintf('\n=== SEQUENTIAL DEPENDENCY (trial-level LME) ===\n');
formula_full = 'correct ~ prev_correct * condition + (1 | subj)';
formula_add  = 'correct ~ prev_correct + condition + (1 | subj)';
mdl_full = fitglme(trl_seq, formula_full, ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
mdl_add  = fitglme(trl_seq, formula_add,  ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
% likelihood ratio test: interaction term
% with this:
lrt_tbl = compare(mdl_add, mdl_full);
lrt_p   = lrt_tbl.pValue(2);
% BIC-based BF (Wagenmakers 2007): BF10 ≈ exp((BIC_null - BIC_full)/2)
bf_int = exp((mdl_add.ModelCriterion.BIC - mdl_full.ModelCriterion.BIC)/2);
fprintf('  Interaction (prev_correct x condition) LRT: p=%.4f\n', lrt_p);
if bf_int < 1
    fprintf('  BF01=%.2f — no condition-dependent sequential effect\n', 1/bf_int);
else
    fprintf('  BF10=%.2f — condition-dependent sequential effect present\n', bf_int);
end
% main effect of prev_correct: is there any dependency at all?
coef     = mdl_full.Coefficients;
prev_row = strcmp(coef.Name,'prev_correct');
fprintf('  prev_correct main effect: beta=%.3f, SE=%.3f, t=%.2f, p=%.4f\n', ...
    coef.Estimate(prev_row), coef.SE(prev_row), ...
    coef.tStat(prev_row), coef.pValue(prev_row));
stat_results.seq.lrt_p  = lrt_p;
stat_results.seq.bf_int = bf_int;
stat_results.seq.coef   = coef;

save(fullfile(res_dir,'stat_results.mat'),'stat_results');

%%%%%%%%%%%%%%%%%%%%%%%
%% visualization
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== GENERATING FIGURES ===\n');

%% fig 1: 1-back
figure('color','w','Position',[100 100 1200 400]);
subplot(1,3,1);
raincloud([b1_acc_sam, b1_acc_sim, b1_acc_new], {c_same,c_sim,c_new}, ...
    {'Hit(Same)','Hit(Sim)','CR(New)'}, 'Accuracy', '1-Back Accuracy', [0,1]);
add_sig([b1_acc_sam, b1_acc_sim, b1_acc_new], [1 2; 2 3; 1 3]);
subplot(1,3,2);
raincloud([b1_rt_sam, b1_rt_sim], {c_same,c_sim}, ...
    {'Hit(Same)','Hit(Sim)'}, 'RT (s)', '1-Back RT');
add_sig([b1_rt_sam, b1_rt_sim], [1 2]);
subplot(1,3,3); hold on;
fields_1b = {'acc_same','err_same_as_sim','err_same_as_new';
             'err_sim_as_same','acc_sim','err_sim_as_new';
             'err_new_as_same','err_new_as_sim','acc_new'};
mat_1b_mean = cellfun(@(f) mean(get_v('one',f)), fields_1b);
mat_1b_se   = cellfun(@(f) std(get_v('one',f))/sqrt(length(subj_ids)), fields_1b);
draw_matrix(mat_1b_mean, mat_1b_se, {c_same,c_sim,c_new}, ...
    {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
title('Response Matrix');
set(gcf,'PaperPositionMode','auto');
print(gcf, fullfile(fig_dir,'fig_1back.pdf'),'-dpdf','-vector');

%% fig 2: 2-back primary + controls
figure('color','w','Position',[50 50 1400 850]);

subplot(2,4,1);
raincloud([b2_ldi_c, b2_ldi_i, b2_ldi_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'LDI', 'Lure Discrimination Index', [-0.2,1]);
set(gca,'YTick',[0 0.5 1]);
add_sig([b2_ldi_c, b2_ldi_i, b2_ldi_n], [1 2; 2 3; 1 3]);

subplot(2,4,2);
raincloud([b2_hr_c, b2_hr_i, b2_hr_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'P(sim|A-B)', 'Hit Rate', [0,1]);
add_sig([b2_hr_c, b2_hr_i, b2_hr_n], [1 2; 2 3; 1 3]);

subplot(2,4,3);
raincloud([b2_far_c, b2_far_i, b2_far_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'P(sim|A-N)', 'False Alarm Rate', [0,0.15]);
add_sig([b2_far_c, b2_far_i, b2_far_n], [1 2; 2 3; 1 3]);

subplot(2,4,4);
raincloud([b2_dp_c, b2_dp_i, b2_dp_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'd''', 'Target Detection (selectivity control)');
add_sig([b2_dp_c, b2_dp_i, b2_dp_n], [1 2; 2 3; 1 3]);

subplot(2,4,5);
raincloud([b2_rt_l_c, b2_rt_l_i, b2_rt_l_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'RT (s)', 'RT: Lure Trials', [0,1.5]);
set(gca,'YTick',[0 0.5 1 1.5]);
add_sig([b2_rt_l_c, b2_rt_l_i, b2_rt_l_n], [1 2; 2 3; 1 3]);

subplot(2,4,6);
raincloud([b2_rt_t_c, b2_rt_t_i, b2_rt_t_n], {c_comp,c_iso,c_nov}, ...
    {'Compared','Isolated','Novel'}, 'RT (s)', 'RT: Target Trials (selectivity control)', [0,1.5]);
set(gca,'YTick',[0 0.5 1 1.5]);
add_sig([b2_rt_t_c, b2_rt_t_i, b2_rt_t_n], [1 2; 2 3; 1 3]);

subplot(2,4,7);
scatter_with_fit(b1_acc_sim, ldi_adv, c_comp, ...
    '1-Back Sim Accuracy', 'LDI Advantage (Comp - Iso)', ...
    sprintf('Encoding quality → WM benefit\nr=%.2f, p=%.3f', r_enc, p_enc));

set(gcf,'Color','w','Renderer','painters','PaperPositionMode','auto');
print(gcf, fullfile(fig_dir,'fig_2back.pdf'),'-dpdf','-vector');

%% fig 3: confusion matrices
figure('color','w','Position',[50 50 1200 400]);
conds_str   = {'comp','iso','nov'};
cond_labels = {'Compared','Isolated','Novel'};
for ci = 1:3
    subplot(1,3,ci); hold on; c = conds_str{ci};
    fields_2b = {sprintf('acc_AA_%s',c),      sprintf('err_AA_%s_as_k',c), sprintf('err_AA_%s_as_n',c);
                 sprintf('err_AB_%s_as_j',c), sprintf('acc_AB_%s',c),      sprintf('err_AB_%s_as_n',c);
                 sprintf('err_AN_%s_as_j',c), sprintf('err_AN_%s_as_k',c), sprintf('acc_AN_%s',c)};
    mat_mean = cellfun(@(f) mean(get_v('two',f)), fields_2b);
    mat_se   = cellfun(@(f) std(get_v('two',f))/sqrt(length(subj_ids)), fields_2b);
    draw_matrix(mat_mean, mat_se, {c_same,c_sim,c_new}, ...
        {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
    title(cond_labels{ci},'FontSize',14);
end
sgtitle('Response Matrices','FontSize',16);
set(gcf,'PaperPositionMode','auto');
print(gcf, fullfile(fig_dir,'fig_confmat.pdf'),'-dpdf','-vector');

%% fig 4: recognition
figure('color','w','Position',[100 100 600 500]);
data_rec = [d_tot, rec_d_c, rec_d_i]; hold on;
fill([-1 4 4 -1],[-1 -1 0 0],[0.92 0.92 0.92],'EdgeColor','none');
yline(0,'r--','Chance','LineWidth',2,'LabelHorizontalAlignment','left');
raincloud(data_rec,{[0.3 0.3 0.3],c_comp,c_iso}, ...
    {'Overall','Compared','Isolated'},'d''','Recognition Memory');
add_sig(data_rec,[1 0; 2 3]);
set(gcf,'PaperPositionMode','auto');
print(gcf, fullfile(fig_dir,'fig_recognition.pdf'),'-dpdf','-vector');

fprintf('\nDone. All results saved to %s\n', res_dir);

%% ============================================================
%%                        FUNCTIONS
%% ============================================================

function [res, posthoc] = run_rm_anova(lbl, d1, d2, d3, within)
    t   = table(d1,d2,d3,'VariableNames',{'compared','isolated','novel'});
    rm  = fitrm(t,'compared-novel~1','WithinDesign',within);
    tbl = ranova(rm); eps_s = epsilon(rm); m = mauchly(rm);
    res.F      = tbl.F(1);
    res.df1    = tbl.DF(1);
    res.df2_gg = tbl.DF(2)*eps_s.GreenhouseGeisser;
    res.p_gg   = tbl.pValueGG(1);
    res.eps    = eps_s.GreenhouseGeisser;
    res.eta2p  = tbl.SumSq(1)/(tbl.SumSq(1)+tbl.SumSq(2));
    fprintf('\n%s: F(%d,%.1f)=%.2f, p=%.4f (GG), eps=%.3f, eta2p=%.3f, mauchly p=%.3f\n', ...
        lbl, res.df1, res.df2_gg, res.F, res.p_gg, res.eps, res.eta2p, m.pValue);
    % BF via BIC approximation (Wagenmakers 2007)
    n = length(d1); subj = (1:n)';
    bf_tbl  = table(subj,d1,d2,d3,'VariableNames',{'subj','c1','c2','c3'});
    bf_long = stack(bf_tbl,{'c1','c2','c3'},'NewDataVariableName','y','IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj); bf_long.cond = categorical(bf_long.cond);
    try
        mdl_eff  = fitlme(bf_long,'y~cond+(1|subj)','FitMethod','ML');
        mdl_null = fitlme(bf_long,'y~1+(1|subj)',   'FitMethod','ML');
        bf10 = exp((mdl_null.ModelCriterion.BIC - mdl_eff.ModelCriterion.BIC)/2);
        res.bf10 = bf10;
        if bf10<1, fprintf('  BF01=%.2f (evidence for null)\n',1/bf10);
        else,      fprintf('  BF10=%.2f (evidence for effect)\n',bf10); end
    catch
        res.bf10 = NaN;
        fprintf('  BF: could not compute\n');
    end
    % post-hoc with Holm-Bonferroni and Cohen's d
    [~,p1,~,s1]=ttest(d1,d2); [~,p2,~,s2]=ttest(d2,d3); [~,p3,~,s3]=ttest(d1,d3);
    adj = holm_bonf([p1 p2 p3]); df = s1.df;
    posthoc.comp_iso = struct('t',s1.tstat,'df',df,'p',p1,'p_adj',adj(1),'d',s1.tstat/sqrt(df+1));
    posthoc.iso_nov  = struct('t',s2.tstat,'df',df,'p',p2,'p_adj',adj(2),'d',s2.tstat/sqrt(df+1));
    posthoc.comp_nov = struct('t',s3.tstat,'df',df,'p',p3,'p_adj',adj(3),'d',s3.tstat/sqrt(df+1));
    fprintf('  post-hoc (Holm-Bonf):\n');
    fprintf('    comp>iso:  t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df,s1.tstat,posthoc.comp_iso.d,adj(1),stars(adj(1)));
    fprintf('    iso>nov:   t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df,s2.tstat,posthoc.iso_nov.d, adj(2),stars(adj(2)));
    fprintf('    comp>nov:  t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df,s3.tstat,posthoc.comp_nov.d,adj(3),stars(adj(3)));
end

function res = report_t(s, p, lbl)
    d = s.tstat/sqrt(s.df+1);
    fprintf('  %s: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', lbl, s.df, s.tstat, d, p, stars(p));
    res.t = s.tstat; res.df = s.df; res.p = p; res.d = d;
end

function adj = holm_bonf(p)
    n = length(p); [ps,idx] = sort(p); as = min(1,ps.*(n:-1:1));
    for i = 2:n, as(i) = max(as(i),as(i-1)); end
    adj(idx) = as;
end

function s = stars(p)
    s = repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001));
end

function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); hold on;
    d_v = mat(:); d_v = d_v(~isnan(d_v));
    mn = min(d_v); mx = max(d_v); rng = mx-mn; if rng==0, rng=1; end
    auto_lim = [mn-(rng*0.15), mx+(rng*0.15)];
    jit = -0.15 - rand(size(mat))*0.20;
    x_c = repmat(1:n_grps,n_rows,1) + jit;
    plot(x_c',mat','-','Color',[0.7 0.7 0.7 0.4],'LineWidth',0.5);
    for i = 1:n_grps
        d = mat(:,i); d = d(~isnan(d)); if isempty(d), continue; end
        [f,xi] = ksdensity(d); f = f/max(f)*0.4;
        patch([i+f, i*ones(1,length(f))],[xi,fliplr(xi)],cols{i},'EdgeColor','none','FaceAlpha',0.5);
        scatter(x_c(:,i),mat(:,i),25,cols{i},'filled','MarkerFaceAlpha',0.6);
        q = quantile(d,[0.25,0.5,0.75]);
        rectangle('Position',[i-0.06,q(1),0.12,q(3)-q(1)],'FaceColor',cols{i},'EdgeColor','k','LineWidth',1.2);
        plot([i-0.06,i+0.06],[q(2),q(2)],'k-','LineWidth',2);
    end
    set(gca,'XTick',1:n_grps,'XTickLabel',xlbls,'FontSize',14);
    if ~isempty(ttl), title(ttl,'FontSize',14); end
    ylabel(ylbl,'FontSize',14,'FontWeight','bold'); xlim([0.2,n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); else, ylim(auto_lim); end
    grid off; box off; hold off;
end

function add_sig(data, pairs)
    cl = ylim; rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng*0.08; line_lvl = 0; hold on;
    base = max(cl(2), max(data(:))+step*0.5);
    for i = 1:size(pairs,1)
        c1 = pairs(i,1); c2 = pairs(i,2);
        if c2==0, [~,p]=ttest(data(:,c1)); paired=false;
        else,     [~,p]=ttest(data(:,c1),data(:,c2)); paired=true; end
        if p < 0.05
            line_lvl = line_lvl+1; y_p = base+line_lvl*step;
            if p<0.001, txt='***'; elseif p<0.01, txt='**'; else, txt='*'; end
            if paired
                plot([c1,c1,c2,c2],[y_p-step*0.3,y_p,y_p,y_p-step*0.3],'k-','LineWidth',1.2);
                text(mean([c1 c2]),y_p+step*0.1,txt,'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            else
                text(c1,y_p,txt,'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
            end
        end
    end
    if line_lvl>0, ylim([cl(1),base+line_lvl*step+step]); end
    hold off;
end

function draw_matrix(mat_mean, mat_se, cols, ylbl, xlbl)
    hold on;
    for r = 1:3
        for c = 1:3
            v = mat_mean(r,c); se = mat_se(r,c);
            t_c = (v>0.5)*[1 1 1]+(v<=0.5)*[0.2 0.2 0.2];
            patch([c-0.48 c+0.48 c+0.48 c-0.48],[r-0.48 r-0.48 r+0.48 r+0.48], ...
                cols{r},'EdgeColor','none','FaceAlpha',v);
            text(c, r-0.12, sprintf('%.2f',v), 'HorizontalAlignment','center', ...
                'Color',t_c,'FontWeight','bold','FontSize',10);
            text(c, r+0.18, sprintf('±%.2f',se), 'HorizontalAlignment','center', ...
                'Color',t_c,'FontSize',8);
            if r==1, text(c,0.4,xlbl{c},'HorizontalAlignment','center', ...
                'Color',cols{c},'FontWeight','bold','FontSize',10); end
        end
        text(0.4,r,ylbl{r},'HorizontalAlignment','right', ...
            'Color',cols{r},'FontWeight','bold','FontSize',10);
    end
    axis ij equal; xlim([0.2 3.8]); ylim([0.3 3.5]);
    set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none'); hold off;
end

function scatter_with_fit(x, y, col, xlbl, ylbl, ttl)
    valid = ~isnan(x) & ~isnan(y); x = x(valid); y = y(valid); hold on;
    scatter(x,y,55,col,'filled','MarkerFaceAlpha',0.65);
    p_fit = polyfit(x,y,1); x_fit = linspace(min(x),max(x),100);
    plot(x_fit,polyval(p_fit,x_fit),'Color',col,'LineWidth',2);
    xlabel(xlbl,'FontSize',12); ylabel(ylbl,'FontSize',12);
    title(ttl,'FontSize',11); grid off; box off; hold off;
end