% clear; clc; close all;
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %% setup
% %%%%%%%%%%%%%%%%%%%%%%%
% subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617];
% base_dir = '..'; 
% res_dir = fullfile(base_dir, 'results');
% fig_dir = fullfile(base_dir, 'figures');
% min_rt = 0.150;
% % colors
% c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
% c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %% load & compute subject-level stats
% %%%%%%%%%%%%%%%%%%%%%%%
% for s = 1:length(subj_ids)
%     curr_id = subj_ids(s);
%     fprintf('processing subject %d...\n', curr_id);
%     subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', curr_id));
%     load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', curr_id)), 'final_data_output');
%     r1 = final_data_output.results_1_back_all;
%     r2 = final_data_output.results_2_back_all;
%     stats = [];
% 
%     %% 1-back stats
%     r1.resp_key = cellstr(r1.resp_key); r1.resp_key(strcmp(r1.resp_key,'NA')) = {'none'};
%     r1.correct = strcmp(cellstr(r1.corr_resp), r1.resp_key);
%     idx_sim = r1.condition=="compared" & r1.identity=="B";
%     idx_sam = r1.condition=="repeat" & strcmp(r1.corr_resp,'j');
%     idx_cr = (r1.condition=="compared"&r1.identity=="A")|r1.condition=="isolated"|(r1.condition=="repeat"&strcmp(r1.corr_resp,'none'));
%     stats.one.acc_sim = mean(r1.correct(idx_sim));
%     stats.one.acc_same = mean(r1.correct(idx_sam));
%     stats.one.acc_new = mean(r1.correct(idx_cr));
%     v_rt = r1.rt > min_rt;
%     stats.one.rt_sim = median(r1.rt(idx_sim & r1.correct & v_rt), 'omitnan');
%     stats.one.rt_same = median(r1.rt(idx_sam & r1.correct & v_rt), 'omitnan');
%     n_sim = sum(idx_sim); n_sam = sum(idx_sam); n_cr = sum(idx_cr);
%     stats.one.err_sim_as_same = sum(strcmp(r1.resp_key(idx_sim),'j'))/n_sim;
%     stats.one.err_sim_as_new = sum(strcmp(r1.resp_key(idx_sim),'none'))/n_sim;
%     stats.one.err_same_as_sim = sum(strcmp(r1.resp_key(idx_sam),'k'))/n_sam;
%     stats.one.err_same_as_new = sum(strcmp(r1.resp_key(idx_sam),'none'))/n_sam;
%     stats.one.err_new_as_same = sum(strcmp(r1.resp_key(idx_cr),'j'))/n_cr;
%     stats.one.err_new_as_sim = sum(strcmp(r1.resp_key(idx_cr),'k'))/n_cr;
% 
%     %% 2-back stats
%     r2.resp_key = cellstr(r2.resp_key); r2.resp_key(strcmp(r2.resp_key,'NA')) = {'none'};
%     r2.correct = strcmp(cellstr(r2.corr_resp), r2.resp_key);
%     pan = zeros(height(r2),1); 
%     for i=1:height(r2)-2, if strcmp(r2.goal(i),'A-N'), pan(i+2)=1; end; end 
%     real = ~contains(r2.goal, "JUNK"); v_rt = r2.rt > min_rt;
%     aa_idx = real & strcmp(r2.goal,'A-A'); ab_idx = real & strcmp(r2.goal,'A-B'); an_idx = (pan==1);
%     comp_idx = real & strcmp(r2.condition,'compared'); iso_idx = real & strcmp(r2.condition,'isolated'); nov_idx = real & strcmp(r2.condition,'novel');
%     j_idx = real & strcmp(r2.corr_resp,'j'); k_idx = real & strcmp(r2.corr_resp,'k'); n_idx = real & strcmp(r2.corr_resp, 'none');
%     calc_d = @(h,f,nh,nf) norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - norminv(max(1/(2*nf), min(1-1/(2*nf), f)));
%     stats.two.acc_AA_comp = mean(r2.correct(aa_idx & comp_idx & j_idx));
%     stats.two.acc_AA_iso = mean(r2.correct(aa_idx & iso_idx & j_idx));
%     stats.two.acc_AA_nov = mean(r2.correct(aa_idx & nov_idx & j_idx));
%     stats.two.acc_AB_comp = mean(r2.correct(ab_idx & comp_idx & k_idx));
%     stats.two.acc_AB_iso = mean(r2.correct(ab_idx & iso_idx & k_idx));
%     stats.two.acc_AB_nov = mean(r2.correct(ab_idx & nov_idx & k_idx));
%     stats.two.acc_AN_comp = mean(r2.correct(an_idx & comp_idx & n_idx));
%     stats.two.acc_AN_iso = mean(r2.correct(an_idx & iso_idx & n_idx));
%     stats.two.acc_AN_nov = mean(r2.correct(an_idx & nov_idx & n_idx));
%     stats.two.rt_AA_comp = median(r2.rt(aa_idx & comp_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
%     stats.two.rt_AA_iso = median(r2.rt(aa_idx & iso_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
%     stats.two.rt_AA_nov = median(r2.rt(aa_idx & nov_idx & j_idx & r2.correct==1 & v_rt), 'omitnan');
%     stats.two.rt_AB_comp = median(r2.rt(ab_idx & comp_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');
%     stats.two.rt_AB_iso = median(r2.rt(ab_idx & iso_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');
%     stats.two.rt_AB_nov = median(r2.rt(ab_idx & nov_idx & k_idx & r2.correct==1 & v_rt), 'omitnan');
%     n_AA_comp = sum(aa_idx & comp_idx & j_idx); n_AA_iso = sum(aa_idx & iso_idx & j_idx); n_AA_nov = sum(aa_idx & nov_idx & j_idx);
%     n_AB_comp = sum(ab_idx & comp_idx & k_idx); n_AB_iso = sum(ab_idx & iso_idx & k_idx); n_AB_nov = sum(ab_idx & nov_idx & k_idx);
%     n_AN_comp = sum(an_idx & comp_idx & n_idx); n_AN_iso = sum(an_idx & iso_idx & n_idx); n_AN_nov = sum(an_idx & nov_idx & n_idx);
%     stats.two.err_AA_comp_as_k = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'k')) / n_AA_comp;
%     stats.two.err_AA_comp_as_n = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'none')) / n_AA_comp;
%     stats.two.err_AB_comp_as_j = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'j')) / n_AB_comp;
%     stats.two.err_AB_comp_as_n = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'none')) / n_AB_comp;
%     stats.two.err_AN_comp_as_j = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'j')) / n_AN_comp;
%     stats.two.err_AN_comp_as_k = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'k')) / n_AN_comp;
%     stats.two.err_AA_iso_as_k = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'k')) / n_AA_iso;
%     stats.two.err_AA_iso_as_n = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'none')) / n_AA_iso;
%     stats.two.err_AB_iso_as_j = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'j')) / n_AB_iso;
%     stats.two.err_AB_iso_as_n = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'none')) / n_AB_iso;
%     stats.two.err_AN_iso_as_j = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'j')) / n_AN_iso;
%     stats.two.err_AN_iso_as_k = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'k')) / n_AN_iso;
%     stats.two.err_AA_nov_as_k = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'k')) / n_AA_nov;
%     stats.two.err_AA_nov_as_n = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'none')) / n_AA_nov;
%     stats.two.err_AB_nov_as_j = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'j')) / n_AB_nov;
%     stats.two.err_AB_nov_as_n = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'none')) / n_AB_nov;
%     stats.two.err_AN_nov_as_j = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'j')) / n_AN_nov;
%     stats.two.err_AN_nov_as_k = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'k')) / n_AN_nov;
%     stats.two.ldi_comp = stats.two.acc_AB_comp - stats.two.err_AN_comp_as_k;
%     stats.two.ldi_iso = stats.two.acc_AB_iso - stats.two.err_AN_iso_as_k;
%     stats.two.ldi_nov = stats.two.acc_AB_nov - stats.two.err_AN_nov_as_k;
%     stats.two.dprime_comp = calc_d(stats.two.acc_AA_comp, stats.two.err_AN_comp_as_j, n_AA_comp, n_AN_comp);
%     stats.two.dprime_iso = calc_d(stats.two.acc_AA_iso, stats.two.err_AN_iso_as_j, n_AA_iso, n_AN_iso);
%     stats.two.dprime_nov = calc_d(stats.two.acc_AA_nov, stats.two.err_AN_nov_as_j, n_AA_nov, n_AN_nov);
% 
%     %% recognition stats
%     if isfield(final_data_output, 'results_recognition')
%         rec = final_data_output.results_recognition;
%         rec.correct = strcmp(cellstr(rec.corr_resp), cellstr(rec.resp_key));
%         old = rec(rec.trial_type=="old",:); new = rec(rec.trial_type~="old",:);
%         n_new = height(new); n_fa = sum(strcmp(cellstr(new.resp_key),'j') & new.rt>min_rt);
%         far = max(1/(2*n_new), min(1-1/(2*n_new), n_fa/n_new));
%         tc = old(old.condition=="compared",:); nc = height(tc);
%         hc = sum(tc.correct & tc.rt>min_rt)/nc;
%         stats.rec.d_comp = calc_d(hc, far, nc, n_new);
%         ti = old(old.condition=="isolated",:); ni = height(ti);
%         hi = sum(ti.correct & ti.rt>min_rt)/ni;
%         stats.rec.d_iso = calc_d(hi, far, ni, n_new);
%     else
%         fprintf('  recognition data missing for subject %d, filling with NaNs.\n', curr_id);
%         stats.rec.d_comp = NaN; stats.rec.d_iso = NaN;
%     end
%     all_subjs(s).id = curr_id; all_subjs(s).stats = stats;
% end
% save(fullfile(res_dir, 'all_subjs_stats.mat'), 'all_subjs');
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %% extract group-level vectors
% %%%%%%%%%%%%%%%%%%%%%%%
% get_v = @(f1, f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs);
% b1_acc_sam = get_v('one','acc_same')'; b1_acc_sim = get_v('one','acc_sim')'; b1_acc_new = get_v('one','acc_new')';
% b1_rt_sam = get_v('one','rt_same')'; b1_rt_sim = get_v('one','rt_sim')';
% b2_ldi_c = get_v('two','ldi_comp')'; b2_ldi_i = get_v('two','ldi_iso')'; b2_ldi_n = get_v('two','ldi_nov')';
% b2_dp_c = get_v('two','dprime_comp')'; b2_dp_i = get_v('two','dprime_iso')'; b2_dp_n = get_v('two','dprime_nov')';
% b2_rt_l_c = get_v('two','rt_AB_comp')'; b2_rt_l_i = get_v('two','rt_AB_iso')'; b2_rt_l_n = get_v('two','rt_AB_nov')';
% b2_rt_t_c = get_v('two','rt_AA_comp')'; b2_rt_t_i = get_v('two','rt_AA_iso')'; b2_rt_t_n = get_v('two','rt_AA_nov')';
% rec_d_c = get_v('rec','d_comp'); rec_d_i = get_v('rec','d_iso');
% load(fullfile(res_dir, 'all_trials_preprocessed.mat'), 'all_preprocessed');
% [pup, pup_mean, t_pup, pup_subjs] = extract_pupil_subj(all_preprocessed);
% fprintf('pupil: %d subjs, %d samples (%.2fs)\n', length(pup_subjs), length(t_pup), t_pup(end));
% load(fullfile(res_dir, 'gaze_reinstat_res_m.mat'));
% [reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso, match_bb_comp, match_bb_iso, match_ba_comp, match_ba_iso, ...
%  baseline_bb_comp, baseline_bb_iso, baseline_ba_comp, baseline_ba_iso] = extract_gaze_subj(reinstat_res, subj_ids);
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %% statistical tests
% %%%%%%%%%%%%%%%%%%%%%%%
% within = table(categorical({'compared';'isolated';'novel'}), 'VariableNames', {'Condition'});
% stat_results = struct();
% 
% %% 1-back tests
% fprintf('\n=== 1-BACK STATS ===\n');
% [~, stat_results.b1_rt_sam_v_sim.p, ~, stat_results.b1_rt_sam_v_sim.stats] = ttest(b1_rt_sam, b1_rt_sim);
% [~, stat_results.b1_acc_sam_v_sim.p, ~, stat_results.b1_acc_sam_v_sim.stats] = ttest(b1_acc_sam, b1_acc_sim);
% [~, stat_results.b1_acc_sim_v_new.p, ~, stat_results.b1_acc_sim_v_new.stats] = ttest(b1_acc_sim, b1_acc_new);
% fprintf('RT: same v sim: t(%d)=%.2f, p=%.4f\n', stat_results.b1_rt_sam_v_sim.stats.df, stat_results.b1_rt_sam_v_sim.stats.tstat, stat_results.b1_rt_sam_v_sim.p);
% fprintf('Acc: same v sim: t(%d)=%.2f, p=%.4f\n', stat_results.b1_acc_sam_v_sim.stats.df, stat_results.b1_acc_sam_v_sim.stats.tstat, stat_results.b1_acc_sam_v_sim.p);
% fprintf('Acc: sim v new: t(%d)=%.2f, p=%.4f\n', stat_results.b1_acc_sim_v_new.stats.df, stat_results.b1_acc_sim_v_new.stats.tstat, stat_results.b1_acc_sim_v_new.p);
% 
% %% 2-back tests
% fprintf('\n=== 2-BACK STATS ===\n');
% [stat_results.b2_ldi, stat_results.b2_ldi_posthoc] = run_rm_anova('LDI', b2_ldi_c, b2_ldi_i, b2_ldi_n, within);
% [stat_results.b2_dprime, stat_results.b2_dprime_posthoc] = run_rm_anova('d-prime', b2_dp_c, b2_dp_i, b2_dp_n, within);
% [stat_results.b2_rt_lure, stat_results.b2_rt_lure_posthoc] = run_rm_anova('RT lure', b2_rt_l_c, b2_rt_l_i, b2_rt_l_n, within);
% [stat_results.b2_rt_target, stat_results.b2_rt_target_posthoc] = run_rm_anova('RT target', b2_rt_t_c, b2_rt_t_i, b2_rt_t_n, within);
% 
% %% recognition test
% fprintf('\n=== RECOGNITION STATS ===\n');
% d_tot = (rec_d_c + rec_d_i)/2;
% [~, stat_results.rec_overall.p, ~, stat_results.rec_overall.stats] = ttest(d_tot);
% [~, stat_results.rec_comp_v_iso.p, ~, stat_results.rec_comp_v_iso.stats] = ttest(rec_d_c, rec_d_i);
% fprintf('Overall d'': t(%d)=%.2f, p=%.4f\n', stat_results.rec_overall.stats.df, stat_results.rec_overall.stats.tstat, stat_results.rec_overall.p);
% fprintf('Comp v Iso: t(%d)=%.2f, p=%.4f\n', stat_results.rec_comp_v_iso.stats.df, stat_results.rec_comp_v_iso.stats.tstat, stat_results.rec_comp_v_iso.p);
% 
% %% pupil tests
% fprintf('\n=== PUPIL STATS ===\n');
% [stat_results.pup_lure, stat_results.pup_lure_posthoc] = run_rm_anova('Pupil A-B', pup_mean.ab_com, pup_mean.ab_iso, pup_mean.ab_nov, within);
% [stat_results.pup_target, stat_results.pup_target_posthoc] = run_rm_anova('Pupil A-A', pup_mean.aa_com, pup_mean.aa_iso, pup_mean.aa_nov, within);
% 
% fprintf('\n=== PUPIL TIME SERIES (CLUSTER PERM) ===\n');
% t_idx = t_pup >= 0 & t_pup <= 1.5; t_clust = t_pup(t_idx); n_perm = 1000;
% fprintf('A-B trials:\n');
% [stat_results.pup_ts_ab_ci.clusters, stat_results.pup_ts_ab_ci.p, stat_results.pup_ts_ab_ci.t] = cluster_perm_ttest(pup.ab_com(:,t_idx), pup.ab_iso(:,t_idx), n_perm);
% print_clusters('  comp vs iso', stat_results.pup_ts_ab_ci.clusters, stat_results.pup_ts_ab_ci.p, t_clust);
% [stat_results.pup_ts_ab_in.clusters, stat_results.pup_ts_ab_in.p, stat_results.pup_ts_ab_in.t] = cluster_perm_ttest(pup.ab_iso(:,t_idx), pup.ab_nov(:,t_idx), n_perm);
% print_clusters('  iso vs nov', stat_results.pup_ts_ab_in.clusters, stat_results.pup_ts_ab_in.p, t_clust);
% [stat_results.pup_ts_ab_cn.clusters, stat_results.pup_ts_ab_cn.p, stat_results.pup_ts_ab_cn.t] = cluster_perm_ttest(pup.ab_com(:,t_idx), pup.ab_nov(:,t_idx), n_perm);
% print_clusters('  comp vs nov', stat_results.pup_ts_ab_cn.clusters, stat_results.pup_ts_ab_cn.p, t_clust);
% fprintf('A-A trials:\n');
% [stat_results.pup_ts_aa_ci.clusters, stat_results.pup_ts_aa_ci.p, stat_results.pup_ts_aa_ci.t] = cluster_perm_ttest(pup.aa_com(:,t_idx), pup.aa_iso(:,t_idx), n_perm);
% print_clusters('  comp vs iso', stat_results.pup_ts_aa_ci.clusters, stat_results.pup_ts_aa_ci.p, t_clust);
% [stat_results.pup_ts_aa_in.clusters, stat_results.pup_ts_aa_in.p, stat_results.pup_ts_aa_in.t] = cluster_perm_ttest(pup.aa_iso(:,t_idx), pup.aa_nov(:,t_idx), n_perm);
% print_clusters('  iso vs nov', stat_results.pup_ts_aa_in.clusters, stat_results.pup_ts_aa_in.p, t_clust);
% [stat_results.pup_ts_aa_cn.clusters, stat_results.pup_ts_aa_cn.p, stat_results.pup_ts_aa_cn.t] = cluster_perm_ttest(pup.aa_com(:,t_idx), pup.aa_nov(:,t_idx), n_perm);
% print_clusters('  comp vs nov', stat_results.pup_ts_aa_cn.clusters, stat_results.pup_ts_aa_cn.p, t_clust);
% 
% %% gaze reinstatement tests
% fprintf('\n=== GAZE REINSTATEMENT STATS ===\n');
% stat_results.gaze_2x2 = run_2x2_anova(reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso);
% 
% fprintf('\n=== GAZE PERMUTATION (MATCH VS MISMATCH) ===\n');
% [stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_comp.obs] = run_permutation(match_bb_comp, baseline_bb_comp, n_perm);
% [stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_bb_iso.obs] = run_permutation(match_bb_iso, baseline_bb_iso, n_perm);
% [stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_comp.obs] = run_permutation(match_ba_comp, baseline_ba_comp, n_perm);
% [stat_results.gaze_match_ba_iso.p, stat_results.gaze_match_ba_iso.obs] = run_permutation(match_ba_iso, baseline_ba_iso, n_perm);
% fprintf('BB comp: p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_comp.obs);
% fprintf('BB iso:  p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_bb_iso.obs);
% fprintf('BA comp: p=%.4f, obs=%.3f\n', stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_comp.obs);
% fprintf('BA iso:  p=%.4f, obs=%.3f\n', stat_results.gaze_match_ba_iso.p, stat_results.gaze_match_ba_iso.obs);
% 
% %% gaze-behavior correlation
% fprintf('\n=== GAZE-BEHAVIOR CORRELATIONS ===\n');
% valid = ~isnan(reinst_bb_comp) & ~isnan(b2_ldi_c) & b2_ldi_c <= 0.8;
% [stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p] = corr(reinst_bb_comp(valid), b2_ldi_c(valid));
% stat_results.corr_gr_ldi_comp.n = sum(valid);
% fprintf('Comp: r=%.3f, p=%.4f, n=%d\n', stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p, stat_results.corr_gr_ldi_comp.n);
% valid = ~isnan(reinst_bb_iso) & ~isnan(b2_ldi_i) & b2_ldi_i <= 0.8;
% [stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p] = corr(reinst_bb_iso(valid), b2_ldi_i(valid));
% stat_results.corr_gr_ldi_iso.n = sum(valid);
% fprintf('Iso:  r=%.3f, p=%.4f, n=%d\n', stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p, stat_results.corr_gr_ldi_iso.n);

%%%%%%%%%%%%%%%%%%%%%%%
%% visualization
%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('\n=== GENERATING FIGURES ===\n');
% 
% %% fig: 1-back
% figure('color','w','Position',[100 100 1200 400]);
% subplot(1,3,1); 
% data = [b1_acc_sam, b1_acc_sim, b1_acc_new];
% raincloud(data, {c_same, c_sim, c_new}, {'Hit(Same)','Hit(Sim)','CR(New)'}, 'Accuracy', '', [0,1]);
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(1,3,2); 
% data = [b1_rt_sam, b1_rt_sim];
% raincloud(data, {c_same, c_sim}, {'Hit(Same)','Hit(Sim)'}, 'RT (s)', '');
% add_sig(data, [1 2]);
% subplot(1,3,3); hold on;
% mat_1_back = [mean(get_v('one','acc_same')), mean(get_v('one','err_same_as_sim')), mean(get_v('one','err_same_as_new'));
%               mean(get_v('one','err_sim_as_same')), mean(get_v('one','acc_sim')), mean(get_v('one','err_sim_as_new'));
%               mean(get_v('one','err_new_as_same')), mean(get_v('one','err_new_as_sim')), mean(get_v('one','acc_new'))];
% draw_matrix(mat_1_back, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('Response matrix', 'FontSize', 12);
% sgtitle('1-Back Task Performance', 'FontSize', 16);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'behav_1back.pdf'), '-dpdf', '-vector');
% 
% %% fig: 2-back
% figure('color','w','Position',[50 50 1200 900]);
% subplot(3,3,1);
% data = [b2_ldi_c, b2_ldi_i, b2_ldi_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'LDI', 'Lure Discrimination', [0,1]);
% set(gca, 'YTick', [0 0.5 1]); add_sig(data, [1 2; 2 3; 1 3]);
% subplot(3,3,4);
% data = [b2_dp_c, b2_dp_i, b2_dp_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'd''', 'Recognition');
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(3,3,2);
% data = [b2_rt_l_c, b2_rt_l_i, b2_rt_l_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Lure Discrimination)', [0,1.5]);
% add_sig(data, [1 2; 2 3; 1 3]); set(gca, 'YTick', [0 0.5 1 1.5])
% subplot(3,3,5);
% data = [b2_rt_t_c, b2_rt_t_i, b2_rt_t_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Recognition)');
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(3,3,7); hold on;
% mat_comp = [mean(get_v('two','acc_AA_comp')), mean(get_v('two','err_AA_comp_as_k')), mean(get_v('two','err_AA_comp_as_n'));
%             mean(get_v('two','err_AB_comp_as_j')), mean(get_v('two','acc_AB_comp')), mean(get_v('two','err_AB_comp_as_n'));
%             mean(get_v('two','err_AN_comp_as_j')), mean(get_v('two','err_AN_comp_as_k')), mean(get_v('two','acc_AN_comp'))];
% draw_matrix(mat_comp, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('compared', 'FontSize', 14);
% subplot(3,3,8); hold on;
% mat_iso = [mean(get_v('two','acc_AA_iso')), mean(get_v('two','err_AA_iso_as_k')), mean(get_v('two','err_AA_iso_as_n'));
%            mean(get_v('two','err_AB_iso_as_j')), mean(get_v('two','acc_AB_iso')), mean(get_v('two','err_AB_iso_as_n'));
%            mean(get_v('two','err_AN_iso_as_j')), mean(get_v('two','err_AN_iso_as_k')), mean(get_v('two','acc_AN_iso'))];
% draw_matrix(mat_iso, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('isolated', 'FontSize', 14);
% subplot(3,3,9); hold on;
% mat_nov = [mean(get_v('two','acc_AA_nov')), mean(get_v('two','err_AA_nov_as_k')), mean(get_v('two','err_AA_nov_as_n'));
%            mean(get_v('two','err_AB_nov_as_j')), mean(get_v('two','acc_AB_nov')), mean(get_v('two','err_AB_nov_as_n'));
%            mean(get_v('two','err_AN_nov_as_j')), mean(get_v('two','err_AN_nov_as_k')), mean(get_v('two','acc_AN_nov'))];
% draw_matrix(mat_nov, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('novel', 'FontSize', 14);
% sgtitle('2-Back Task Performance', 'FontSize', 16);
% set(gcf, 'Color', 'w'); set(gcf, 'Renderer', 'painters'); set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'behav_2back.pdf'), '-dpdf', '-vector');
% 
% %% fig: recognition
% figure('color','w','Position',[100 100 600 500]);
% data = [d_tot', rec_d_c', rec_d_i'];
% hold on;
% fill([-2, 5, 5, -2], [-2, -2, 0, 0], [0.92 0.92 0.92], 'EdgeColor', 'none'); 
% yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
% raincloud(data, {[0.3 0.3 0.3], c_comp, c_iso}, {'overall','compared','isolated'}, 'd''', 'Validation of Episodic Encoding');
% add_sig(data, [1 0; 2 3]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'behav_recog.pdf'), '-dpdf', '-vector');
% 
% %% fig: gaze match vs mismatch
% figure('color','w','position',[50 50 1000 800]);
% subplot(2,2,1);
% data = [match_bb_comp, baseline_bb_comp];
% raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'B-B compared', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_comp.p);
% subplot(2,2,2);
% data = [match_bb_iso, baseline_bb_iso];
% raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'B-B isolated', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_iso.p);
% subplot(2,2,3);
% data = [match_ba_comp, baseline_ba_comp];
% raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'B-A compared (predictive)', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_ba_comp.p);
% subplot(2,2,4);
% data = [match_ba_iso, baseline_ba_iso];
% raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'B-A isolated (predictive)', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_ba_iso.p);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reins_mmis.pdf'), '-dpdf', '-vector');
% 
% %% fig: gaze permutation distributions
% figure('color','w','position',[50 50 1000 800]);
% perm_titles = {'B-B Compared', 'B-B Isolated', 'B-A Compared', 'B-A Isolated'};
% perm_data_A = {match_bb_comp, match_bb_iso, match_ba_comp, match_ba_iso};
% perm_data_B = {baseline_bb_comp, baseline_bb_iso, baseline_ba_comp, baseline_ba_iso};
% perm_colors = {c_comp, c_iso, c_comp, c_iso};
% for p = 1:4
%     subplot(2,2,p);
%     [p_val, obs_diff, null_dist] = run_permutation_with_dist(perm_data_A{p}, perm_data_B{p}, 1000);
%     histogram(null_dist, 30, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'w', 'FaceAlpha', 0.7); hold on;
%     yl = [0 100]; set(gca, 'YTick', [0 25 50 75 100]);
%     line([obs_diff obs_diff], [0 yl(2)], 'Color', perm_colors{p}, 'LineWidth', 2.5);
%     title(perm_titles{p}, 'FontSize', 14); xlabel('Gaze Reinstatement Index'); ylabel('Frequency');
%     legend('null distribution', 'ground-truth data', 'Location', 'northwest');
%     text(obs_diff, yl(2)*0.8, sprintf(' p = %.3f', p_val), 'Color', perm_colors{p}, 'FontWeight', 'bold');
%     grid off; box off;
% end
% sgtitle('Permutation Tests', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reins_perm.pdf'), '-dpdf', '-vector');

%% fig: gaze reinstatement index
figure('color','w','position',[50 50 800 600]);
data_matrix = [reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso];
labels = {'comp B-B', 'iso B-B', 'comp B-A', 'iso B-A'};
colors = {c_comp, c_iso, c_comp, c_iso}; hold on;
yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
raincloud(data_matrix, colors, labels, 'Gaze Reinstatement Index', 'Gaze Reinstatement', [-0.1 0.4]);
set(gca, 'YTick', [-0.1 0 0.1 0.2 0.3 0.4]);
p_chance_all = [stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_iso.p];
for i = 1:4
    p = p_chance_all(i);
    if p < 0.05
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
        text(i, 0.45, txt, 'HorizontalAlignment', 'center', 'FontSize', 18);
    end
end
pairs_to_test = [1 2; 3 4; 1 3]; 
pvals = [stat_results.gaze_2x2.posthoc.bb_ci.p_adj; stat_results.gaze_2x2.posthoc.ba_ci.p_adj; stat_results.gaze_2x2.posthoc.comp_bb_ba.p_adj];
add_sig_perm(data_matrix, pairs_to_test, pvals);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'gaze_reins_idx.svg'), '-dsvg', '-vector');

% %% fig: gaze-behavior correlation
% figure('color','w','position',[50 50 800 400]);
% subplot(1,2,1); hold on;
% valid = ~isnan(reinst_bb_comp) & ~isnan(b2_ldi_c) & b2_ldi_c <= 0.8;
% x = reinst_bb_comp(valid); y = b2_ldi_c(valid);
% scatter(x, y, 60, c_comp, 'filled', 'MarkerFaceAlpha', 0.6);
% p_fit = polyfit(x, y, 1); x_fit = linspace(min(x), max(x), 100);
% plot(x_fit, polyval(p_fit, x_fit), 'Color', c_comp, 'LineWidth', 2);
% xlabel('Gaze Reinstatement Index', 'FontSize', 12); ylabel('LDI', 'FontSize', 12);
% title(sprintf('compared: r=%.2f, p=%.3f (n=%d)', stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p, stat_results.corr_gr_ldi_comp.n), 'FontSize', 14);
% grid off; box off;
% subplot(1,2,2); hold on;
% valid = ~isnan(reinst_bb_iso) & ~isnan(b2_ldi_i) & b2_ldi_i <= 0.8;
% x = reinst_bb_iso(valid); y = b2_ldi_i(valid);
% scatter(x, y, 60, c_iso, 'filled', 'MarkerFaceAlpha', 0.6);
% p_fit = polyfit(x, y, 1); x_fit = linspace(min(x), max(x), 100);
% plot(x_fit, polyval(p_fit, x_fit), 'Color', c_iso, 'LineWidth', 2);
% xlabel('Gaze Reinstatement Index', 'FontSize', 12); ylabel('LDI', 'FontSize', 12);
% title(sprintf('isolated: r=%.2f, p=%.3f (n=%d)', stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p, stat_results.corr_gr_ldi_iso.n), 'FontSize', 14);
% grid off; box off;
% sgtitle('Predict LDI from Gaze Reinstatement', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reinstat_ldi.pdf'), '-dpdf', '-vector');
% 
% %% fig: pupil timeseries
% figure('color','w','position',[50 50 1200 500]);
% subplot(1,2,1); hold on;
% plot_pup_ts(t_pup, pup.ab_com, c_comp, 'compared');
% plot_pup_ts(t_pup, pup.ab_iso, c_iso, 'isolated');
% plot_pup_ts(t_pup, pup.ab_nov, c_nov, 'novel');
% xlabel('Time from stimulus onset (s)', 'FontSize', 14); ylabel('Pupil size change (a.u.)', 'FontSize', 14);
% title('A-B Lure Discrimination', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
% subplot(1,2,2); hold on;
% plot_pup_ts(t_pup, pup.aa_com, c_comp, 'compared');
% plot_pup_ts(t_pup, pup.aa_iso, c_iso, 'isolated');
% plot_pup_ts(t_pup, pup.aa_nov, c_nov, 'novel');
% xlabel('Time from stimulus onset (s)', 'FontSize', 14); ylabel('Pupil size change (a.u.)', 'FontSize', 14);
% title('A-A Target Detection', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
% sgtitle('Pupil Responses over Time', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'pupil_ts.pdf'), '-dpdf', '-vector');
% 
% %% fig: pupil mean
% figure('color','w','position',[50 50 1000 400]);
% subplot(1,2,1);
% data = [pup_mean.ab_com, pup_mean.ab_iso, pup_mean.ab_nov];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'Mean Pupil (a.u.)', 'A-B Lure Discrimination');
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(1,2,2);
% data = [pup_mean.aa_com, pup_mean.aa_iso, pup_mean.aa_nov];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'Mean Pupil (a.u.)', 'A-A Target Detection');
% add_sig(data, [1 2; 2 3; 1 3]);
% sgtitle('Pupil Dilation by Condition', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'pupil_mean.pdf'), '-dpdf', '-vector');

%%%%%%%%%%%%%%%%%%%%%%%
%% save comprehensive results
%%%%%%%%%%%%%%%%%%%%%%%
% save(fullfile(res_dir, 'all_subj_results.mat'), ...
%     'subj_ids', 'all_subjs', ...
%     'b1_acc_sam', 'b1_acc_sim', 'b1_acc_new', 'b1_rt_sam', 'b1_rt_sim', ...
%     'b2_ldi_c', 'b2_ldi_i', 'b2_ldi_n', 'b2_dp_c', 'b2_dp_i', 'b2_dp_n', ...
%     'b2_rt_l_c', 'b2_rt_l_i', 'b2_rt_l_n', 'b2_rt_t_c', 'b2_rt_t_i', 'b2_rt_t_n', ...
%     'rec_d_c', 'rec_d_i', 'd_tot', ...
%     'reinst_bb_comp', 'reinst_bb_iso', 'reinst_ba_comp', 'reinst_ba_iso', ...
%     'match_bb_comp', 'match_bb_iso', 'match_ba_comp', 'match_ba_iso', ...
%     'baseline_bb_comp', 'baseline_bb_iso', 'baseline_ba_comp', 'baseline_ba_iso', ...
%     'pup', 'pup_mean', 't_pup', 'pup_subjs', 't_clust', ...
%     'stat_results');
% fprintf('\nALL RESULTS SAVED: %s\n', fullfile(res_dir, 'all_subj_results.mat'));

%%%%%%%%%%%%%%%%%%%%%%%
%% functions
%%%%%%%%%%%%%%%%%%%%%%%
function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); hold on;
    d_v = mat(:); d_v = d_v(~isnan(d_v)); 
    mn = min(d_v); mx = max(d_v); rng = mx-mn; if rng==0, rng=1; end
    auto_lim = [mn-(rng*0.15), mx+(rng*0.15)];
    jit = -0.15 - (rand(size(mat)) * 0.20); x_c = repmat(1:n_grps, n_rows, 1) + jit;
    plot(x_c', mat', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);
    for i = 1:n_grps
        d = mat(:,i); d = d(~isnan(d)); if isempty(d), continue; end
        [f, xi] = ksdensity(d); f = f/max(f)*0.4;
        patch([i+f, i*ones(1,length(f))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',0.5);
        scatter(x_c(:,i), mat(:,i), 20, cols{i}, 'filled', 'MarkerFaceAlpha',0.6);
        q = quantile(d, [0.25, 0.5, 0.75]);
        rectangle('Position', [i-0.06, q(1), 0.12, q(3)-q(1)], 'FaceColor', cols{i}, 'EdgeColor','k', 'LineWidth',1.2);
        plot([i-0.06, i+0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 12);
    if ~isempty(ttl), title(ttl, 'FontSize', 14); end
    ylabel(ylbl,'FontSize',14,'FontWeight','bold'); xlim([0.2, n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); else, ylim(auto_lim); end
    grid off; set(gca,'GridAlpha',0.1); box off; hold off;
end

function add_sig(data, pairs)
    [~, n_grps] = size(data);
    cl = ylim; y_top = cl(2); rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng * 0.08; line_lvl = 0; hold on;
    real_mx = max(data(:)); if isnan(real_mx), real_mx=y_top; end
    base = max(y_top, real_mx + step*0.5);
    for i = 1:size(pairs,1)
        c1 = pairs(i,1); c2 = pairs(i,2);
        if c2==0, [~, p]=ttest(data(:,c1)); paired=false; else, [~, p]=ttest(data(:,c1),data(:,c2)); paired=true; end
        if p < 0.05
            line_lvl=line_lvl+1; y_p = base + (line_lvl*step);
            if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
            if paired, plot([c1,c1,c2,c2],[y_p-step*0.3,y_p,y_p,y_p-step*0.3],'k-','LineWidth',1.2);
               text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            else, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold'); end
        end
    end
    if line_lvl>0, ylim([cl(1), base+(line_lvl*step)+step]); end; hold off;
end

function draw_matrix(mat, cols, ylbl, xlbl)
    hold on;
    for r=1:3
        for c=1:3
            v = mat(r,c); t_c = (v>0.5)*[1 1 1] + (v<=0.5)*[0.2 0.2 0.2];           
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], cols{r}, 'EdgeColor','none','FaceAlpha',v);
            text(c, r, sprintf('%.2f',v), 'HorizontalAlignment','center', 'Color',t_c, 'FontWeight','bold', 'FontSize',10);
            if r==1, text(c, 0.4, xlbl{c}, 'HorizontalAlignment','center', 'Color',cols{c}, 'FontWeight','bold', 'FontSize',10); end
        end
        text(0.4, r, ylbl{r}, 'HorizontalAlignment','right', 'Color',cols{r}, 'FontWeight','bold', 'FontSize',10);
    end 
    axis ij equal; xlim([0.2 3.8]); ylim([0.3 3.5]); set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); hold off;
end

function [p_val, obs_diff] = run_permutation(cond_A, cond_B, n_perms)
    valid = ~isnan(cond_A) & ~isnan(cond_B);
    diffs = cond_A(valid) - cond_B(valid);
    obs_diff = mean(diffs); null_dist = zeros(n_perms, 1);
    for i = 1:n_perms
        signs = sign(rand(size(diffs)) - 0.5);
        null_dist(i) = mean(diffs .* signs);
    end
    p_val = mean(abs(null_dist) >= abs(obs_diff));
end

function add_sig_perm(data, pairs, pvals)
    cl = ylim; rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng * 0.08; line_lvl = 0; hold on;
    base = max(cl(2), max(data(:)) + step*0.5);
    for i = 1:size(pairs,1)
        p = pvals(i); c1 = pairs(i,1); c2 = pairs(i,2);
        if p >= 0.05, continue; end
        line_lvl = line_lvl+1; y_p = base + (line_lvl*step);
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
        if c2 == 0, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
        else, plot([c1 c1 c2 c2], [y_p-step*0.3 y_p y_p y_p-step*0.3], 'k-', 'LineWidth', 1.2);
              text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold'); end
    end
    if line_lvl>0, ylim([cl(1), base+(line_lvl*step)+step]); end; hold off;
end

function [p_val, obs_diff, null_dist] = run_permutation_with_dist(cond_A, cond_B, n_perms)
    valid = ~isnan(cond_A) & ~isnan(cond_B);
    diffs = cond_A(valid) - cond_B(valid);
    obs_diff = mean(diffs); null_dist = zeros(n_perms, 1);
    for i = 1:n_perms
        signs = sign(rand(size(diffs)) - 0.5);
        null_dist(i) = mean(diffs .* signs);
    end
    p_val = mean(abs(null_dist) >= abs(obs_diff));
end

function a = avg_pup_cond(sv, cond, goal, maxs)
    tr = sv(strcmp(sv.condition,cond) & strcmp(sv.goal,goal), :);
    if height(tr)==0, a=nan(1,maxs); return; end
    m = nan(height(tr), maxs);
    for i = 1:height(tr)
        p = tr.pupil_preprocessed{i} - mean(tr.baseline_pupil_preprocessed{i},'omitnan');
        m(i,1:length(p)) = p';
    end
    a = mean(m,1,'omitnan');
end

function plot_pup_ts(t, d, col, lbl)
    mu = mean(d,1,'omitnan'); nv = sum(~isnan(d(:,1)));
    if nv==0, return; end
    se = std(d,0,1,'omitnan')/sqrt(nv);
    plot(t, mu, 'Color', col, 'LineWidth', 2.5, 'DisplayName', lbl);
    vi = ~isnan(mu);
    if sum(vi)>1
        fill([t(vi) fliplr(t(vi))], [mu(vi)+se(vi) fliplr(mu(vi)-se(vi))], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function [res, posthoc] = run_rm_anova(lbl, d1, d2, d3, within)
    t = table(d1, d2, d3, 'VariableNames', {'compared','isolated','novel'});
    rm = fitrm(t, 'compared-novel~1', 'WithinDesign', within);
    tbl = ranova(rm); eps = epsilon(rm); m = mauchly(rm);
    res.F = tbl.F(1); res.df1 = tbl.DF(1); res.df2_gg = tbl.DF(2)*eps.GreenhouseGeisser;
    res.p_gg = tbl.pValueGG(1); res.mauchly_p = m.pValue; res.eps = eps.GreenhouseGeisser;
    fprintf('\n%s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f\n', lbl, res.df1, res.df2_gg, res.F, res.p_gg, res.mauchly_p, res.eps);
    n = length(d1); subj = (1:n)';
    bf_tbl = table(subj, d1, d2, d3, 'VariableNames', {'subj','c1','c2','c3'});
    bf_long = stack(bf_tbl, {'c1','c2','c3'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj); bf_long.cond = categorical(bf_long.cond);
    res.bf10 = bf.anova(bf_long, 'y~cond', 'treatAsRandom', {'subj'}, 'verbose', false);
    if res.bf10<1, fprintf('  BF01=%.2f (Evidence for Null)\n', 1/res.bf10); else, fprintf('  BF10=%.2f (Evidence for Effect)\n', res.bf10); end
    [~,p1,~,s1] = ttest(d1,d2); [~,p2,~,s2] = ttest(d2,d3); [~,p3,~,s3] = ttest(d1,d3);
    pvals = [p1 p2 p3]; adj = holm_bonf(pvals); df = s1.df;
    posthoc.comp_iso.t = s1.tstat; posthoc.comp_iso.df = df; posthoc.comp_iso.p = p1; posthoc.comp_iso.p_adj = adj(1);
    posthoc.iso_nov.t = s2.tstat; posthoc.iso_nov.df = df; posthoc.iso_nov.p = p2; posthoc.iso_nov.p_adj = adj(2);
    posthoc.comp_nov.t = s3.tstat; posthoc.comp_nov.df = df; posthoc.comp_nov.p = p3; posthoc.comp_nov.p_adj = adj(3);
    fprintf('  post-hoc (holm-bonf):\n');
    fprintf('    comp-iso: t(%d)=%.2f, p=%.4f, p_adj=%.4f %s\n', df, s1.tstat, p1, adj(1), repmat('*',1,(adj(1)<0.05)+(adj(1)<0.01)+(adj(1)<0.001)));
    fprintf('    iso-nov: t(%d)=%.2f, p=%.4f, p_adj=%.4f %s\n', df, s2.tstat, p2, adj(2), repmat('*',1,(adj(2)<0.05)+(adj(2)<0.01)+(adj(2)<0.001)));
    fprintf('    comp-nov: t(%d)=%.2f, p=%.4f, p_adj=%.4f %s\n', df, s3.tstat, p3, adj(3), repmat('*',1,(adj(3)<0.05)+(adj(3)<0.01)+(adj(3)<0.001)));
end

function adj = holm_bonf(p)
    n = length(p); [ps, idx] = sort(p); as = min(1, ps .* (n:-1:1));
    for i = 2:n, as(i) = max(as(i), as(i-1)); end
    adj(idx) = as;
end

function [reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso, match_bb_comp, match_bb_iso, match_ba_comp, match_ba_iso, ...
          baseline_bb_comp, baseline_bb_iso, baseline_ba_comp, baseline_ba_iso] = extract_gaze_subj(reinstat_res, subj_ids)
    n = length(subj_ids);
    reinst_bb_comp = nan(n,1); reinst_bb_iso = nan(n,1); reinst_ba_comp = nan(n,1); reinst_ba_iso = nan(n,1);
    match_bb_comp = nan(n,1); match_bb_iso = nan(n,1); match_ba_comp = nan(n,1); match_ba_iso = nan(n,1);
    baseline_bb_comp = nan(n,1); baseline_bb_iso = nan(n,1); baseline_ba_comp = nan(n,1); baseline_ba_iso = nan(n,1);
    for i = 1:n
        sid = subj_ids(i);
        reinst_bb_comp(i) = mean(reinstat_res.bb_compared.reinst_index(reinstat_res.bb_compared.subj_id==sid),'omitnan');
        reinst_bb_iso(i) = mean(reinstat_res.bb_isolated.reinst_index(reinstat_res.bb_isolated.subj_id==sid),'omitnan');
        reinst_ba_comp(i) = mean(reinstat_res.ba_compared.reinst_index(reinstat_res.ba_compared.subj_id==sid),'omitnan');
        reinst_ba_iso(i) = mean(reinstat_res.ba_isolated.reinst_index(reinstat_res.ba_isolated.subj_id==sid),'omitnan');
        match_bb_comp(i) = mean(reinstat_res.bb_compared.match_score(reinstat_res.bb_compared.subj_id==sid),'omitnan');
        match_bb_iso(i) = mean(reinstat_res.bb_isolated.match_score(reinstat_res.bb_isolated.subj_id==sid),'omitnan');
        match_ba_comp(i) = mean(reinstat_res.ba_compared.match_score(reinstat_res.ba_compared.subj_id==sid),'omitnan');
        match_ba_iso(i) = mean(reinstat_res.ba_isolated.match_score(reinstat_res.ba_isolated.subj_id==sid),'omitnan');
        baseline_bb_comp(i) = mean(reinstat_res.bb_compared.baseline_score(reinstat_res.bb_compared.subj_id==sid),'omitnan');
        baseline_bb_iso(i) = mean(reinstat_res.bb_isolated.baseline_score(reinstat_res.bb_isolated.subj_id==sid),'omitnan');
        baseline_ba_comp(i) = mean(reinstat_res.ba_compared.baseline_score(reinstat_res.ba_compared.subj_id==sid),'omitnan');
        baseline_ba_iso(i) = mean(reinstat_res.ba_isolated.baseline_score(reinstat_res.ba_isolated.subj_id==sid),'omitnan');
    end
end

function res = run_2x2_anova(r_bb_c, r_bb_i, r_ba_c, r_ba_i)
    t = table(r_bb_c, r_bb_i, r_ba_c, r_ba_i, 'VariableNames', {'BBC','BBI','BAC','BAI'});
    within = table({'BB';'BB';'BA';'BA'}, {'Comp';'Iso';'Comp';'Iso'}, 'VariableNames', {'PT','Cond'});
    rm = fitrm(t, 'BBC-BAI~1', 'WithinDesign', within); tbl = ranova(rm, 'WithinModel', 'PT*Cond');
    fprintf('\n=== 2x2 ANOVA (Gaze) ===\n');
    rows = {'PT', 'Cond', 'PT:Cond'}; term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p = tbl.pValue(idx); s_ef = tbl.SumSq(idx); 
            if idx<height(tbl), s_er = tbl.SumSq(idx+1); else, s_er = NaN; end
            eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F=%.2f, p=%.4f, eta=%.2f\n', rows{i}, f, p, eta);
            if i==1, res.main_pt.F=f; res.main_pt.p=p; res.main_pt.eta=eta; end
            if i==2, res.main_cond.F=f; res.main_cond.p=p; res.main_cond.eta=eta; end
            if i==3, res.interact.F=f; res.interact.p=p; res.interact.eta=eta; end
        end
    end
    fprintf('  -- Post-Hoc (Holm-Bonferroni) --\n');
    pairs = {r_bb_c, r_bb_i, 'BB: Comp-Iso'; r_ba_c, r_ba_i, 'BA: Comp-Iso'; r_bb_c, r_ba_c, 'Comp: BB-BA'; r_bb_i, r_ba_i, 'Iso:  BB-BA'};
    p_raw = zeros(1,4); stats = cell(1,4);
    for k=1:4, [~, p_raw(k), ~, stats{k}] = ttest(pairs{k,1}, pairs{k,2}); end
    p_adj = holm_bonf(p_raw);
    res.posthoc.bb_ci.t = stats{1}.tstat; res.posthoc.bb_ci.df = stats{1}.df; res.posthoc.bb_ci.p = p_raw(1); res.posthoc.bb_ci.p_adj = p_adj(1);
    res.posthoc.ba_ci.t = stats{2}.tstat; res.posthoc.ba_ci.df = stats{2}.df; res.posthoc.ba_ci.p = p_raw(2); res.posthoc.ba_ci.p_adj = p_adj(2);
    res.posthoc.comp_bb_ba.t = stats{3}.tstat; res.posthoc.comp_bb_ba.df = stats{3}.df; res.posthoc.comp_bb_ba.p = p_raw(3); res.posthoc.comp_bb_ba.p_adj = p_adj(3);
    res.posthoc.iso_bb_ba.t = stats{4}.tstat; res.posthoc.iso_bb_ba.df = stats{4}.df; res.posthoc.iso_bb_ba.p = p_raw(4); res.posthoc.iso_bb_ba.p_adj = p_adj(4);
    for k=1:4
        s = stats{k}; sig = repmat('*',1, p_adj(k)<0.05);
        fprintf('    %s: t(%d)=%.2f, p=%.4f, p_adj=%.4f %s\n', pairs{k,3}, s.df, s.tstat, p_raw(k), p_adj(k), sig);
    end
end

function [pup, pup_mean, t_pup, pup_subjs] = extract_pupil_subj(all_preprocessed)
    valid = all_preprocessed(all_preprocessed.preprocess_success, :);
    max_samp = max(cellfun(@length, valid.pupil_preprocessed));
    sr = valid.sample_rate(1); t_pup = (0:max_samp-1)/sr;
    pup_subjs = unique(valid.subj_id); n = length(pup_subjs);
    conds = {'compared','isolated','novel'}; goals = {'A-B','A-A'}; gtag = {'ab','aa'};
    pup = struct(); pup_mean = struct();
    for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn) = nan(n, max_samp); end, end
    for s = 1:n
        sv = valid(valid.subj_id == pup_subjs(s), :);
        for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn)(s,:) = avg_pup_cond(sv, conds{k}, goals{g}, max_samp); end, end
    end
    t_win = t_pup >= 0.5 & t_pup <= 1.5;
    for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup_mean.(fn) = mean(pup.(fn)(:,t_win), 2, 'omitnan'); end, end
end

function [clusters, p_vals, t_obs] = cluster_perm_ttest(d1, d2, n_perm)
    [n_subj, n_t] = size(d1);
    t_obs = zeros(1, n_t);
    for t = 1:n_t, [~,~,~,s] = ttest(d1(:,t), d2(:,t)); t_obs(t) = s.tstat; end
    t_thresh = tinv(0.975, n_subj-1);
    [clusters, cl_stats] = find_clusters(abs(t_obs) > t_thresh, t_obs);
    if isempty(clusters), p_vals = []; return; end
    null_max = zeros(n_perm, 1);
    for p = 1:n_perm
        signs = (rand(n_subj,1) > 0.5)*2 - 1;
        d1p = d1; d2p = d2;
        for s = 1:n_subj, if signs(s)<0, tmp = d1p(s,:); d1p(s,:) = d2p(s,:); d2p(s,:) = tmp; end, end
        t_perm = zeros(1, n_t);
        for t = 1:n_t, [~,~,~,s] = ttest(d1p(:,t), d2p(:,t)); t_perm(t) = s.tstat; end
        [~, perm_stats] = find_clusters(abs(t_perm) > t_thresh, t_perm);
        if ~isempty(perm_stats), null_max(p) = max(perm_stats); end
    end
    p_vals = arrayfun(@(i) mean(null_max >= cl_stats(i)), 1:length(clusters));
end

function [clusters, cl_stats] = find_clusters(above_thresh, t_vals)
    clusters = {}; cl_stats = []; in_cl = false; cl_start = 0;
    for t = 1:length(above_thresh)
        if above_thresh(t) && ~in_cl, cl_start = t; in_cl = true;
        elseif ~above_thresh(t) && in_cl, clusters{end+1} = cl_start:(t-1); cl_stats(end+1) = sum(abs(t_vals(cl_start:(t-1)))); in_cl = false; end
    end
    if in_cl, clusters{end+1} = cl_start:length(above_thresh); cl_stats(end+1) = sum(abs(t_vals(cl_start:end))); end
end

function print_clusters(lbl, clusters, p_vals, t_vec)
    if isempty(clusters), fprintf('%s: no sig clusters\n', lbl); return; end
    for i = 1:length(clusters)
        if p_vals(i) < 0.05
            fprintf('%s: %.3f-%.3fs, p=%.3f %s\n', lbl, t_vec(clusters{i}(1)), t_vec(clusters{i}(end)), p_vals(i), repmat('*',1,(p_vals(i)<0.05)+(p_vals(i)<0.01)+(p_vals(i)<0.001)));
        end
    end
    if all(p_vals >= 0.05), fprintf('%s: no sig clusters\n', lbl); end
end