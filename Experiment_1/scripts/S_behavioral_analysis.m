clear; clc; close all;

% %%%%%%%%%%%%%%%%%%%%%%%
% %% setup
% %%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631];
addpath(genpath('/Users/bai/Documents/GitHub/MEOW/toolbox/bayesFactor-master'));
base_dir = '..'; 
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');
min_rt = 0.150;
% colors
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
% load(fullfile(res_dir, 'all_subj_results.mat'));
% load(fullfile(res_dir, 'spatial_entropy_results.mat'));
% load(fullfile(res_dir, 'gaze_cumu_res.mat'));
% load(fullfile(res_dir, 'all_trials_pupil.mat'), 'all_preprocessed');
% [pup, pup_mean, t_pup, pup_subjs] = extract_pupil_subj(all_preprocessed);
% fprintf('pupil: %d subjs, %d samples (%.2fs)\n', length(pup_subjs), length(t_pup), t_pup(end));
load(fullfile(res_dir, 'gaze_reinstat_res_full.mat'));
[reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso, match_bb_comp, match_bb_iso, match_ba_comp, match_ba_iso, ...
 baseline_bb_comp, baseline_bb_iso, baseline_ba_comp, baseline_ba_iso] = extract_gaze_subj(reinstat_res, subj_ids);
load(fullfile(res_dir, 'gaze_a2b2_vs_a1b1.mat'));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% math
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% get_v = @(f1, f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs);
% b1_acc_sam = get_v('one','acc_same')'; b1_acc_sim = get_v('one','acc_sim')'; b1_acc_new = get_v('one','acc_new')';
% b1_rt_sam = get_v('one','rt_same')'; b1_rt_sim = get_v('one','rt_sim')';
% b2_ldi_c = get_v('two','ldi_comp')'; b2_ldi_i = get_v('two','ldi_iso')'; b2_ldi_n = get_v('two','ldi_nov')';
% b2_dp_c = get_v('two','dprime_comp')'; b2_dp_i = get_v('two','dprime_iso')'; b2_dp_n = get_v('two','dprime_nov')';
% b2_rt_l_c = get_v('two','rt_AB_comp')'; b2_rt_l_i = get_v('two','rt_AB_iso')'; b2_rt_l_n = get_v('two','rt_AB_nov')';
% b2_rt_t_c = get_v('two','rt_AA_comp')'; b2_rt_t_i = get_v('two','rt_AA_iso')'; b2_rt_t_n = get_v('two','rt_AA_nov')';
% rec_d_c = get_v('rec','d_comp'); rec_d_i = get_v('rec','d_iso');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % B2B1 gaze consistency
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaze_bb_corr = nan(length(subj_ids), 2); gaze_bb_incorr = nan(length(subj_ids), 2); gaze_bb_overall = nan(length(subj_ids), 2);
for c = 1:2
    if c==1, bb_data = reinstat_res.bb_compared; else, bb_data = reinstat_res.bb_isolated; end
    for s = 1:length(subj_ids)
        sid = subj_ids(s); bb_subj = bb_data(bb_data.subj_id == sid, :); if height(bb_subj) == 0, continue; end
        gaze_bb_overall(s,c) = mean(bb_subj.reinst_index, 'omitnan');
        corr_trials = bb_subj.correct == 1; incorr_trials = bb_subj.correct == 0;
        if sum(corr_trials) > 0, gaze_bb_corr(s,c) = mean(bb_subj.reinst_index(corr_trials), 'omitnan'); end
        if sum(incorr_trials) > 0, gaze_bb_incorr(s,c) = mean(bb_subj.reinst_index(incorr_trials), 'omitnan'); end
    end
end
exc_bb = [609, 606, 608]; remaining_subj_ids_bb = subj_ids;
for e = exc_bb, idx = find(remaining_subj_ids_bb == e); if ~isempty(idx), gaze_bb_corr(idx,:) = []; gaze_bb_incorr(idx,:) = []; gaze_bb_overall(idx,:) = []; remaining_subj_ids_bb(idx) = []; end; end
fprintf('B2-B1 final n=%d\n', length(remaining_subj_ids_bb));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A2B1 gaze reinstatement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaze_ba_corr = nan(length(subj_ids), 2); gaze_ba_incorr = nan(length(subj_ids), 2); gaze_ba_overall = nan(length(subj_ids), 2);
for c = 1:2
    if c==1, ba_data = reinstat_res.ba_compared; bb_data = reinstat_res.bb_compared; else, ba_data = reinstat_res.ba_isolated; bb_data = reinstat_res.bb_isolated; end
    for s = 1:length(subj_ids)
        sid = subj_ids(s); ba_subj = ba_data(ba_data.subj_id == sid, :); if height(ba_subj) == 0, continue; end
        ba_subj.b2_correct = nan(height(ba_subj), 1);
        for i = 1:height(ba_subj), matching_b2 = bb_data(bb_data.subj_id == sid & bb_data.tr_1b_b == ba_subj.tr_1b_b(i), :); if height(matching_b2) == 1, ba_subj.b2_correct(i) = matching_b2.correct; end; end
        corr_trials = ba_subj.b2_correct == 1; incorr_trials = ba_subj.b2_correct == 0; gaze_ba_overall(s,c) = mean(ba_subj.reinst_index, 'omitnan');
        if sum(corr_trials) > 0, gaze_ba_corr(s,c) = mean(ba_subj.reinst_index(corr_trials), 'omitnan'); end
        if sum(incorr_trials) > 0, gaze_ba_incorr(s,c) = mean(ba_subj.reinst_index(incorr_trials), 'omitnan'); end
    end
end
exc_ba = [609, 606, 608, 618]; remaining_subj_ids_ba = subj_ids;
for e = exc_ba, idx = find(remaining_subj_ids_ba == e); if ~isempty(idx), gaze_ba_corr(idx,:) = []; gaze_ba_incorr(idx,:) = []; gaze_ba_overall(idx,:) = []; remaining_subj_ids_ba(idx) = []; end; end
fprintf('A2-B1 final n=%d\n', length(remaining_subj_ids_ba));
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%
% % % extract cumulative gaze reinstatement
% % % %%%%%%%%%%%%%%%%%%%%%%%
% % subj_traj_comp = cumulative_res.subj_traj_comp;
% % subj_traj_iso = cumulative_res.subj_traj_iso;
% % n_fix_to_plot = cumulative_res.n_fix_to_plot;
% % early_range = cumulative_res.early_range;
% % late_range = cumulative_res.late_range;
% % 
% % mean_comp = mean(subj_traj_comp, 1, 'omitnan');
% % sem_comp = std(subj_traj_comp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_comp), 1));
% % mean_iso = mean(subj_traj_iso, 1, 'omitnan');
% % sem_iso = std(subj_traj_iso, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_iso), 1));
% % fix_numbers = 1:n_fix_to_plot;
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%
% % % % extract pupil time series
% % % %%%%%%%%%%%%%%%%%%%%%%%
% % valid = all_preprocessed(all_preprocessed.preprocess_success, :);
% % max_samp = max(cellfun(@length, valid.pupil_preprocessed));
% % sr = valid.sample_rate(1);
% % t_pup = (0:max_samp-1)/sr;
% % pup_subjs = unique(valid.subj_id);
% % n_subjs = length(pup_subjs);
% % conds = {'compared', 'isolated', 'novel'};
% % goals = {'A-B', 'A-A'};
% % gtag = {'ab', 'aa'};
% % t_win = t_pup >= 1.0 & t_pup <= 1.5;
% % 
% % pup = struct(); pup_mean = struct();
% % for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn) = nan(n_subjs, max_samp); end, end
% % pup_corr = struct(); pup_incorr = struct();
% % for k = 1:3
% %     fn = sprintf('ab_%s', conds{k}(1:3));
% %     pup_corr.(fn) = nan(n_subjs, max_samp);
% %     pup_incorr.(fn) = nan(n_subjs, max_samp);
% % end
% % 
% % for s = 1:n_subjs
% %     sid = pup_subjs(s); sv = valid(valid.subj_id == sid, :);
% %     for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn)(s,:) = avg_pup_cond(sv, conds{k}, goals{g}, max_samp); end, end
% %     for k = 1:3
% %         fn = sprintf('ab_%s', conds{k}(1:3));
% %         sv_corr = sv(strcmp(sv.condition, conds{k}) & strcmp(sv.goal, 'A-B') & sv.correct == 1, :);
% %         pup_corr.(fn)(s,:) = avg_pup_cond_correctness(sv_corr, max_samp);
% %         sv_incorr = sv(strcmp(sv.condition, conds{k}) & strcmp(sv.goal, 'A-B') & sv.correct == 0, :);
% %         pup_incorr.(fn)(s,:) = avg_pup_cond_correctness(sv_incorr, max_samp);
% %     end
% % end
% % for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup_mean.(fn) = mean(pup.(fn)(:, t_win), 2, 'omitnan'); end, end
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% statistical tests
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% within = table(categorical({'compared';'isolated';'novel'}), 'VariableNames', {'Condition'});
% stat_results = struct();

%% 1-back tests
% fprintf('\n=== 1-BACK STATS ===\n');
% [~, stat_results.b1_rt_sam_v_sim.p, ~, stat_results.b1_rt_sam_v_sim.stats] = ttest(b1_rt_sam, b1_rt_sim);
% [~, stat_results.b1_acc_sam_v_sim.p, ~, stat_results.b1_acc_sam_v_sim.stats] = ttest(b1_acc_sam, b1_acc_sim);
% [~, stat_results.b1_acc_sim_v_new.p, ~, stat_results.b1_acc_sim_v_new.stats] = ttest(b1_acc_sim, b1_acc_new);
% fprintf('RT: same v sim: t(%d)=%.2f, p=%.4f\n', stat_results.b1_rt_sam_v_sim.stats.df, stat_results.b1_rt_sam_v_sim.stats.tstat, stat_results.b1_rt_sam_v_sim.p);
% fprintf('Acc: same v sim: t(%d)=%.2f, p=%.4f\n', stat_results.b1_acc_sam_v_sim.stats.df, stat_results.b1_acc_sam_v_sim.stats.tstat, stat_results.b1_acc_sam_v_sim.p);
% fprintf('Acc: sim v new: t(%d)=%.2f, p=%.4f\n', stat_results.b1_acc_sim_v_new.stats.df, stat_results.b1_acc_sim_v_new.stats.tstat, stat_results.b1_acc_sim_v_new.p);
% 
%% 2-back tests
% fprintf('\n=== 2-BACK STATS ===\n');
% [stat_results.b2_ldi, stat_results.b2_ldi_posthoc] = run_rm_anova('LDI', b2_ldi_c, b2_ldi_i, b2_ldi_n, within);
% [stat_results.b2_dprime, stat_results.b2_dprime_posthoc] = run_rm_anova('d-prime', b2_dp_c, b2_dp_i, b2_dp_n, within);
% [stat_results.b2_rt_lure, stat_results.b2_rt_lure_posthoc] = run_rm_anova('RT lure', b2_rt_l_c, b2_rt_l_i, b2_rt_l_n, within);
% [stat_results.b2_rt_target, stat_results.b2_rt_target_posthoc] = run_rm_anova('RT target', b2_rt_t_c, b2_rt_t_i, b2_rt_t_n, within);
% data_2back = [b2_ldi_c, b2_ldi_i, b2_ldi_n, b2_dp_c, b2_dp_i, b2_dp_n];
% [stat_results.b2_task_cond, stat_results.b2_task_cond_posthoc] = ...
%     run_2factor_rm_anova('2-back: Task × Condition', data_2back, ...
%                          'task', 'condition', ...
%                          {'ldi', 'dprime'}, {'compared', 'isolated', 'novel'});
% data_2back_rt = [b2_rt_l_c, b2_rt_l_i, b2_rt_l_n, b2_rt_t_c, b2_rt_t_i, b2_rt_t_n];
% [stat_results.b2_rt_type_cond, stat_results.b2_rt_type_cond_posthoc] = ...
%     run_2factor_rm_anova('2-back RT: Trial Type × Condition', data_2back_rt, ...
%                          'trial_type', 'condition', ...
%                          {'lure', 'target'}, {'compared', 'isolated', 'novel'});
% %% recognition test
% fprintf('\n=== RECOGNITION STATS ===\n');
% d_tot = (rec_d_c + rec_d_i)/2;
% [~, stat_results.rec_overall.p, ~, stat_results.rec_overall.stats] = ttest(d_tot);
% [~, stat_results.rec_comp_v_iso.p, ~, stat_results.rec_comp_v_iso.stats] = ttest(rec_d_c, rec_d_i);
% fprintf('Overall d'': t(%d)=%.2f, p=%.4f\n', stat_results.rec_overall.stats.df, stat_results.rec_overall.stats.tstat, stat_results.rec_overall.p);
% fprintf('Comp v Iso: t(%d)=%.2f, p=%.4f\n', stat_results.rec_comp_v_iso.stats.df, stat_results.rec_comp_v_iso.stats.tstat, stat_results.rec_comp_v_iso.p);
% 
% [stat_results.gaze_bb_2x2, stat_results.gaze_bb_2x2_ph] = run_2x2_anova_correctness('B2B1 consistency', gaze_bb_corr(:,1), gaze_bb_incorr(:,1), gaze_bb_corr(:,2), gaze_bb_incorr(:,2), gaze_bb_overall(:,1), gaze_bb_overall(:,2));
% [stat_results.gaze_ba_2x2, stat_results.gaze_ba_2x2_ph] = run_2x2_anova_correctness('A2-B1 gaze', gaze_ba_corr(:,1), gaze_ba_incorr(:,1), gaze_ba_corr(:,2), gaze_ba_incorr(:,2), gaze_ba_overall(:,1), gaze_ba_overall(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A2B2 vs A1B1 gaze overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stat_results.gaze_a2b2_a1b1_2x2, stat_results.gaze_a2b2_a1b1_2x2_ph] = run_2x2_anova_time_cond('A2B2 vs A1B1 gaze', gaze_a1b1_comp, gaze_a2b2_comp, gaze_a1b1_iso, gaze_a2b2_iso);


% n_perm=1000;
% fprintf('\ngaze permutation match vs mismatch\n');
% [stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_comp.obs] = run_permutation(match_bb_comp, baseline_bb_comp, n_perm);
% [stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_bb_iso.obs] = run_permutation(match_bb_iso, baseline_bb_iso, n_perm);
% [stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_comp.obs] = run_permutation(match_ba_comp, baseline_ba_comp, n_perm);
% [stat_results.gaze_match_ba_iso.p, stat_results.gaze_match_ba_iso.obs] = run_permutation(match_ba_iso, baseline_ba_iso, n_perm);
% fprintf('BB comp: p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_comp.obs);
% fprintf('BB iso:  p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_bb_iso.obs);
% fprintf('BA comp: p=%.4f, obs=%.3f\n', stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_comp.obs);
% fprintf('BA iso:  p=%.4f, obs=%.3f\n', stat_results.gaze_match_ba_iso.p, stat_results.gaze_match_ba_iso.obs);
% 
% %% spatial entropy
% [stat_results.entropy_2x2, stat_results.entropy_2x2_ph] = run_2x2_anova_item_cond('Spatial Entropy 2×2', entropy_a_comp_subj, entropy_a_iso_subj, entropy_b_comp_subj, entropy_b_iso_subj);
% 
% %% gaze cumu by fixation
% pvals_fixation = nan(1, n_fix_to_plot);
% tvals_fixation = nan(1, n_fix_to_plot);
% dvals_fixation = nan(1, n_fix_to_plot);
% for f = 1:n_fix_to_plot
%     comp_fix = subj_traj_comp(:, f);
%     iso_fix = subj_traj_iso(:, f);
%     valid = ~isnan(comp_fix) & ~isnan(iso_fix);
% 
%     if sum(valid) > 2
%         [~, pvals_fixation(f), ~, stats] = ttest(comp_fix(valid), iso_fix(valid));
%         tvals_fixation(f) = stats.tstat;
%         dvals_fixation(f) = stats.tstat / sqrt(stats.df + 1);
% 
%         fprintf('fixation %d: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', ...
%             f, stats.df, tvals_fixation(f), dvals_fixation(f), pvals_fixation(f), ...
%             repmat('*', 1, (pvals_fixation(f)<0.05)+(pvals_fixation(f)<0.01)+(pvals_fixation(f)<0.001)));
%     end
% end
% sig_fixations = find(pvals_fixation < 0.05);
% if ~isempty(sig_fixations)
%     divergence_point = sig_fixations(1);
%     fprintf('\ndivergence point: Fixation %d (p=%.4f, d=%.2f)\n', ...
%         divergence_point, pvals_fixation(divergence_point), dvals_fixation(divergence_point));
% else
%     fprintf('\nno significant divergence found\n');
%     divergence_point = [];
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% % extract gaze reinstatement by correctness
% %%%%%%%%%%%%%%%%%%%%%%%
% clear gaze_corr gaze_incorr gaze_overall extracted_subj_ids
% gaze_datasets = {reinstat_res.bb_compared, reinstat_res.bb_isolated};
% gaze_corr = nan(length(subj_ids), 2);
% gaze_incorr = nan(length(subj_ids), 2);
% gaze_overall = nan(length(subj_ids), 2);
% extracted_subj_ids = subj_ids;
% 
% for c = 1:2
%     cond_data = gaze_datasets{c};
%     for s = 1:length(extracted_subj_ids)
%         sid = extracted_subj_ids(s);
%         subj_idx = cond_data.subj_id == sid;
%         if sum(subj_idx) > 0 && ismember('correct', cond_data.Properties.VariableNames)
%             subj_data = cond_data(subj_idx, :);
%             corr_trials = subj_data.correct == 1;
%             incorr_trials = subj_data.correct == 0;
%             gaze_overall(s,c) = mean(subj_data.reinst_index, 'omitnan');
%             if sum(corr_trials) > 0, gaze_corr(s,c) = mean(subj_data.reinst_index(corr_trials), 'omitnan'); end
%             if sum(incorr_trials) > 0, gaze_incorr(s,c) = mean(subj_data.reinst_index(incorr_trials), 'omitnan'); end
%         end
%     end
% end
% 
% exc=[609,606, 608]; for e=exc, idx=find(extracted_subj_ids==e); gaze_corr(idx,:)=[]; gaze_incorr(idx,:)=[]; gaze_overall(idx,:)=[]; extracted_subj_ids(idx)=[]; end; fprintf('Excluded: %s, Final n=%d, Diff SD=%.3f\n', mat2str(exc), length(extracted_subj_ids), std(gaze_corr(:,1)-gaze_incorr(:,1),'omitnan'));
% 
% % % gaze-behavior correlation
% % fprintf('\ngaze behavior correlation\n');
% % valid = ~isnan(reinst_bb_comp) & ~isnan(b2_ldi_c) & b2_ldi_c <= 0.8;
% % [stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p] = corr(reinst_bb_comp(valid), b2_ldi_c(valid));
% % stat_results.corr_gr_ldi_comp.n = sum(valid);
% % fprintf('Comp: r=%.3f, p=%.4f, n=%d\n', stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p, stat_results.corr_gr_ldi_comp.n);
% % valid = ~isnan(reinst_bb_iso) & ~isnan(b2_ldi_i) & b2_ldi_i <= 0.8;
% % [stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p] = corr(reinst_bb_iso(valid), b2_ldi_i(valid));
% % stat_results.corr_gr_ldi_iso.n = sum(valid);
% % fprintf('Iso:  r=%.3f, p=%.4f, n=%d\n', stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p, stat_results.corr_gr_ldi_iso.n);
% % 
% % % A2B1 gaze reinstatement x LDI
% % fprintf('\nA2B1 gaze reinstatement x LDI:\n');
% % valid_comp = ~isnan(gaze_ba_comp) & ~isnan(ldi_comp_ba) & ldi_comp_ba <= 0.8;
% % valid_iso = ~isnan(gaze_ba_iso) & ~isnan(ldi_iso_ba) & ldi_iso_ba <= 0.8;
% % [r_ba_comp, p_ba_comp] = corr(gaze_ba_comp(valid_comp), ldi_comp_ba(valid_comp));
% % [r_ba_iso, p_ba_iso] = corr(gaze_ba_iso(valid_iso), ldi_iso_ba(valid_iso));
% % fprintf('  Compared: r=%.3f, p=%.4f, n=%d\n', r_ba_comp, p_ba_comp, sum(valid_comp));
% % fprintf('  Isolated: r=%.3f, p=%.4f, n=%d\n', r_ba_iso, p_ba_iso, sum(valid_iso));
% % 
% % figure('color','w','position',[50 50 500 450]); hold on;
% % scatter(gaze_ba_comp(valid_comp), ldi_comp_ba(valid_comp), 100, c_comp, 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', sprintf('Compared (r=%.2f, p=%.3f)', r_ba_comp, p_ba_comp));
% % scatter(gaze_ba_iso(valid_iso), ldi_iso_ba(valid_iso), 100, c_iso, 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', sprintf('Isolated (r=%.2f, p=%.3f)', r_ba_iso, p_ba_iso));
% % p_comp = polyfit(gaze_ba_comp(valid_comp), ldi_comp_ba(valid_comp), 1);
% % p_iso = polyfit(gaze_ba_iso(valid_iso), ldi_iso_ba(valid_iso), 1);
% % x_comp = linspace(min(gaze_ba_comp(valid_comp)), max(gaze_ba_comp(valid_comp)), 100);
% % x_iso = linspace(min(gaze_ba_iso(valid_iso)), max(gaze_ba_iso(valid_iso)), 100);
% % plot(x_comp, polyval(p_comp, x_comp), '-', 'Color', c_comp, 'LineWidth', 2.5, 'HandleVisibility', 'off');
% % plot(x_iso, polyval(p_iso, x_iso), '-', 'Color', c_iso, 'LineWidth', 2.5, 'HandleVisibility', 'off');
% % xlabel('A2-B1 gaze reinstatement (r)', 'FontSize', 13); ylabel('LDI', 'FontSize', 13);
% % title('A2-B1 Gaze Reinstatement × LDI', 'FontSize', 14);
% % legend('Location', 'best', 'FontSize', 11); set(gca, 'FontSize', 12, 'Box', 'off', 'LineWidth', 1.5);
% % print(gcf, fullfile(fig_dir, 'corr_a2b1_reinstatement_ldi.pdf'), '-dpdf', '-vector');
% 
% %% pupil anova (condition x correctness)
% pup_corr_win_com = mean(pup_corr.ab_com(:, t_win), 2, 'omitnan');
% pup_incorr_win_com = mean(pup_incorr.ab_com(:, t_win), 2, 'omitnan');
% pup_corr_win_iso = mean(pup_corr.ab_iso(:, t_win), 2, 'omitnan');
% pup_incorr_win_iso = mean(pup_incorr.ab_iso(:, t_win), 2, 'omitnan');
% pup_corr_win_nov = mean(pup_corr.ab_nov(:, t_win), 2, 'omitnan');
% pup_incorr_win_nov = mean(pup_incorr.ab_nov(:, t_win), 2, 'omitnan');
% cond_comp = mean(pup.ab_com(:, t_win), 2, 'omitnan');
% cond_iso = mean(pup.ab_iso(:, t_win), 2, 'omitnan');
% cond_nov = mean(pup.ab_nov(:, t_win), 2, 'omitnan');
% pup_corr_win_com(pup_subjs==624)=[]; pup_incorr_win_com(pup_subjs==624)=[]; pup_corr_win_iso(pup_subjs==624)=[]; pup_incorr_win_iso(pup_subjs==624)=[]; pup_corr_win_nov(pup_subjs==624)=[]; pup_incorr_win_nov(pup_subjs==624)=[]; cond_comp(pup_subjs==624)=[]; cond_iso(pup_subjs==624)=[]; cond_nov(pup_subjs==624)=[];
% [stat_results.pup_3x2, stat_results.pup_3x2_ph] = run_3x2_anova_correctness('pupil 3x2', ...
%     pup_corr_win_com, pup_incorr_win_com, pup_corr_win_iso, pup_incorr_win_iso, pup_corr_win_nov, pup_incorr_win_nov, ...
%     cond_comp, cond_iso, cond_nov);
% 
% 
% % fprintf('\npupil timeseries perm tests\n');
% % t_idx = t_pup >= 0 & t_pup <= 1.5; t_clust = t_pup(t_idx); n_perm = 1000;
% % fprintf('Compared: correct vs incorrect\n');
% % [stat_results.pup_ts_corr_com.clusters, stat_results.pup_ts_corr_com.p, stat_results.pup_ts_corr_com.t] = ...
% %     cluster_perm_ttest(pup_corr.ab_com(:,t_idx), pup_incorr.ab_com(:,t_idx), n_perm);
% % print_clusters('  ', stat_results.pup_ts_corr_com.clusters, stat_results.pup_ts_corr_com.p, t_clust);
% % 
% % fprintf('Isolated: correct vs incorrect\n');
% % [stat_results.pup_ts_corr_iso.clusters, stat_results.pup_ts_corr_iso.p, stat_results.pup_ts_corr_iso.t] = ...
% %     cluster_perm_ttest(pup_corr.ab_iso(:,t_idx), pup_incorr.ab_iso(:,t_idx), n_perm);
% % print_clusters('  ', stat_results.pup_ts_corr_iso.clusters, stat_results.pup_ts_corr_iso.p, t_clust);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% fig: 1-back
% figure('color','w','Position',[100 100 1200 400]);
% subplot(1,3,1);
% data = [b1_acc_sam, b1_acc_sim, b1_acc_new];
% raincloud(data, {c_same, c_sim, c_new}, {'same','similar','different'}, 'Accuracy', '', [0,1]);
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(1,3,2);
% data = [b1_rt_sam, b1_rt_sim];
% raincloud(data, {c_same, c_sim}, {'same','similar'}, 'RT (s)', '');
% add_sig(data, [1 2]);
% subplot(1,3,3); hold on;
% mat_1_back = [mean(get_v('one','acc_same')), mean(get_v('one','err_same_as_sim')), mean(get_v('one','err_same_as_new'));
%     mean(get_v('one','err_sim_as_same')), mean(get_v('one','acc_sim')), mean(get_v('one','err_sim_as_new'));
%     mean(get_v('one','err_new_as_same')), mean(get_v('one','err_new_as_sim')), mean(get_v('one','acc_new'))];
% draw_matrix(mat_1_back, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('Response matrix', 'FontSize', 12);
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
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'd''', 'Target Detection');
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(3,3,2);
% data = [b2_rt_l_c, b2_rt_l_i, b2_rt_l_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Lure Discrimination)', [0,1.5]);
% add_sig(data, [1 2; 2 3; 1 3]); set(gca, 'YTick', [0 0.5 1 1.5])
% subplot(3,3,5);
% data = [b2_rt_t_c, b2_rt_t_i, b2_rt_t_n];
% raincloud(data, {c_comp, c_iso, c_nov}, {'compared','isolated','novel'}, 'RT (s)', 'RT (Target Detection)');
% add_sig(data, [1 2; 2 3; 1 3]);
% subplot(3,3,7); hold on;
% mat_comp = [mean(get_v('two','acc_AA_comp')), mean(get_v('two','err_AA_comp_as_k')), mean(get_v('two','err_AA_comp_as_n'));
%     mean(get_v('two','err_AB_comp_as_j')), mean(get_v('two','acc_AB_comp')), mean(get_v('two','err_AB_comp_as_n'));
%     mean(get_v('two','err_AN_comp_as_j')), mean(get_v('two','err_AN_comp_as_k')), mean(get_v('two','acc_AN_comp'))];
% draw_matrix(mat_comp, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('compared', 'FontSize', 14);
% subplot(3,3,8); hold on;
% mat_iso = [mean(get_v('two','acc_AA_iso')), mean(get_v('two','err_AA_iso_as_k')), mean(get_v('two','err_AA_iso_as_n'));
%     mean(get_v('two','err_AB_iso_as_j')), mean(get_v('two','acc_AB_iso')), mean(get_v('two','err_AB_iso_as_n'));
%     mean(get_v('two','err_AN_iso_as_j')), mean(get_v('two','err_AN_iso_as_k')), mean(get_v('two','acc_AN_iso'))];
% draw_matrix(mat_iso, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('isolated', 'FontSize', 14);
% subplot(3,3,9); hold on;
% mat_nov = [mean(get_v('two','acc_AA_nov')), mean(get_v('two','err_AA_nov_as_k')), mean(get_v('two','err_AA_nov_as_n'));
%     mean(get_v('two','err_AB_nov_as_j')), mean(get_v('two','acc_AB_nov')), mean(get_v('two','err_AB_nov_as_n'));
%     mean(get_v('two','err_AN_nov_as_j')), mean(get_v('two','err_AN_nov_as_k')), mean(get_v('two','acc_AN_nov'))];
% draw_matrix(mat_nov, {c_same, c_sim, c_new}, {'Exp Same','Exp Sim','Exp New'}, {'Resp Same','Resp Sim','Resp New'});
% title('novel', 'FontSize', 14);
% set(gcf, 'Color', 'w'); set(gcf, 'Renderer', 'painters'); set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'behav_2back.pdf'), '-dpdf', '-vector');
% 
% %% fig: recognition
% figure('color','w','Position',[100 100 600 500]);
% data = [d_tot', rec_d_c', rec_d_i'];
% hold on;
% fill([-2, 5, 5, -2], [-2, -2, 0, 0], [0.92 0.92 0.92], 'EdgeColor', 'none');
% yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
% raincloud(data, {[0.3 0.3 0.3], c_comp, c_iso}, {'compared+isolated','compared','isolated'}, 'd''', 'Validation of Episodic Encoding');
% add_sig(data, [1 0; 2 3]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'behav_recog.pdf'), '-dpdf', '-vector');
% 
% %% fig: gaze match vs mismatch
% figure('color','w','position',[50 50 1000 800]);
% subplot(2,2,1);
% data = [match_bb_comp, baseline_bb_comp];
% raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'compared B2B1', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_comp.p);
% subplot(2,2,2);
% data = [match_bb_iso, baseline_bb_iso];
% raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'isolated B2B1', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_iso.p);
% subplot(2,2,3);
% data = [match_ba_comp, baseline_ba_comp];
% raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'compared A2B1', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_ba_comp.p);
% subplot(2,2,4);
% data = [match_ba_iso, baseline_ba_iso];
% raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'isolated A2B1', [0,1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_ba_iso.p);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reins_mmis.pdf'), '-dpdf', '-vector');
% 
% %% fig: gaze permutation distributions
% figure('color','w','position',[50 50 1000 800]);
% perm_titles = {'compared B2B1', 'isolated B2B1', 'compared A2B1', 'isolated A2B1'};
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
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reins_perm.pdf'), '-dpdf', '-vector');
% 
% %% fig: gaze reinstatement score
% figure('color','w','position',[50 50 800 600]);
% data_matrix = [reinst_bb_comp, reinst_bb_iso, reinst_ba_comp, reinst_ba_iso];
% labels = {'compared B2B1', 'isolated B2B1', 'compared A2B1', 'isolated A2B1'};
% colors = {c_comp, c_iso, c_comp, c_iso}; hold on;
% yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
% raincloud(data_matrix, colors, labels, 'Gaze Reinstatement Score', 'Gaze Reinstatement', [-0.1 0.4]);
% set(gca, 'YTick', [-0.1 0 0.1 0.2 0.3 0.4]);
% p_chance_all = [stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_iso.p];
% for i = 1:4
%     p = p_chance_all(i);
%     if p < 0.05
%         if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
%         text(i, 0.45, txt, 'HorizontalAlignment', 'center', 'FontSize', 18);
%     end
% end
% pairs_to_test = [1 2; 3 4; 1 3]; 
% pvals = [stat_results.gaze_match_bb_comp.p, stat_results.gaze_match_bb_iso.p, stat_results.gaze_match_ba_comp.p, stat_results.gaze_match_ba_iso.p];
% add_sig_perm(data_matrix, pairs_to_test, pvals);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reins_idx.pdf'), '-dpdf', '-vector');
% 
%% fig. gaze by condition and correctness
figure('color','w','position',[50 50 900 500]);
subplot(1,2,1); raincloud([gaze_bb_corr(:,1), gaze_bb_incorr(:,1)], {c_comp, c_comp*0.5+[1 1 1]*0.5}, {'B2 correct', 'B2 incorrect'}, 'Gaze Consistency (r)', 'compared', []); add_sig([gaze_bb_corr(:,1), gaze_bb_incorr(:,1)], [1 2]); set(gca, 'YTick', [-0.2 0 0.2 0.4], 'YLim', [-0.2 0.5], 'XTick', [1 2]);
subplot(1,2,2); raincloud([gaze_bb_corr(:,2), gaze_bb_incorr(:,2)], {c_iso, c_iso*0.5+[1 1 1]*0.5}, {'B2 correct', 'B2 incorrect'}, 'Gaze Consistency (r)', 'isolated', []); add_sig([gaze_bb_corr(:,2), gaze_bb_incorr(:,2)], [1 2]); set(gca, 'YTick', [-0.2 0 0.2 0.4], 'YLim', [-0.2 0.5], 'XTick', [1 2]);
print(gcf, fullfile(fig_dir, 'gaze_b2b1_cor.pdf'), '-dpdf', '-vector');

figure('color','w','position',[50 50 900 500]);
subplot(1,2,1); raincloud([gaze_ba_corr(:,1), gaze_ba_incorr(:,1)], {c_comp, c_comp*0.5+[1 1 1]*0.5}, {'B2 correct', 'B2 incorrect'}, 'Gaze reinstatement (r)', 'compared', []); add_sig([gaze_ba_corr(:,1), gaze_ba_incorr(:,1)], [1 2]); set(gca, 'YTick', [-0.2 0 0.2 0.4], 'YLim', [-0.2 0.5], 'XTick', [1 2]);
subplot(1,2,2); raincloud([gaze_ba_corr(:,2), gaze_ba_incorr(:,2)], {c_iso, c_iso*0.5+[1 1 1]*0.5}, {'B2 correct', 'B2 incorrect'}, 'Gaze reinstatement (r)', 'isolated', []); add_sig([gaze_ba_corr(:,2), gaze_ba_incorr(:,2)], [1 2]); set(gca, 'YTick', [-0.2 0 0.2 0.4], 'YLim', [-0.2 0.5], 'XTick', [1 2]);
print(gcf, fullfile(fig_dir, 'gaze_a2b1_cor.pdf'), '-dpdf', '-vector');

%% fig. spatial alignment
figure('color','w','position',[50 50 900 500]);
subplot(1,2,1); raincloud([gaze_a2b2_comp, gaze_a1b1_comp], {c_comp, c_comp*0.5+[1 1 1]*0.5}, {'A2B2', 'A1B1'}, 'Gaze Overlap (r)', 'compared', []); add_sig([gaze_a2b2_comp, gaze_a1b1_comp], [1 2]); set(gca, 'YTick', [0.4 0.5 0.6 0.7 0.8 0.9 1.0], 'YLim', [0.4 1.0], 'XTick', [1 2]);
subplot(1,2,2); raincloud([gaze_a2b2_iso, gaze_a1b1_iso], {c_iso, c_iso*0.5+[1 1 1]*0.5}, {'A2B2', 'A1B1'}, 'Gaze Overlap (r)', 'isolated', []); add_sig([gaze_a2b2_iso, gaze_a1b1_iso], [1 2]); set(gca, 'YTick', [0.4 0.5 0.6 0.7 0.8 0.9 1.0], 'YLim', [0.4 1.0], 'XTick', [1 2]);
print(gcf, fullfile(fig_dir, 'gaze_a2b2_vs_a1b1.pdf'), '-dpdf', '-vector');

% % fig. encoding spatial entropy: A B items (Compared vs Isolated)
% figure('Position', [100, 100, 1000, 500]);
% subplot(1,2,1);
% data_a = [entropy_a_comp_subj, entropy_a_iso_subj];
% raincloud(data_a, {c_comp, c_iso}, {'compared', 'isolated'}, 'Spatial Entropy (bits)', 'Encoding A', [6 6.8]);
% set(gca, 'YTick', [6 6.2 6.4 6.6 6.8]); add_sig(data_a, [1 2]);
% subplot(1,2,2);
% data_b = [entropy_b_comp_subj, entropy_b_iso_subj];
% raincloud(data_b, {c_comp, c_iso}, {'compared', 'isolated'}, 'Spatial Entropy (bits)', 'Encoding B', [6 6.8]);
% set(gca, 'YTick', [6 6.2 6.4 6.6 6.8]); add_sig(data_b, [1 2]);
% print(gcf, fullfile(fig_dir, 'spatial_entropy_item.pdf'), '-dpdf', '-vector');
% 
% % fig. encoding AB similarity
% figure('Position', [100, 100, 500, 500]);
% valid = ~isnan(subj_overlap_comp) & ~isnan(subj_overlap_iso);
% comp_col = subj_overlap_comp(valid);
% iso_col = subj_overlap_iso(valid);
% if isrow(comp_col), comp_col = comp_col'; end
% if isrow(iso_col), iso_col = iso_col'; end
% data_overlap = [comp_col, iso_col];
% raincloud(data_overlap, {c_comp, c_iso}, {'compared', 'isolated'}, 'Spatial Similarity (r)', 'A-B Spatial Similarity', [0 1]);
% set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig(data_overlap, [1 2]);
% print(gcf, fullfile(fig_dir, 'spatial_entropy_relational.pdf'), '-dpdf', '-vector');
% 
% % fig. gaze cumu by fixation
% figure('color', 'w', 'position', [50 50 600 600]);
% subplot('Position', [0.12, 0.15, 0.8, 0.75]); hold on;
% fix_numbers_plot = 1:n_fix_to_plot;
% fill([fix_numbers_plot, fliplr(fix_numbers_plot)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], c_comp, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([fix_numbers_plot, fliplr(fix_numbers_plot)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], c_iso, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(fix_numbers_plot, mean_comp, 'o-', 'Color', c_comp, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', c_comp);
% plot(fix_numbers_plot, mean_iso, 's-', 'Color', c_iso, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', c_iso);
% yline(0, 'k--', 'LineWidth', 1.5);
% y_top = max([mean_comp + sem_comp, mean_iso + sem_iso]) + 0.02;
% ylims = ylim;
% y_sig = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
% for f = 1:n_fix_to_plot
%     if pvals_fixation(f) < 0.05
%         if pvals_fixation(f) < 0.001, txt = '***';
%         elseif pvals_fixation(f) < 0.01, txt = '**';
%         else, txt = '*'; end
%         text(f, y_sig, txt, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
%     end
% end
% if ~isempty(divergence_point)
%     ylims = ylim;
%     y_range = ylims(2) - ylims(1);
%     line_bottom = ylims(1) + y_range * 0.1;
%     line_top = ylims(1) + y_range * 0.7;
%     plot([divergence_point, divergence_point], [line_bottom, line_top], '--', 'Color', [1 0.4 0.6], 'LineWidth', 2);
%     text(divergence_point, line_top + y_range*0.05, sprintf('Start to diverge'), ...
%         'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', [1 0.4 0.6], 'FontWeight', 'bold');
% end
% title('Fixation-by-Fixation Reinstatement Buildup', 'FontSize', 16, 'FontWeight', 'bold');
% xlabel('Cumulative Fixation Number', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Gaze Reinstatement Score', 'FontSize', 16, 'FontWeight', 'bold');
% legend({'', '', 'compared', 'isolated'}, 'Location', 'northwest', 'FontSize', 16, 'Box', 'off');
% xlim([0, n_fix_to_plot + 0.5]);
% set(gca, 'XTick', 0:n_fix_to_plot, 'FontSize', 16, 'LineWidth', 1.5);
% box off; grid off;
% print(gcf, fullfile(fig_dir, 'gaze_cumu.pdf'), '-dpdf', '-vector');
% 
% 
% % EM-WM relationships
% fprintf('\n=== WM-EM RELATIONSHIP ===\n');
% wm = (b2_dp_c + b2_dp_i + b2_dp_n)/3;
% lc = b2_ldi_c;
% li = b2_ldi_i;
% ln = b2_ldi_n;
% em_rel = lc - li;
% em_comp = lc - ln;
% em_iso = li - ln;
% em_avg = (lc + li)/2 - ln;
% valid = ~isnan(wm) & ~isnan(lc) & ~isnan(li) & ~isnan(ln);
% wm = wm(valid);
% lc = lc(valid);
% li = li(valid);
% ln = ln(valid);
% em_rel = lc - li;
% em_comp = lc - ln;
% em_iso = li - ln;
% em_avg = (lc + li)/2 - ln;
% fprintf('\n--- LDI Correlations ---\n');
% [r_cn,p_cn] = corr(lc,ln); fprintf('compared vs novel: r=%.3f, p=%.4f\n',r_cn,p_cn);
% [r_ci,p_ci] = corr(lc,li); fprintf('compared vs isolated: r=%.3f, p=%.4f\n',r_ci,p_ci);
% [r_in,p_in] = corr(li,ln); fprintf('isolated vs novel: r=%.3f, p=%.4f\n',r_in,p_in);
% fprintf('\n--- WM predicts EM benefit ---\n');
% [r_wm_rel,p_wm_rel] = corr(wm,em_rel); fprintf('compared-isolated: r=%.3f, p=%.4f\n',r_wm_rel,p_wm_rel);
% [r_wm_comp,p_wm_comp] = corr(wm,em_comp); fprintf('compared-novel: r=%.3f, p=%.4f\n',r_wm_comp,p_wm_comp);
% [r_wm_iso,p_wm_iso] = corr(wm,em_iso); fprintf('isolated-novel: r=%.3f, p=%.4f\n',r_wm_iso,p_wm_iso);
% [r_wm_avg,p_wm_avg] = corr(wm,em_avg); fprintf('(compared+isolated)/2-novel: r=%.3f, p=%.4f\n',r_wm_avg,p_wm_avg);
% figure('color','w','position',[50 50 1600 400]);
% subplot(1,5,1); scatter(ln,lc,80,'filled','MarkerFaceAlpha',0.6); lsline;
% xlabel('LDI novel','FontSize',14,'FontWeight','bold'); ylabel('LDI compared','FontSize',14,'FontWeight','bold');
% title(sprintf('r=%.3f, p=%.4f',r_cn,p_cn),'FontSize',12); axis square; grid off; box off; set(gca,'FontSize',12,'LineWidth',1.5);
% subplot(1,5,2); scatter(wm,em_rel,80,'filled','MarkerFaceAlpha',0.6); lsline;
% xlabel('WM (d'' total)','FontSize',14,'FontWeight','bold'); ylabel('compared-isolated','FontSize',14,'FontWeight','bold');
% title(sprintf('r=%.3f, p=%.4f',r_wm_rel,p_wm_rel),'FontSize',12); axis square; grid off; box off; set(gca,'FontSize',12,'LineWidth',1.5);
% subplot(1,5,3); scatter(wm,em_comp,80,'filled','MarkerFaceAlpha',0.6); lsline;
% xlabel('WM (d'' total)','FontSize',14,'FontWeight','bold'); ylabel('compared-novel','FontSize',14,'FontWeight','bold');
% title(sprintf('r=%.3f, p=%.4f',r_wm_comp,p_wm_comp),'FontSize',12); axis square; grid off; box off; set(gca,'FontSize',12,'LineWidth',1.5);
% subplot(1,5,4); scatter(wm,em_iso,80,'filled','MarkerFaceAlpha',0.6); lsline;
% xlabel('WM (d'' total)','FontSize',14,'FontWeight','bold'); ylabel('isolated-novel','FontSize',14,'FontWeight','bold');
% title(sprintf('r=%.3f, p=%.4f',r_wm_iso,p_wm_iso),'FontSize',12); axis square; grid off; box off; set(gca,'FontSize',12,'LineWidth',1.5);
% subplot(1,5,5); scatter(wm,em_avg,80,'filled','MarkerFaceAlpha',0.6); lsline;
% xlabel('WM (d'' total)','FontSize',14,'FontWeight','bold'); ylabel('(compared+isolated)/2-novel','FontSize',14,'FontWeight','bold');
% title(sprintf('r=%.3f, p=%.4f',r_wm_avg,p_wm_avg),'FontSize',12); axis square; grid off; box off; set(gca,'FontSize',12,'LineWidth',1.5);
% sgtitle('WM-EM Relationship','FontSize',16,'FontWeight','bold');
% print(gcf,fullfile(fig_dir,'wm_em_relationship.pdf'),'-dpdf','-vector','-fillpage');
% stat_results.wm_em.n=sum(valid);
% stat_results.wm_em.r_cn=r_cn; stat_results.wm_em.p_cn=p_cn;
% stat_results.wm_em.r_ci=r_ci; stat_results.wm_em.p_ci=p_ci;
% stat_results.wm_em.r_in=r_in; stat_results.wm_em.p_in=p_in;
% stat_results.wm_em.r_wm_rel=r_wm_rel; stat_results.wm_em.p_wm_rel=p_wm_rel;
% stat_results.wm_em.r_wm_comp=r_wm_comp; stat_results.wm_em.p_wm_comp=p_wm_comp;
% stat_results.wm_em.r_wm_iso=r_wm_iso; stat_results.wm_em.p_wm_iso=p_wm_iso;
% stat_results.wm_em.r_wm_avg=r_wm_avg; stat_results.wm_em.p_wm_avg=p_wm_avg;
% 
% %% fig: gaze-behavior correlation
% figure('color','w','position',[50 50 800 400]);
% subplot(1,2,1); hold on;
% valid = ~isnan(reinst_bb_comp) & ~isnan(b2_ldi_c) & b2_ldi_c <= 0.8;
% x = reinst_bb_comp(valid); y = b2_ldi_c(valid);
% scatter(x, y, 60, c_comp, 'filled', 'MarkerFaceAlpha', 0.6);
% p_fit = polyfit(x, y, 1); x_fit = linspace(min(x), max(x), 100);
% plot(x_fit, polyval(p_fit, x_fit), 'Color', c_comp, 'LineWidth', 2);
% xlabel('Gaze Reinstatement Score', 'FontSize', 12); ylabel('LDI', 'FontSize', 12);
% title(sprintf('compared: r=%.2f, p=%.3f (n=%d)', stat_results.corr_gr_ldi_comp.r, stat_results.corr_gr_ldi_comp.p, stat_results.corr_gr_ldi_comp.n), 'FontSize', 14);
% grid off; box off;
% subplot(1,2,2); hold on;
% valid = ~isnan(reinst_bb_iso) & ~isnan(b2_ldi_i) & b2_ldi_i <= 0.8;
% x = reinst_bb_iso(valid); y = b2_ldi_i(valid);
% scatter(x, y, 60, c_iso, 'filled', 'MarkerFaceAlpha', 0.6);
% plot(x_fit, polyval(p_fit, x_fit), 'Color', c_iso, 'LineWidth', 2);
% xlabel('Gaze Reinstatement Score', 'FontSize', 12); ylabel('LDI', 'FontSize', 12);
% title(sprintf('isolated: r=%.2f, p=%.3f (n=%d)', stat_results.corr_gr_ldi_iso.r, stat_results.corr_gr_ldi_iso.p, stat_results.corr_gr_ldi_iso.n), 'FontSize', 14);
% grid off; box off;
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'gaze_reinstat_ldi.pdf'), '-dpdf', '-vector');
% 
% %% fig. pupil timeseries
% figure('color','w','position',[50 50 1200 500]);
% subplot(1,2,1); hold on;
% plot_pup_ts(t_pup, pup.ab_com, c_comp, 'compared');
% plot_pup_ts(t_pup, pup.ab_iso, c_iso, 'isolated');
% plot_pup_ts(t_pup, pup.ab_nov, c_nov, 'novel');
% xlabel('Time (s)', 'FontSize', 14); ylabel('Baseline-corrected pupil size change', 'FontSize', 14);
% title('Lure discrimination', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
% subplot(1,2,2); hold on;
% plot_pup_ts(t_pup, pup.aa_com, c_comp, 'compared');
% plot_pup_ts(t_pup, pup.aa_iso, c_iso, 'isolated');
% plot_pup_ts(t_pup, pup.aa_nov, c_nov, 'novel');
% xlabel('Time (s)', 'FontSize', 14); ylabel('Baseline-corrected pupil size change', 'FontSize', 14);
% title('Target detection', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'pupil_condition.pdf'), '-dpdf', '-vector');
% 
% %% fig. pupil by condition and correctness
% figure('color','w','position',[50 50 1400 500]);
% subplot(1,3,1);
% mat_comp = [pup_corr_win_com, pup_incorr_win_com];
% raincloud(mat_comp, {c_comp, c_comp*0.5+[1 1 1]*0.5}, {'correct', 'incorrect'}, 'Pupil change (win 1-1.5s)', 'compared', []);
% add_sig(mat_comp, [1 2]);
% subplot(1,3,2);
% mat_iso = [pup_corr_win_iso, pup_incorr_win_iso];
% raincloud(mat_iso, {c_iso, c_iso*0.5+[1 1 1]*0.5}, {'correct', 'incorrect'}, 'Pupil change (win 1-1.5s)', 'isolated', []);
% add_sig(mat_iso, [1 2]);
% subplot(1,3,3);
% mat_nov = [pup_corr_win_nov, pup_incorr_win_nov];
% raincloud(mat_nov, {c_nov, c_nov*0.5+[1 1 1]*0.5}, {'correct', 'incorrect'}, 'Pupil change (win 1-1.5s)', 'novel', []);
% add_sig(mat_nov, [1 2]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'pupil_correctness_raincloud.pdf'), '-dpdf', '-vector');



%%%%%%%%%%%%%%%%%%%%%%
% save comprehensive results
%%%%%%%%%%%%%%%%%%%%%%
% save(fullfile(res_dir, 'all_subj_results.mat'), ...
%     'subj_ids', 'extracted_subj_ids', ...
%     'b1_acc_sam', 'b1_acc_sim', 'b1_acc_new', 'b1_rt_sam', 'b1_rt_sim', ...
%     'b2_ldi_c', 'b2_ldi_i', 'b2_ldi_n', 'b2_dp_c', 'b2_dp_i', 'b2_dp_n', ...
%     'b2_rt_l_c', 'b2_rt_l_i', 'b2_rt_l_n', 'b2_rt_t_c', 'b2_rt_t_i', 'b2_rt_t_n', ...
%     'rec_d_c', 'rec_d_i', 'd_tot', ...
%     'reinst_bb_comp', 'reinst_bb_iso', 'reinst_ba_comp', 'reinst_ba_iso', ...
%     'match_bb_comp', 'match_bb_iso', 'match_ba_comp', 'match_ba_iso', ...
%     'baseline_bb_comp', 'baseline_bb_iso', 'baseline_ba_comp', 'baseline_ba_iso', ...
%     'gaze_corr', 'gaze_incorr', 'gaze_overall', ...
%     'pup_corr_win_com', 'pup_incorr_win_com', 'pup_corr_win_iso', 'pup_incorr_win_iso', 'pup_corr_win_nov', 'pup_incorr_win_nov', ...
%     'cond_comp', 'cond_iso', 'cond_nov', ...
%     'pup', 'pup_mean', 'pup_corr', 'pup_incorr', 't_pup', 't_win', 'pup_subjs', ...
%     'stat_results', ...
%     'c_comp', 'c_iso', 'c_nov');
% fprintf('\nresults saved: %s\n', fullfile(res_dir, 'all_subj_results.mat'));
% save(fullfile(res_dir, 'comprehensive_behavioral_results.mat'), ...
%     'subj_ids', 'all_subj_stats', ...
%     'acc_AB_comp', 'acc_AB_iso', 'acc_AB_nov', ...
%     'acc_AA_comp', 'acc_AA_iso', 'acc_AA_nov', ...
%     'rt_AB_comp', 'rt_AB_iso', 'rt_AB_nov', ...
%     'rt_AA_comp', 'rt_AA_iso', 'rt_AA_nov', ...
%     'ldi_comp', 'ldi_iso', 'ldi_nov', ...
%     'dprime_comp', 'dprime_iso', 'dprime_nov', ...
%     'reinst_bb_comp', 'reinst_bb_iso', 'reinst_ba_comp', 'reinst_ba_iso', ...
%     'match_bb_comp', 'match_bb_iso', 'match_ba_comp', 'match_ba_iso', ...
%     'baseline_bb_comp', 'baseline_bb_iso', 'baseline_ba_comp', 'baseline_ba_iso', ...
%     'a2b2_comp_subj', 'a2b2_iso_subj', 'a1b1_comp_subj', 'a1b1_iso_subj', ...
%     'pup', 'pup_mean', 't_pup', 'pup_subjs', ...
%     'entropy_comp', 'entropy_iso', 'entropy_nov', ... % add actual variable names from spatial_entropy_results
%     'min_rt', 'c_comp', 'c_iso', 'c_nov', 'c_sim', 'c_same', 'c_new', ...
%     '-v7.3');



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
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 18);
    if ~isempty(ttl), title(ttl, 'FontSize', 20); end
    ylabel(ylbl,'FontSize',18,'FontWeight','bold'); xlim([0.2, n_grps+0.8]);
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
               text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',20,'FontWeight','bold');
            else, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
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

function [res, posthoc] = run_2x2_anova_correctness(lbl, cc, ci, ic, ii, cond_comp, cond_iso)
    fprintf('\n%s:\n', lbl);
    t = table(cc, ci, ic, ii, 'VariableNames', {'cc','ci','ic','ii'});
    within = table({'comp';'comp';'iso';'iso'}, {'corr';'incorr';'corr';'incorr'}, 'VariableNames', {'cond','correctness'});
    rm = fitrm(t, 'cc-ii~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'cond*correctness'); eps = epsilon(rm); m = mauchly(rm);
    rows = {'(Intercept):cond', '(Intercept):correctness', '(Intercept):cond:correctness'};
    labels = {'condition', 'correctness', 'interaction'};
    term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            mauch_p = m.pValue(1); ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f, eta=%.2f %s\n', ...
                labels{i}, df1, df2_gg, f, p_gg, mauch_p, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.cond.F=f; res.cond.p=p_gg; res.cond.eta=eta; res.cond.mauchly_p=mauch_p; res.cond.eps=ep; end
            if i==2, res.corr.F=f; res.corr.p=p_gg; res.corr.eta=eta; res.corr.mauchly_p=mauch_p; res.corr.eps=ep; end
            if i==3, res.int.F=f; res.int.p=p_gg; res.int.eta=eta; res.int.mauchly_p=mauch_p; res.int.eps=ep; end
        end
    end
    fprintf('  condition post-hoc:\n');
    [~,pc,~,sc] = ttest(cond_comp, cond_iso);
    posthoc.cond_comp_iso.t=sc.tstat; posthoc.cond_comp_iso.p=pc; posthoc.cond_comp_iso.d=sc.tstat/sqrt(sc.df+1);
    fprintf('    comp>iso: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', sc.df, sc.tstat, posthoc.cond_comp_iso.d, pc, repmat('*',1,(pc<0.05)+(pc<0.01)+(pc<0.001)));
    fprintf('  correctness post-hoc:\n');
    [~,p1,~,s1] = ttest(cc, ci); [~,p2,~,s2] = ttest(ic, ii); [~,p3,~,s3] = ttest(cc-ci, ic-ii);
    df = s1.df;
    posthoc.comp_corr_incorr.t=s1.tstat; posthoc.comp_corr_incorr.p=p1; posthoc.comp_corr_incorr.d=s1.tstat/sqrt(df+1);
    posthoc.iso_corr_incorr.t=s2.tstat; posthoc.iso_corr_incorr.p=p2; posthoc.iso_corr_incorr.d=s2.tstat/sqrt(df+1);
    posthoc.comp_vs_iso.t=s3.tstat; posthoc.comp_vs_iso.p=p3; posthoc.comp_vs_iso.d=s3.tstat/sqrt(df+1);
    fprintf('    comp corr>incorr: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s1.tstat, posthoc.comp_corr_incorr.d, p1, repmat('*',1,(p1<0.05)+(p1<0.01)+(p1<0.001)));
    fprintf('    iso corr>incorr: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s2.tstat, posthoc.iso_corr_incorr.d, p2, repmat('*',1,(p2<0.05)+(p2<0.01)+(p2<0.01)+(p2<0.001)));
    fprintf('    (comp-iso) diff: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s3.tstat, posthoc.comp_vs_iso.d, p3, repmat('*',1,(p3<0.05)+(p3<0.01)+(p3<0.001)));
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

function [res, posthoc] = run_3x2_anova_correctness(lbl, cc, ci, ic, ii, nc, ni, cond_comp, cond_iso, cond_nov)
    t = table(cc, ci, ic, ii, nc, ni, 'VariableNames', {'cc','ci','ic','ii','nc','ni'});
    within = table({'comp';'comp';'iso';'iso';'nov';'nov'}, {'corr';'incorr';'corr';'incorr';'corr';'incorr'}, 'VariableNames', {'cond','correctness'});
    rm = fitrm(t, 'cc-ni~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'cond*correctness'); eps = epsilon(rm); m = mauchly(rm);
    rows = {'(Intercept):cond', '(Intercept):correctness', '(Intercept):cond:correctness'};
    labels = {'condition', 'correctness', 'interaction'};
    term_names = string(tbl.Properties.RowNames);
    fprintf('\n%s:\n', lbl);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            mauch_p = m.pValue(1); ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f, eta=%.2f %s\n', ...
                labels{i}, df1, df2_gg, f, p_gg, mauch_p, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.cond.F=f; res.cond.p=p_gg; res.cond.eta=eta; res.cond.mauchly_p=mauch_p; res.cond.eps=ep; end
            if i==2, res.corr.F=f; res.corr.p=p_gg; res.corr.eta=eta; res.corr.mauchly_p=mauch_p; res.corr.eps=ep; end
            if i==3, res.int.F=f; res.int.p=p_gg; res.int.eta=eta; res.int.mauchly_p=mauch_p; res.int.eps=ep; end
        end
    end
    n = length(cc); subj = (1:n)';
    bf_tbl = table(subj, cc, ci, ic, ii, nc, ni, 'VariableNames', {'subj','cc','ci','ic','ii','nc','ni'});
    bf_long = stack(bf_tbl, {'cc','ci','ic','ii','nc','ni'}, 'NewDataVariableName','y', 'IndexVariableName','cell');
    bf_long.subj = categorical(bf_long.subj); bf_long.cell = categorical(bf_long.cell);
    res.bf10 = bf.anova(bf_long, 'y~cell', 'treatAsRandom', {'subj'}, 'verbose', false);
    if res.bf10<1, fprintf('  BF01=%.2f (evidence for null)\n', 1/res.bf10); else, fprintf('  BF10=%.2f (evidence for effect)\n', res.bf10); end
    fprintf('  condition post-hoc (holm-bonf):\n');
    [~,pc1,~,sc1] = ttest(cond_comp, cond_iso); 
    [~,pc2,~,sc2] = ttest(cond_iso, cond_nov); 
    [~,pc3,~,sc3] = ttest(cond_comp, cond_nov);
    pvals_cond = [pc1 pc2 pc3]; adj_cond = holm_bonf(pvals_cond); df = sc1.df;
    posthoc.cond_comp_iso.t=sc1.tstat; posthoc.cond_comp_iso.p=pc1; posthoc.cond_comp_iso.p_adj=adj_cond(1); posthoc.cond_comp_iso.d=sc1.tstat/sqrt(df+1);
    posthoc.cond_iso_nov.t=sc2.tstat; posthoc.cond_iso_nov.p=pc2; posthoc.cond_iso_nov.p_adj=adj_cond(2); posthoc.cond_iso_nov.d=sc2.tstat/sqrt(df+1);
    posthoc.cond_comp_nov.t=sc3.tstat; posthoc.cond_comp_nov.p=pc3; posthoc.cond_comp_nov.p_adj=adj_cond(3); posthoc.cond_comp_nov.d=sc3.tstat/sqrt(df+1);
    fprintf('    comp>iso: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, sc1.tstat, posthoc.cond_comp_iso.d, adj_cond(1), repmat('*',1,(adj_cond(1)<0.05)+(adj_cond(1)<0.01)+(adj_cond(1)<0.001)));
    fprintf('    iso>nov: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, sc2.tstat, posthoc.cond_iso_nov.d, adj_cond(2), repmat('*',1,(adj_cond(2)<0.05)+(adj_cond(2)<0.01)+(adj_cond(2)<0.001)));
    fprintf('    comp>nov: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, sc3.tstat, posthoc.cond_comp_nov.d, adj_cond(3), repmat('*',1,(adj_cond(3)<0.05)+(adj_cond(3)<0.01)+(adj_cond(3)<0.001)));
    
    fprintf('  correctness post-hoc (holm-bonf):\n');
    [~,p1,~,s1] = ttest(cc, ci); [~,p2,~,s2] = ttest(ic, ii); [~,p3,~,s3] = ttest(nc, ni);
    [~,p4,~,s4] = ttest(cc-ci, ic-ii); [~,p5,~,s5] = ttest(cc-ci, nc-ni); [~,p6,~,s6] = ttest(ic-ii, nc-ni);
    pvals = [p1 p2 p3 p4 p5 p6]; adj = holm_bonf(pvals);
    posthoc.comp_corr_incorr.t=s1.tstat; posthoc.comp_corr_incorr.p=p1; posthoc.comp_corr_incorr.p_adj=adj(1); posthoc.comp_corr_incorr.d=s1.tstat/sqrt(df+1);
    posthoc.iso_corr_incorr.t=s2.tstat; posthoc.iso_corr_incorr.p=p2; posthoc.iso_corr_incorr.p_adj=adj(2); posthoc.iso_corr_incorr.d=s2.tstat/sqrt(df+1);
    posthoc.nov_corr_incorr.t=s3.tstat; posthoc.nov_corr_incorr.p=p3; posthoc.nov_corr_incorr.p_adj=adj(3); posthoc.nov_corr_incorr.d=s3.tstat/sqrt(df+1);
    posthoc.comp_vs_iso.t=s4.tstat; posthoc.comp_vs_iso.p=p4; posthoc.comp_vs_iso.p_adj=adj(4); posthoc.comp_vs_iso.d=s4.tstat/sqrt(df+1);
    posthoc.comp_vs_nov.t=s5.tstat; posthoc.comp_vs_nov.p=p5; posthoc.comp_vs_nov.p_adj=adj(5); posthoc.comp_vs_nov.d=s5.tstat/sqrt(df+1);
    posthoc.iso_vs_nov.t=s6.tstat; posthoc.iso_vs_nov.p=p6; posthoc.iso_vs_nov.p_adj=adj(6); posthoc.iso_vs_nov.d=s6.tstat/sqrt(df+1);
    fprintf('    comp corr>incorr: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s1.tstat, posthoc.comp_corr_incorr.d, adj(1), repmat('*',1,(adj(1)<0.05)+(adj(1)<0.01)+(adj(1)<0.001)));
    fprintf('    iso corr>incorr: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s2.tstat, posthoc.iso_corr_incorr.d, adj(2), repmat('*',1,(adj(2)<0.05)+(adj(2)<0.01)+(adj(2)<0.001)));
    fprintf('    nov corr>incorr: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s3.tstat, posthoc.nov_corr_incorr.d, adj(3), repmat('*',1,(adj(3)<0.05)+(adj(3)<0.01)+(adj(3)<0.001)));
    fprintf('    (comp-iso) diff: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s4.tstat, posthoc.comp_vs_iso.d, adj(4), repmat('*',1,(adj(4)<0.05)+(adj(4)<0.01)+(adj(4)<0.001)));
    fprintf('    (comp-nov) diff: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s5.tstat, posthoc.comp_vs_nov.d, adj(5), repmat('*',1,(adj(5)<0.05)+(adj(5)<0.01)+(adj(5)<0.001)));
    fprintf('    (iso-nov) diff: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s6.tstat, posthoc.iso_vs_nov.d, adj(6), repmat('*',1,(adj(6)<0.05)+(adj(6)<0.01)+(adj(6)<0.001)));
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

function a = avg_pup_cond_correctness(sv, maxs)
    if height(sv)==0, a=nan(1,maxs); return; end
    m = nan(height(sv), maxs);
    for i = 1:height(sv)
        p = sv.pupil_preprocessed{i} - mean(sv.baseline_pupil_preprocessed{i},'omitnan');
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
function [res, posthoc] = run_2x2_anova_item_cond(lbl, a_comp, a_iso, b_comp, b_iso)
    valid = ~isnan(a_comp) & ~isnan(a_iso) & ~isnan(b_comp) & ~isnan(b_iso);
    a_comp = a_comp(valid); a_iso = a_iso(valid); b_comp = b_comp(valid); b_iso = b_iso(valid);
    fprintf('\n%s (n=%d):\n', lbl, sum(valid));
    t = table(a_comp, a_iso, b_comp, b_iso, 'VariableNames', {'ac','ai','bc','bi'});
    within = table({'a';'a';'b';'b'}, {'comp';'iso';'comp';'iso'}, 'VariableNames', {'item','cond'});
    rm = fitrm(t, 'ac-bi~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'item*cond'); eps = epsilon(rm); m = mauchly(rm);
    rows = {'(Intercept):item', '(Intercept):cond', '(Intercept):item:cond'};
    labels = {'item', 'condition', 'interaction'};
    term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            mauch_p = m.pValue(1); ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f, eta=%.2f %s\n', labels{i}, df1, df2_gg, f, p_gg, mauch_p, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.item.F=f; res.item.p=p_gg; res.item.eta=eta; res.item.mauchly_p=mauch_p; res.item.eps=ep; end
            if i==2, res.cond.F=f; res.cond.p=p_gg; res.cond.eta=eta; res.cond.mauchly_p=mauch_p; res.cond.eps=ep; end
            if i==3, res.int.F=f; res.int.p=p_gg; res.int.eta=eta; res.int.mauchly_p=mauch_p; res.int.eps=ep; end
        end
    end
    
    fprintf('  item post-hoc (holm-bonf):\n');
    item_a = (a_comp + a_iso) / 2;
    item_b = (b_comp + b_iso) / 2;
    [~,pi,~,si] = ttest(item_a, item_b);
    posthoc.item_a_b.t=si.tstat; posthoc.item_a_b.p=pi; posthoc.item_a_b.p_adj=pi; posthoc.item_a_b.d=si.tstat/sqrt(si.df+1);
    fprintf('    a>b: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', si.df, si.tstat, posthoc.item_a_b.d, pi, repmat('*',1,(pi<0.05)+(pi<0.01)+(pi<0.001)));
    
    fprintf('  condition post-hoc (holm-bonf):\n');
    cond_comp_all = (a_comp + b_comp) / 2;
    cond_iso_all = (a_iso + b_iso) / 2;
    [~,pc,~,sc] = ttest(cond_comp_all, cond_iso_all);
    posthoc.cond_comp_iso.t=sc.tstat; posthoc.cond_comp_iso.p=pc; posthoc.cond_comp_iso.p_adj=pc; posthoc.cond_comp_iso.d=sc.tstat/sqrt(sc.df+1);
    fprintf('    comp>iso: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', sc.df, sc.tstat, posthoc.cond_comp_iso.d, pc, repmat('*',1,(pc<0.05)+(pc<0.01)+(pc<0.001)));
    
    fprintf('  simple effects (holm-bonf):\n');
    [~,p1,~,s1] = ttest(a_comp, a_iso); [~,p2,~,s2] = ttest(b_comp, b_iso);
    [~,p3,~,s3] = ttest(a_comp, b_comp); [~,p4,~,s4] = ttest(a_iso, b_iso);
    [~,p5,~,s5] = ttest(a_comp-a_iso, b_comp-b_iso);
    pvals = [p1 p2 p3 p4 p5]; adj = holm_bonf(pvals); df = s1.df;
    posthoc.a_comp_iso.t=s1.tstat; posthoc.a_comp_iso.p=p1; posthoc.a_comp_iso.p_adj=adj(1); posthoc.a_comp_iso.d=s1.tstat/sqrt(df+1);
    posthoc.b_comp_iso.t=s2.tstat; posthoc.b_comp_iso.p=p2; posthoc.b_comp_iso.p_adj=adj(2); posthoc.b_comp_iso.d=s2.tstat/sqrt(df+1);
    posthoc.comp_a_b.t=s3.tstat; posthoc.comp_a_b.p=p3; posthoc.comp_a_b.p_adj=adj(3); posthoc.comp_a_b.d=s3.tstat/sqrt(df+1);
    posthoc.iso_a_b.t=s4.tstat; posthoc.iso_a_b.p=p4; posthoc.iso_a_b.p_adj=adj(4); posthoc.iso_a_b.d=s4.tstat/sqrt(df+1);
    posthoc.a_vs_b_diff.t=s5.tstat; posthoc.a_vs_b_diff.p=p5; posthoc.a_vs_b_diff.p_adj=adj(5); posthoc.a_vs_b_diff.d=s5.tstat/sqrt(df+1);
    fprintf('    a: comp>iso: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s1.tstat, posthoc.a_comp_iso.d, adj(1), repmat('*',1,(adj(1)<0.05)+(adj(1)<0.01)+(adj(1)<0.001)));
    fprintf('    b: comp>iso: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s2.tstat, posthoc.b_comp_iso.d, adj(2), repmat('*',1,(adj(2)<0.05)+(adj(2)<0.01)+(adj(2)<0.001)));
    fprintf('    comp: a>b: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s3.tstat, posthoc.comp_a_b.d, adj(3), repmat('*',1,(adj(3)<0.05)+(adj(3)<0.01)+(adj(3)<0.001)));
    fprintf('    iso: a>b: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s4.tstat, posthoc.iso_a_b.d, adj(4), repmat('*',1,(adj(4)<0.05)+(adj(4)<0.01)+(adj(4)<0.001)));
    fprintf('    (a-b) diff: t(%d)=%.2f, d=%.2f, p_adj=%.4f %s\n', df, s5.tstat, posthoc.a_vs_b_diff.d, adj(5), repmat('*',1,(adj(5)<0.05)+(adj(5)<0.01)+(adj(5)<0.001)));
end
function [clusters, p_vals, t_obs] = cluster_perm_ttest(d1, d2, n_perm)
    [n_subj, n_t] = size(d1); t_obs = zeros(1, n_t);
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

function print_clusters(prefix, clusters, p_vals, t_vec)
    if isempty(clusters), fprintf('%sno sig clusters\n', prefix); return; end
    for i = 1:length(clusters)
        if p_vals(i) < 0.05
            fprintf('%s%.3f-%.3fs, p=%.3f %s\n', prefix, t_vec(clusters{i}(1)), t_vec(clusters{i}(end)), p_vals(i), repmat('*',1,(p_vals(i)<0.05)+(p_vals(i)<0.01)+(p_vals(i)<0.001)));
        end
    end
    if all(p_vals >= 0.05), fprintf('%sno sig clusters\n', prefix); end
end

function [res, posthoc] = run_2x2_anova_time_cond(lbl, a1_comp, a2_comp, a1_iso, a2_iso)
    valid = ~isnan(a1_comp) & ~isnan(a2_comp) & ~isnan(a1_iso) & ~isnan(a2_iso);
    a1_comp = a1_comp(valid); a2_comp = a2_comp(valid); a1_iso = a1_iso(valid); a2_iso = a2_iso(valid);
    fprintf('\n%s (n=%d):\n', lbl, sum(valid));
    
    fprintf('  Means:\n');
    fprintf('    A1B1 compared: M=%.3f, SD=%.3f\n', mean(a1_comp), std(a1_comp));
    fprintf('    A2B2 compared: M=%.3f, SD=%.3f\n', mean(a2_comp), std(a2_comp));
    fprintf('    A1B1 isolated: M=%.3f, SD=%.3f\n', mean(a1_iso), std(a1_iso));
    fprintf('    A2B2 isolated: M=%.3f, SD=%.3f\n', mean(a2_iso), std(a2_iso));
    
    t = table(a1_comp, a2_comp, a1_iso, a2_iso, 'VariableNames', {'a1c','a2c','a1i','a2i'});
    within = table({'a1';'a2';'a1';'a2'}, {'comp';'comp';'iso';'iso'}, 'VariableNames', {'time','cond'});
    rm = fitrm(t, 'a1c-a2i~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'time*cond'); eps = epsilon(rm); m = mauchly(rm);
    rows = {'(Intercept):time', '(Intercept):cond', '(Intercept):time:cond'};
    labels = {'time', 'condition', 'interaction'};
    term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            mauch_p = m.pValue(1); ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f, eta=%.2f %s\n', labels{i}, df1, df2_gg, f, p_gg, mauch_p, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.time.F=f; res.time.p=p_gg; res.time.eta=eta; res.time.mauchly_p=mauch_p; res.time.eps=ep; end
            if i==2, res.cond.F=f; res.cond.p=p_gg; res.cond.eta=eta; res.cond.mauchly_p=mauch_p; res.cond.eps=ep; end
            if i==3, res.int.F=f; res.int.p=p_gg; res.int.eta=eta; res.int.mauchly_p=mauch_p; res.int.eps=ep; end
        end
    end
    
    fprintf('  simple effects (raw p):\n');
    [~,p1,~,s1] = ttest(a2_comp, a1_comp); [~,p2,~,s2] = ttest(a2_iso, a1_iso);
    [~,p3,~,s3] = ttest(a1_comp, a1_iso); [~,p4,~,s4] = ttest(a2_comp, a2_iso);
    df = s1.df;
    posthoc.comp_a2_a1.t=s1.tstat; posthoc.comp_a2_a1.p=p1; posthoc.comp_a2_a1.d=s1.tstat/sqrt(df+1);
    posthoc.iso_a2_a1.t=s2.tstat; posthoc.iso_a2_a1.p=p2; posthoc.iso_a2_a1.d=s2.tstat/sqrt(df+1);
    posthoc.a1_comp_iso.t=s3.tstat; posthoc.a1_comp_iso.p=p3; posthoc.a1_comp_iso.d=s3.tstat/sqrt(df+1);
    posthoc.a2_comp_iso.t=s4.tstat; posthoc.a2_comp_iso.p=p4; posthoc.a2_comp_iso.d=s4.tstat/sqrt(df+1);
    fprintf('    comp: a2b2>a1b1: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s1.tstat, posthoc.comp_a2_a1.d, p1, repmat('*',1,(p1<0.05)+(p1<0.01)+(p1<0.001)));
    fprintf('    iso: a2b2>a1b1: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s2.tstat, posthoc.iso_a2_a1.d, p2, repmat('*',1,(p2<0.05)+(p2<0.01)+(p2<0.001)));
    fprintf('    a1b1: comp>iso: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s3.tstat, posthoc.a1_comp_iso.d, p3, repmat('*',1,(p3<0.05)+(p3<0.01)+(p3<0.001)));
    fprintf('    a2b2: comp>iso: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s4.tstat, posthoc.a2_comp_iso.d, p4, repmat('*',1,(p4<0.05)+(p4<0.01)+(p4<0.001)));
end

function [res, posthoc] = run_2factor_rm_anova(lbl, data_matrix, factor1_name, factor2_name, factor1_levels, factor2_levels)
    fprintf('\n=== %s ===\n', lbl);
    valid = ~any(isnan(data_matrix), 2); data_clean = data_matrix(valid, :); n_subj = size(data_clean, 1);
    fprintf('n=%d (excluded %d)\n', n_subj, sum(~valid));
    n_cols = size(data_clean, 2); var_names = cell(1, n_cols); idx = 1;
    for f1 = 1:length(factor1_levels)
        for f2 = 1:length(factor2_levels)
            var_names{idx} = sprintf('%s_%s', factor1_levels{f1}, factor2_levels{f2}); idx = idx + 1;
        end
    end
    t = array2table(data_clean, 'VariableNames', var_names);
    f1 = repmat(factor1_levels(:), length(factor2_levels), 1); f2 = repelem(factor2_levels(:), length(factor1_levels));
    within = table(f1, f2, 'VariableNames', {factor1_name, factor2_name});
    model_formula = sprintf('%s-%s~1', var_names{1}, var_names{end}); rm = fitrm(t, model_formula, 'WithinDesign', within);
    within_model = sprintf('%s*%s', factor1_name, factor2_name); tbl = ranova(rm, 'WithinModel', within_model);
    eps = epsilon(rm); m = mauchly(rm);
    rows = {sprintf('(Intercept):%s', factor1_name), sprintf('(Intercept):%s', factor2_name), sprintf('(Intercept):%s:%s', factor1_name, factor2_name)};
    labels = {factor1_name, factor2_name, 'interaction'}; term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            if i <= height(m), mauch_p = m.pValue(i); else, mauch_p = 1.0; end
            ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            mauch_str = ''; if mauch_p < 1, mauch_str = sprintf(', mauchly p=%.3f', mauch_p); end
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG)%s, eps=%.3f, eta=%.2f %s\n', ...
                labels{i}, df1, df2_gg, f, p_gg, mauch_str, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.factor1.F=f; res.factor1.p=p_gg; res.factor1.eta=eta; res.factor1.mauchly_p=mauch_p; res.factor1.eps=ep;
            elseif i==2, res.factor2.F=f; res.factor2.p=p_gg; res.factor2.eta=eta; res.factor2.mauchly_p=mauch_p; res.factor2.eps=ep;
            else, res.interaction.F=f; res.interaction.p=p_gg; res.interaction.eta=eta; res.interaction.mauchly_p=mauch_p; res.interaction.eps=ep; end
        end
    end
    fprintf('\n  Simple effects by %s:\n', factor1_name); posthoc = struct();
    for f1 = 1:length(factor1_levels)
        cols_f1 = (f1-1)*length(factor2_levels) + (1:length(factor2_levels)); data_f1 = data_clean(:, cols_f1);
        fprintf('    %s:\n', factor1_levels{f1});
        for c1 = 1:length(factor2_levels)-1
            for c2 = c1+1:length(factor2_levels)
                [~, p, ~, s] = ttest(data_f1(:,c1), data_f1(:,c2)); d = s.tstat / sqrt(n_subj);
                comp_name = sprintf('%s_%s_vs_%s', factor1_levels{f1}, factor2_levels{c1}, factor2_levels{c2});
                posthoc.(comp_name).t = s.tstat; posthoc.(comp_name).p = p; posthoc.(comp_name).d = d;
                fprintf('      %s vs %s: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', factor2_levels{c1}, factor2_levels{c2}, s.df, s.tstat, d, p, ...
                    repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001)));
            end
        end
    end
end