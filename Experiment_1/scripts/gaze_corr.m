%%%%%%%%%%%%%%%%%%%%%%%
%% extract correctness-split data
%%%%%%%%%%%%%%%%%%%%%%%
% gaze by correctness
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617];
base_dir = '..'; 
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');
min_rt = 0.150;
% colors
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
load(fullfile(res_dir, 'gaze_reinstat_res_m.mat'));
[reinst_bb_comp_corr, reinst_bb_comp_incorr, reinst_bb_iso_corr, reinst_bb_iso_incorr, ...
 match_bb_comp_corr, match_bb_comp_incorr, match_bb_iso_corr, match_bb_iso_incorr, ...
 baseline_bb_comp_corr, baseline_bb_comp_incorr, baseline_bb_iso_corr, baseline_bb_iso_incorr, ...
 n_bb_comp_corr, n_bb_comp_incorr, n_bb_iso_corr, n_bb_iso_incorr] = extract_gaze_by_correctness(reinstat_res, subj_ids);

% pupil by correctness
[pup_corr, pup_incorr, pup_mean_corr, pup_mean_incorr, pup_subjs_corr, ...
 n_ab_comp_corr, n_ab_comp_incorr, n_ab_iso_corr, n_ab_iso_incorr] = extract_pupil_by_correctness(all_preprocessed);

%% gaze by correctness
fprintf('\n=== GAZE BY CORRECTNESS ===\n');
fprintf('BB Compared:\n');
[stat_results.gaze_corr_bb_comp.p, stat_results.gaze_corr_bb_comp.stats] = do_ttest(reinst_bb_comp_corr, reinst_bb_comp_incorr);
fprintf('  t(%d)=%.2f, p=%.4f, n_corr=%d, n_incorr=%d\n', stat_results.gaze_corr_bb_comp.stats.df, stat_results.gaze_corr_bb_comp.stats.tstat, stat_results.gaze_corr_bb_comp.p, sum(~isnan(reinst_bb_comp_corr)), sum(~isnan(reinst_bb_comp_incorr)));
fprintf('BB Isolated:\n');
[stat_results.gaze_corr_bb_iso.p, stat_results.gaze_corr_bb_iso.stats] = do_ttest(reinst_bb_iso_corr, reinst_bb_iso_incorr);
fprintf('  t(%d)=%.2f, p=%.4f, n_corr=%d, n_incorr=%d\n', stat_results.gaze_corr_bb_iso.stats.df, stat_results.gaze_corr_bb_iso.stats.tstat, stat_results.gaze_corr_bb_iso.p, sum(~isnan(reinst_bb_iso_corr)), sum(~isnan(reinst_bb_iso_incorr)));

fprintf('\n=== GAZE PERMUTATION BY CORRECTNESS ===\n');
[stat_results.gaze_match_bb_comp_corr.p, stat_results.gaze_match_bb_comp_corr.obs] = run_permutation(match_bb_comp_corr, baseline_bb_comp_corr, n_perm);
[stat_results.gaze_match_bb_comp_incorr.p, stat_results.gaze_match_bb_comp_incorr.obs] = run_permutation(match_bb_comp_incorr, baseline_bb_comp_incorr, n_perm);
[stat_results.gaze_match_bb_iso_corr.p, stat_results.gaze_match_bb_iso_corr.obs] = run_permutation(match_bb_iso_corr, baseline_bb_iso_corr, n_perm);
[stat_results.gaze_match_bb_iso_incorr.p, stat_results.gaze_match_bb_iso_incorr.obs] = run_permutation(match_bb_iso_incorr, baseline_bb_iso_incorr, n_perm);
fprintf('BB comp correct:   p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_comp_corr.p, stat_results.gaze_match_bb_comp_corr.obs);
fprintf('BB comp incorrect: p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_comp_incorr.p, stat_results.gaze_match_bb_comp_incorr.obs);
fprintf('BB iso correct:    p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_iso_corr.p, stat_results.gaze_match_bb_iso_corr.obs);
fprintf('BB iso incorrect:  p=%.4f, obs=%.3f\n', stat_results.gaze_match_bb_iso_incorr.p, stat_results.gaze_match_bb_iso_incorr.obs);

%% pupil by correctness
fprintf('\n=== PUPIL BY CORRECTNESS ===\n');
fprintf('A-B Compared:\n');
[stat_results.pup_corr_ab_comp.p, stat_results.pup_corr_ab_comp.stats] = do_ttest(pup_mean_corr.ab_com, pup_mean_incorr.ab_com);
fprintf('  t(%d)=%.2f, p=%.4f, n_corr=%d, n_incorr=%d\n', stat_results.pup_corr_ab_comp.stats.df, stat_results.pup_corr_ab_comp.stats.tstat, stat_results.pup_corr_ab_comp.p, sum(~isnan(pup_mean_corr.ab_com)), sum(~isnan(pup_mean_incorr.ab_com)));
fprintf('A-B Isolated:\n');
[stat_results.pup_corr_ab_iso.p, stat_results.pup_corr_ab_iso.stats] = do_ttest(pup_mean_corr.ab_iso, pup_mean_incorr.ab_iso);
fprintf('  t(%d)=%.2f, p=%.4f, n_corr=%d, n_incorr=%d\n', stat_results.pup_corr_ab_iso.stats.df, stat_results.pup_corr_ab_iso.stats.tstat, stat_results.pup_corr_ab_iso.p, sum(~isnan(pup_mean_corr.ab_iso)), sum(~isnan(pup_mean_incorr.ab_iso)));

fprintf('\n=== PUPIL TIME SERIES BY CORRECTNESS (CLUSTER PERM) ===\n');
fprintf('A-B Compared:\n');
[stat_results.pup_ts_corr_ab_comp.clusters, stat_results.pup_ts_corr_ab_comp.p, stat_results.pup_ts_corr_ab_comp.t] = cluster_perm_ttest(pup_corr.ab_com(:,t_idx), pup_incorr.ab_com(:,t_idx), n_perm);
print_clusters('  corr vs incorr', stat_results.pup_ts_corr_ab_comp.clusters, stat_results.pup_ts_corr_ab_comp.p, t_clust);
fprintf('A-B Isolated:\n');
[stat_results.pup_ts_corr_ab_iso.clusters, stat_results.pup_ts_corr_ab_iso.p, stat_results.pup_ts_corr_ab_iso.t] = cluster_perm_ttest(pup_corr.ab_iso(:,t_idx), pup_incorr.ab_iso(:,t_idx), n_perm);
print_clusters('  corr vs incorr', stat_results.pup_ts_corr_ab_iso.clusters, stat_results.pup_ts_corr_ab_iso.p, t_clust);

%% fig: gaze by correctness
figure('color','w','position',[50 50 1000 800]);
subplot(2,2,1);
data = [match_bb_comp_corr, baseline_bb_comp_corr];
raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'BB Comp Correct', [0,1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_comp_corr.p);
subplot(2,2,2);
data = [match_bb_comp_incorr, baseline_bb_comp_incorr];
raincloud(data, {c_comp, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'BB Comp Incorrect', [0,1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_comp_incorr.p);
subplot(2,2,3);
data = [match_bb_iso_corr, baseline_bb_iso_corr];
raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'BB Iso Correct', [0,1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_iso_corr.p);
subplot(2,2,4);
data = [match_bb_iso_incorr, baseline_bb_iso_incorr];
raincloud(data, {c_iso, [200 200 200]/255}, {'match','mismatch'}, 'Spatial Similarity', 'BB Iso Incorrect', [0,1]);
set(gca, 'YTick', [0 0.25 0.5 0.75 1]); add_sig_perm(data, [1 2], stat_results.gaze_match_bb_iso_incorr.p);
sgtitle('Gaze Reinstatement by Correctness', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'gaze_reins_correctness.pdf'), '-dpdf', '-vector');

%% fig: gaze reinstatement index by correctness
figure('color','w','position',[50 50 800 600]);
data_matrix = [reinst_bb_comp_corr, reinst_bb_comp_incorr, reinst_bb_iso_corr, reinst_bb_iso_incorr];
labels = {'comp corr', 'comp incorr', 'iso corr', 'iso incorr'};
colors = {c_comp, [c_comp*0.5 + [1 1 1]*0.5], c_iso, [c_iso*0.5 + [1 1 1]*0.5]}; hold on;
yline(0, 'r--', 'Chance', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
raincloud(data_matrix, colors, labels, 'Gaze Reinstatement Index', 'Gaze by Correctness', [-0.1 0.4]);
set(gca, 'YTick', [-0.1 0 0.1 0.2 0.3 0.4]);
p_chance_all = [stat_results.gaze_match_bb_comp_corr.p, stat_results.gaze_match_bb_comp_incorr.p, stat_results.gaze_match_bb_iso_corr.p, stat_results.gaze_match_bb_iso_incorr.p];
for i = 1:4
    p = p_chance_all(i);
    if p < 0.05
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
        text(i, 0.45, txt, 'HorizontalAlignment', 'center', 'FontSize', 18);
    end
end
pairs_to_test = [1 2; 3 4];
pvals = [stat_results.gaze_corr_bb_comp.p; stat_results.gaze_corr_bb_iso.p];
add_sig_perm(data_matrix, pairs_to_test, pvals);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'gaze_reins_idx_correctness.pdf'), '-dpdf', '-vector');

%% fig: pupil by correctness timeseries
figure('color','w','position',[50 50 1200 500]);
subplot(1,2,1); hold on;
plot_pup_ts(t_pup, pup_corr.ab_com, c_comp, 'correct');
plot_pup_ts(t_pup, pup_incorr.ab_com, [c_comp*0.5 + [1 1 1]*0.5], 'incorrect');
xlabel('Time from stimulus onset (s)', 'FontSize', 14); ylabel('Pupil size change (a.u.)', 'FontSize', 14);
title('A-B Compared by Correctness', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
subplot(1,2,2); hold on;
plot_pup_ts(t_pup, pup_corr.ab_iso, c_iso, 'correct');
plot_pup_ts(t_pup, pup_incorr.ab_iso, [c_iso*0.5 + [1 1 1]*0.5], 'incorrect');
xlabel('Time from stimulus onset (s)', 'FontSize', 14); ylabel('Pupil size change (a.u.)', 'FontSize', 14);
title('A-B Isolated by Correctness', 'FontSize', 14); legend('Location', 'best'); grid off; box off; xlim([0 1.5]);
sgtitle('Pupil by Correctness', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'pupil_ts_correctness.pdf'), '-dpdf', '-vector');

%% fig: pupil mean by correctness
figure('color','w','position',[50 50 1000 400]);
subplot(1,2,1);
data = [pup_mean_corr.ab_com, pup_mean_incorr.ab_com];
raincloud(data, {c_comp, [c_comp*0.5 + [1 1 1]*0.5]}, {'correct','incorrect'}, 'Mean Pupil (a.u.)', 'A-B Compared');
add_sig(data, [1 2]);
subplot(1,2,2);
data = [pup_mean_corr.ab_iso, pup_mean_incorr.ab_iso];
raincloud(data, {c_iso, [c_iso*0.5 + [1 1 1]*0.5]}, {'correct','incorrect'}, 'Mean Pupil (a.u.)', 'A-B Isolated');
add_sig(data, [1 2]);
sgtitle('Pupil by Correctness', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(fig_dir, 'pupil_mean_correctness.pdf'), '-dpdf', '-vector');

function [p, stats] = do_ttest(d1, d2)
    valid = ~isnan(d1) & ~isnan(d2);
    [~, p, ~, stats] = ttest(d1(valid), d2(valid));
end

function [reinst_bb_comp_corr, reinst_bb_comp_incorr, reinst_bb_iso_corr, reinst_bb_iso_incorr, ...
          match_bb_comp_corr, match_bb_comp_incorr, match_bb_iso_corr, match_bb_iso_incorr, ...
          baseline_bb_comp_corr, baseline_bb_comp_incorr, baseline_bb_iso_corr, baseline_bb_iso_incorr, ...
          n_bb_comp_corr, n_bb_comp_incorr, n_bb_iso_corr, n_bb_iso_incorr] = extract_gaze_by_correctness(reinstat_res, subj_ids)
    n = length(subj_ids);
    reinst_bb_comp_corr = nan(n,1); reinst_bb_comp_incorr = nan(n,1); reinst_bb_iso_corr = nan(n,1); reinst_bb_iso_incorr = nan(n,1);
    match_bb_comp_corr = nan(n,1); match_bb_comp_incorr = nan(n,1); match_bb_iso_corr = nan(n,1); match_bb_iso_incorr = nan(n,1);
    baseline_bb_comp_corr = nan(n,1); baseline_bb_comp_incorr = nan(n,1); baseline_bb_iso_corr = nan(n,1); baseline_bb_iso_incorr = nan(n,1);
    n_bb_comp_corr = nan(n,1); n_bb_comp_incorr = nan(n,1); n_bb_iso_corr = nan(n,1); n_bb_iso_incorr = nan(n,1);
    for i = 1:n
        sid = subj_ids(i);
        comp_corr = reinstat_res.bb_compared.subj_id==sid & reinstat_res.bb_compared.correct==1;
        comp_incorr = reinstat_res.bb_compared.subj_id==sid & reinstat_res.bb_compared.correct==0;
        iso_corr = reinstat_res.bb_isolated.subj_id==sid & reinstat_res.bb_isolated.correct==1;
        iso_incorr = reinstat_res.bb_isolated.subj_id==sid & reinstat_res.bb_isolated.correct==0;
        reinst_bb_comp_corr(i) = mean(reinstat_res.bb_compared.reinst_index(comp_corr),'omitnan');
        reinst_bb_comp_incorr(i) = mean(reinstat_res.bb_compared.reinst_index(comp_incorr),'omitnan');
        reinst_bb_iso_corr(i) = mean(reinstat_res.bb_isolated.reinst_index(iso_corr),'omitnan');
        reinst_bb_iso_incorr(i) = mean(reinstat_res.bb_isolated.reinst_index(iso_incorr),'omitnan');
        match_bb_comp_corr(i) = mean(reinstat_res.bb_compared.match_score(comp_corr),'omitnan');
        match_bb_comp_incorr(i) = mean(reinstat_res.bb_compared.match_score(comp_incorr),'omitnan');
        match_bb_iso_corr(i) = mean(reinstat_res.bb_isolated.match_score(iso_corr),'omitnan');
        match_bb_iso_incorr(i) = mean(reinstat_res.bb_isolated.match_score(iso_incorr),'omitnan');
        baseline_bb_comp_corr(i) = mean(reinstat_res.bb_compared.baseline_score(comp_corr),'omitnan');
        baseline_bb_comp_incorr(i) = mean(reinstat_res.bb_compared.baseline_score(comp_incorr),'omitnan');
        baseline_bb_iso_corr(i) = mean(reinstat_res.bb_isolated.baseline_score(iso_corr),'omitnan');
        baseline_bb_iso_incorr(i) = mean(reinstat_res.bb_isolated.baseline_score(iso_incorr),'omitnan');
        n_bb_comp_corr(i) = sum(comp_corr); n_bb_comp_incorr(i) = sum(comp_incorr);
        n_bb_iso_corr(i) = sum(iso_corr); n_bb_iso_incorr(i) = sum(iso_incorr);
    end
end

function [pup_corr, pup_incorr, pup_mean_corr, pup_mean_incorr, pup_subjs_corr, ...
          n_ab_comp_corr, n_ab_comp_incorr, n_ab_iso_corr, n_ab_iso_incorr] = extract_pupil_by_correctness(all_preprocessed)
    valid = all_preprocessed(all_preprocessed.preprocess_success, :);
    max_samp = max(cellfun(@length, valid.pupil_preprocessed));
    sr = valid.sample_rate(1); t_pup = (0:max_samp-1)/sr;
    pup_subjs_corr = unique(valid.subj_id); n = length(pup_subjs_corr);
    pup_corr = struct(); pup_incorr = struct(); pup_mean_corr = struct(); pup_mean_incorr = struct();
    n_ab_comp_corr = nan(n,1); n_ab_comp_incorr = nan(n,1); n_ab_iso_corr = nan(n,1); n_ab_iso_incorr = nan(n,1);
    for s = 1:n
        sv = valid(valid.subj_id == pup_subjs_corr(s), :);
        pup_corr.ab_com(s,:) = avg_pup_cond_corr(sv, 'compared', 'A-B', max_samp, 1);
        pup_incorr.ab_com(s,:) = avg_pup_cond_corr(sv, 'compared', 'A-B', max_samp, 0);
        pup_corr.ab_iso(s,:) = avg_pup_cond_corr(sv, 'isolated', 'A-B', max_samp, 1);
        pup_incorr.ab_iso(s,:) = avg_pup_cond_corr(sv, 'isolated', 'A-B', max_samp, 0);
        n_ab_comp_corr(s) = sum(strcmp(sv.condition,'compared') & strcmp(sv.goal,'A-B') & sv.correct==1);
        n_ab_comp_incorr(s) = sum(strcmp(sv.condition,'compared') & strcmp(sv.goal,'A-B') & sv.correct==0);
        n_ab_iso_corr(s) = sum(strcmp(sv.condition,'isolated') & strcmp(sv.goal,'A-B') & sv.correct==1);
        n_ab_iso_incorr(s) = sum(strcmp(sv.condition,'isolated') & strcmp(sv.goal,'A-B') & sv.correct==0);
    end
    t_win = t_pup >= 0.5 & t_pup <= 1.5;
    pup_mean_corr.ab_com = mean(pup_corr.ab_com(:,t_win), 2, 'omitnan');
    pup_mean_incorr.ab_com = mean(pup_incorr.ab_com(:,t_win), 2, 'omitnan');
    pup_mean_corr.ab_iso = mean(pup_corr.ab_iso(:,t_win), 2, 'omitnan');
    pup_mean_incorr.ab_iso = mean(pup_incorr.ab_iso(:,t_win), 2, 'omitnan');
end

function a = avg_pup_cond_corr(sv, cond, goal, maxs, correct_flag)
    tr = sv(strcmp(sv.condition,cond) & strcmp(sv.goal,goal) & sv.correct==correct_flag, :);
    if height(tr)==0, a=nan(1,maxs); return; end
    m = nan(height(tr), maxs);
    for i = 1:height(tr)
        p = tr.pupil_preprocessed{i} - mean(tr.baseline_pupil_preprocessed{i},'omitnan');
        m(i,1:length(p)) = p';
    end
    a = mean(m,1,'omitnan');
end

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