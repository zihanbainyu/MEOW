clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631];
addpath(genpath('/Users/bai/Documents/GitHub/MEOW/toolbox/bayesFactor-master'));
base_dir = '..'; 
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');
min_rt = 0.150;
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
cohend = @(x,y) mean(x-y,'omitnan') / std(x-y,'omitnan');
load(fullfile(res_dir, 'all_subjs_stats.mat'), 'all_subjs');
% load(fullfile(res_dir, 'spatial_entropy_results.mat'));
load(fullfile(res_dir, 'all_trials_pupil.mat'), 'all_preprocessed');
[pup, pup_mean, t_pup, pup_subjs] = extract_pupil_subj(all_preprocessed);
fprintf('pupil: %d subjs, %d samples (%.2fs)\n', length(pup_subjs), length(t_pup), t_pup(end));

get_v = @(f1, f2) arrayfun(@(x) x.stats.(f1).(f2), all_subjs);
b1_acc_sam = get_v('one','acc_same')'; b1_acc_sim = get_v('one','acc_sim')'; b1_acc_new = get_v('one','acc_new')';
b1_rt_sam = get_v('one','rt_same')'; b1_rt_sim = get_v('one','rt_sim')';
b2_ldi_c = get_v('two','ldi_comp')'; b2_ldi_i = get_v('two','ldi_iso')'; b2_ldi_n = get_v('two','ldi_nov')';
b2_dp_c = get_v('two','dprime_comp')'; b2_dp_i = get_v('two','dprime_iso')'; b2_dp_n = get_v('two','dprime_nov')';
b2_rt_l_c = get_v('two','rt_AB_comp')'; b2_rt_l_i = get_v('two','rt_AB_iso')'; b2_rt_l_n = get_v('two','rt_AB_nov')';
b2_rt_t_c = get_v('two','rt_AA_comp')'; b2_rt_t_i = get_v('two','rt_AA_iso')'; b2_rt_t_n = get_v('two','rt_AA_nov')';
rec_d_c = get_v('rec','d_comp'); rec_d_i = get_v('rec','d_iso');


valid = all_preprocessed(all_preprocessed.preprocess_success, :);
max_samp = max(cellfun(@length, valid.pupil_preprocessed));
sr = valid.sample_rate(1);
t_pup = (0:max_samp-1)/sr;
pup_subjs = unique(valid.subj_id);
n_subjs = length(pup_subjs);
conds = {'compared','isolated','novel'}; goals = {'A-B','A-A'}; gtag = {'ab','aa'};
t_win = t_pup >= 1.0 & t_pup <= 1.5;
fprintf('pupil: %d subjs, %d samples (%.2fs)\n', n_subjs, max_samp, t_pup(end));

pup = struct(); pup_corr = struct(); pup_incorr = struct();
for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn) = nan(n_subjs, max_samp); end, end
for k = 1:3, fn = sprintf('ab_%s', conds{k}(1:3)); pup_corr.(fn) = nan(n_subjs, max_samp); pup_incorr.(fn) = nan(n_subjs, max_samp); end

for s = 1:n_subjs
    sid = pup_subjs(s); sv = valid(valid.subj_id == sid, :);
    for g = 1:2, for k = 1:3, fn = sprintf('%s_%s', gtag{g}, conds{k}(1:3)); pup.(fn)(s,:) = avg_pup_cond(sv, conds{k}, goals{g}, max_samp); end, end
    % for k = 1:3
    %     fn = sprintf('ab_%s', conds{k}(1:3));
    %     pup_corr.(fn)(s,:) = avg_pup_cond_correctness(sv(strcmp(sv.condition,conds{k}) & strcmp(sv.goal,'A-B') & sv.correct==1, :), max_samp);
    %     pup_incorr.(fn)(s,:) = avg_pup_cond_correctness(sv(strcmp(sv.condition,conds{k}) & strcmp(sv.goal,'A-B') & sv.correct==0, :), max_samp);
    % end
    for k = 1:3
    fn = sprintf('ab_%s', conds{k}(1:3));
    sv_lure = sv(strcmp(sv.condition, conds{k}) & strcmp(sv.goal, 'A-B') & strcmp(cellstr(sv.corr_resp), 'k'), :);
    pup_corr.(fn)(s,:) = avg_pup_cond_correctness(sv_lure(sv_lure.correct == 1, :), max_samp);
    pup_incorr.(fn)(s,:) = avg_pup_cond_correctness(sv_lure(sv_lure.correct == 0, :), max_samp);
    end
end

exc_pup = false(size(pup_subjs));
pup_lure_corr  = [mean(pup_corr.ab_com(~exc_pup, t_win),2,'omitnan'), mean(pup_corr.ab_iso(~exc_pup, t_win),2,'omitnan'), mean(pup_corr.ab_nov(~exc_pup, t_win),2,'omitnan')];
pup_lure_incorr = [mean(pup_incorr.ab_com(~exc_pup, t_win),2,'omitnan'), mean(pup_incorr.ab_iso(~exc_pup, t_win),2,'omitnan'), mean(pup_incorr.ab_nov(~exc_pup, t_win),2,'omitnan')];
pup_lure_cond  = [mean(pup.ab_com(~exc_pup, t_win),2,'omitnan'), mean(pup.ab_iso(~exc_pup, t_win),2,'omitnan'), mean(pup.ab_nov(~exc_pup, t_win),2,'omitnan')];
pup_targ_cond  = [mean(pup.aa_com(~exc_pup, t_win),2,'omitnan'), mean(pup.aa_iso(~exc_pup, t_win),2,'omitnan'), mean(pup.aa_nov(~exc_pup, t_win),2,'omitnan')];

fprintf('\n--- Pupil: Lure trials (condition × correctness) ---\n');
[stat_results.pup_3x2, stat_results.pup_3x2_ph] = run_3x2_anova_correctness('pupil lure 3x2', ...
    pup_lure_corr(:,1), pup_lure_incorr(:,1), pup_lure_corr(:,2), pup_lure_incorr(:,2), pup_lure_corr(:,3), pup_lure_incorr(:,3), ...
    pup_lure_cond(:,1), pup_lure_cond(:,2), pup_lure_cond(:,3));

fprintf('\n--- Pupil: Target trials (condition main effect) ---\n');
within_3 = table(categorical({'compared';'isolated';'novel'}), 'VariableNames', {'Condition'});
[stat_results.pup_targ_cond, stat_results.pup_targ_cond_ph] = run_rm_anova('Pupil target', pup_targ_cond(:,1), pup_targ_cond(:,2), pup_targ_cond(:,3), within_3);



n_perm   = 1000;
alpha    = 0.05;
f_crit   = finv(1 - 0.05, 2, (n_subjs-1)*2);
t_crit   = tinv(1 - 0.025, n_subjs-1);   % two-tailed

trial_types = {'lure', 'target'};
lure_data   = cat(3, pup.ab_com(~exc_pup,:), pup.ab_iso(~exc_pup,:), pup.ab_nov(~exc_pup,:));
target_data = cat(3, pup.aa_com(~exc_pup,:), pup.aa_iso(~exc_pup,:), pup.aa_nov(~exc_pup,:));
all_data    = {lure_data, target_data};

cond_names  = {'compared','isolated','novel'};
pairs       = [1 2; 1 3; 2 3];
pair_names  = {'compared vs isolated', 'compared vs novel', 'isolated vs novel'};
pair_short_names = {'comp vs iso', 'comp vs nov', 'iso vs nov'};
cluster_results = struct();

for tt = 1:2
    D   = all_data{tt};   % [n_subjs x n_time x 3]
    n_t = size(D, 2);
    cluster_results(tt).trial_type = trial_types{tt};
    cluster_results(tt).condition_data = D;

    %% --- ANOVA ---
    F_obs = zeros(1, n_t);
    for ti = 1:n_t, F_obs(ti) = rm_anova_F(D(:,ti,:)); end

    [obs_clusts, obs_clust_stats] = get_clusters_F(F_obs, f_crit);

    rng(42);
    max_null_F = zeros(n_perm, 1);
    for perm = 1:n_perm
        D_perm = D;
        for s = 1:n_subjs, D_perm(s,:,:) = D(s,:,randperm(3)); end
        F_perm = zeros(1, n_t);
        for ti = 1:n_t, F_perm(ti) = rm_anova_F(D_perm(:,ti,:)); end
        [~, ps] = get_clusters_F(F_perm, f_crit);
        if ~isempty(ps), max_null_F(perm) = max(ps); end
    end

    sig_anova_mask = false(1, n_t);
    anova_clusters = struct('idx', {}, 'stat', {}, 'pval', {}, 'sig', {});
    for ci = 1:length(obs_clusts)
        pval = mean(max_null_F >= obs_clust_stats(ci));
        sig_anova_mask(obs_clusts{ci}) = pval < alpha;
        anova_clusters(ci).idx = obs_clusts{ci};
        anova_clusters(ci).stat = obs_clust_stats(ci);
        anova_clusters(ci).pval = pval;
        anova_clusters(ci).sig = pval < alpha;
    end

    cluster_results(tt).F_obs = F_obs;
    cluster_results(tt).f_crit = f_crit;
    cluster_results(tt).anova_mask = sig_anova_mask;
    cluster_results(tt).anova_clusters = anova_clusters;
    cluster_results(tt).max_null_F = max_null_F;

    fprintf('\n=== %s trials: rm-ANOVA ===\n', trial_types{tt});
    fprintf('F_crit = %.3f\n', f_crit);
    if isempty(obs_clusts)
        fprintf('  no significant clusters\n');
    end
    for ci = 1:length(obs_clusts)
        idx   = obs_clusts{ci};
        cstat = obs_clust_stats(ci);
        pval  = mean(max_null_F >= cstat);
        fprintf('  cluster %d: %.3f–%.3fs  F_sum=%.2f  p=%.4f%s\n', ...
            ci, t_pup(idx(1)), t_pup(idx(end)), cstat, pval, repmat('  ***',1,pval<alpha));
    end

    %% --- pairwise post-hoc cluster permutation ---
    fprintf('\n  -- post-hoc pairwise (sign-flip permutation) --\n');
    cluster_results(tt).pairwise = struct('name', {}, 't_obs', {}, 'clusters', {}, 'max_null_T', {});

    for p = 1:3
        c1 = pairs(p,1); c2 = pairs(p,2);
        diff_obs = D(:,:,c1) - D(:,:,c2);   % [n_subjs x n_time]

        % observed t at every timepoint
        t_obs = mean(diff_obs,1,'omitnan') ./ (std(diff_obs,0,1,'omitnan') / sqrt(n_subjs));

        [obs_clusts_t, obs_clust_stats_t] = get_clusters_T(t_obs, t_crit);
        clusts_t = struct('idx', {}, 'stat', {}, 'pval', {}, 'sig', {}, 'direction', {}, 'overlap', {});

        % sign-flip null
        rng(42);
        max_null_T = zeros(n_perm, 1);
        for perm = 1:n_perm
            signs      = (randi(2, n_subjs, 1)*2 - 3);
            diff_perm  = diff_obs .* signs;
            t_perm     = mean(diff_perm,1,'omitnan') ./ (std(diff_perm,0,1,'omitnan') / sqrt(n_subjs));
            [~, ps]    = get_clusters_T(t_perm, t_crit);
            if ~isempty(ps), max_null_T(perm) = max(abs(ps)); end
        end

        fprintf('\n  %s:\n', pair_names{p});
        if isempty(obs_clusts_t)
            fprintf('    no clusters above threshold\n');
        end
        for ci = 1:length(obs_clusts_t)
            idx   = obs_clusts_t{ci};
            cstat = obs_clust_stats_t(ci);
            pval  = mean(max_null_T >= abs(cstat));
            % direction: which condition is larger
            dir_mean = mean(mean(diff_obs(:,idx),2,'omitnan'),'omitnan');
            if dir_mean > 0
                dir_str = sprintf('%s > %s', cond_names{c1}, cond_names{c2});
            else
                dir_str = sprintf('%s > %s', cond_names{c2}, cond_names{c1});
            end
            overlap = sum(sig_anova_mask(idx)) / length(idx);
            sig_ph  = pval < alpha && overlap > 0;
            clusts_t(ci).idx = idx;
            clusts_t(ci).stat = cstat;
            clusts_t(ci).pval = pval;
            clusts_t(ci).sig = sig_ph;
            clusts_t(ci).direction = dir_str;
            clusts_t(ci).overlap = overlap;
            fprintf('    cluster %d: %.3f–%.3fs  T_sum=%.2f  p=%.4f  [%s]%s\n', ...
                ci, t_pup(idx(1)), t_pup(idx(end)), cstat, pval, dir_str, repmat('  ***',1,sig_ph));
        end
        if isempty(obs_clusts_t)
            clusts_t = struct('idx', {}, 'stat', {}, 'pval', {}, 'sig', {}, 'direction', {}, 'overlap', {});
        end
        cluster_results(tt).pairwise(p).name = pair_names{p};
        cluster_results(tt).pairwise(p).short_name = pair_short_names{p};
        cluster_results(tt).pairwise(p).t_obs = t_obs;
        cluster_results(tt).pairwise(p).clusters = clusts_t;
        cluster_results(tt).pairwise(p).max_null_T = max_null_T;
    end
end

%% Publication-ready figures
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
c_clust  = [108 108 196]/255;   % purple — significant cluster fill
c_ns     = [180 180 180]/255;   % gray — non-sig cluster fill
c_win    = [224 124 60]/255;    % orange — analysis window
c_thresh = [100 100 100]/255;   % gray dashed — threshold line
alpha_fill = 0.20;

win_t = [1.0 1.5];

main_fig = plot_pupil_main_figure(t_pup, cluster_results, cond_names, c_comp, c_iso, c_nov, win_t);
export_pub_figure(main_fig, fig_dir, 'main_pupil_timecourses');

supp_fig = plot_pupil_cluster_supplement(t_pup, cluster_results, win_t, t_crit, ...
    c_clust, c_ns, c_win, c_thresh, alpha_fill);
export_pub_figure(supp_fig, fig_dir, 'supp_pupil_cluster_permutation');

dist_fig = plot_pupil_null_distributions(cluster_results, c_clust, c_thresh);
export_pub_figure(dist_fig, fig_dir, 'supp_pupil_null_distributions');
fprintf('Publication-ready pupil figures saved to %s\n', fig_dir);

%% ── local helpers ────────────────────────────────────────────────────────

function fig = plot_pupil_main_figure(t_pup, cluster_results, cond_names, c_comp, c_iso, c_nov, win_t)
cols = [c_comp; c_iso; c_nov];
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 12]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

for tt = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    D = cluster_results(tt).condition_data;
    yl = get_data_ylim(D);
    patch(ax, [win_t(1) win_t(2) win_t(2) win_t(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.88 0.88 0.88], 'EdgeColor','none', 'FaceAlpha',0.45, 'HandleVisibility','off');
    for c = 1:3
        plot_mean_sem(ax, t_pup, D(:,:,c), cols(c,:), titlecase(cond_names{c}));
    end
    yline(ax, 0, '-', 'Color',[0.82 0.82 0.82], 'LineWidth',0.6, 'HandleVisibility','off');
    add_sig_bar(ax, t_pup, cluster_results(tt).anova_mask, yl);
    format_pub_axes(ax);
    ylim(ax, yl);
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Baseline-corrected pupil');
    title(ax, task_label(cluster_results(tt).trial_type), ...
        'FontName','Helvetica', 'FontSize',9, 'FontWeight','bold');
    if tt == 1
        legend(ax, 'Location','northwest', 'Box','off', 'FontSize',7);
    end
end
title(tl, 'pupil response by condition', 'FontName','Helvetica', 'FontSize',10, 'FontWeight','bold');
end

function fig = plot_pupil_cluster_supplement(t_pup, cluster_results, win_t, t_crit, c_clust, c_ns, c_win, c_thresh, alpha_fill)
fig = figure('Color','w','Units','centimeters','Position',[2 2 30 32]);
tl = tiledlayout(fig, 4, 2, 'TileSpacing','tight', 'Padding','compact');

for tt = 1:2
    ax = nexttile(tl, tt); hold(ax,'on');
    res = cluster_results(tt);
    f_ylim = [0 max(max(res.F_obs, [], 'omitnan') * 1.12, eps)];
    ylim(ax, f_ylim);
    fill_win(ax, win_t, c_win, alpha_fill*0.6);
    shade_clusters_F(ax, t_pup, res.F_obs, res.anova_mask, c_clust, alpha_fill);
    yline(ax, res.f_crit, '--', 'Color', c_thresh, 'LineWidth',2.0, 'HandleVisibility','off');
    plot(ax, t_pup, res.F_obs, 'Color', c_clust, 'LineWidth',4.0);
    label_clusters(ax, t_pup, res.F_obs, res.anova_clusters, c_clust, c_ns);
    format_panel(ax, 'F', false);
    set(ax, 'YTick', unique(round([0 res.f_crit f_ylim(2)], 1)));
    title(ax, sprintf('%s: rm-anova', task_label(res.trial_type)), ...
        'FontSize',20, 'FontWeight','normal', 'FontName','Helvetica');
end

for p = 1:3
    for tt = 1:2
        ax = nexttile(tl, 2 + (p-1)*2 + tt); hold(ax,'on');
        res = cluster_results(tt).pairwise(p);
        tmax = max(abs(res.t_obs), [], 'omitnan');
        t_ylim = [-1 1] * max(tmax * 1.12, t_crit * 1.20);
        ylim(ax, t_ylim);
        fill_win(ax, win_t, c_win, alpha_fill*0.6);
        for ci = 1:length(res.clusters)
            cl = res.clusters(ci);
            mask = false(size(t_pup));
            mask(cl.idx) = true;
            col = c_ns;
            if cl.sig, col = c_clust; end
            shade_clusters_T(ax, t_pup, res.t_obs, mask, col, alpha_fill);
        end
        yline(ax, t_crit, '--', 'Color', c_thresh, 'LineWidth',1.8, 'HandleVisibility','off');
        yline(ax, -t_crit, '--', 'Color', c_thresh, 'LineWidth',1.8, 'HandleVisibility','off');
        yline(ax, 0, '-', 'Color',[0.84 0.84 0.84], 'LineWidth',1.2, 'HandleVisibility','off');
        plot(ax, t_pup, res.t_obs, 'Color', c_clust, 'LineWidth',3.2);
        label_clusters(ax, t_pup, res.t_obs, res.clusters, c_clust, c_ns);
        format_panel(ax, 't', p == 3);
        set(ax, 'YTick', round([-t_crit 0 t_crit], 1));
        title(ax, sprintf('%s: %s', task_label(cluster_results(tt).trial_type), res.name), ...
            'FontSize',20, 'FontWeight','normal', 'FontName','Helvetica');
    end
end
xlabel(tl, 'Time (s)', 'FontName','Helvetica', 'FontSize',20);
title(tl, 'cluster permutation tests', 'FontName','Helvetica', 'FontSize',22, 'FontWeight','bold');
end

function fig = plot_pupil_null_distributions(cluster_results, c_obs, c_null)
fig = figure('Color','w','Units','centimeters','Position',[2 2 30 32]);
tl = tiledlayout(fig, 4, 2, 'TileSpacing','tight', 'Padding','compact');

for tt = 1:2
    ax = nexttile(tl, tt); hold(ax, 'on');
    obs = max_cluster_stat(cluster_results(tt).anova_clusters, false);
    plot_null_panel(ax, cluster_results(tt).max_null_F, obs, c_obs, c_null);
    title(ax, sprintf('%s: rm-anova', task_label(cluster_results(tt).trial_type)), ...
        'FontSize',20, 'FontWeight','normal', 'FontName','Helvetica');
    ylabel(ax, 'count', 'FontSize',20, 'FontName','Helvetica');
end

for p = 1:3
    for tt = 1:2
        ax = nexttile(tl, 2 + (p-1)*2 + tt); hold(ax, 'on');
        res = cluster_results(tt).pairwise(p);
        obs = max_cluster_stat(res.clusters, true);
        plot_null_panel(ax, res.max_null_T, obs, c_obs, c_null);
        title(ax, sprintf('%s: %s', task_label(cluster_results(tt).trial_type), res.name), ...
            'FontSize',20, 'FontWeight','normal', 'FontName','Helvetica');
        if tt == 1, ylabel(ax, 'count', 'FontSize',20, 'FontName','Helvetica'); end
        if p == 3, xlabel(ax, 'max cluster mass', 'FontSize',20, 'FontName','Helvetica'); end
    end
end
title(tl, 'observed cluster mass vs null distribution', ...
    'FontName','Helvetica', 'FontSize',22, 'FontWeight','bold');
end

function export_pub_figure(fig, fig_dir, base_name)
set(fig, 'Renderer','painters', 'InvertHardcopy','off');
pdf_file = fullfile(fig_dir, sprintf('%s.pdf', base_name));
tiff_file = fullfile(fig_dir, sprintf('%s.tiff', base_name));
fig_file = fullfile(fig_dir, sprintf('%s.fig', base_name));
try
    exportgraphics(fig, pdf_file, 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(fig, tiff_file, 'Resolution',600, 'BackgroundColor','white');
catch
    print(fig, pdf_file, '-dpdf', '-painters', '-r600');
    print(fig, tiff_file, '-dtiff', '-r600');
end
savefig(fig, fig_file);
end

function plot_null_panel(ax, null_vals, obs_val, c_obs, c_null)
histogram(ax, null_vals, 24, 'FaceColor',[0.78 0.78 0.78], ...
    'EdgeColor','none', 'FaceAlpha',1);
xline(ax, obs_val, '-', 'Color',c_obs, 'LineWidth',4.0);
set(ax, 'FontSize',20, 'FontName','Helvetica', 'TickDir','out', ...
    'Box','off', 'LineWidth',1.5, 'YTick', []);
ax.XAxis.Exponent = 0;
xl = xlim(ax);
xlim(ax, [min([xl(1), obs_val * 0.95, 0]), max([xl(2), obs_val * 1.05, eps])]);
text(ax, 0.98, 0.90, 'null', 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize',20, 'Color',c_null, 'FontName','Helvetica');
text(ax, 0.98, 0.76, observed_label(obs_val), 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize',20, 'Color',c_obs, 'FontName','Helvetica');
end

function val = max_cluster_stat(clusters, use_abs)
val = 0;
if isempty(clusters), return; end
stats = [clusters.stat];
if isempty(stats), return; end
if use_abs, stats = abs(stats); end
val = max(stats, [], 'omitnan');
end

function lbl = observed_label(obs_val)
if obs_val == 0
    lbl = 'observed = 0';
else
    lbl = 'observed';
end
end

function plot_mean_sem(ax, t, d, col, lbl)
mu = mean(d, 1, 'omitnan');
n = sum(~isnan(d), 1);
se = std(d, 0, 1, 'omitnan') ./ sqrt(n);
vi = isfinite(mu) & isfinite(se);
fill(ax, [t(vi) fliplr(t(vi))], [mu(vi)+se(vi) fliplr(mu(vi)-se(vi))], ...
    col, 'FaceAlpha',0.22, 'EdgeColor','none', 'HandleVisibility','off');
plot(ax, t, mu, 'Color', col, 'LineWidth',1.8, 'DisplayName',lbl);
end

function yl = get_data_ylim(D)
mu = squeeze(mean(D, 1, 'omitnan'));
lo = min(mu(:), [], 'omitnan');
hi = max(mu(:), [], 'omitnan');
pad = max((hi - lo) * 0.18, eps);
yl = [lo - pad, hi + pad];
end

function add_sig_bar(ax, t, mask, yl)
d = diff([0 mask 0]);
starts = find(d == 1);
ends = find(d == -1) - 1;
y = yl(2) - 0.06 * range(yl);
for c = 1:length(starts)
    plot(ax, [t(starts(c)) t(ends(c))], [y y], 'k-', 'LineWidth',2.0, 'HandleVisibility','off');
end
end

function label_clusters(ax, t, yvec, clusts, c_sig, c_ns)
yl = ylim(ax);
for ci = 1:length(clusts)
    cl = clusts(ci);
    if isempty(cl.idx), continue; end
    [~, pk] = max(abs(yvec(cl.idx)));
    pk_t = t(cl.idx(pk));
    pk_v = yvec(cl.idx(pk));
    lbl = sprintf('p = %.3f', cl.pval);
    col = c_sig;
    if ~cl.sig
        lbl = sprintf('p = %.3f ns', cl.pval);
        col = c_ns * 0.55;
    end
    yoff = 0.08 * range(yl) * sign(max(pk_v, eps));
    text(ax, pk_t, pk_v + yoff, lbl, 'FontSize',20, ...
        'HorizontalAlignment','center', 'Color',col, 'FontName','Helvetica');
end
end

function format_pub_axes(ax)
set(ax, 'FontSize',7.5, 'FontName','Helvetica', 'TickDir','out', 'Box','off', ...
    'TickLength',[0.015 0.01], 'LineWidth',0.7);
xlim(ax, [0, ax.XLim(2)]);
end

function s = titlecase(s)
parts = split(string(s));
for i = 1:numel(parts)
    p = char(parts(i));
    if ~isempty(p), parts(i) = string([upper(p(1)) lower(p(2:end))]); end
end
s = char(strjoin(parts, ' '));
end

function lbl = task_label(trial_type)
switch lower(trial_type)
    case 'lure'
        lbl = 'mnemonic discrimination';
    case 'target'
        lbl = 'same item detection';
    otherwise
        lbl = lower(trial_type);
end
end

function fill_win(ax, win_t, col, alph)
yl = ylim(ax);
x = [win_t(1) win_t(2) win_t(2) win_t(1)];
y = [yl(1) yl(1) yl(2) yl(2)];
fill(ax, x, y, col, 'EdgeColor','none','FaceAlpha',alph,'HandleVisibility','off');
% dashed vertical borders
plot(ax, [win_t(1) win_t(1)], yl, '--', 'Color',col, 'LineWidth',0.9, 'HandleVisibility','off');
plot(ax, [win_t(2) win_t(2)], yl, '--', 'Color',col, 'LineWidth',0.9, 'HandleVisibility','off');
end

function shade_clusters_F(ax, t, Fvec, mask, col, alph)
d = diff([0 mask 0]);
starts = find(d == 1); ends = find(d == -1) - 1;
for c = 1:length(starts)
    idx = starts(c):ends(c);
    xt = [t(idx(1)) t(idx(end)) t(idx(end)) t(idx(1))];
    yt = [0 0 max(Fvec(idx))*1.15 max(Fvec(idx))*1.15];
    fill(ax, xt, yt, col, 'EdgeColor','none','FaceAlpha',alph,'HandleVisibility','off');
end
end

function shade_clusters_T(ax, t, tvec, mask, col, alph)
d = diff([0 mask 0]);
starts = find(d == 1); ends = find(d == -1) - 1;
yl = ylim(ax);
for c = 1:length(starts)
    idx = starts(c):ends(c);
    xt = [t(idx(1)) t(idx(end)) t(idx(end)) t(idx(1))];
    yt = [yl(1) yl(1) yl(2) yl(2)];
    fill(ax, xt, yt, col, 'EdgeColor','none','FaceAlpha',alph,'HandleVisibility','off');
end
end

function format_panel(ax, ylabel_str, show_xlabel)
set(ax, 'FontSize',20,'FontName','Helvetica','TickDir','out','Box','off', ...
    'XLim',[0 1.5],'XTick',[0 1 1.5],'TickLength',[0.015 0.01],'LineWidth',1.5);
ylabel(ax, ylabel_str, 'FontSize',20,'FontName','Helvetica');
if ~show_xlabel
    set(ax,'XTickLabel',[]);
end
end
%% helpers

function F = rm_anova_F(x)
% x: [n_subjs x 1 x 3], one timepoint
% repeated-measures one-way ANOVA, 3 conditions
x = squeeze(x);                         % [n_subjs x 3]
x = x(~any(isnan(x),2), :);            % drop subjects with any nan at this tp
[n, k] = size(x);
grand   = mean(x(:));
cond_m  = mean(x, 1);
subj_m  = mean(x, 2);

SS_cond = n * sum((cond_m - grand).^2);
SS_subj = k * sum((subj_m - grand).^2);
SS_tot  = sum((x - grand).^2, 'all');
SS_err  = SS_tot - SS_cond - SS_subj;

df_cond = k - 1;
df_err  = (k-1)*(n-1);
F = (SS_cond/df_cond) / (SS_err/df_err);
end

function [clusts, clust_stats] = get_clusters(Fvec, f_crit)
above = Fvec > f_crit;
d = diff([0 above 0]);
starts = find(d ==  1);
ends   = find(d == -1) - 1;
clusts = {}; clust_stats = [];
for c = 1:length(starts)
    idx = starts(c):ends(c);
    clusts{end+1}      = idx;           %#ok<AGROW>
    clust_stats(end+1) = sum(Fvec(idx)); %#ok<AGROW>
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

function [res, posthoc] = run_3x2_anova_correctness(lbl, cc, ci, ic, ii, nc, ni, cond_comp, cond_iso, cond_nov)
    cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    sig = @(p) repmat('*',1,(p<.05)+(p<.01)+(p<.001));
    n = size(cc,1);
    data_tbl = table(cc, ci, ic, ii, nc, ni, ...
        'VariableNames', {'comp_corr', 'comp_incorr', 'iso_corr', 'iso_incorr', 'nov_corr', 'nov_incorr'});
    
    Condition = categorical({'comp'; 'comp'; 'iso'; 'iso'; 'nov'; 'nov'});
    Correctness = categorical({'corr'; 'incorr'; 'corr'; 'incorr'; 'corr'; 'incorr'});
    within_design = table(Condition, Correctness, 'VariableNames', {'Condition', 'Correctness'});
    
    rm = fitrm(data_tbl, 'comp_corr-nov_incorr ~ 1', 'WithinDesign', within_design);
    
    atbl = ranova(rm, 'WithinModel', 'Condition*Correctness');
    
    fprintf('\n%s:\n', lbl);
    eff_names = {'(Intercept):Condition', '(Intercept):Correctness', '(Intercept):Condition:Correctness'};
    err_names = {'Error(Condition)', 'Error(Correctness)', 'Error(Condition:Correctness)'};
    labels = {'condition', 'correctness', 'interaction'};
    fnames = {'cond', 'corr', 'int'};
    
    for i = 1:3
        F = atbl{eff_names{i}, 'F'};
        p = atbl{eff_names{i}, 'pValue'};
        df1 = atbl{eff_names{i}, 'DF'};
        df_err = atbl{err_names{i}, 'DF'};
        ss = atbl{eff_names{i}, 'SumSq'};
        ss_err = atbl{err_names{i}, 'SumSq'};
        
        eta = ss / (ss + ss_err); 
        
        fprintf('  %s: F(%d,%d) = %.2f, p = %.4f, eta2_p = %.3f %s\n', labels{i}, df1, df_err, F, p, eta, sig(p));
        res.(fnames{i}) = struct('F',F,'p',p,'eta2_p',eta,'df1',df1,'df2',df_err);
    end

    subj_col = (1:n)';
    bf_cond = table(subj_col, cond_comp, cond_iso, cond_nov, 'VariableNames',{'subj','comp','iso','nov'});
    bf_cond_long = stack(bf_cond, {'comp','iso','nov'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_cond_long.subj = categorical(bf_cond_long.subj); bf_cond_long.cond = categorical(bf_cond_long.cond);
    
    res.bf10_cond = bf.anova(bf_cond_long, 'y~cond', 'treatAsRandom',{'subj'}, 'verbose',false);
    if res.bf10_cond < 1
        fprintf('  condition BF01 = %.2f\n', 1/res.bf10_cond); 
    else
        fprintf('  condition BF10 = %.2f\n', res.bf10_cond); 
    end
    
    all_corr = mean([cc ic nc],2); all_incorr = mean([ci ii ni],2);
    bf_corr = table(subj_col, all_corr, all_incorr, 'VariableNames',{'subj','corr','incorr'});
    bf_corr_long = stack(bf_corr, {'corr','incorr'}, 'NewDataVariableName','y', 'IndexVariableName','acc');
    bf_corr_long.subj = categorical(bf_corr_long.subj); bf_corr_long.acc = categorical(bf_corr_long.acc);
    res.bf10_corr = bf.anova(bf_corr_long, 'y~acc', 'treatAsRandom',{'subj'}, 'verbose',false);
    
    if res.bf10_corr < 1
        fprintf('  correctness BF01 = %.2f\n', 1/res.bf10_corr); 
    else
        fprintf('  correctness BF10 = %.2f\n', res.bf10_corr); 
    end
    
    fprintf('  condition post-hocs:\n');
    [~,p1,~,s1]=ttest(cond_comp,cond_iso); [~,p2,~,s2]=ttest(cond_iso,cond_nov); [~,p3,~,s3]=ttest(cond_comp,cond_nov);
    posthoc.cond_ci=struct('t',s1.tstat,'df',s1.df,'d',cd(cond_comp,cond_iso),'p',p1);
    posthoc.cond_in=struct('t',s2.tstat,'df',s2.df,'d',cd(cond_iso,cond_nov),'p',p2);
    posthoc.cond_cn=struct('t',s3.tstat,'df',s3.df,'d',cd(cond_comp,cond_nov),'p',p3);
    fprintf('    comp v iso: t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s1.df,s1.tstat,cd(cond_comp,cond_iso),p1,sig(p1));
    fprintf('    iso v nov:  t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s2.df,s2.tstat,cd(cond_iso,cond_nov),p2,sig(p2));
    fprintf('    comp v nov: t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s3.df,s3.tstat,cd(cond_comp,cond_nov),p3,sig(p3));
    
    fprintf('  correctness within condition:\n');
    [~,p1,~,s1]=ttest(cc,ci,'Tail','right'); [~,p2,~,s2]=ttest(ic,ii,'Tail','right'); [~,p3,~,s3]=ttest(nc,ni,'Tail','right');
    posthoc.corr_comp=struct('t',s1.tstat,'df',s1.df,'d',cd(cc,ci),'p',p1);
    posthoc.corr_iso=struct('t',s2.tstat,'df',s2.df,'d',cd(ic,ii),'p',p2);
    posthoc.corr_nov=struct('t',s3.tstat,'df',s3.df,'d',cd(nc,ni),'p',p3);
    fprintf('    comp corr v incorr: t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s1.df,s1.tstat,cd(cc,ci),p1,sig(p1));
    fprintf('    iso corr v incorr:  t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s2.df,s2.tstat,cd(ic,ii),p2,sig(p2));
    fprintf('    nov corr v incorr:  t(%d) = %.2f, d = %.2f, p = %.4f %s\n',s3.df,s3.tstat,cd(nc,ni),p3,sig(p3));

    try
        bf_comp = bf.ttest(cc-ci); posthoc.corr_comp.bf10 = bf_comp;
        bf_iso  = bf.ttest(ic-ii); posthoc.corr_iso.bf10  = bf_iso;
        bf_nov  = bf.ttest(nc-ni); posthoc.corr_nov.bf10  = bf_nov;
        fprintf('    BF10 comp corr v incorr: %.3f\n', bf_comp);
        fprintf('    BF10 iso  corr v incorr: %.3f (BF01 = %.3f)\n', bf_iso, 1/bf_iso);
        fprintf('    BF10 nov  corr v incorr: %.3f (BF01 = %.3f)\n', bf_nov, 1/bf_nov);
    catch ME
        fprintf('    BF computation failed: %s\n', ME.message);
    end
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

function [res, posthoc] = run_rm_anova_2x2(lbl, ac, ai, bc, bi)
    within = table({'A2A1'; 'A2A1'; 'B1A2'; 'B1A2'}, {'Comp'; 'Iso'; 'Comp'; 'Iso'}, ...
        'VariableNames', {'PairType', 'TaskCond'});
    t = table(ac, ai, bc, bi, 'VariableNames', {'ac','ai','bc','bi'});
    rm = fitrm(t, 'ac-bi~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'PairType*TaskCond');
    res.df1 = tbl.DF(7); res.df2 = tbl.DF(8);
    res.F_pair = tbl.F(3); res.p_pair = tbl.pValue(3); res.eta2_pair = tbl.SumSq(3)/(tbl.SumSq(3)+tbl.SumSq(4));
    res.F_cond = tbl.F(5); res.p_cond = tbl.pValue(5); res.eta2_cond = tbl.SumSq(5)/(tbl.SumSq(5)+tbl.SumSq(6));
    res.F_int  = tbl.F(7); res.p_int  = tbl.pValue(7); res.eta2_int  = tbl.SumSq(7)/(tbl.SumSq(7)+tbl.SumSq(8));
    
    fprintf('\n%s:\n', lbl);
    fprintf('  Main PairType: F(%d,%d) = %.2f, p = %.4f, eta2_p = %.3f\n', res.df1, res.df2, res.F_pair, res.p_pair, res.eta2_pair);
    fprintf('  Main TaskCond: F(%d,%d) = %.2f, p = %.4f, eta2_p = %.3f\n', res.df1, res.df2, res.F_cond, res.p_cond, res.eta2_cond);
    fprintf('  Interaction:   F(%d,%d) = %.2f, p = %.4f, eta2_p = %.3f\n', res.df1, res.df2, res.F_int, res.p_int, res.eta2_int);
    
    n = length(ac); subj = (1:n)';
    bf_tbl = table(subj, ac, ai, bc, bi, 'VariableNames', {'subj','ac','ai','bc','bi'});
    bf_long = stack(bf_tbl, {'ac','ai','bc','bi'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj);
    bf_long.PairType = categorical(contains(string(bf_long.cond), 'a')); % true if ac/ai
    bf_long.TaskCond = categorical(contains(string(bf_long.cond), 'c')); % true if ac/bc
    
    try
        res.bf10 = bf.anova(bf_long, 'y~PairType*TaskCond', 'treatAsRandom', {'subj'}, 'verbose', false);
        if res.bf10 < 1, fprintf('  BF01 = %.2f\n', 1/res.bf10); else, fprintf('  BF10 = %.2f\n', res.bf10); end
    catch
        fprintf('  (Bayes Factor skipped: ensure toolbox supports 2x2 formulation)\n');
    end
        cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    [~,p1,~,s1] = ttest(ac, ai); [~,p2,~,s2] = ttest(bc, bi); 
    [~,p3,~,s3] = ttest(ac, bc); [~,p4,~,s4] = ttest(ai, bi);
    
    df = s1.df;
    posthoc.a2a1_comp_iso = struct('t',s1.tstat,'df',df,'d',cd(ac,ai),'p',p1);
    posthoc.b1a2_comp_iso = struct('t',s2.tstat,'df',df,'d',cd(bc,bi),'p',p2);
    posthoc.comp_a2a1_b1a2 = struct('t',s3.tstat,'df',df,'d',cd(ac,bc),'p',p3);
    posthoc.iso_a2a1_b1a2  = struct('t',s4.tstat,'df',df,'d',cd(ai,bi),'p',p4);
    
    fprintf('  post-hoc:\n');
    fprintf('    A2A1 (comp-iso):  t(%d) = %.2f, d = %.2f, p = %.4f %s\n', df, s1.tstat, cd(ac,ai), p1, repmat('*',1,(p1<.05)+(p1<.01)+(p1<.001)));
    fprintf('    B1A2 (comp-iso):  t(%d) = %.2f, d = %.2f, p = %.4f %s\n', df, s2.tstat, cd(bc,bi), p2, repmat('*',1,(p2<.05)+(p2<.01)+(p2<.001)));
    fprintf('    comp (A2A1-B1A2): t(%d) = %.2f, d = %.2f, p = %.4f %s\n', df, s3.tstat, cd(ac,bc), p3, repmat('*',1,(p3<.05)+(p3<.01)+(p3<.001)));
    fprintf('    iso  (A2A1-B1A2): t(%d) = %.2f, d = %.2f, p = %.4f %s\n', df, s4.tstat, cd(ai,bi), p4, repmat('*',1,(p4<.05)+(p4<.01)+(p4<.001)));
end
function [res, posthoc] = run_rm_anova(lbl, d1, d2, d3, within)
    t = table(d1, d2, d3, 'VariableNames', {'compared','isolated','novel'});
    rm = fitrm(t, 'compared-novel~1', 'WithinDesign', within);
    tbl = ranova(rm); eps = epsilon(rm); m = mauchly(rm);
    res.F = tbl.F(1); res.df1 = tbl.DF(1); res.df2_gg = tbl.DF(2)*eps.GreenhouseGeisser;
    res.p_gg = tbl.pValueGG(1); res.eps = eps.GreenhouseGeisser;
    res.eta2_p = tbl.SumSq(1) / (tbl.SumSq(1) + tbl.SumSq(2));
    fprintf('\n%s: F(%d,%.1f) = %.2f, p = %.4f, eta2_p = %.3f, eps = %.3f\n', ...
        lbl, res.df1, res.df2_gg, res.F, res.p_gg, res.eta2_p, res.eps);
    n = length(d1); subj = (1:n)';
    bf_tbl = table(subj, d1, d2, d3, 'VariableNames', {'subj','c1','c2','c3'});
    bf_long = stack(bf_tbl, {'c1','c2','c3'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj); bf_long.cond = categorical(bf_long.cond);
    res.bf10 = bf.anova(bf_long, 'y~cond', 'treatAsRandom', {'subj'}, 'verbose', false);
    if res.bf10 < 1, fprintf('  BF01 = %.2f\n', 1/res.bf10); else, fprintf('  BF10 = %.2f\n', res.bf10); end
    cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    [~,p1,~,s1] = ttest(d1,d2); [~,p2,~,s2] = ttest(d2,d3); [~,p3,~,s3] = ttest(d1,d3);
    adj = holm_bonf([p1 p2 p3]); df = s1.df;
    posthoc.comp_iso = struct('t',s1.tstat,'df',df,'d',cd(d1,d2),'p',p1,'p_adj',adj(1));
    posthoc.iso_nov  = struct('t',s2.tstat,'df',df,'d',cd(d2,d3),'p',p2,'p_adj',adj(2));
    posthoc.comp_nov = struct('t',s3.tstat,'df',df,'d',cd(d1,d3),'p',p3,'p_adj',adj(3));
    fprintf('  post-hoc:\n');
    fprintf('    comp-iso: t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s1.tstat, cd(d1,d2), adj(1), repmat('*',1,(adj(1)<.05)+(adj(1)<.01)+(adj(1)<.001)));
    fprintf('    iso-nov:  t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s2.tstat, cd(d2,d3), adj(2), repmat('*',1,(adj(2)<.05)+(adj(2)<.01)+(adj(2)<.001)));
    fprintf('    comp-nov: t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s3.tstat, cd(d1,d3), adj(3), repmat('*',1,(adj(3)<.05)+(adj(3)<.01)+(adj(3)<.001)));
end

function adj = holm_bonf(p)
    n = length(p); [ps, idx] = sort(p); as = min(1, ps .* (n:-1:1));
    for i = 2:n, as(i) = max(as(i), as(i-1)); end
    adj(idx) = as;
end


function [clusts, clust_stats] = get_clusters_F(Fvec, f_crit)
above = Fvec > f_crit;
d = diff([0 above 0]);
starts = find(d ==  1); ends = find(d == -1) - 1;
clusts = {}; clust_stats = [];
for c = 1:length(starts)
    idx = starts(c):ends(c);
    clusts{end+1}      = idx;            %#ok<AGROW>
    clust_stats(end+1) = sum(Fvec(idx)); %#ok<AGROW>
end
end

function [clusts, clust_stats] = get_clusters_T(tvec, t_crit)
% separate positive and negative clusters
clusts = {}; clust_stats = [];
for sign_dir = [1 -1]
    above = (sign_dir * tvec) > t_crit;
    d = diff([0 above 0]);
    starts = find(d ==  1); ends = find(d == -1) - 1;
    for c = 1:length(starts)
        idx = starts(c):ends(c);
        clusts{end+1}      = idx;            %#ok<AGROW>
        clust_stats(end+1) = sum(tvec(idx)); %#ok<AGROW>
    end
end
end
