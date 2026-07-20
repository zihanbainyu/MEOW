%% Subsequent memory: A1-B2 reinstatement predicting B2 accuracy
clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');

c_comp = [180 174 211]/255; c_iso = [176 230 255]/255;

%% Load baseline-corrected A1-B2 results
load(fullfile(res_dir, 'gaze_reinstat_res_a1b2.mat'), 'reinstat_res_a1b2');

results_comp = reinstat_res_a1b2.a1b2_compared;
results_iso  = reinstat_res_a1b2.a1b2_isolated;

fprintf('A1-B2 pairs loaded: comp=%d, iso=%d\n', height(results_comp), height(results_iso));
fprintf('B2 correct rate: comp=%.1f%%, iso=%.1f%%\n', ...
    100*mean(results_comp.correct, 'omitnan'), ...
    100*mean(results_iso.correct, 'omitnan'));

%% Subject-level means by condition x accuracy
subj_ids = unique([results_comp.subj_id; results_iso.subj_id]);
n_subj = length(subj_ids);
sim_comp_corr = nan(n_subj, 1); sim_comp_incr = nan(n_subj, 1);
sim_iso_corr  = nan(n_subj, 1); sim_iso_incr  = nan(n_subj, 1);
sim_comp_all  = nan(n_subj, 1); sim_iso_all   = nan(n_subj, 1);

for s = 1:n_subj
    sid = subj_ids(s);
    cc = results_comp.reinst_index(results_comp.subj_id == sid & results_comp.correct == 1);
    ci = results_comp.reinst_index(results_comp.subj_id == sid & results_comp.correct == 0);
    ic = results_iso.reinst_index(results_iso.subj_id == sid & results_iso.correct == 1);
    ii_t = results_iso.reinst_index(results_iso.subj_id == sid & results_iso.correct == 0);
    ca = results_comp.reinst_index(results_comp.subj_id == sid);
    ia = results_iso.reinst_index(results_iso.subj_id == sid);
    if ~isempty(cc), sim_comp_corr(s) = mean(cc, 'omitnan'); end
    if ~isempty(ci), sim_comp_incr(s) = mean(ci, 'omitnan'); end
    if ~isempty(ic), sim_iso_corr(s)  = mean(ic, 'omitnan'); end
    if ~isempty(ii_t), sim_iso_incr(s) = mean(ii_t, 'omitnan'); end
    if ~isempty(ca), sim_comp_all(s) = mean(ca, 'omitnan'); end
    if ~isempty(ia), sim_iso_all(s)  = mean(ia, 'omitnan'); end
end

%% 2x2 ANOVA (condition x correctness)
v = ~isnan(sim_comp_corr) & ~isnan(sim_comp_incr) & ~isnan(sim_iso_corr) & ~isnan(sim_iso_incr);
fprintf('\nSubjects with all 4 cells: %d\n', sum(v));

fprintf('\n=== A1-B2 reinstatement (baseline-corrected) by B2 accuracy ===\n');
t = table(sim_comp_corr(v), sim_comp_incr(v), sim_iso_corr(v), sim_iso_incr(v), ...
    'VariableNames', {'cc','ci','ic','ii'});
within = table({'comp';'comp';'iso';'iso'}, {'corr';'incorr';'corr';'incorr'}, ...
    'VariableNames', {'cond','correctness'});
rm = fitrm(t, 'cc-ii~1', 'WithinDesign', within);
tbl = ranova(rm, 'WithinModel', 'cond*correctness');

rows = {'(Intercept):cond', '(Intercept):correctness', '(Intercept):cond:correctness'};
labels = {'condition', 'correctness', 'interaction'};
term_names = string(tbl.Properties.RowNames);
for i = 1:3
    idx = find(contains(term_names, rows{i}), 1);
    if ~isempty(idx)
        f = tbl.F(idx); p_val = tbl.pValue(idx);
        df1 = tbl.DF(idx); df2 = tbl.DF(idx+1);
        s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1);
        eta = s_ef / (s_ef + s_er);
        fprintf('  %s: F(%d,%d) = %.2f, p = %.4f, eta2 = %.3f %s\n', ...
            labels{i}, df1, df2, f, p_val, eta, ...
            repmat('*', 1, (p_val<0.05) + (p_val<0.01) + (p_val<0.001)));
    end
end

%% Post-hoc paired t-tests
fprintf('\n  Post-hoc tests:\n');
[~, p1, ~, s1] = ttest(sim_comp_corr(v), sim_comp_incr(v));
[~, p2, ~, s2] = ttest(sim_iso_corr(v),  sim_iso_incr(v));
[~, pc, ~, sc] = ttest(sim_comp_all(v),   sim_iso_all(v));
[~, p3, ~, s3] = ttest(sim_comp_corr(v)-sim_comp_incr(v), sim_iso_corr(v)-sim_iso_incr(v));
n_v = sum(v);
fprintf('    comp corr > incorr: t(%d) = %.2f, p = %.4f, d = %.2f %s\n', ...
    s1.df, s1.tstat, p1, s1.tstat/sqrt(n_v), repmat('*', 1, (p1<0.05)+(p1<0.01)+(p1<0.001)));
fprintf('    iso corr > incorr:  t(%d) = %.2f, p = %.4f, d = %.2f %s\n', ...
    s2.df, s2.tstat, p2, s2.tstat/sqrt(n_v), repmat('*', 1, (p2<0.05)+(p2<0.01)+(p2<0.001)));
fprintf('    comp > iso (all):   t(%d) = %.2f, p = %.4f, d = %.2f %s\n', ...
    sc.df, sc.tstat, pc, sc.tstat/sqrt(n_v), repmat('*', 1, (pc<0.05)+(pc<0.01)+(pc<0.001)));
fprintf('    (comp-iso) diff:    t(%d) = %.2f, p = %.4f, d = %.2f %s\n', ...
    s3.df, s3.tstat, p3, s3.tstat/sqrt(n_v), repmat('*', 1, (p3<0.05)+(p3<0.01)+(p3<0.001)));

%% Cell means
fprintf('\n  Cell means (baseline-corrected):\n');
fprintf('    comp correct:   M = %.4f (SD = %.4f)\n', mean(sim_comp_corr(v)), std(sim_comp_corr(v)));
fprintf('    comp incorrect: M = %.4f (SD = %.4f)\n', mean(sim_comp_incr(v)), std(sim_comp_incr(v)));
fprintf('    iso correct:    M = %.4f (SD = %.4f)\n', mean(sim_iso_corr(v)),  std(sim_iso_corr(v)));
fprintf('    iso incorrect:  M = %.4f (SD = %.4f)\n', mean(sim_iso_incr(v)), std(sim_iso_incr(v)));

%% Raincloud plot
data_sme = [sim_comp_corr(v), sim_comp_incr(v), sim_iso_corr(v), sim_iso_incr(v)];
cols = {c_comp, c_comp*0.5+0.5, c_iso, c_iso*0.5+0.5};
xlbls = {'Comp correct', 'Comp incorrect', 'Iso correct', 'Iso incorrect'};

figure('Color', 'w', 'Position', [50 50 900 600]);
hold on;

[n_rows, n_grps] = size(data_sme);
d_v = data_sme(:); d_v = d_v(~isnan(d_v));
mn = min(d_v); mx = max(d_v); rng_v = mx - mn; if rng_v == 0, rng_v = 1; end

jit = -0.15 - (rand(size(data_sme)) * 0.20);
x_c = repmat(1:n_grps, n_rows, 1) + jit;
plot(x_c', data_sme', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);

for i = 1:n_grps
    d = data_sme(:, i); d_clean = d(~isnan(d));
    if isempty(d_clean), continue; end
    [f, xi] = ksdensity(d_clean);
    f = f / max(f) * 0.4;
    patch([i + f, i * ones(1, length(f))], [xi, fliplr(xi)], cols{i}, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    scatter(x_c(:, i), data_sme(:, i), 20, cols{i}, 'filled', 'MarkerFaceAlpha', 0.6);
    q = quantile(d_clean, [0.25, 0.5, 0.75]);
    rectangle('Position', [i - 0.06, q(1), 0.12, q(3) - q(1)], ...
        'FaceColor', cols{i}, 'EdgeColor', 'k', 'LineWidth', 1.2);
    plot([i - 0.06, i + 0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
end

yline(0, 'r--', 'LineWidth', 2);

% Significance brackets
pairs = [1 2; 3 4; 1 3; 2 4];
cl = ylim; y_top = cl(2); step = rng_v * 0.08;
real_mx = max(data_sme(:)); base = max(y_top, real_mx + step * 0.5);
line_lvl = 0;
for i = 1:size(pairs, 1)
    c1 = pairs(i, 1); c2 = pairs(i, 2);
    [~, p] = ttest(data_sme(:, c1), data_sme(:, c2));
    if p < 0.05
        line_lvl = line_lvl + 1;
        y_p = base + (line_lvl * step);
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
        plot([c1, c1, c2, c2], [y_p - step*0.3, y_p, y_p, y_p - step*0.3], 'k-', 'LineWidth', 1.2);
        text(mean([c1, c2]), y_p + step*0.1, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 18, 'FontWeight', 'bold');
    end
end
if line_lvl > 0, ylim([cl(1), base + (line_lvl * step) + step]); end

set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 18);
ylabel('A1-B2 Gaze Similarity', 'FontSize', 20);
xlim([0.2, n_grps + 0.8]);
box off;
hold off;

print(gcf, fullfile(fig_dir, 'fig_a1b2.pdf'), '-dpdf', '-vector');
fprintf('\nFigure saved.\n');