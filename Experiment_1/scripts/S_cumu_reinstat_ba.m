clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');
load(fullfile(res_dir, 'cumu_reinstat_ba.mat'));

fix_cols_comp = cumulative_results_comp{:, 4:3+n_fix_to_plot};
fix_cols_iso = cumulative_results_iso{:, 4:3+n_fix_to_plot};
x = (1:n_fix_to_plot)';
x_c = x - mean(x);

slope_comp = nan(height(cumulative_results_comp), 1);
for i = 1:height(cumulative_results_comp)
    y = fix_cols_comp(i, :)';
    valid = ~isnan(y);
    if sum(valid) >= 2
        b = x_c(valid) \ y(valid);
        slope_comp(i) = b;
    end
end

slope_iso = nan(height(cumulative_results_iso), 1);
for i = 1:height(cumulative_results_iso)
    y = fix_cols_iso(i, :)';
    valid = ~isnan(y);
    if sum(valid) >= 2
        b = x_c(valid) \ y(valid);
        slope_iso(i) = b;
    end
end

T_comp = table(cumulative_results_comp.subj_id, slope_comp, 'VariableNames', {'subj_id','slope'});
T_iso = table(cumulative_results_iso.subj_id, slope_iso, 'VariableNames', {'subj_id','slope'});
T_comp = T_comp(~isnan(T_comp.slope), :);
T_iso = T_iso(~isnan(T_iso.slope), :);

subj_slope_comp = grpstats(T_comp, 'subj_id', 'mean', 'DataVars', 'slope');
subj_slope_iso = grpstats(T_iso, 'subj_id', 'mean', 'DataVars', 'slope');
subj_slope_comp.Properties.VariableNames{'mean_slope'} = 'slope';
subj_slope_iso.Properties.VariableNames{'mean_slope'} = 'slope';

common_subjs = intersect(subj_slope_comp.subj_id, subj_slope_iso.subj_id);
[~, ic] = ismember(common_subjs, subj_slope_comp.subj_id);
[~, ii] = ismember(common_subjs, subj_slope_iso.subj_id);
slopes_c = subj_slope_comp.slope(ic);
slopes_i = subj_slope_iso.slope(ii);
n_subj = length(common_subjs);

fprintf('=== Per-Trial Slope Analysis ===\n');
fprintf('N subjects = %d\n\n', n_subj);

fprintf('--- One-sample t-tests (slope > 0) ---\n');
[~, p_c, ~, st_c] = ttest(slopes_c, 0, 'Tail', 'right');
[~, p_i, ~, st_i] = ttest(slopes_i, 0, 'Tail', 'right');
d_c = mean(slopes_c) / std(slopes_c);
d_i = mean(slopes_i) / std(slopes_i);
fprintf('Compared:  M=%.4f (SD=%.4f), t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    mean(slopes_c), std(slopes_c), st_c.df, st_c.tstat, p_c, d_c);
fprintf('Isolated:  M=%.4f (SD=%.4f), t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    mean(slopes_i), std(slopes_i), st_i.df, st_i.tstat, p_i, d_i);

fprintf('\n--- Paired t-test (compared > isolated slopes) ---\n');
[~, p_diff, ci_diff, st_diff] = ttest(slopes_c, slopes_i, 'Tail', 'right');
d_diff = mean(slopes_c - slopes_i) / std(slopes_c - slopes_i);
fprintf('Difference: M=%.4f (SD=%.4f), t(%d)=%.3f, p=%.4f, d=%.3f\n', ...
    mean(slopes_c - slopes_i), std(slopes_c - slopes_i), st_diff.df, st_diff.tstat, p_diff, d_diff);
fprintf('95%% CI: [%.4f, %.4f]\n', ci_diff(1), ci_diff(2));

[p_wilcox, ~, wstats] = signrank(slopes_c, slopes_i);
fprintf('\nWilcoxon signed-rank: z=%.3f, p=%.4f\n', wstats.zval, p_wilcox);

%% Plot
valid_all = ~any(isnan(subj_traj_comp), 2) & ~any(isnan(subj_traj_iso), 2);
data_comp = subj_traj_comp(valid_all, :);
data_iso = subj_traj_iso(valid_all, :);

mean_comp = mean(data_comp, 1);
sem_comp = std(data_comp, 0, 1) / sqrt(n_subj);
mean_iso = mean(data_iso, 1);
sem_iso = std(data_iso, 0, 1) / sqrt(n_subj);
fix_numbers = 1:n_fix_to_plot;
x_fit = linspace(1, n_fix_to_plot, 50);
x_fit_c = x_fit - mean(x);

comp_color = [0.58, 0.58, 0.78];
iso_color = [0.62, 0.86, 0.96];
comp_fit = mean(slopes_c) * x_fit_c + mean(mean_comp);
iso_fit = mean(slopes_i) * x_fit_c + mean(mean_iso);

p_c_str = regexprep(sprintf('%.2f', p_c), '^0', '');
p_i_str = regexprep(sprintf('%.2f', p_i), '^0', '');

figure('Position', [100, 100, 250, 320], 'Color', 'w');
hold on;
fill([fix_numbers, fliplr(fix_numbers)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], ...
    comp_color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([fix_numbers, fliplr(fix_numbers)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], ...
    iso_color, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_comp = plot(x_fit, comp_fit, '-', 'Color', comp_color, 'LineWidth', 2.5);
h_iso = plot(x_fit, iso_fit, '-', 'Color', iso_color, 'LineWidth', 2.5);
plot([0.9, n_fix_to_plot + 0.35], [0, 0], '--', 'Color', [0.45, 0.45, 0.45], ...
    'LineWidth', 1, 'HandleVisibility', 'off');

xlim([0.95, n_fix_to_plot + 0.55]);
ylim([-0.018, max([0.08, mean_comp + sem_comp, mean_iso + sem_iso]) + 0.005]);
set(gca, 'XTick', 1:n_fix_to_plot, 'YTick', 0:0.02:0.06, 'FontSize', 10, ...
    'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
xlabel('Cumulative fixation number', 'FontSize', 10);
ylabel('');

legend([h_comp, h_iso], {['compared: p = ' p_c_str], ['isolated:    p = ' p_i_str]}, ...
    'Location', 'southwest', 'Box', 'off', 'FontSize', 9);

x_bracket = n_fix_to_plot + 0.22;
y_bracket_low = min(comp_fit(end), iso_fit(end));
y_bracket_high = max(comp_fit(end), iso_fit(end));
plot([x_bracket, x_bracket], [y_bracket_low, y_bracket_high], 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
plot([x_bracket - 0.06, x_bracket], [y_bracket_low, y_bracket_low], 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
plot([x_bracket - 0.06, x_bracket], [y_bracket_high, y_bracket_high], 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
text(x_bracket + 0.08, mean([y_bracket_low, y_bracket_high]), '*', ...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

set(gcf, 'PaperPositionMode', 'auto');

print(gcf, fullfile(res_dir, 'cumulative_reinstatement_ba_slopes.pdf'), '-dpdf', '-vector');

save(fullfile(res_dir, 'cumulative_reinstatement_ba_slope_stats.mat'), ...
    'slopes_c', 'slopes_i', 'common_subjs', 'p_diff', 'st_diff', 'd_diff');

fprintf('\nDone.\n');