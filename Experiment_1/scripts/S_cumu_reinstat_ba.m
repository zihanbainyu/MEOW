clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results', 'cumulative_reinstatement_ba');
load(fullfile(res_dir, 'cumulative_reinstatement_ba_results.mat'));

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

fprintf('\n--- Paired t-test (compared vs isolated slopes) ---\n');
[~, p_diff, ci_diff, st_diff] = ttest(slopes_c, slopes_i);
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

figure('Position', [100, 100, 900, 450]);

subplot(1, 2, 1); hold on;
fill([fix_numbers, fliplr(fix_numbers)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], ...
    [0.2, 0.6, 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([fix_numbers, fliplr(fix_numbers)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], ...
    [0.8, 0.4, 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(fix_numbers, mean_comp, 'o', 'Color', [0.2, 0.6, 0.8], 'MarkerSize', 9, 'MarkerFaceColor', [0.2, 0.6, 0.8], 'LineWidth', 1.5);
plot(fix_numbers, mean_iso, 's', 'Color', [0.8, 0.4, 0.4], 'MarkerSize', 9, 'MarkerFaceColor', [0.8, 0.4, 0.4], 'LineWidth', 1.5);
plot(x_fit, mean(slopes_c) * x_fit_c + mean(mean_comp), '-', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2.5);
plot(x_fit, mean(slopes_i) * x_fit_c + mean(mean_iso), '-', 'Color', [0.8, 0.4, 0.4], 'LineWidth', 2.5);
plot([0.5, n_fix_to_plot + 0.5], [0, 0], 'k--', 'LineWidth', 1);
xlabel('Cumulative Fixation Number', 'FontSize', 11);
ylabel('Cross-Item Reinstatement', 'FontSize', 11);
title('B\rightarrowA Reinstatement Buildup', 'FontSize', 13);
legend({'', '', 'Compared', 'Isolated'}, 'Location', 'best', 'FontSize', 10);
xlim([0.5, n_fix_to_plot + 0.5]); box off; grid on;

subplot(1, 2, 2); hold on;
bar_data = [mean(slopes_c), mean(slopes_i)];
bar_err = [std(slopes_c)/sqrt(n_subj), std(slopes_i)/sqrt(n_subj)];
b = bar(bar_data, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.6);
b.CData(1,:) = [0.2, 0.6, 0.8]; b.CData(2,:) = [0.8, 0.4, 0.4];
errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
scatter(ones(n_subj, 1) * 0.85 + rand(n_subj, 1) * 0.3, slopes_c, 20, [0.2, 0.6, 0.8], 'filled', 'MarkerFaceAlpha', 0.4);
scatter(ones(n_subj, 1) * 1.85 + rand(n_subj, 1) * 0.3, slopes_i, 20, [0.8, 0.4, 0.4], 'filled', 'MarkerFaceAlpha', 0.4);
for s = 1:n_subj
    plot([1, 2], [slopes_c(s), slopes_i(s)], '-', 'Color', [0.5 0.5 0.5 0.15], 'LineWidth', 0.5);
end
set(gca, 'XTick', 1:2, 'XTickLabel', {'Compared', 'Isolated'});
ylabel('Reinstatement Buildup Rate (slope)', 'FontSize', 11);
title('Per-Subject Slopes', 'FontSize', 13);
if p_diff < 0.05
    y_sig = max(bar_data) + max(bar_err) + 0.01;
    plot([1, 2], [y_sig, y_sig], 'k-', 'LineWidth', 1.5);
    if p_diff < 0.001, p_str = '***';
    elseif p_diff < 0.01, p_str = '**';
    else, p_str = '*';
    end
    text(1.5, y_sig + 0.003, p_str, 'HorizontalAlignment', 'center', 'FontSize', 14);
end
box off;

print(gcf, fullfile(res_dir, 'cumulative_reinstatement_ba_slopes.pdf'), '-dpdf', '-vector');

save(fullfile(res_dir, 'cumulative_reinstatement_ba_slope_stats.mat'), ...
    'slopes_c', 'slopes_i', 'common_subjs', 'p_diff', 'st_diff', 'd_diff');

fprintf('\nDone.\n');