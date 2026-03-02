%% Load and plot cumulative reinstatement results
clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('Loading cumulative reinstatement results...\n');
load(fullfile(res_dir, 'gaze_entropy_temporal/cumulative_reinstatement_results.mat'), 'cumulative_res');

cumulative_results_comp = cumulative_res.compared;
cumulative_results_iso = cumulative_res.isolated;
subj_traj_comp = cumulative_res.subj_traj_comp;
subj_traj_iso = cumulative_res.subj_traj_iso;
n_fix_to_plot = cumulative_res.n_fix_to_plot;
early_range = cumulative_res.early_range;
late_range = cumulative_res.late_range;

mean_comp = mean(subj_traj_comp, 1, 'omitnan');
sem_comp = std(subj_traj_comp, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_comp), 1));
mean_iso = mean(subj_traj_iso, 1, 'omitnan');
sem_iso = std(subj_traj_iso, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(subj_traj_iso), 1));
fix_numbers = 1:n_fix_to_plot;

%% Main trajectory plot
figure('Position', [100, 100, 800, 500]); hold on;
fill([fix_numbers, fliplr(fix_numbers)], [mean_comp + sem_comp, fliplr(mean_comp - sem_comp)], [0.2, 0.6, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([fix_numbers, fliplr(fix_numbers)], [mean_iso + sem_iso, fliplr(mean_iso - sem_iso)], [0.8, 0.4, 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(fix_numbers, mean_comp, 'o-', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.6, 0.8]);
plot(fix_numbers, mean_iso, 's-', 'Color', [0.8, 0.4, 0.4], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.8, 0.4, 0.4]);
plot([0.5, n_fix_to_plot + 0.5], [0, 0], 'k--', 'LineWidth', 1);
xlabel('Cumulative Fixation Number', 'FontSize', 12);
ylabel('Gaze Reinstatement (baseline-corrected)', 'FontSize', 12);
title('Fixation-by-Fixation Reinstatement Buildup', 'FontSize', 14);
legend({'', '', 'Compared', 'Isolated'}, 'Location', 'best', 'FontSize', 11);
xlim([0.5, n_fix_to_plot + 0.5]); box off; grid on;
print(gcf, fullfile(fig_dir, 'cumulative_reinstatement.pdf'), '-dpdf', '-vector');

%% Early vs Late comparison
early_comp = mean(subj_traj_comp(:, early_range), 2, 'omitnan');
late_comp = mean(subj_traj_comp(:, late_range), 2, 'omitnan');
early_iso = mean(subj_traj_iso(:, early_range), 2, 'omitnan');
late_iso = mean(subj_traj_iso(:, late_range), 2, 'omitnan');

valid = ~isnan(early_comp) & ~isnan(early_iso);
[~, p_early] = ttest(early_comp(valid), early_iso(valid));
fprintf('\nEarly (Fix %d-%d): Compared M=%.3f, Isolated M=%.3f, p=%.4f\n', early_range(1), early_range(end), mean(early_comp(valid)), mean(early_iso(valid)), p_early);

valid = ~isnan(late_comp) & ~isnan(late_iso);
[~, p_late] = ttest(late_comp(valid), late_iso(valid));
fprintf('Late (Fix %d-%d): Compared M=%.3f, Isolated M=%.3f, p=%.4f\n', late_range(1), late_range(end), mean(late_comp(valid)), mean(late_iso(valid)), p_late);

figure('Position', [100, 100, 900, 400]);
subplot(1, 2, 1); hold on;
valid = ~isnan(early_comp) & ~isnan(early_iso);
bar_data = [mean(early_comp(valid)), mean(early_iso(valid))];
bar_err = [std(early_comp(valid))/sqrt(sum(valid)), std(early_iso(valid))/sqrt(sum(valid))];
b = bar(bar_data, 'FaceColor', 'flat'); b.CData(1,:) = [0.2, 0.6, 0.8]; b.CData(2,:) = [0.8, 0.4, 0.4];
errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTickLabel', {'Compared', 'Isolated'}); ylabel('Mean Reinstatement (r)', 'FontSize', 12);
title(sprintf('Early (Fix %d-%d)', early_range(1), early_range(end)), 'FontSize', 12);
if p_early < 0.05
    y_sig = max(bar_data) + max(bar_err) + 0.05; plot([1, 2], [y_sig, y_sig], 'k-', 'LineWidth', 1);
    text(1.5, y_sig + 0.02, sprintf('p=%.3f', p_early), 'HorizontalAlignment', 'center', 'FontSize', 10);
end
box off;

subplot(1, 2, 2); hold on;
valid = ~isnan(late_comp) & ~isnan(late_iso);
bar_data = [mean(late_comp(valid)), mean(late_iso(valid))];
bar_err = [std(late_comp(valid))/sqrt(sum(valid)), std(late_iso(valid))/sqrt(sum(valid))];
b = bar(bar_data, 'FaceColor', 'flat'); b.CData(1,:) = [0.2, 0.6, 0.8]; b.CData(2,:) = [0.8, 0.4, 0.4];
errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTickLabel', {'Compared', 'Isolated'}); ylabel('Mean Reinstatement (r)', 'FontSize', 12);
title(sprintf('Late (Fix %d-%d)', late_range(1), late_range(end)), 'FontSize', 12);
if p_late < 0.05
    y_sig = max(bar_data) + max(bar_err) + 0.05; plot([1, 2], [y_sig, y_sig], 'k-', 'LineWidth', 1);
    text(1.5, y_sig + 0.02, sprintf('p=%.3f', p_late), 'HorizontalAlignment', 'center', 'FontSize', 10);
end
box off;

sgtitle('Cumulative Reinstatement: Early vs Late', 'FontSize', 14, 'FontWeight', 'bold');
print(gcf, fullfile(fig_dir, 'cumulative_reinstatement_early_late.pdf'), '-dpdf', '-vector');

fprintf('\nDone!\n');