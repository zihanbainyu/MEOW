%% Load parsed data
% fprintf('=== Loading parsed trial data ===\n');
% base_dir = '..';
% res_dir = fullfile(base_dir, 'results');

% load(fullfile(res_dir, 'all_trials_raw.mat'), 'all_trials_table');
% fprintf('Loaded %d trials\n', height(all_trials_table));
% 
% %% Preprocess pupil data (blink interpolation + filtering per trial)
% fprintf('\n=== Preprocessing Pupil Data ===\n');
% 
% % Add PCDM helper functions to path
% pcdm_dir = fullfile(base_dir, 'scripts/PCDM-main');
% addpath(genpath(pcdm_dir));
% 
% all_preprocessed = all_trials_table;
% all_preprocessed.pupil_preprocessed = cell(height(all_preprocessed), 1);
% all_preprocessed.baseline_pupil_preprocessed = cell(height(all_preprocessed), 1);
% all_preprocessed.preprocess_success = false(height(all_preprocessed), 1);
% 
% for tr_idx = 1:height(all_preprocessed)
%     if mod(tr_idx, 100) == 0
%         fprintf('Processing trial %d/%d...\n', tr_idx, height(all_preprocessed));
%     end
% 
%     trial = all_preprocessed(tr_idx, :);
%     pupil_raw = trial.pupil{1};
%     baseline_raw = trial.baseline_pupil{1};
%     sample_rate = trial.sample_rate;
% 
%     % Check if trial has enough valid data
%     valid_pct = sum(~isnan(pupil_raw) & pupil_raw > 0) / length(pupil_raw);
% 
%     if valid_pct < 0.5  % Less than 50% valid data
%         all_preprocessed.pupil_preprocessed{tr_idx} = nan(size(pupil_raw));
%         all_preprocessed.baseline_pupil_preprocessed{tr_idx} = nan(size(baseline_raw));
%         all_preprocessed.preprocess_success(tr_idx) = false;
%         continue;
%     end
% 
%     try
%         % Concatenate baseline + trial for this trial
%         concat_pupil = [baseline_raw; pupil_raw];
%         baseline_len = length(baseline_raw);
%         trial_len = length(pupil_raw);
% 
%         % Step 1: Blink interpolation
%         pupil_interp = concat_pupil;
%         pupil_interp(isnan(pupil_interp)) = 0;
% 
%         % Detect blinks by velocity threshold
%         velocity_thresh = sample_rate / 10;
%         pupil_diff = abs(diff(pupil_interp));
%         aboveT = find(pupil_diff > velocity_thresh);
%         pupil_interp(aboveT) = 0;
% 
%         % Interpolate blinks
%         pupil_interp = blinkinterp(pupil_interp', sample_rate, 5, 4, 50, 75);
% 
%         % Handle non-finite values from interpolation
%         if any(~isfinite(pupil_interp))
%             pupil_interp(~isfinite(pupil_interp)) = NaN;
%             pupil_interp = fillmissing(pupil_interp, 'linear', 'EndValues', 'nearest');
%         end
% 
%         % Final check for finite values
%         if ~all(isfinite(pupil_interp))
%             error('Still contains non-finite values after fillmissing');
%         end
% 
%         % Step 2: Bandpass filter (0.03 - 10 Hz)
%         baseline = nanmean(pupil_interp);
%         pupil_filtered = myBWfilter(pupil_interp, [0.03, 10], sample_rate, 'bandpass');
% 
%         % Step 3: Split back into baseline and trial
%         all_preprocessed.baseline_pupil_preprocessed{tr_idx} = pupil_filtered(1:baseline_len);
%         all_preprocessed.pupil_preprocessed{tr_idx} = pupil_filtered(baseline_len+1:end);
%         all_preprocessed.preprocess_success(tr_idx) = true;
% 
%     catch ME
%         % If preprocessing fails, store NaN and mark as failed
%         all_preprocessed.pupil_preprocessed{tr_idx} = nan(size(pupil_raw));
%         all_preprocessed.baseline_pupil_preprocessed{tr_idx} = nan(size(baseline_raw));
%         all_preprocessed.preprocess_success(tr_idx) = false;
%     end
% end
% 
% fprintf('\n=== Preprocessing Complete ===\n');
% fprintf('Total trials: %d\n', height(all_preprocessed));
% 
% % Statistics
% n_success = sum(all_preprocessed.preprocess_success);
% fprintf('Successfully preprocessed: %d (%.1f%%)\n', n_success, 100*n_success/height(all_preprocessed));
% 
% fprintf('\nPreprocessing success by condition:\n');
% for cond = ["compared", "isolated", "novel"]
%     for goal = ["A-A", "A-B"]
%         idx = strcmp(all_preprocessed.condition, cond) & strcmp(all_preprocessed.goal, goal);
%         n_total = sum(idx);
%         n_success = sum(idx & all_preprocessed.preprocess_success);
%         fprintf('  %s %s: %d/%d (%.1f%%)\n', goal, cond, n_success, n_total, 100*n_success/n_total);
%     end
% end
% 
% % Save complete dataset with preprocessed pupil
% save(fullfile(res_dir, 'all_trials_preprocessed.mat'), 'all_preprocessed', '-v7.3');
% fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_preprocessed.mat'));

%% Load preprocessed data
% fprintf('=== Loading preprocessed data ===\n');
% base_dir = '..';
% res_dir = fullfile(base_dir, 'results');
% 
% load(fullfile(res_dir, 'all_trials_preprocessed.mat'), 'all_preprocessed');
% fprintf('Loaded %d trials\n', height(all_preprocessed));

%% Visualize raw vs preprocessed pupil data
% fprintf('\n=== Visualizing Raw vs Preprocessed Pupil ===\n');
% 
% % Get one example trial from each condition
% conditions = ["compared", "isolated", "novel"];
% colors = {[180 174 211]/255, [176 230 255]/255, [183 210 205]/255};
% 
% % A-B trials
% figure('color','w','position',[50 50 1400 800]);
% 
% for c_idx = 1:length(conditions)
%     cond = conditions(c_idx);
% 
%     % Get one trial for this condition (only successful preprocessing)
%     trials_cond = all_preprocessed(strcmp(all_preprocessed.condition, cond) & ...
%                                    strcmp(all_preprocessed.goal, 'A-B') & ...
%                                    all_preprocessed.preprocess_success, :);
% 
%     if height(trials_cond) > 0
%         % Take first trial
%         trial = trials_cond(1, :);
%         time_vec = (0:length(trial.pupil{1})-1) / trial.sample_rate;
% 
%         % Raw data
%         subplot(2, 3, c_idx); hold on;
%         plot(time_vec, trial.pupil{1}, 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
%         title(sprintf('%s (raw)\nSubj %d, Trial %d', cond, trial.subj_id, trial.trial_id_in_run), 'FontSize', 11);
%         xlabel('Time (s)', 'FontSize', 10);
%         ylabel('Pupil size', 'FontSize', 10);
%         xlim([0 1.5]); ylim([5000 7200]);
%         grid on; box off;
% 
%         % Preprocessed data
%         subplot(2, 3, c_idx + 3); hold on;
%         plot(time_vec, trial.pupil_preprocessed{1}, 'Color', colors{c_idx}, 'LineWidth', 2);
%         title(sprintf('%s (preprocessed)\nSubj %d, Trial %d', cond, trial.subj_id, trial.trial_id_in_run), 'FontSize', 11);
%         xlabel('Time (s)', 'FontSize', 10);
%         ylabel('Baseline-Corrected Pupil size', 'FontSize', 10);
%         xlim([0 1.5]); ylim([-1000 300]);
%         grid on; box off;
%     end
% end
% 
% sgtitle('Representative Pupil Responses (Raw vs.Preprocessed)', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(res_dir, 'rep_pupil_ts_raw_prep.pdf'), '-dpdf', '-vector');

% %% Group-level visualization of preprocessed pupil data
fprintf('\n=== Group-Level Pupil Analysis ===\n');

% load(fullfile(res_dir, 'all_trials_preprocessed.mat'));
% Colors
c_comp = [180 174 211]/255; 
c_iso = [176 230 255]/255; 
c_nov = [183 210 205]/255;

% Filter for successful preprocessing only
valid_trials = all_preprocessed(all_preprocessed.preprocess_success, :);

fprintf('Valid trials for analysis: %d/%d (%.1f%%)\n', ...
    height(valid_trials), height(all_preprocessed), 100*height(valid_trials)/height(all_preprocessed));

% Find maximum trial length
max_samples = max(cellfun(@length, valid_trials.pupil_preprocessed));
sample_rate = valid_trials.sample_rate(1);
time_vec = (0:max_samples-1) / sample_rate;

fprintf('Max trial length: %d samples (%.2f s)\n', max_samples, max_samples/sample_rate);

%% Extract trial-averaged pupil by subject and condition

unique_subjs = unique(valid_trials.subj_id);
n_subj = length(unique_subjs);

% Initialize storage (subjects x timepoints)
pupil_ab_comp = nan(n_subj, max_samples);
pupil_ab_iso = nan(n_subj, max_samples);
pupil_ab_nov = nan(n_subj, max_samples);
pupil_aa_comp = nan(n_subj, max_samples);
pupil_aa_iso = nan(n_subj, max_samples);
pupil_aa_nov = nan(n_subj, max_samples);

for s_idx = 1:n_subj
    subj_id = unique_subjs(s_idx);

    % Get trials for this subject
    subj_trials = valid_trials(valid_trials.subj_id == subj_id, :);

    % Helper function to average trials by condition
    avg_condition = @(cond, goal) average_trials_by_condition(subj_trials, cond, goal, max_samples);

    pupil_ab_comp(s_idx, :) = avg_condition('compared', 'A-B');
    pupil_ab_iso(s_idx, :) = avg_condition('isolated', 'A-B');
    pupil_ab_nov(s_idx, :) = avg_condition('novel', 'A-B');
    pupil_aa_comp(s_idx, :) = avg_condition('compared', 'A-A');
    pupil_aa_iso(s_idx, :) = avg_condition('isolated', 'A-A');
    pupil_aa_nov(s_idx, :) = avg_condition('novel', 'A-A');
end

fprintf('\nSubjects with valid data per condition:\n');
fprintf('  A-B compared: %d/%d\n', sum(~isnan(pupil_ab_comp(:,1))), n_subj);
fprintf('  A-B isolated: %d/%d\n', sum(~isnan(pupil_ab_iso(:,1))), n_subj);
fprintf('  A-B novel: %d/%d\n', sum(~isnan(pupil_ab_nov(:,1))), n_subj);
fprintf('  A-A compared: %d/%d\n', sum(~isnan(pupil_aa_comp(:,1))), n_subj);
fprintf('  A-A isolated: %d/%d\n', sum(~isnan(pupil_aa_iso(:,1))), n_subj);
fprintf('  A-A novel: %d/%d\n', sum(~isnan(pupil_aa_nov(:,1))), n_subj);

%% Visualization - Group average
% figure('color','w','position',[50 50 1200 500]);
% 
% % A-B trials
% subplot(1,2,1); hold on;
% plot_pupil_timeseries(time_vec, pupil_ab_comp, c_comp, 'compared');
% plot_pupil_timeseries(time_vec, pupil_ab_iso, c_iso, 'isolated');
% plot_pupil_timeseries(time_vec, pupil_ab_nov, c_nov, 'novel');
% xlabel('Time from stimulus onset (s)', 'FontSize', 16);
% ylabel('Pupil size change (a.u.)', 'FontSize', 16);
% title('A-B Lure Discrimination', 'FontSize', 16);
% legend('Location', 'best');
% grid on; box off;
% xlim([0 1.5]);
% % yline(0, 'k--', 'LineWidth', 1);
% 
% % A-A trials
% subplot(1,2,2); hold on;
% plot_pupil_timeseries(time_vec, pupil_aa_comp, c_comp, 'compared');
% plot_pupil_timeseries(time_vec, pupil_aa_iso, c_iso, 'isolated');
% plot_pupil_timeseries(time_vec, pupil_aa_nov, c_nov, 'novel');
% xlabel('Time from stimulus onset (s)', 'FontSize', 16);
% ylabel('Pupil size change (a.u.)', 'FontSize', 16);
% title('A-A Target Detection', 'FontSize', 16);
% legend('Location', 'best');
% grid on; box off;
% xlim([0 1.5]);
% % yline(0, 'k--', 'LineWidth', 1);
% 
% sgtitle('Pupil Responses over Time', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(res_dir, 'pupil_ts.pdf'), '-dpdf', '-vector');

% 
% %% Individual subjects - all 3 conditions
% n_cols = 4;
% n_rows = ceil(n_subj / n_cols);
% 
% % A-B trials
% figure('color','w','position',[50 50 1600 1200]);
% for s_idx = 1:n_subj
%     subplot(n_rows, n_cols, s_idx); hold on;
% 
%     if ~isnan(pupil_ab_comp(s_idx, 1))
%         plot(time_vec, pupil_ab_comp(s_idx, :), 'Color', c_comp, 'LineWidth', 2);
%     end
%     if ~isnan(pupil_ab_iso(s_idx, 1))
%         plot(time_vec, pupil_ab_iso(s_idx, :), 'Color', c_iso, 'LineWidth', 2);
%     end
%     if ~isnan(pupil_ab_nov(s_idx, 1))
%         plot(time_vec, pupil_ab_nov(s_idx, :), 'Color', c_nov, 'LineWidth', 2);
%     end
% 
%     title(sprintf('Sub %d', unique_subjs(s_idx)), 'FontSize', 10);
%     xlabel('Time (s)', 'FontSize', 8);
%     ylabel('Pupil', 'FontSize', 8);
%     xlim([0 1.5]);
%     yline(0, 'k--', 'LineWidth', 0.5);
%     grid on; box off;
% 
%     if s_idx == 1
%         legend('compared', 'isolated', 'novel', 'Location', 'best', 'FontSize', 8);
%     end
% end
% sgtitle('A-B Trials: Individual Subjects', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(res_dir, 'Individual_AB_Pupil_Timeseries.pdf'), '-dpdf', '-vector');
% 
% % A-A trials
% figure('color','w','position',[50 50 1600 1200]);
% for s_idx = 1:n_subj
%     subplot(n_rows, n_cols, s_idx); hold on;
% 
%     if ~isnan(pupil_aa_comp(s_idx, 1))
%         plot(time_vec, pupil_aa_comp(s_idx, :), 'Color', c_comp, 'LineWidth', 2);
%     end
%     if ~isnan(pupil_aa_iso(s_idx, 1))
%         plot(time_vec, pupil_aa_iso(s_idx, :), 'Color', c_iso, 'LineWidth', 2);
%     end
%     if ~isnan(pupil_aa_nov(s_idx, 1))
%         plot(time_vec, pupil_aa_nov(s_idx, :), 'Color', c_nov, 'LineWidth', 2);
%     end
% 
%     title(sprintf('Sub %d', unique_subjs(s_idx)), 'FontSize', 10);
%     xlabel('Time (s)', 'FontSize', 8);
%     ylabel('Pupil', 'FontSize', 8);
%     xlim([0 1.5]);
%     yline(0, 'k--', 'LineWidth', 0.5);
%     grid on; box off;
% 
%     if s_idx == 1
%         legend('compared', 'isolated', 'novel', 'Location', 'best', 'FontSize', 8);
%     end
% end
% sgtitle('A-A Trials: Individual Subjects', 'FontSize', 16, 'FontWeight', 'bold');
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(res_dir, 'Individual_AA_Pupil_Timeseries.pdf'), '-dpdf', '-vector');
% 
% fprintf('\n=== Visualization Complete ===\n');


%% functions
function t = compute_t(x,y)
    [~,~,~,stats] = ttest(x,y); t = stats.tstat;
end
function s = sig_star(p)
    if p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s='n.s.'; end
end
function write_cluster_results(fid, clusters, p_values, time_vec)
    if isempty(clusters)
        fprintf(fid, '  No significant clusters found (no timepoints exceeded threshold)\n');
        return;
    end
    
    sig_found = false;
    for i = 1:length(clusters)
        cluster_idx = clusters{i};
        t_start = time_vec(cluster_idx(1));
        t_end = time_vec(cluster_idx(end));
        p = p_values(i);
        
        if p < 0.05
            fprintf(fid, '  Cluster %d: %.3f-%.3f s, p=%.3f %s\n', ...
                i, t_start, t_end, p, sig_star(p));
            sig_found = true;
        end
    end
    
    if ~sig_found
        fprintf(fid, '  No significant clusters (p<0.05)\n');
        fprintf(fid, '  (%d clusters tested but none survived correction)\n', length(clusters));
    end
end
function [clusters, p_values, obs_t] = cluster_permutation_test(data1, data2, time_vec, n_perm)
    n_subj = size(data1, 1); n_time = size(data1, 2);
    fprintf('  Computing observed t-statistics...\n');
    obs_t = zeros(1, n_time);
    for t = 1:n_time, [~, ~, ~, stats] = ttest(data1(:,t), data2(:,t)); obs_t(t) = stats.tstat; end
    t_thresh = tinv(0.975, n_subj-1); above_thresh = abs(obs_t) > t_thresh;
    fprintf('  Identifying observed clusters...\n');
    clusters = {}; cluster_stats = []; in_cluster = false; cluster_start = 0;
    for t = 1:n_time
        if above_thresh(t) && ~in_cluster, cluster_start = t; in_cluster = true;
        elseif ~above_thresh(t) && in_cluster, clusters{end+1} = cluster_start:(t-1); cluster_stats(end+1) = sum(abs(obs_t(cluster_start:(t-1)))); in_cluster = false; end
    end
    if in_cluster, clusters{end+1} = cluster_start:n_time; cluster_stats(end+1) = sum(abs(obs_t(cluster_start:n_time))); end
    fprintf('  Found %d observed clusters\n', length(clusters));
    if isempty(clusters), p_values = []; return; end
    fprintf('  Running %d permutations: ', n_perm); null_max_cluster = zeros(n_perm, 1);
    for perm = 1:n_perm
        if mod(perm, 100) == 0, fprintf('%d...', perm); end
        signs = (rand(n_subj, 1) > 0.5) * 2 - 1; perm_data1 = data1; perm_data2 = data2;
        for s = 1:n_subj, if signs(s) < 0, temp = perm_data1(s,:); perm_data1(s,:) = perm_data2(s,:); perm_data2(s,:) = temp; end, end
        perm_t = zeros(1, n_time); for t = 1:n_time, [~, ~, ~, stats] = ttest(perm_data1(:,t), perm_data2(:,t)); perm_t(t) = stats.tstat; end
        perm_above = abs(perm_t) > t_thresh; perm_cluster_stats = []; in_cluster = false; cluster_start = 0;
        for t = 1:n_time
            if perm_above(t) && ~in_cluster, cluster_start = t; in_cluster = true;
            elseif ~perm_above(t) && in_cluster, perm_cluster_stats(end+1) = sum(abs(perm_t(cluster_start:(t-1)))); in_cluster = false; end
        end
        if in_cluster, perm_cluster_stats(end+1) = sum(abs(perm_t(cluster_start:n_time))); end
        if ~isempty(perm_cluster_stats), null_max_cluster(perm) = max(perm_cluster_stats); end
    end
    fprintf(' Done!\n');
    p_values = zeros(1, length(clusters)); for c = 1:length(clusters), p_values(c) = mean(null_max_cluster >= cluster_stats(c)); end
end
function display_clusters(clusters, p_values, time_vec)
    if isempty(clusters), fprintf('    No significant clusters found\n'); return; end
    sig_found = false;
    for i = 1:length(clusters)
        if p_values(i) < 0.05
            fprintf('    Cluster %d: %.3f-%.3f s, p=%.3f %s\n', i, time_vec(clusters{i}(1)), time_vec(clusters{i}(end)), p_values(i), sig_star(p_values(i)));
            sig_found = true;
        end
    end
    if ~sig_found, fprintf('    No significant clusters (p<0.05)\n'); end
end
function plot_with_sig(time_vec, data1, data2, color1, color2, labels, clusters, p_values, obs_t)
    plot_pupil_timeseries(time_vec, data1, color1, labels{1});
    plot_pupil_timeseries(time_vec, data2, color2, labels{2});
    if ~isempty(clusters)
        yl = ylim;
        for i = 1:length(clusters)
            if p_values(i) < 0.05
                cluster_idx = clusters{i}; t_start = time_vec(cluster_idx(1)); t_end = time_vec(cluster_idx(end));
                mean_t = mean(obs_t(cluster_idx)); shade_color = color1; if mean_t < 0, shade_color = color2; end
                fill([t_start t_end t_end t_start], [yl(1) yl(1) yl(2) yl(2)], shade_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            end
        end
    end
    legend('Location', 'best');
end
function avg_pupil = average_trials_by_condition(subj_trials, cond, goal, max_samples)
    % Filter trials
    trials_filt = subj_trials(strcmp(subj_trials.condition, cond) & ...
                              strcmp(subj_trials.goal, goal), :);
    
    if height(trials_filt) == 0
        avg_pupil = nan(1, max_samples);
        return;
    end
    
    % Accumulate trials (pad to max length)
    trial_matrix = nan(height(trials_filt), max_samples);
    for tr = 1:height(trials_filt)
        pupil_data = trials_filt.pupil_preprocessed{tr};
        baseline_pupil = trials_filt.baseline_pupil_preprocessed{tr};
        
        % Baseline correct: subtract mean of the 200-sample baseline period
        baseline = mean(baseline_pupil, 'omitnan');
        pupil_data = pupil_data - baseline;
        
        trial_matrix(tr, 1:length(pupil_data)) = pupil_data';
    end
    
    % Average across trials
    avg_pupil = mean(trial_matrix, 1, 'omitnan');
end
function plot_pupil_timeseries(time_vec, data, color, label)
    mean_pupil = mean(data, 1, 'omitnan');
    n_valid = sum(~isnan(data(:,1)));
    
    if n_valid == 0
        warning('No valid data for condition: %s', label);
        return;
    end
    
    sem_pupil = std(data, 0, 1, 'omitnan') / sqrt(n_valid);
    
    % Plot mean line
    plot(time_vec, mean_pupil, 'Color', color, 'LineWidth', 2.5, 'DisplayName', label);
    
    % Plot shaded error region (SEM)
    valid_idx = ~isnan(mean_pupil);
    if sum(valid_idx) > 1
        fill([time_vec(valid_idx), fliplr(time_vec(valid_idx))], ...
             [mean_pupil(valid_idx) + sem_pupil(valid_idx), fliplr(mean_pupil(valid_idx) - sem_pupil(valid_idx))], ...
             color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end