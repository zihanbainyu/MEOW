clear; clc; close all;
fprintf('=== Parse All Trials from All Subjects ===\n');

base_dir = '..';
subj_ids = [501, 601, 602, 603, 604, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 617];
res_dir = fullfile(base_dir, 'results');

% C = struct('W_CM',60, 'W_PX',1920, 'H_PX',1080, 'D_CM',100, 'SR',1000, 'BASELINE_DUR',0.5);
% 
% % Master dataset - all trials from all subjects
% all_trials_data = [];
% trial_counter = 0;
% 
% for s_idx = 1:length(subj_ids)
%     subj_id = subj_ids(s_idx);
%     fprintf('\n--- Subject %d (%d/%d) ---\n', subj_id, s_idx, length(subj_ids));
% 
%     r_dir = fullfile(base_dir, 'data', sprintf('sub%03d', subj_id));
% 
%     % Load behavioral data
%     concat_file = fullfile(r_dir, sprintf('sub%03d_concat.mat', subj_id));
%     if ~isfile(concat_file)
%         fprintf('  Behavioral data not found, skipping\n');
%         continue;
%     end
%     load(concat_file, 'final_data_output');
% 
%     for run = 1:4
%         fprintf('  Run %d: ', run);
% 
%         % Load/parse ASC
%         asc_f = fullfile(r_dir, sprintf('%03d_2_%01d.asc', subj_id, run));
%         mat_f = fullfile(r_dir, sprintf('%03d_2_%01d_raw_v2.mat', subj_id, run));
% 
%         if ~isfile(asc_f)
%             fprintf('ASC not found\n');
%             continue;
%         end
% 
%         if isfile(mat_f)
%             load(mat_f, 'raw');
%         else
%             raw = parse_asc(asc_f);
%             save(mat_f, 'raw');
%         end
% 
%         % Get behavioral data for this run
%         behav = final_data_output.results_2_back_all(final_data_output.results_2_back_all.block==run, :);
%         n_trials = height(behav);
% 
%         % Convert to DVA
%         [x_dva, y_dva] = px2dva(raw.xp, raw.yp, C);
% 
%         % Extract each trial
%         ev_msg = raw.ev.msg;
%         ev_ts = raw.ev.ts;
% 
%         valid_count = 0;
%         for tr = 1:n_trials
%             % Find trial markers
%             tid = find(contains(ev_msg, sprintf('TRIALID %d', tr)), 1);
%             if isempty(tid), continue; end
% 
%             onset_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'STIM_ONSET'), 1);
%             if isempty(onset_idx), continue; end
%             t_onset = ev_ts(tid + onset_idx - 1);
% 
%             result_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'TRIAL_RESULT'), 1);
%             if isempty(result_idx), continue; end
%             t_trial_end = ev_ts(tid + result_idx - 1);
% 
%             % Baseline period
%             t_baseline_start = t_onset - (C.BASELINE_DUR * 1000);
% 
%             % Find sample indices
%             idx_base_s = find(raw.ts >= t_baseline_start, 1);
%             idx_base_e = find(raw.ts < t_onset, 1, 'last');
%             idx_trial_s = find(raw.ts >= t_onset, 1);
%             idx_trial_e = find(raw.ts <= t_trial_end, 1, 'last');
% 
%             if isempty(idx_base_s) || isempty(idx_base_e) || idx_base_e <= idx_base_s
%                 continue;
%             end
%             if isempty(idx_trial_s) || isempty(idx_trial_e) || idx_trial_e <= idx_trial_s
%                 continue;
%             end
% 
%             trial_counter = trial_counter + 1;
% 
%             % Store trial data
%             trial_data = struct();
%             trial_data.global_trial_id = trial_counter;
%             trial_data.subj_id = subj_id;
%             trial_data.run = run;
%             trial_data.trial_id_in_run = tr;
% 
%             % Eye tracking data
%             trial_data.x_dva = x_dva(idx_trial_s:idx_trial_e);
%             trial_data.y_dva = y_dva(idx_trial_s:idx_trial_e);
%             trial_data.pupil = raw.pp(idx_trial_s:idx_trial_e);
%             trial_data.timestamps = raw.ts(idx_trial_s:idx_trial_e);
% 
%             trial_data.baseline_x_dva = x_dva(idx_base_s:idx_base_e);
%             trial_data.baseline_y_dva = y_dva(idx_base_s:idx_base_e);
%             trial_data.baseline_pupil = raw.pp(idx_base_s:idx_base_e);
%             trial_data.baseline_timestamps = raw.ts(idx_base_s:idx_base_e);
% 
%             trial_data.sample_rate = C.SR;
%             trial_data.trial_length = idx_trial_e - idx_trial_s + 1;
%             trial_data.trial_duration_s = trial_data.trial_length / C.SR;
% 
%             % Behavioral data
%             trial_data.condition = char(behav.condition(tr));
%             trial_data.goal = char(behav.goal(tr));
%             trial_data.identity = char(behav.identity(tr));
%             trial_data.stim_id = char(behav.stim_id(tr));
%             trial_data.corr_resp = char(behav.corr_resp(tr));
%             trial_data.resp_key = char(behav.resp_key(tr));
%             if strcmp(trial_data.resp_key, 'NA')
%                 trial_data.resp_key = 'none';
%             end
%             trial_data.rt = behav.rt(tr);
%             trial_data.correct = strcmp(trial_data.corr_resp, trial_data.resp_key);
% 
%             all_trials_data = [all_trials_data; trial_data];
%             valid_count = valid_count + 1;
%         end
% 
%         fprintf('%d/%d trials\n', valid_count, n_trials);
%     end
% end
% 
% fprintf('\n=== Parsing Complete ===\n');
% fprintf('Total trials extracted: %d\n', length(all_trials_data));
% 
% % Convert to table
% all_trials_table = struct2table(all_trials_data);
% 
% % Summary statistics
% fprintf('\nTrial counts by condition:\n');
% for cond = ["compared", "isolated", "novel"]
%     for goal = ["A-A", "A-B"]
%         n_total = sum(strcmp(all_trials_table.condition, cond) & strcmp(all_trials_table.goal, goal));
%         n_correct = sum(strcmp(all_trials_table.condition, cond) & strcmp(all_trials_table.goal, goal) & all_trials_table.correct);
%         fprintf('  %s %s: %d total (%d correct, %.1f%%)\n', goal, cond, n_total, n_correct, 100*n_correct/n_total);
%     end
% end
% 
% fprintf('\nTrials per subject:\n');
% for s_idx = 1:length(subj_ids)
%     n = sum(all_trials_table.subj_id == subj_ids(s_idx));
%     fprintf('  Subject %d: %d trials\n', subj_ids(s_idx), n);
% end
% 
% % Save complete dataset
% save(fullfile(res_dir, 'all_trials_raw.mat'), 'all_trials_data', 'all_trials_table', '-v7.3');
% fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_raw.mat'));
% 
% fprintf('\n=== Ready for PCDM preprocessing ===\n');
% fprintf('Next: Filter trials and preprocess with PCDM\n');
% 
% %% Preprocess with PCDM
% fprintf('\n=== Preprocessing with PCDM ===\n');
% 
% % Add PCDM to path
% pcdm_dir = fullfile(base_dir, 'scripts/PCDM-main');
% addpath(genpath(pcdm_dir));
% 
% % Prepare data structure for PCDM preprocessing
% % Group trials by subject for preprocessing
% unique_subjs = unique(all_trials_table.subj_id);
% n_subj = length(unique_subjs);
% 
% all_preprocessed = [];
% 
% for s_idx = 1:n_subj
%     subj_id = unique_subjs(s_idx);
%     fprintf('Preprocessing subject %d (%d/%d)...\n', subj_id, s_idx, n_subj);
% 
%     % Get all trials for this subject
%     subj_trials = all_trials_table(all_trials_table.subj_id == subj_id, :);
% 
%     % Prepare input structure for PCDM
%     in = struct();
%     in.xPos = {[]};
%     in.yPos = {[]};
%     in.pupilArea = {[]};
%     in.sampleRate = {C.SR};
%     in.startInds = {[]};
%     in.trialTypes = {[]};
% 
%     % Concatenate all trials into continuous timeseries
%     concat_x = [];
%     concat_y = [];
%     concat_pupil = [];
%     trial_starts = [];
%     trial_types = [];
% 
%     current_sample = 1;
%     for tr = 1:height(subj_trials)
%         trial = subj_trials(tr, :);
% 
%         % Concatenate data
%         concat_x = [concat_x; trial.x_dva{1}];
%         concat_y = [concat_y; trial.y_dva{1}];
%         concat_pupil = [concat_pupil; trial.pupil{1}];
% 
%         % Record trial start/end
%         trial_end = current_sample + length(trial.pupil{1}) - 1;
%         trial_starts = [trial_starts; current_sample, trial_end];
% 
%         % Assign trial type (we'll use 1 for now, will filter later)
%         trial_types = [trial_types; 1];
% 
%         current_sample = trial_end + 1;
%     end
% 
%     in.xPos{1} = concat_x;
%     in.yPos{1} = concat_y;
%     in.pupilArea{1} = concat_pupil;
%     in.startInds{1} = trial_starts;
%     in.trialTypes{1} = trial_types;
% 
%     % Run PCDM preprocessing (dataAnalysis only, no model fitting)
%     in.putativeIRFdur = 4;
%     in.predictionWindow = 1.5;
% 
%     try
%         % Call dataAnalysis directly (not fitModel)
%         op.interpolateBlinks = 1;
%         op.downsampleRate = 10;
%         op.fitTimeseries = 0;
% 
%         d = dataAnalysis(in, op);
% 
%         % Extract preprocessed pupil for each trial
%         pupil_preprocessed = d.pupilTS{1};  % Preprocessed continuous timeseries
% 
%         for tr = 1:height(subj_trials)
%             trial_start = trial_starts(tr, 1);
%             trial_end = trial_starts(tr, 2);
% 
%             % Add preprocessed pupil to trial
%             preprocessed_trial = subj_trials(tr, :);
%             preprocessed_trial.pupil_preprocessed = {pupil_preprocessed(trial_start:trial_end)};
% 
%             all_preprocessed = [all_preprocessed; preprocessed_trial];
%         end
% 
%         fprintf('  Success: %d trials preprocessed\n', height(subj_trials));
% 
%     catch ME
%         fprintf('  ERROR: %s\n', ME.message);
%         % If preprocessing fails, keep original data without preprocessed pupil
%         for tr = 1:height(subj_trials)
%             preprocessed_trial = subj_trials(tr, :);
%             preprocessed_trial.pupil_preprocessed = {nan(size(subj_trials.pupil{tr}))};
%             all_preprocessed = [all_preprocessed; preprocessed_trial];
%         end
%     end
% end
% 
% fprintf('\n=== Preprocessing Complete ===\n');
% fprintf('Total trials with preprocessed pupil: %d\n', height(all_preprocessed));
% 
% % Check how many trials have valid preprocessed data
% n_valid = sum(cellfun(@(x) ~all(isnan(x)), all_preprocessed.pupil_preprocessed));
% fprintf('Valid preprocessed trials: %d (%.1f%%)\n', n_valid, 100*n_valid/height(all_preprocessed));
% 
% % Save complete dataset with preprocessed pupil
% save(fullfile(res_dir, 'all_trials_preprocessed.mat'), 'all_preprocessed', '-v7.3');
% fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_preprocessed.mat'));
% 
% fprintf('\nDataset ready for analysis!\n');
% fprintf('Columns include:\n');
% fprintf('  - global_trial_id, subj_id, run, trial_id_in_run\n');
% fprintf('  - condition, goal, identity, correct, rt\n');
% fprintf('  - pupil (raw), pupil_preprocessed (blink-interpolated, filtered)\n');
% fprintf('  - x_dva, y_dva, timestamps\n');



%% Group-level visualization of preprocessed pupil data
fprintf('\n=== Group-Level Pupil Analysis ===\n');
load(fullfile(res_dir, 'all_trials_preprocessed.mat'));
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
figure('color','w','position',[50 50 1200 500]);

% A-B trials
subplot(1,2,1); hold on;
plot_pupil_timeseries(time_vec, pupil_ab_comp, c_comp, 'compared');
plot_pupil_timeseries(time_vec, pupil_ab_iso, c_iso, 'isolated');
plot_pupil_timeseries(time_vec, pupil_ab_nov, c_nov, 'novel');
xlabel('Time from stimulus onset (s)', 'FontSize', 12);
ylabel('Pupil size (preprocessed)', 'FontSize', 12);
title('A-B Trials (Lure Discrimination)', 'FontSize', 14);
legend('Location', 'best');
grid on; box off;
xlim([0 1.5]);

% A-A trials
subplot(1,2,2); hold on;
plot_pupil_timeseries(time_vec, pupil_aa_comp, c_comp, 'compared');
plot_pupil_timeseries(time_vec, pupil_aa_iso, c_iso, 'isolated');
plot_pupil_timeseries(time_vec, pupil_aa_nov, c_nov, 'novel');
xlabel('Time from stimulus onset (s)', 'FontSize', 12);
ylabel('Pupil size (preprocessed)', 'FontSize', 12);
title('A-A Trials (Recognition)', 'FontSize', 14);
legend('Location', 'best');
grid on; box off;
xlim([0 1.5]);

sgtitle('Group-Level Pupil Timeseries by Condition', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(res_dir, 'Group_Pupil_Timeseries.pdf'), '-dpdf', '-vector');

%% Individual subjects - all 3 conditions
n_cols = 4;
n_rows = ceil(n_subj / n_cols);

% A-B trials
figure('color','w','position',[50 50 1600 1200]);
for s_idx = 1:n_subj
    subplot(n_rows, n_cols, s_idx); hold on;
    
    if ~isnan(pupil_ab_comp(s_idx, 1))
        plot(time_vec, pupil_ab_comp(s_idx, :), 'Color', c_comp, 'LineWidth', 2);
    end
    if ~isnan(pupil_ab_iso(s_idx, 1))
        plot(time_vec, pupil_ab_iso(s_idx, :), 'Color', c_iso, 'LineWidth', 2);
    end
    if ~isnan(pupil_ab_nov(s_idx, 1))
        plot(time_vec, pupil_ab_nov(s_idx, :), 'Color', c_nov, 'LineWidth', 2);
    end
    
    title(sprintf('Sub %d', unique_subjs(s_idx)), 'FontSize', 10);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Pupil', 'FontSize', 8);
    xlim([0 1.5]);
    grid on; box off;
    
    if s_idx == 1
        legend('compared', 'isolated', 'novel', 'Location', 'best', 'FontSize', 8);
    end
end
sgtitle('A-B Trials: Individual Subjects', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(res_dir, 'Individual_AB_Pupil_Timeseries.pdf'), '-dpdf', '-vector');

% A-A trials
figure('color','w','position',[50 50 1600 1200]);
for s_idx = 1:n_subj
    subplot(n_rows, n_cols, s_idx); hold on;
    
    if ~isnan(pupil_aa_comp(s_idx, 1))
        plot(time_vec, pupil_aa_comp(s_idx, :), 'Color', c_comp, 'LineWidth', 2);
    end
    if ~isnan(pupil_aa_iso(s_idx, 1))
        plot(time_vec, pupil_aa_iso(s_idx, :), 'Color', c_iso, 'LineWidth', 2);
    end
    if ~isnan(pupil_aa_nov(s_idx, 1))
        plot(time_vec, pupil_aa_nov(s_idx, :), 'Color', c_nov, 'LineWidth', 2);
    end
    
    title(sprintf('Sub %d', unique_subjs(s_idx)), 'FontSize', 10);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Pupil', 'FontSize', 8);
    xlim([0 1.5]);
    grid on; box off;
    
    if s_idx == 1
        legend('compared', 'isolated', 'novel', 'Location', 'best', 'FontSize', 8);
    end
end
sgtitle('A-A Trials: Individual Subjects', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, fullfile(res_dir, 'Individual_AA_Pupil_Timeseries.pdf'), '-dpdf', '-vector');

fprintf('\n=== Visualization Complete ===\n');

%% Helper functions
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

%% Helper functions
function r = parse_asc(f)
    fid = fopen(f); 
    N=5e6; 
    ts=zeros(N,1); 
    xp=nan(N,1); 
    yp=nan(N,1); 
    pp=nan(N,1);
    ev_ts=zeros(5000,1); 
    ev_msg=cell(5000,1); 
    c=0; 
    ec=0;
    
    while ~feof(fid)
        l = fgetl(fid); 
        if isempty(l), continue; end
        
        if l(1)>='0' && l(1)<='9'
            c=c+1; 
            d = sscanf(l, '%d %f %f %d');
            if numel(d)==4
                ts(c)=d(1); 
                xp(c)=d(2); 
                yp(c)=d(3); 
                pp(c)=d(4);
            else
                ts(c)=sscanf(l,'%d',1);
            end
        elseif startsWith(l, 'MSG')
            ec=ec+1; 
            dat = textscan(l, 'MSG %d %s', 'Delimiter', '');
            ev_ts(ec)=dat{1}; 
            ev_msg{ec}=l;
        end
    end
    fclose(fid);
    
    r.ts=ts(1:c); 
    r.xp=xp(1:c); 
    r.yp=yp(1:c); 
    r.pp=pp(1:c);
    r.ev.ts=ev_ts(1:ec); 
    r.ev.msg=ev_msg(1:ec);
end

function [xd, yd] = px2dva(x, y, C)
    px2cm = C.W_PX / C.W_CM; 
    cx = C.W_PX/2; 
    cy = C.H_PX/2;
    xc = nan(size(x)); 
    yc = nan(size(y));
    ok = (x>0 & x<C.W_PX & y>0 & y<C.H_PX);
    xc(ok) = (x(ok)-cx)/px2cm; 
    yc(ok) = (y(ok)-cy)/px2cm;
    xd = atan(xc/C.D_CM) * (180/pi); 
    yd = atan(yc/C.D_CM) * (180/pi);
end