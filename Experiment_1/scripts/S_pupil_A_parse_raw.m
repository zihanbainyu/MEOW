clear; clc; close all;
fprintf('=== Parse All Trials from All Subjects ===\n');

base_dir = '..';
subj_ids = [618, 619, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631];
res_dir = fullfile(base_dir, 'results');

C = struct('W_CM',60, 'W_PX',1920, 'H_PX',1080, 'D_CM',100, 'SR',1000, 'BASELINE_DUR',0.5);

% Master dataset - all trials from all subjects
all_trials_data = [];
trial_counter = 0;

for s_idx = 1:length(subj_ids)
    subj_id = subj_ids(s_idx);
    fprintf('\n--- Subject %d (%d/%d) ---\n', subj_id, s_idx, length(subj_ids));

    r_dir = fullfile(base_dir, 'data', sprintf('sub%03d', subj_id));

    % Load behavioral data
    concat_file = fullfile(r_dir, sprintf('sub%03d_concat.mat', subj_id));
    if ~isfile(concat_file)
        fprintf('  Behavioral data not found, skipping\n');
        continue;
    end
    load(concat_file, 'final_data_output');

    for run = 1:4
        fprintf('  Run %d: ', run);

        % Load/parse ASC
        asc_f = fullfile(r_dir, sprintf('%03d_2_%01d.asc', subj_id, run));
      
        raw = parse_asc(asc_f);

        % Get behavioral data for this run
        behav = final_data_output.results_2_back_all(final_data_output.results_2_back_all.block==run, :);
        n_trials = height(behav);

        % Convert to DVA
        [x_dva, y_dva] = px2dva(raw.xp, raw.yp, C);

        % Extract each trial
        ev_msg = raw.ev.msg;
        ev_ts = raw.ev.ts;

        valid_count = 0;
        for tr = 1:n_trials
            % Find trial markers
            tid = find(contains(ev_msg, sprintf('TRIALID %d', tr)), 1);
            if isempty(tid), continue; end

            onset_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'STIM_ONSET'), 1);
            if isempty(onset_idx), continue; end
            t_onset = ev_ts(tid + onset_idx - 1);

            result_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'TRIAL_RESULT'), 1);
            if isempty(result_idx), continue; end
            t_trial_end = ev_ts(tid + result_idx - 1);

            % Baseline period
            t_baseline_start = t_onset - (C.BASELINE_DUR * 1000);

            % Find sample indices
            idx_base_s = find(raw.ts >= t_baseline_start, 1);
            idx_base_e = find(raw.ts < t_onset, 1, 'last');
            idx_trial_s = find(raw.ts >= t_onset, 1);
            idx_trial_e = find(raw.ts <= t_trial_end, 1, 'last');

            if isempty(idx_base_s) || isempty(idx_base_e) || idx_base_e <= idx_base_s
                continue;
            end
            if isempty(idx_trial_s) || isempty(idx_trial_e) || idx_trial_e <= idx_trial_s
                continue;
            end

            trial_counter = trial_counter + 1;

            % Store trial data
            trial_data = struct();
            trial_data.global_trial_id = trial_counter;
            trial_data.subj_id = subj_id;
            trial_data.run = run;
            trial_data.trial_id_in_run = tr;

            % Eye tracking data
            trial_data.x_dva = x_dva(idx_trial_s:idx_trial_e);
            trial_data.y_dva = y_dva(idx_trial_s:idx_trial_e);
            trial_data.pupil = raw.pp(idx_trial_s:idx_trial_e);
            trial_data.timestamps = raw.ts(idx_trial_s:idx_trial_e);

            trial_data.baseline_x_dva = x_dva(idx_base_s:idx_base_e);
            trial_data.baseline_y_dva = y_dva(idx_base_s:idx_base_e);
            trial_data.baseline_pupil = raw.pp(idx_base_s:idx_base_e);
            trial_data.baseline_timestamps = raw.ts(idx_base_s:idx_base_e);

            trial_data.sample_rate = C.SR;
            trial_data.trial_length = idx_trial_e - idx_trial_s + 1;
            trial_data.trial_duration_s = trial_data.trial_length / C.SR;

            % Behavioral data
            trial_data.condition = char(behav.condition(tr));
            trial_data.goal = char(behav.goal(tr));
            trial_data.identity = char(behav.identity(tr));
            trial_data.stim_id = char(behav.stim_id(tr));
            trial_data.corr_resp = char(behav.corr_resp(tr));
            trial_data.resp_key = char(behav.resp_key(tr));
            if strcmp(trial_data.resp_key, 'NA')
                trial_data.resp_key = 'none';
            end
            trial_data.rt = behav.rt(tr);
            trial_data.correct = strcmp(trial_data.corr_resp, trial_data.resp_key);

            all_trials_data = [all_trials_data; trial_data];
            valid_count = valid_count + 1;
        end

        fprintf('%d/%d trials\n', valid_count, n_trials);
    end
end

fprintf('\n=== Parsing Complete ===\n');
fprintf('Total trials extracted: %d\n', length(all_trials_data));

% Convert to table
all_trials_table = struct2table(all_trials_data);

% Summary statistics
fprintf('\nTrial counts by condition:\n');
for cond = ["compared", "isolated", "novel"]
    for goal = ["A-A", "A-B"]
        n_total = sum(strcmp(all_trials_table.condition, cond) & strcmp(all_trials_table.goal, goal));
        n_correct = sum(strcmp(all_trials_table.condition, cond) & strcmp(all_trials_table.goal, goal) & all_trials_table.correct);
        fprintf('  %s %s: %d total (%d correct, %.1f%%)\n', goal, cond, n_total, n_correct, 100*n_correct/n_total);
    end
end

fprintf('\nTrials per subject:\n');
for s_idx = 1:length(subj_ids)
    n = sum(all_trials_table.subj_id == subj_ids(s_idx));
    fprintf('  Subject %d: %d trials\n', subj_ids(s_idx), n);
end
% Save complete dataset
save(fullfile(res_dir, 'all_trials_raw_half.mat'), 'all_trials_data', 'all_trials_table', '-v7.3');
fprintf('\nSaved to: %s\n', fullfile(res_dir, 'all_trials_raw_half.mat'));

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