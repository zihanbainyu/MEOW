clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj_ids = [501, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617];
gaze_subj_ids = [501, 601, 602, 603, 604, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 617];

base_dir = '..'; 
out_dir = fullfile(base_dir, 'data', 'eye_movement_data');
res_dir = fullfile(base_dir, 'results');
min_rt = 0.150;

% Colors
c_comp = [180 174 211]/255; 
c_iso = [176 230 255]/255; 
c_nov = [183 210 205]/255; 
c_sim = [255 191 205]/255; 
c_same = [97 125 184]/255; 
c_new = [219 219 219]/255;

% Gaze reinstatement parameters
spatial_params = struct('xres', 1920, 'yres', 1080, 'roi_x', [760, 1160], ...
                        'roi_y', [340, 740], 'grid_x', 20, 'grid_y', 20, 'sigma', 2);
n_baseline = 20;
n_perm = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: MATH - BEHAVIORAL METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('==================================================================\n');
fprintf('COMPUTING BEHAVIORAL METRICS\n');
fprintf('==================================================================\n');

all_subjs = struct();

for s = 1:length(subj_ids)
    curr_id = subj_ids(s);
    fprintf('Processing subject %d...\n', curr_id);
    
    subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', curr_id));
    load(fullfile(subj_dir, sprintf('sub%03d_concat.mat', curr_id)), 'final_data_output');
    
    r1 = final_data_output.results_1_back_all;
    r2 = final_data_output.results_2_back_all;
    stats = struct();
    
    % --- 1-back metrics ---
    r1.resp_key = cellstr(r1.resp_key); 
    r1.resp_key(strcmp(r1.resp_key, 'NA')) = {'none'};
    r1.correct = strcmp(cellstr(r1.corr_resp), r1.resp_key);
    
    idx_sim = r1.condition == "compared" & r1.identity == "B";
    idx_sam = r1.condition == "repeat" & strcmp(r1.corr_resp, 'j');
    idx_cr = (r1.condition == "compared" & r1.identity == "A") | ...
             r1.condition == "isolated" | ...
             (r1.condition == "repeat" & strcmp(r1.corr_resp, 'none'));
    
    stats.one.acc_sim = mean(r1.correct(idx_sim));
    stats.one.acc_same = mean(r1.correct(idx_sam));
    stats.one.acc_new = mean(r1.correct(idx_cr));
    
    v_rt = r1.rt > min_rt;
    stats.one.rt_sim = median(r1.rt(idx_sim & r1.correct & v_rt), 'omitnan');
    stats.one.rt_same = median(r1.rt(idx_sam & r1.correct & v_rt), 'omitnan');
    
    n_sim = sum(idx_sim); n_sam = sum(idx_sam); n_cr = sum(idx_cr);
    stats.one.err_sim_as_same = sum(strcmp(r1.resp_key(idx_sim), 'j')) / n_sim;
    stats.one.err_sim_as_new = sum(strcmp(r1.resp_key(idx_sim), 'none')) / n_sim;
    stats.one.err_same_as_sim = sum(strcmp(r1.resp_key(idx_sam), 'k')) / n_sam;
    stats.one.err_same_as_new = sum(strcmp(r1.resp_key(idx_sam), 'none')) / n_sam;
    stats.one.err_new_as_same = sum(strcmp(r1.resp_key(idx_cr), 'j')) / n_cr;
    stats.one.err_new_as_sim = sum(strcmp(r1.resp_key(idx_cr), 'k')) / n_cr;
    
    % --- 2-back metrics ---
    r2.resp_key = cellstr(r2.resp_key); 
    r2.resp_key(strcmp(r2.resp_key, 'NA')) = {'none'};
    r2.correct = strcmp(cellstr(r2.corr_resp), r2.resp_key);
    
    pan = zeros(height(r2), 1); 
    for i = 1:height(r2)-2
        if strcmp(r2.goal(i), 'A-N')
            pan(i+2) = 1; 
        end
    end 
    
    real = ~contains(r2.goal, "JUNK");
    v_rt = r2.rt > min_rt;
    aa_idx = real & strcmp(r2.goal, 'A-A'); 
    ab_idx = real & strcmp(r2.goal, 'A-B'); 
    an_idx = (pan == 1);
    comp_idx = real & strcmp(r2.condition, 'compared');
    iso_idx = real & strcmp(r2.condition, 'isolated');
    nov_idx = real & strcmp(r2.condition, 'novel');
    j_idx = real & strcmp(r2.corr_resp, 'j'); 
    k_idx = real & strcmp(r2.corr_resp, 'k');
    n_idx = real & strcmp(r2.corr_resp, 'none');
    
    calc_d = @(h, f, nh, nf) norminv(max(1/(2*nh), min(1-1/(2*nh), h))) - ...
                             norminv(max(1/(2*nf), min(1-1/(2*nf), f)));
    
    stats.two.acc_AA_comp = mean(r2.correct(aa_idx & comp_idx & j_idx));
    stats.two.acc_AA_iso = mean(r2.correct(aa_idx & iso_idx & j_idx));
    stats.two.acc_AA_nov = mean(r2.correct(aa_idx & nov_idx & j_idx));
    stats.two.acc_AB_comp = mean(r2.correct(ab_idx & comp_idx & k_idx));
    stats.two.acc_AB_iso = mean(r2.correct(ab_idx & iso_idx & k_idx));
    stats.two.acc_AB_nov = mean(r2.correct(ab_idx & nov_idx & k_idx));
    stats.two.acc_AN_comp = mean(r2.correct(an_idx & comp_idx & n_idx));
    stats.two.acc_AN_iso = mean(r2.correct(an_idx & iso_idx & n_idx));
    stats.two.acc_AN_nov = mean(r2.correct(an_idx & nov_idx & n_idx));
    
    stats.two.rt_AA_comp = median(r2.rt(aa_idx & comp_idx & j_idx & r2.correct == 1 & v_rt), 'omitnan');
    stats.two.rt_AA_iso = median(r2.rt(aa_idx & iso_idx & j_idx & r2.correct == 1 & v_rt), 'omitnan');
    stats.two.rt_AA_nov = median(r2.rt(aa_idx & nov_idx & j_idx & r2.correct == 1 & v_rt), 'omitnan');
    stats.two.rt_AB_comp = median(r2.rt(ab_idx & comp_idx & k_idx & r2.correct == 1 & v_rt), 'omitnan');
    stats.two.rt_AB_iso = median(r2.rt(ab_idx & iso_idx & k_idx & r2.correct == 1 & v_rt), 'omitnan');
    stats.two.rt_AB_nov = median(r2.rt(ab_idx & nov_idx & k_idx & r2.correct == 1 & v_rt), 'omitnan');
    
    n_AA_comp = sum(aa_idx & comp_idx & j_idx);
    n_AA_iso = sum(aa_idx & iso_idx & j_idx);
    n_AA_nov = sum(aa_idx & nov_idx & j_idx);
    n_AB_comp = sum(ab_idx & comp_idx & k_idx);
    n_AB_iso = sum(ab_idx & iso_idx & k_idx);
    n_AB_nov = sum(ab_idx & nov_idx & k_idx);
    n_AN_comp = sum(an_idx & comp_idx & n_idx);
    n_AN_iso = sum(an_idx & iso_idx & n_idx);
    n_AN_nov = sum(an_idx & nov_idx & n_idx);
    
    stats.two.err_AA_comp_as_k = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'k')) / n_AA_comp;
    stats.two.err_AA_comp_as_n = sum(strcmp(r2.resp_key(aa_idx & comp_idx & j_idx), 'none')) / n_AA_comp;
    stats.two.err_AB_comp_as_j = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'j')) / n_AB_comp;
    stats.two.err_AB_comp_as_n = sum(strcmp(r2.resp_key(ab_idx & comp_idx & k_idx), 'none')) / n_AB_comp;
    stats.two.err_AN_comp_as_j = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'j')) / n_AN_comp;
    stats.two.err_AN_comp_as_k = sum(strcmp(r2.resp_key(an_idx & comp_idx & n_idx), 'k')) / n_AN_comp;
    
    stats.two.err_AA_iso_as_k = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'k')) / n_AA_iso;
    stats.two.err_AA_iso_as_n = sum(strcmp(r2.resp_key(aa_idx & iso_idx & j_idx), 'none')) / n_AA_iso;
    stats.two.err_AB_iso_as_j = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'j')) / n_AB_iso;
    stats.two.err_AB_iso_as_n = sum(strcmp(r2.resp_key(ab_idx & iso_idx & k_idx), 'none')) / n_AB_iso;
    stats.two.err_AN_iso_as_j = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'j')) / n_AN_iso;
    stats.two.err_AN_iso_as_k = sum(strcmp(r2.resp_key(an_idx & iso_idx & n_idx), 'k')) / n_AN_iso;
    
    stats.two.err_AA_nov_as_k = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'k')) / n_AA_nov;
    stats.two.err_AA_nov_as_n = sum(strcmp(r2.resp_key(aa_idx & nov_idx & j_idx), 'none')) / n_AA_nov;
    stats.two.err_AB_nov_as_j = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'j')) / n_AB_nov;
    stats.two.err_AB_nov_as_n = sum(strcmp(r2.resp_key(ab_idx & nov_idx & k_idx), 'none')) / n_AB_nov;
    stats.two.err_AN_nov_as_j = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'j')) / n_AN_nov;
    stats.two.err_AN_nov_as_k = sum(strcmp(r2.resp_key(an_idx & nov_idx & n_idx), 'k')) / n_AN_nov;
    
    stats.two.ldi_comp = stats.two.acc_AB_comp - stats.two.err_AN_comp_as_k;
    stats.two.ldi_iso = stats.two.acc_AB_iso - stats.two.err_AN_iso_as_k;
    stats.two.ldi_nov = stats.two.acc_AB_nov - stats.two.err_AN_nov_as_k;
    
    stats.two.dprime_comp = calc_d(stats.two.acc_AA_comp, stats.two.err_AN_comp_as_j, n_AA_comp, n_AN_comp);
    stats.two.dprime_iso = calc_d(stats.two.acc_AA_iso, stats.two.err_AN_iso_as_j, n_AA_iso, n_AN_iso);
    stats.two.dprime_nov = calc_d(stats.two.acc_AA_nov, stats.two.err_AN_nov_as_j, n_AA_nov, n_AN_nov);
    
    % --- Recognition metrics ---
    if isfield(final_data_output, 'results_recognition')
        rec = final_data_output.results_recognition;
        rec.correct = strcmp(cellstr(rec.corr_resp), cellstr(rec.resp_key));
        old = rec(rec.trial_type == "old", :); 
        new = rec(rec.trial_type ~= "old", :);
        n_new = height(new); 
        n_fa = sum(strcmp(cellstr(new.resp_key), 'j') & new.rt > min_rt);
        far = max(1/(2*n_new), min(1-1/(2*n_new), n_fa/n_new));
        
        tc = old(old.condition == "compared", :); 
        nc = height(tc);
        hc = sum(tc.correct & tc.rt > min_rt) / nc;
        stats.rec.d_comp = calc_d(hc, far, nc, n_new);
        
        ti = old(old.condition == "isolated", :); 
        ni = height(ti);
        hi = sum(ti.correct & ti.rt > min_rt) / ni;
        stats.rec.d_iso = calc_d(hi, far, ni, n_new);
    else
        fprintf('  Recognition data missing for subject %d\n', curr_id);
        stats.rec.d_comp = NaN;
        stats.rec.d_iso = NaN;
    end
    
    all_subjs(s).id = curr_id;
    all_subjs(s).stats = stats;
end

fprintf('\nBehavioral metrics computed for %d subjects\n', length(subj_ids));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: MATH - GAZE REINSTATEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==================================================================\n');
fprintf('COMPUTING GAZE REINSTATEMENT METRICS\n');
fprintf('==================================================================\n');

load(fullfile(out_dir, 'group_eye_movement_m.mat'));

% Define trial indices
comp_idx = strcmp(Mw.condition, 'compared'); 
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B'); 
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back'); 
task_2b_idx = strcmp(Mw.task, '2_back');

one_b_comp_idx = task_1b_idx & b_idx & comp_idx; 
one_b_iso_idx = task_1b_idx & b_idx & iso_idx;
two_b_comp_idx = task_2b_idx & ab_idx & b_idx & comp_idx; 
two_b_iso_idx = task_2b_idx & ab_idx & b_idx & iso_idx;

tr_one_b_comp = unique(Mw(one_b_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_one_b_iso = unique(Mw(one_b_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_comp = unique(Mw(two_b_comp_idx, {'subj_id', 'trial_id', 'stim_id'}));
tr_two_b_iso = unique(Mw(two_b_iso_idx, {'subj_id', 'trial_id', 'stim_id'}));

% Match trials across tasks
pairs_comp = innerjoin(tr_one_b_comp, tr_two_b_comp, 'Keys', {'subj_id', 'stim_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id'}, 'RightVariables', {'trial_id'});
pairs_comp.Properties.VariableNames = {'subj_id', 'tr_one_b', 'stim_id', 'tr_two_b'};

pairs_iso = innerjoin(tr_one_b_iso, tr_two_b_iso, 'Keys', {'subj_id', 'stim_id'}, ...
    'LeftVariables', {'subj_id', 'trial_id', 'stim_id'}, 'RightVariables', {'trial_id'});
pairs_iso.Properties.VariableNames = {'subj_id', 'tr_one_b', 'stim_id', 'tr_two_b'};

fprintf('1-back trials: compared=%d, isolated=%d\n', height(tr_one_b_comp), height(tr_one_b_iso));
fprintf('2-back trials: compared=%d, isolated=%d\n', height(tr_two_b_comp), height(tr_two_b_iso));
fprintf('Matched pairs: compared=%d, isolated=%d\n', height(pairs_comp), height(pairs_iso));

% Create baseline pools per subject
unique_subjs = unique(pairs_comp.subj_id);
subj_baseline_comp = struct(); 
subj_baseline_iso = struct();

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx);
    pool_comp = tr_one_b_comp(tr_one_b_comp.subj_id == sid, :);
    pool_iso = tr_one_b_iso(tr_one_b_iso.subj_id == sid, :);
    
    if height(pool_comp) >= n_baseline
        rng(sid); 
        subj_baseline_comp.(sprintf('s%d', sid)) = pool_comp.trial_id(randperm(height(pool_comp), n_baseline));
    end
    
    if height(pool_iso) >= n_baseline
        rng(sid + 10000);
        subj_baseline_iso.(sprintf('s%d', sid)) = pool_iso.trial_id(randperm(height(pool_iso), n_baseline));
    end
end

% Compute compared condition reinstatement
results_comp = table('Size', [height(pairs_comp), 6], ...
    'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'subj_id', 'tr_one_b', 'tr_two_b', 'match_score', 'baseline_score', 'reinst_index'});

n_comp = 0;
for i = 1:height(pairs_comp)
    pair = pairs_comp(i, :); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_comp, skey)
        continue; 
    end
    
    % Encoding map
    fix_enc = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_one_b, :);
    map_enc = create_fixation_map(fix_enc.x, fix_enc.y, fix_enc.dur, spatial_params);
    
    % Retrieval map
    fix_ret = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_two_b, :);
    map_ret = create_fixation_map(fix_ret.x, fix_ret.y, fix_ret.dur, spatial_params);
    
    % Match score
    match_score = corr(map_enc(:), map_ret(:));
    
    % Baseline scores
    baseline_trials = subj_baseline_comp.(skey);
    baseline_scores = zeros(length(baseline_trials), 1);
    for b = 1:length(baseline_trials)
        fix_base = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_trials(b), :);
        map_base = create_fixation_map(fix_base.x, fix_base.y, fix_base.dur, spatial_params);
        baseline_scores(b) = corr(map_base(:), map_ret(:));
    end
    baseline_mean = mean(baseline_scores);
    
    n_comp = n_comp + 1;
    results_comp.subj_id(n_comp) = pair.subj_id;
    results_comp.tr_one_b(n_comp) = pair.tr_one_b;
    results_comp.tr_two_b(n_comp) = pair.tr_two_b;
    results_comp.match_score(n_comp) = match_score;
    results_comp.baseline_score(n_comp) = baseline_mean;
    results_comp.reinst_index(n_comp) = match_score - baseline_mean;
end
results_comp = results_comp(1:n_comp, :);

% Compute isolated condition reinstatement
results_iso = table('Size', [height(pairs_iso), 6], ...
    'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'subj_id', 'tr_one_b', 'tr_two_b', 'match_score', 'baseline_score', 'reinst_index'});

n_iso = 0;
for i = 1:height(pairs_iso)
    pair = pairs_iso(i, :); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_iso, skey)
        continue; 
    end
    
    % Encoding map
    fix_enc = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_one_b, :);
    map_enc = create_fixation_map(fix_enc.x, fix_enc.y, fix_enc.dur, spatial_params);
    
    % Retrieval map
    fix_ret = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_two_b, :);
    map_ret = create_fixation_map(fix_ret.x, fix_ret.y, fix_ret.dur, spatial_params);
    
    % Match score
    match_score = corr(map_enc(:), map_ret(:));
    
    % Baseline scores
    baseline_trials = subj_baseline_iso.(skey);
    baseline_scores = zeros(length(baseline_trials), 1);
    for b = 1:length(baseline_trials)
        fix_base = Mw(Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_trials(b), :);
        map_base = create_fixation_map(fix_base.x, fix_base.y, fix_base.dur, spatial_params);
        baseline_scores(b) = corr(map_base(:), map_ret(:));
    end
    baseline_mean = mean(baseline_scores);
    
    n_iso = n_iso + 1;
    results_iso.subj_id(n_iso) = pair.subj_id;
    results_iso.tr_one_b(n_iso) = pair.tr_one_b;
    results_iso.tr_two_b(n_iso) = pair.tr_two_b;
    results_iso.match_score(n_iso) = match_score;
    results_iso.baseline_score(n_iso) = baseline_mean;
    results_iso.reinst_index(n_iso) = match_score - baseline_mean;
end
results_iso = results_iso(1:n_iso, :);

fprintf('\nGaze reinstatement computed:\n');
fprintf('  Compared: %d trial pairs\n', height(results_comp));
fprintf('  Isolated: %d trial pairs\n', height(results_iso));

% Aggregate to subject level
n_gaze_subjs = length(gaze_subj_ids);
ldi_comp_matched = nan(n_gaze_subjs, 1); 
ldi_iso_matched = nan(n_gaze_subjs, 1);
reinst_comp_matched = nan(n_gaze_subjs, 1); 
reinst_iso_matched = nan(n_gaze_subjs, 1);

for i = 1:n_gaze_subjs
    sid = gaze_subj_ids(i);
    subj_idx = find([all_subjs.id] == sid);
    if ~isempty(subj_idx)
        ldi_comp_matched(i) = all_subjs(subj_idx).stats.two.ldi_comp;
        ldi_iso_matched(i) = all_subjs(subj_idx).stats.two.ldi_iso;
    end
    reinst_comp_matched(i) = mean(results_comp.reinst_index(results_comp.subj_id == sid), 'omitnan');
    reinst_iso_matched(i) = mean(results_iso.reinst_index(results_iso.subj_id == sid), 'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 3: VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==================================================================\n');
fprintf('GENERATING VISUALIZATIONS\n');
fprintf('==================================================================\n');

% Helper function for extracting values
get_v = @(task, field) arrayfun(@(x) x.stats.(task).(field), all_subjs);
fmt_desc = @(v) sprintf('M=%.2f SD=%.2f', mean(v, 'omitnan'), std(v, 'omitnan'));

% --- Figure 1: 1-back behavioral performance ---
figure('Color', 'w', 'Position', [100 100 1000 400]);

subplot(1, 2, 1);
b1_acc_sam = get_v('one', 'acc_same')'; 
b1_acc_sim = get_v('one', 'acc_sim')'; 
b1_acc_new = get_v('one', 'acc_new')';
raincloud([b1_acc_sam, b1_acc_sim, b1_acc_new], {c_same, c_sim, c_new}, ...
    {'Same', 'Similar', 'New'}, 'Accuracy', '1-back Accuracy', []);
add_sig([b1_acc_sam, b1_acc_sim, b1_acc_new], [1 2; 2 3]);

subplot(1, 2, 2);
b1_rt_sam = get_v('one', 'rt_same')'; 
b1_rt_sim = get_v('one', 'rt_sim')';
raincloud([b1_rt_sam, b1_rt_sim], {c_same, c_sim}, ...
    {'Same', 'Similar'}, 'RT (s)', '1-back RT', []);
add_sig([b1_rt_sam, b1_rt_sim], [1 2]);

% --- Figure 2: 2-back LDI performance ---
figure('Color', 'w', 'Position', [100 100 1000 400]);

subplot(1, 2, 1);
b2_ldi_c = get_v('two', 'ldi_comp')'; 
b2_ldi_i = get_v('two', 'ldi_iso')'; 
b2_ldi_n = get_v('two', 'ldi_nov')';
raincloud([b2_ldi_c, b2_ldi_i, b2_ldi_n], {c_comp, c_iso, c_nov}, ...
    {'Compared', 'Isolated', 'Novel'}, 'LDI Score', '2-back Lure Discrimination', []);
add_sig([b2_ldi_c, b2_ldi_i, b2_ldi_n], [1 2; 2 3; 1 3]);

subplot(1, 2, 2);
b2_dp_c = get_v('two', 'dprime_comp')'; 
b2_dp_i = get_v('two', 'dprime_iso')'; 
b2_dp_n = get_v('two', 'dprime_nov')';
raincloud([b2_dp_c, b2_dp_i, b2_dp_n], {c_comp, c_iso, c_nov}, ...
    {'Compared', 'Isolated', 'Novel'}, 'd-prime', '2-back Target Detection', []);
add_sig([b2_dp_c, b2_dp_i, b2_dp_n], [1 2; 2 3; 1 3]);

% --- Figure 3: 2-back RT performance ---
figure('Color', 'w', 'Position', [100 100 1000 400]);

subplot(1, 2, 1);
b2_rt_l_c = get_v('two', 'rt_AB_comp')'; 
b2_rt_l_i = get_v('two', 'rt_AB_iso')'; 
b2_rt_l_n = get_v('two', 'rt_AB_nov')';
raincloud([b2_rt_l_c, b2_rt_l_i, b2_rt_l_n], {c_comp, c_iso, c_nov}, ...
    {'Compared', 'Isolated', 'Novel'}, 'RT (s)', '2-back Lure RT', []);
add_sig([b2_rt_l_c, b2_rt_l_i, b2_rt_l_n], [1 2; 2 3; 1 3]);

subplot(1, 2, 2);
b2_rt_t_c = get_v('two', 'rt_AA_comp')'; 
b2_rt_t_i = get_v('two', 'rt_AA_iso')'; 
b2_rt_t_n = get_v('two', 'rt_AA_nov')';
raincloud([b2_rt_t_c, b2_rt_t_i, b2_rt_t_n], {c_comp, c_iso, c_nov}, ...
    {'Compared', 'Isolated', 'Novel'}, 'RT (s)', '2-back Target RT', []);
add_sig([b2_rt_t_c, b2_rt_t_i, b2_rt_t_n], [1 2; 2 3; 1 3]);

% --- Figure 4: Gaze reinstatement by subject ---
figure('Color', 'w', 'Position', [100 100 1200 400]);

subplot(1, 2, 1);
bar(1:n_gaze_subjs, reinst_comp_matched, 'FaceColor', c_comp, 'EdgeColor', 'k');
set(gca, 'XTick', 1:n_gaze_subjs, 'XTickLabel', gaze_subj_ids, 'FontSize', 10);
xlabel('Subject ID', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title('Gaze Reinstatement - Compared', 'FontSize', 14, 'FontWeight', 'bold'); 
grid on;

subplot(1, 2, 2);
bar(1:n_gaze_subjs, reinst_iso_matched, 'FaceColor', c_iso, 'EdgeColor', 'k');
set(gca, 'XTick', 1:n_gaze_subjs, 'XTickLabel', gaze_subj_ids, 'FontSize', 10);
xlabel('Subject ID', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title('Gaze Reinstatement - Isolated', 'FontSize', 14, 'FontWeight', 'bold'); 
grid on;

% --- Figure 5: Brain-behavior correlations ---
figure('Color', 'w', 'Position', [100 100 1200 500]);

valid_comp = ~isnan(ldi_comp_matched) & ~isnan(reinst_comp_matched);
[r_comp, p_comp] = corr(ldi_comp_matched(valid_comp), reinst_comp_matched(valid_comp));

subplot(1, 2, 1); 
hold on;
scatter(ldi_comp_matched, reinst_comp_matched, 60, c_comp, 'filled', 'MarkerFaceAlpha', 0.6);
p = polyfit(ldi_comp_matched(valid_comp), reinst_comp_matched(valid_comp), 1);
x_fit = linspace(min(ldi_comp_matched(valid_comp)), max(ldi_comp_matched(valid_comp)), 100);
plot(x_fit, polyval(p, x_fit), 'k-', 'LineWidth', 2);
xlabel('LDI Score', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Compared: r=%.3f, p=%.4f', r_comp, p_comp), 'FontSize', 14, 'FontWeight', 'bold');
grid on; 
hold off;

valid_iso = ~isnan(ldi_iso_matched) & ~isnan(reinst_iso_matched);
[r_iso, p_iso] = corr(ldi_iso_matched(valid_iso), reinst_iso_matched(valid_iso));

subplot(1, 2, 2); 
hold on;
scatter(ldi_iso_matched, reinst_iso_matched, 60, c_iso, 'filled', 'MarkerFaceAlpha', 0.6);
p = polyfit(ldi_iso_matched(valid_iso), reinst_iso_matched(valid_iso), 1);
x_fit = linspace(min(ldi_iso_matched(valid_iso)), max(ldi_iso_matched(valid_iso)), 100);
plot(x_fit, polyval(p, x_fit), 'k-', 'LineWidth', 2);
xlabel('LDI Score', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Reinstatement Index', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Isolated: r=%.3f, p=%.4f', r_iso, p_iso), 'FontSize', 14, 'FontWeight', 'bold');
grid on; 
hold off;

fprintf('All visualizations generated\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 4: STATISTICAL SUMMARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==================================================================\n');
fprintf('STATISTICAL SUMMARIES\n');
fprintf('==================================================================\n');

% --- 1-back statistics ---
fprintf('\n--- 1-BACK PERFORMANCE ---\n');
fprintf('Accuracy: Same %s, Similar %s, New %s\n', ...
    fmt_desc(b1_acc_sam), fmt_desc(b1_acc_sim), fmt_desc(b1_acc_new));
do_ttest_print(b1_acc_sam, b1_acc_sim, 'Same vs Similar');
do_ttest_print(b1_acc_sim, b1_acc_new, 'Similar vs New');

fprintf('\nRT: Same %s, Similar %s\n', fmt_desc(b1_rt_sam), fmt_desc(b1_rt_sim));
do_ttest_print(b1_rt_sam, b1_rt_sim, 'Same vs Similar');

% --- 2-back statistics ---
fprintf('\n--- 2-BACK PERFORMANCE ---\n');
fprintf('LDI: Compared %s, Isolated %s, Novel %s\n', ...
    fmt_desc(b2_ldi_c), fmt_desc(b2_ldi_i), fmt_desc(b2_ldi_n));
do_ttest_print(b2_ldi_c, b2_ldi_i, 'Compared vs Isolated');
do_ttest_print(b2_ldi_i, b2_ldi_n, 'Isolated vs Novel');
do_ttest_print(b2_ldi_c, b2_ldi_n, 'Compared vs Novel');

fprintf('\nd-prime: Compared %s, Isolated %s, Novel %s\n', ...
    fmt_desc(b2_dp_c), fmt_desc(b2_dp_i), fmt_desc(b2_dp_n));
do_ttest_print(b2_dp_c, b2_dp_i, 'Compared vs Isolated');
do_ttest_print(b2_dp_i, b2_dp_n, 'Isolated vs Novel');
do_ttest_print(b2_dp_c, b2_dp_n, 'Compared vs Novel');

fprintf('\nLure RT: Compared %s, Isolated %s, Novel %s\n', ...
    fmt_desc(b2_rt_l_c), fmt_desc(b2_rt_l_i), fmt_desc(b2_rt_l_n));
do_ttest_print(b2_rt_l_c, b2_rt_l_i, 'Compared vs Isolated');
do_ttest_print(b2_rt_l_i, b2_rt_l_n, 'Isolated vs Novel');
do_ttest_print(b2_rt_l_c, b2_rt_l_n, 'Compared vs Novel');

fprintf('\nTarget RT: Compared %s, Isolated %s, Novel %s\n', ...
    fmt_desc(b2_rt_t_c), fmt_desc(b2_rt_t_i), fmt_desc(b2_rt_t_n));
do_ttest_print(b2_rt_t_c, b2_rt_t_i, 'Compared vs Isolated');
do_ttest_print(b2_rt_t_i, b2_rt_t_n, 'Isolated vs Novel');
do_ttest_print(b2_rt_t_c, b2_rt_t_n, 'Compared vs Novel');

% --- Gaze reinstatement statistics ---
fprintf('\n--- GAZE REINSTATEMENT ---\n');
fprintf('Compared (n=%d):\n', height(results_comp));
fprintf('  Match: M=%.3f, SD=%.3f\n', ...
    mean(results_comp.match_score, 'omitnan'), std(results_comp.match_score, 'omitnan'));
fprintf('  Baseline: M=%.3f, SD=%.3f\n', ...
    mean(results_comp.baseline_score, 'omitnan'), std(results_comp.baseline_score, 'omitnan'));
fprintf('  Reinstatement Index: M=%.3f, SD=%.3f\n', ...
    mean(results_comp.reinst_index, 'omitnan'), std(results_comp.reinst_index, 'omitnan'));

fprintf('\nIsolated (n=%d):\n', height(results_iso));
fprintf('  Match: M=%.3f, SD=%.3f\n', ...
    mean(results_iso.match_score, 'omitnan'), std(results_iso.match_score, 'omitnan'));
fprintf('  Baseline: M=%.3f, SD=%.3f\n', ...
    mean(results_iso.baseline_score, 'omitnan'), std(results_iso.baseline_score, 'omitnan'));
fprintf('  Reinstatement Index: M=%.3f, SD=%.3f\n', ...
    mean(results_iso.reinst_index, 'omitnan'), std(results_iso.reinst_index, 'omitnan'));

% Permutation tests
fprintf('\n--- PERMUTATION TESTS ---\n');

diffs_c = results_comp.match_score - results_comp.baseline_score; 
diffs_c = diffs_c(~isnan(diffs_c));
obs_c = mean(diffs_c); 
null_c = zeros(n_perm, 1);
for p = 1:n_perm
    null_c(p) = mean(diffs_c .* ((rand(size(diffs_c)) > 0.5) * 2 - 1));
end
p_c = mean(abs(null_c) >= abs(obs_c));
fprintf('Compared reinstatement > 0: M=%.3f, p=%.4f\n', obs_c, p_c);

diffs_i = results_iso.match_score - results_iso.baseline_score; 
diffs_i = diffs_i(~isnan(diffs_i));
obs_i = mean(diffs_i); 
null_i = zeros(n_perm, 1);
for p = 1:n_perm
    null_i(p) = mean(diffs_i .* ((rand(size(diffs_i)) > 0.5) * 2 - 1));
end
p_i = mean(abs(null_i) >= abs(obs_i));
fprintf('Isolated reinstatement > 0: M=%.3f, p=%.4f\n', obs_i, p_i);

obs_diff = obs_c - obs_i; 
pooled = [diffs_c; diffs_i]; 
n_c = length(diffs_c); 
null_diff = zeros(n_perm, 1);
for p = 1:n_perm
    perm_idx = randperm(length(pooled));
    null_diff(p) = mean(pooled(perm_idx(1:n_c))) - mean(pooled(perm_idx(n_c+1:end)));
end
p_diff = mean(abs(null_diff) >= abs(obs_diff));
fprintf('Compared > Isolated: Diff=%.3f, p=%.4f\n', obs_diff, p_diff);

% Brain-behavior correlations
fprintf('\n--- BRAIN-BEHAVIOR CORRELATIONS ---\n');
fprintf('Compared: r=%.3f, p=%.4f\n', r_comp, p_comp);
fprintf('Isolated: r=%.3f, p=%.4f\n', r_iso, p_iso);

% Save results
reinstat_res = struct('compared', results_comp, 'isolated', results_iso, ...
                      'spatial_params', spatial_params);
save(fullfile(res_dir, 'gaze_reinstat_res.mat'), 'reinstat_res');
fprintf('\nResults saved to: %s\n', fullfile(res_dir, 'gaze_reinstat_res.mat'));

fprintf('\n==================================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('==================================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); 
    hold on;
    d_v = mat(:); 
    d_v = d_v(~isnan(d_v)); 
    mn = min(d_v); 
    mx = max(d_v); 
    rng = mx - mn; 
    if rng == 0
        rng = 1; 
    end
    auto_lim = [mn - (rng * 0.15), mx + (rng * 0.15)];
    jit = -0.15 - (rand(size(mat)) * 0.20); 
    x_c = repmat(1:n_grps, n_rows, 1) + jit;
    plot(x_c', mat', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);
    for i = 1:n_grps
        d = mat(:, i); 
        d = d(~isnan(d)); 
        if isempty(d)
            continue; 
        end
        [f, xi] = ksdensity(d); 
        f = f / max(f) * 0.4;
        patch([i+f, i*ones(1, length(f))], [xi, fliplr(xi)], cols{i}, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.5);
        scatter(x_c(:, i), mat(:, i), 20, cols{i}, 'filled', 'MarkerFaceAlpha', 0.6);
        q = quantile(d, [0.25, 0.5, 0.75]);
        rectangle('Position', [i-0.06, q(1), 0.12, q(3)-q(1)], ...
                  'FaceColor', cols{i}, 'EdgeColor', 'k', 'LineWidth', 1.2);
        plot([i-0.06, i+0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 12);
    if ~isempty(ttl)
        title(ttl, 'FontSize', 14); 
    end
    ylabel(ylbl, 'FontSize', 14, 'FontWeight', 'bold'); 
    xlim([0.2, n_grps + 0.8]);
    if nargin > 5 && ~isempty(ylims)
        ylim(ylims); 
    else
        ylim(auto_lim); 
    end
    grid off; 
    set(gca, 'GridAlpha', 0.1); 
    box off; 
    hold off;
end

function add_sig(data, pairs)
    [~, n_grps] = size(data);
    cl = ylim; 
    y_top = cl(2); 
    rng = cl(2) - cl(1); 
    if rng == 0
        rng = 1; 
    end
    step = rng * 0.08; 
    line_lvl = 0; 
    hold on;
    real_mx = max(data(:)); 
    if isnan(real_mx)
        real_mx = y_top; 
    end
    base = max(y_top, real_mx + step * 0.5);
    for i = 1:size(pairs, 1)
        c1 = pairs(i, 1); 
        c2 = pairs(i, 2);
        if c2 == 0
            [~, p] = ttest(data(:, c1)); 
            paired = false; 
        else
            [~, p] = ttest(data(:, c1), data(:, c2)); 
            paired = true; 
        end
        if p < 0.05
            line_lvl = line_lvl + 1; 
            y_p = base + (line_lvl * step);
            if p < 0.001
                txt = '***';
            elseif p < 0.01
                txt = '**';
            else
                txt = '*'; 
            end
            if paired
                plot([c1, c1, c2, c2], [y_p - step*0.3, y_p, y_p, y_p - step*0.3], ...
                     'k-', 'LineWidth', 1.2);
                text(mean([c1 c2]), y_p + step*0.1, txt, ...
                     'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
            else
                text(c1, y_p, txt, 'HorizontalAlignment', 'center', ...
                     'FontSize', 16, 'FontWeight', 'bold'); 
            end
        end
    end
    if line_lvl > 0
        ylim([cl(1), base + (line_lvl * step) + step]); 
    end
    hold off;
end

function do_ttest_print(d1, d2, lbl)
    [~, p, ~, s] = ttest(d1, d2);
    sig = repmat('*', 1, (p < 0.05) + (p < 0.01) + (p < 0.001));
    fprintf('  %s: t(%d)=%.2f, p=%.4f %s\n', lbl, s.df, s.tstat, p, sig);
end

function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & ...
             y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); 
    y = y(in_roi); 
    dur = dur(in_roi);
    
    if isempty(x)
        map = zeros(params.grid_y, params.grid_x); 
        return; 
    end
    
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x + 1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y + 1);
    map = zeros(params.grid_y, params.grid_x);
    
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx)
            map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i); 
        end
    end
    
    if params.sigma > 0
        map = imgaussfilt(map, params.sigma); 
    end
    
    if sum(map(:)) > 0
        map = map / sum(map(:)); 
    end
end