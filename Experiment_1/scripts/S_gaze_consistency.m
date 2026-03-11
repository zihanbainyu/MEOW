clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');
out_dir = fullfile(res_dir, 'gaze_aa_consistency');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

load(fullfile(base_dir, 'data', 'eye_movement_data', 'group_eye_movement_combined.mat'), 'Mw');
load(fullfile(res_dir, 'gaze_reinstat_res_full.mat'), 'reinstat_res');

spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);

comp_idx = strcmp(Mw.condition, 'compared');
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B');
a_idx = strcmp(Mw.identity, 'A');
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back');
task_2b_idx = strcmp(Mw.task, '2_back');

tr_one_b_a_comp = unique(Mw(task_1b_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_a_iso = unique(Mw(task_1b_idx & a_idx & iso_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_comp = unique(Mw(task_1b_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_iso = unique(Mw(task_1b_idx & b_idx & iso_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_comp = unique(Mw(task_2b_idx & ab_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_iso = unique(Mw(task_2b_idx & ab_idx & a_idx & iso_idx, {'subj_id','trial_id','stim_id'}));

pairs_aa_comp = innerjoin(tr_one_b_a_comp, tr_two_b_a_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_comp.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

pairs_aa_iso = innerjoin(tr_one_b_a_iso, tr_two_b_a_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_aa_iso.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id','tr_2b_a'};

tr_two_b_a_iso_tmp = tr_two_b_a_iso;
tr_two_b_a_iso_tmp.base_id = regexprep(tr_two_b_a_iso_tmp.stim_id, '_A_', '_BASE_');
tr_one_b_b_iso_tmp = tr_one_b_b_iso;
tr_one_b_b_iso_tmp.base_id = regexprep(tr_one_b_b_iso_tmp.stim_id, '_B_', '_BASE_');
pairs_ab_iso = innerjoin(tr_two_b_a_iso_tmp, tr_one_b_b_iso_tmp, 'Keys', {'subj_id','base_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
pairs_ab_iso.Properties.VariableNames = {'subj_id','tr_2b_a','stim_id_a','base_id','tr_1b_b','stim_id_b'};

fprintf('A2-A1 compared: %d pairs\nA2-A1 isolated: %d pairs\nA2-B1 isolated: %d pairs\n', ...
    height(pairs_aa_comp), height(pairs_aa_iso), height(pairs_ab_iso));

unique_subjs = unique([pairs_aa_comp.subj_id; pairs_aa_iso.subj_id; pairs_ab_iso.subj_id]);
n_subjs = length(unique_subjs);

subj_pool_b_comp = struct(); subj_pool_b_iso = struct();
for s = 1:n_subjs
    sid = unique_subjs(s); skey = sprintf('s%d', sid);
    sb = tr_one_b_b_comp(tr_one_b_b_comp.subj_id == sid, :);
    if height(sb) > 0, subj_pool_b_comp.(skey) = sb; end
    sb = tr_one_b_b_iso(tr_one_b_b_iso.subj_id == sid, :);
    if height(sb) > 0, subj_pool_b_iso.(skey) = sb; end
end

gaze_aa_comp = nan(n_subjs, 1); gaze_aa_iso = nan(n_subjs, 1); gaze_ab_iso = nan(n_subjs, 1);
gaze_aa_corr = nan(n_subjs, 2); gaze_aa_incorr = nan(n_subjs, 2);

for s = 1:n_subjs
    sid = unique_subjs(s); skey = sprintf('s%d', sid);

    for c = 1:2
        if c == 1, pairs = pairs_aa_comp(pairs_aa_comp.subj_id == sid, :); pool_struct = subj_pool_b_comp; bb = reinstat_res.bb_compared;
        else, pairs = pairs_aa_iso(pairs_aa_iso.subj_id == sid, :); pool_struct = subj_pool_b_iso; bb = reinstat_res.bb_isolated; end
        if height(pairs) == 0 || ~isfield(pool_struct, skey), continue; end
        pool = pool_struct.(skey);
        r_all = []; r_corr = []; r_incorr = [];
        for i = 1:height(pairs)
            fix_2b = Mw(task_2b_idx & Mw.subj_id == sid & Mw.trial_id == pairs.tr_2b_a(i), :);
            fix_1b = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == pairs.tr_1b_a(i), :);
            if height(fix_2b) < 2 || height(fix_1b) < 2, continue; end
            map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
            map_1b = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
            if sum(map_2b(:)) == 0 || sum(map_1b(:)) == 0, continue; end
            r_match = corr(map_2b(:), map_1b(:));
            partner_b = regexprep(pairs.stim_id{i}, '_A_', '_B_');
            bp = pool(~strcmp(pool.stim_id, partner_b), :);
            if height(bp) < 1, continue; end
            mm = nan(height(bp), 1);
            for j = 1:height(bp)
                fx = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == bp.trial_id(j), :);
                if height(fx) < 2, continue; end
                mp = create_fixation_map(fx.x, fx.y, fx.dur, spatial_params);
                if sum(mp(:)) > 0, mm(j) = corr(map_2b(:), mp(:)); end
            end
            r_val = r_match - mean(mm, 'omitnan');
            r_all(end+1) = r_val;
            mb2 = bb(bb.subj_id == sid & strcmp(bb.stim_id, partner_b), :);
            if height(mb2) == 1
                if mb2.correct == 1, r_corr(end+1) = r_val; else, r_incorr(end+1) = r_val; end
            end
        end
        if c == 1, gaze_aa_comp(s) = mean(r_all, 'omitnan'); else, gaze_aa_iso(s) = mean(r_all, 'omitnan'); end
        if ~isempty(r_corr), gaze_aa_corr(s,c) = mean(r_corr); end
        if ~isempty(r_incorr), gaze_aa_incorr(s,c) = mean(r_incorr); end
    end

    sp = pairs_ab_iso(pairs_ab_iso.subj_id == sid, :);
    if height(sp) > 0 && isfield(subj_pool_b_iso, skey)
        pool = subj_pool_b_iso.(skey);
        r_all = [];
        for i = 1:height(sp)
            fix_2b = Mw(task_2b_idx & Mw.subj_id == sid & Mw.trial_id == sp.tr_2b_a(i), :);
            fix_1b = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == sp.tr_1b_b(i), :);
            if height(fix_2b) < 2 || height(fix_1b) < 2, continue; end
            map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
            map_1b = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
            if sum(map_2b(:)) == 0 || sum(map_1b(:)) == 0, continue; end
            r_match = corr(map_2b(:), map_1b(:));
            bp = pool(pool.trial_id ~= sp.tr_1b_b(i), :);
            if height(bp) < 1, continue; end
            mm = nan(height(bp), 1);
            for j = 1:height(bp)
                fx = Mw(task_1b_idx & Mw.subj_id == sid & Mw.trial_id == bp.trial_id(j), :);
                if height(fx) < 2, continue; end
                mp = create_fixation_map(fx.x, fx.y, fx.dur, spatial_params);
                if sum(mp(:)) > 0, mm(j) = corr(map_2b(:), mp(:)); end
            end
            r_all(end+1) = r_match - mean(mm, 'omitnan');
        end
        gaze_ab_iso(s) = mean(r_all, 'omitnan');
    end

    if mod(s, 5) == 0, fprintf('  %d/%d subjects\n', s, n_subjs); end
end

valid = ~isnan(gaze_aa_comp) & ~isnan(gaze_aa_iso) & ~isnan(gaze_ab_iso);
gc = gaze_aa_comp(valid); gi = gaze_aa_iso(valid); gabi = gaze_ab_iso(valid);

fprintf('\n--- Baseline-Corrected Gaze Consistency (N=%d) ---\n', sum(valid));
fprintf('A2-A1 compared: M=%.4f (SD=%.4f)\n', mean(gc), std(gc));
fprintf('A2-A1 isolated: M=%.4f (SD=%.4f)\n', mean(gi), std(gi));
fprintf('A2-B1 isolated: M=%.4f (SD=%.4f)\n', mean(gabi), std(gabi));

fprintf('\nPairwise tests:\n');
[~, p1, ~, s1] = ttest(gc, gi); d1 = mean(gc-gi)/std(gc-gi);
fprintf('  A2A1 comp vs A2A1 iso: t(%d)=%.3f, p=%.4f, d=%.3f\n', s1.df, s1.tstat, p1, d1);
[~, p2, ~, s2] = ttest(gi, gabi); d2 = mean(gi-gabi)/std(gi-gabi);
fprintf('  A2A1 iso vs A2B1 iso:  t(%d)=%.3f, p=%.4f, d=%.3f\n', s2.df, s2.tstat, p2, d2);
[~, p3, ~, s3] = ttest(gc, gabi); d3 = mean(gc-gabi)/std(gc-gabi);
fprintf('  A2A1 comp vs A2B1 iso: t(%d)=%.3f, p=%.4f, d=%.3f\n', s3.df, s3.tstat, p3, d3);

save(fullfile(out_dir, 'gaze_aa_consistency_results_new.mat'), ...
    'gaze_aa_comp', 'gaze_aa_iso', 'gaze_ab_iso', 'gaze_aa_corr', 'gaze_aa_incorr', ...
    'unique_subjs', 'pairs_aa_comp', 'pairs_aa_iso', 'spatial_params');

fprintf('Done.\n');


function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
    if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
    map = zeros(params.grid_y, params.grid_x);
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx), map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i); end
    end
    if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
    if sum(map(:)) > 0, map = map / sum(map(:)); end
end