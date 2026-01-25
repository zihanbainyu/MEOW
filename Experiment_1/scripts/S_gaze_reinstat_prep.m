clear; clc; close all;

sub_list = [501, 601, 602, 603, 604, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 617];
tasks = [1, 2]; task_lbls = {'1_back', '2_back'};
base_dir = '..'; out_dir = fullfile(base_dir, 'data', 'eye_movement_data'); behav_dir = fullfile(base_dir, 'data'); res_dir = fullfile(base_dir, 'results');

% [all_fix, all_sac, all_blk, all_sum, all_behav] = deal({});
% 
% for s_idx = 1:length(sub_list)
%     sid = sub_list(s_idx); fprintf('subj %d\n', sid);
%     s_fldr = sprintf('sub%03d', sid); f_path = fullfile(behav_dir, s_fldr);
% 
%     mat_file = fullfile(f_path, sprintf('sub%03d_concat.mat', sid));
%     if isfile(mat_file)
%         load(mat_file, 'final_data_output');
%         r1 = final_data_output.results_1_back_all;
%         r1.task = repmat({'1_back'}, height(r1), 1); r1.subj_id = repmat(sid, height(r1), 1); r1.trial_id = (1:height(r1))';
%         r2 = final_data_output.results_2_back_all;
%         r2.task = repmat({'2_back'}, height(r2), 1); r2.subj_id = repmat(sid, height(r2), 1); r2.trial_id = (1:height(r2))';
% 
%         f1 = fieldnames(r1); f2 = fieldnames(r2); all_f = unique([f1; f2]);
%         for i = 1:length(all_f)
%             if ~ismember(all_f{i}, f1), r1.(all_f{i}) = repmat({''}, height(r1), 1); end
%             if ~ismember(all_f{i}, f2), r2.(all_f{i}) = repmat({''}, height(r2), 1); end
%         end
%         all_behav = [all_behav; {r1}; {r2}];
%     end
% 
%     acc_tr_t1 = 0; acc_tr_t2 = 0;
%     for b = 1:4
%         for t_idx = 1:length(tasks)
%             tsk = tasks(t_idx); t_lbl = task_lbls{t_idx};
%             fn = fullfile(f_path, sprintf('%03d_%01d_%01d.asc', sid, tsk, b));
%             if ~isfile(fn), continue; end
% 
%             fid = fopen(fn); [f_dat, s_dat, b_dat] = deal(cell(2000,1));
%             [n_f, n_s, n_b, blk_max_raw_tr] = deal(0); [curr_tr, st_t, st_id] = deal(0, -1, 'NA');
% 
%             while ~feof(fid)
%                 ln = fgetl(fid);
%                 if startsWith(ln, 'MSG')
%                     if contains(ln, 'TRIALID')
%                         raw_tr = sscanf(ln, 'MSG %*d TRIALID %d');
%                         if raw_tr > 0
%                             blk_max_raw_tr = max(blk_max_raw_tr, raw_tr);
%                             curr_tr = raw_tr + (tsk==1)*acc_tr_t1 + (tsk==2)*acc_tr_t2;
%                             all_sum{end+1,1} = {sid, b, t_lbl, curr_tr, 'NA', NaN, 0, 0, 0}; st_t = -1;
%                         end
%                     elseif contains(ln, 'STIM_ONSET') && curr_tr > 0
%                         tmp = textscan(ln, '%*s %f %*s %s'); st_t = tmp{1}; st_id = tmp{2}{1};
%                         all_sum{end,1}{5} = st_id; all_sum{end,1}{6} = st_t;
%                     end
%                 elseif curr_tr > 0 && st_t > 0
%                     if startsWith(ln, 'EFIX')
%                         d = sscanf(ln, 'EFIX %c %d %d %d %f %f %d');
%                         if numel(d) == 7
%                             n_f = n_f+1; f_dat{n_f} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7)};
%                             all_sum{end,1}{7} = all_sum{end,1}{7} + 1;
%                         end
%                     elseif startsWith(ln, 'ESACC')
%                         d = sscanf(strrep(ln, '.', 'NaN'), 'ESACC %c %d %d %d %f %f %f %f');
%                         if numel(d) == 8
%                             n_s = n_s+1; s_dat{n_s} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7), d(8)};
%                             all_sum{end,1}{8} = all_sum{end,1}{8} + 1;
%                         end
%                     elseif startsWith(ln, 'EBLINK')
%                         d = sscanf(ln, 'EBLINK %c %d %d %d');
%                         if numel(d) == 4
%                             n_b = n_b+1; b_dat{n_b} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4)};
%                             all_sum{end,1}{9} = all_sum{end,1}{9} + 1;
%                         end
%                     end
%                 end
%             end
%             fclose(fid);
% 
%             if tsk == 1, acc_tr_t1 = acc_tr_t1 + blk_max_raw_tr; else, acc_tr_t2 = acc_tr_t2 + blk_max_raw_tr; end
%             if n_f > 0, all_fix = [all_fix; f_dat(1:n_f)]; end
%             if n_s > 0, all_sac = [all_sac; s_dat(1:n_s)]; end
%             if n_b > 0, all_blk = [all_blk; b_dat(1:n_b)]; end
%         end
%     end
% end
% 
% v_fix = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','x','y','pupil'};
% v_sac = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','sx','sy','ex','ey'};
% v_blk = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur'};
% v_sum = {'subj_id','block','task','trial_id','stim_id','onset_time','n_fix','n_sac','n_blink'};
% 
% fix = cell2table(vertcat(all_fix{:}), 'VariableNames', v_fix);
% sac = cell2table(vertcat(all_sac{:}), 'VariableNames', v_sac);
% blk = cell2table(vertcat(all_blk{:}), 'VariableNames', v_blk);
% sum_tab = cell2table(vertcat(all_sum{:}), 'VariableNames', v_sum);
% 
% all_behav = vertcat(all_behav{:});
% M = outerjoin(fix, all_behav, 'Keys', {'subj_id','task','trial_id','block'}, 'MergeKeys', true, 'Type', 'inner');
% M = removevars(M, 'stim_id_all_behav'); M = renamevars(M, 'stim_id_fix', 'stim_id');
% 
% eye_data = struct(); eye_data.fixations = fix; eye_data.saccades = sac; eye_data.blinks = blk; eye_data.summary = sum_tab;
% eye_data.merged = M; eye_data.merged.resp_key = cellstr(eye_data.merged.resp_key);
% eye_data.merged.resp_key(strcmp(eye_data.merged.resp_key,'NA')) = {'none'};
% eye_data.merged.correct = strcmp(cellstr(eye_data.merged.corr_resp), eye_data.merged.resp_key);
% save(fullfile(out_dir, 'group_eye_movement.mat'), 'eye_data', '-v7.3');


load(fullfile(out_dir, 'group_eye_movement_m.mat'));

comp_idx = strcmp(Mw.condition, 'compared'); 
iso_idx = strcmp(Mw.condition, 'isolated');
ab_idx = strcmp(Mw.goal, 'A-B'); 
a_idx = strcmp(Mw.identity, 'A');
b_idx = strcmp(Mw.identity, 'B');
task_1b_idx = strcmp(Mw.task, '1_back'); 
task_2b_idx = strcmp(Mw.task, '2_back');

one_b_b_comp_idx = task_1b_idx & b_idx & comp_idx; 
one_b_b_iso_idx = task_1b_idx & b_idx & iso_idx;
two_b_b_comp_idx = task_2b_idx & ab_idx & b_idx & comp_idx; 
two_b_b_iso_idx = task_2b_idx & ab_idx & b_idx & iso_idx;
two_b_a_comp_idx = task_2b_idx & ab_idx & a_idx & comp_idx; 
two_b_a_iso_idx = task_2b_idx & ab_idx & a_idx & iso_idx;

tr_one_b_b_comp = unique(Mw(one_b_b_comp_idx, {'subj_id','trial_id','stim_id'}));
tr_one_b_b_iso = unique(Mw(one_b_b_iso_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_b_comp = unique(Mw(two_b_b_comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_b_iso = unique(Mw(two_b_b_iso_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_comp = unique(Mw(two_b_a_comp_idx, {'subj_id','trial_id','stim_id'}));
tr_two_b_a_iso = unique(Mw(two_b_a_iso_idx, {'subj_id','trial_id','stim_id'}));

pairs_bb_comp = innerjoin(tr_one_b_b_comp, tr_two_b_b_comp, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_bb_comp.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id','tr_2b_b'};

pairs_bb_iso = innerjoin(tr_one_b_b_iso, tr_two_b_b_iso, 'Keys', {'subj_id','stim_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id'}, 'RightVariables', {'trial_id'});
pairs_bb_iso.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id','tr_2b_b'};

tr_one_b_b_comp.base_id = regexprep(tr_one_b_b_comp.stim_id, '_B_', '_BASE_');
tr_one_b_b_iso.base_id = regexprep(tr_one_b_b_iso.stim_id, '_B_', '_BASE_');
tr_two_b_a_comp.base_id = regexprep(tr_two_b_a_comp.stim_id, '_A_', '_BASE_');
tr_two_b_a_iso.base_id = regexprep(tr_two_b_a_iso.stim_id, '_A_', '_BASE_');

pairs_ba_comp = innerjoin(tr_one_b_b_comp, tr_two_b_a_comp, 'Keys', {'subj_id','base_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
pairs_ba_comp.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id_b','base_id','tr_2b_a','stim_id_a'};

pairs_ba_iso = innerjoin(tr_one_b_b_iso, tr_two_b_a_iso, 'Keys', {'subj_id','base_id'}, ...
    'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
pairs_ba_iso.Properties.VariableNames = {'subj_id','tr_1b_b','stim_id_b','base_id','tr_2b_a','stim_id_a'};

fprintf('valid fixations: %d\n', height(Mw));
fprintf('1back B: comp=%d, iso=%d\n', height(tr_one_b_b_comp), height(tr_one_b_b_iso));
fprintf('2back B: comp=%d, iso=%d\n', height(tr_two_b_b_comp), height(tr_two_b_b_iso));
fprintf('2back A: comp=%d, iso=%d\n', height(tr_two_b_a_comp), height(tr_two_b_a_iso));
fprintf('matched pairs B-B: comp=%d, iso=%d\n', height(pairs_bb_comp), height(pairs_bb_iso));
fprintf('matched pairs B-A: comp=%d, iso=%d\n\n', height(pairs_ba_comp), height(pairs_ba_iso));

spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);
n_baseline = 20;

unique_subjs = unique(pairs_bb_comp.subj_id);
subj_baseline_b_comp = struct(); 
subj_baseline_b_iso = struct();

for s_idx = 1:length(unique_subjs)
    sid = unique_subjs(s_idx);
    pool_b_comp = tr_one_b_b_comp(tr_one_b_b_comp.subj_id==sid, :);
    pool_b_iso = tr_one_b_b_iso(tr_one_b_b_iso.subj_id==sid, :);
    
    if height(pool_b_comp) >= n_baseline
        rng(sid); 
        subj_baseline_b_comp.(sprintf('s%d', sid)) = pool_b_comp.trial_id(randperm(height(pool_b_comp), n_baseline));
    end
    
    if height(pool_b_iso) >= n_baseline
        rng(sid+10000);
        subj_baseline_b_iso.(sprintf('s%d', sid)) = pool_b_iso.trial_id(randperm(height(pool_b_iso), n_baseline));
    end
end

results_bb_comp = zeros(height(pairs_bb_comp), 6);
n_bb_comp = 0;

for i = 1:height(pairs_bb_comp)
    pair = pairs_bb_comp(i,:); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_comp, skey), continue; end
    
    baseline_trials = subj_baseline_b_comp.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_1b_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_2b_b, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==baseline_trials(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_bb_comp = n_bb_comp + 1;
    results_bb_comp(n_bb_comp,:) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_b, match_score, baseline_score, match_score - baseline_score];
end

results_bb_comp = results_bb_comp(1:n_bb_comp,:);
results_bb_iso = zeros(height(pairs_bb_iso), 6);
n_bb_iso = 0;

for i = 1:height(pairs_bb_iso)
    pair = pairs_bb_iso(i,:); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_iso, skey), continue; end
    
    baseline_trials = subj_baseline_b_iso.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_1b_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_2b_b, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==baseline_trials(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_bb_iso = n_bb_iso + 1;
    results_bb_iso(n_bb_iso,:) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_b, match_score, baseline_score, match_score - baseline_score];
end

results_bb_iso = results_bb_iso(1:n_bb_iso,:);

results_ba_comp = zeros(height(pairs_ba_comp), 6);
n_ba_comp = 0;

for i = 1:height(pairs_ba_comp)
    pair = pairs_ba_comp(i,:); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_comp, skey), continue; end
    
    baseline_trials = subj_baseline_b_comp.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_1b_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_2b_a, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==baseline_trials(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_ba_comp = n_ba_comp + 1;
    results_ba_comp(n_ba_comp,:) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_a, match_score, baseline_score, match_score - baseline_score];
end

results_ba_comp = results_ba_comp(1:n_ba_comp,:);
results_ba_iso = zeros(height(pairs_ba_iso), 6);
n_ba_iso = 0;

for i = 1:height(pairs_ba_iso)
    pair = pairs_ba_iso(i,:); 
    skey = sprintf('s%d', pair.subj_id);
    if ~isfield(subj_baseline_b_iso, skey), continue; end
    
    baseline_trials = subj_baseline_b_iso.(skey);
    baseline_trials = baseline_trials(baseline_trials ~= pair.tr_1b_b);
    if length(baseline_trials) < n_baseline, continue; end
    baseline_trials = baseline_trials(1:n_baseline);
    
    fix_1b_match = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_1b_b, :);
    fix_2b = Mw(task_2b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==pair.tr_2b_a, :);
    if height(fix_1b_match) < 2 || height(fix_2b) < 2, continue; end
    
    map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);
    map_1b_match = create_fixation_map(fix_1b_match.x, fix_1b_match.y, fix_1b_match.dur, spatial_params);
    match_score = corr(map_1b_match(:), map_2b(:));
    
    mismatch_scores = nan(n_baseline, 1);
    for j = 1:n_baseline
        fix_1b_mismatch = Mw(task_1b_idx & Mw.subj_id==pair.subj_id & Mw.trial_id==baseline_trials(j), :);
        if height(fix_1b_mismatch) < 2, continue; end
        map_1b_mismatch = create_fixation_map(fix_1b_mismatch.x, fix_1b_mismatch.y, fix_1b_mismatch.dur, spatial_params);
        mismatch_scores(j) = corr(map_1b_mismatch(:), map_2b(:));
    end
    
    baseline_score = mean(mismatch_scores, 'omitnan');
    n_ba_iso = n_ba_iso + 1;
    results_ba_iso(n_ba_iso,:) = [pair.subj_id, pair.tr_1b_b, pair.tr_2b_a, match_score, baseline_score, match_score - baseline_score];
end

results_ba_iso = results_ba_iso(1:n_ba_iso,:);

results_bb_comp = array2table(results_bb_comp, 'VariableNames', {'subj_id','tr_1b_b','tr_2b_b','match_score','baseline_score','reinst_index'});
results_bb_iso = array2table(results_bb_iso, 'VariableNames', {'subj_id','tr_1b_b','tr_2b_b','match_score','baseline_score','reinst_index'});
results_ba_comp = array2table(results_ba_comp, 'VariableNames', {'subj_id','tr_1b_b','tr_2b_a','match_score','baseline_score','reinst_index'});
results_ba_iso = array2table(results_ba_iso, 'VariableNames', {'subj_id','tr_1b_b','tr_2b_a','match_score','baseline_score','reinst_index'});

reinstat_res = struct();
reinstat_res.bb_compared = results_bb_comp;
reinstat_res.bb_isolated = results_bb_iso;
reinstat_res.ba_compared = results_ba_comp;
reinstat_res.ba_isolated = results_ba_iso;
reinstat_res.spatial_params = spatial_params;

save(fullfile(res_dir, 'gaze_reinstat_res_m.mat'), 'reinstat_res');

function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); 
    y = y(in_roi); 
    dur = dur(in_roi);
    if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
    map = zeros(params.grid_y, params.grid_x);
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx)
            map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i); 
        end
    end
    if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
    if sum(map(:)) > 0, map = map / sum(map(:)); end
end