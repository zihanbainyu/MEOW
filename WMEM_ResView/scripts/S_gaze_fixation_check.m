%==========================================================================
%   WMEM_ResView -- restricted-viewing manipulation check (eye data)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  Config
%  ------------------------------------------------------------------------
sub_list = [601];
tasks = [1, 2]; task_lbls = {'1_back', '2_back'};
base_dir = '..';
behav_dir = fullfile(base_dir, 'data');
out_dir   = fullfile(base_dir, 'data', 'eye_movement_data');
res_dir   = fullfile(base_dir, 'results');
fig_dir   = fullfile(base_dir, 'figures');
for d = {out_dir, res_dir, fig_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

% Enforced radii from scripts/main.m
FIX_GATE_TOL_PX = 100;   % onset-gate window radius (~1.6 deg)
FIX_TOL_PX      = 150;   % break window radius during viewing (~2.4 deg)

% Stimulus extent. The run scripts call Screen('DrawTexture', ..., [], [], 0)
% with an empty destination rect, so each image is drawn at its native size,
% centred. Verified against a real stimulus further down.
IMG_SIZE_PX = 400;

% Condition colours, matching S_behavioral_analysis_indiv.m
c_comp = [180 174 211]/255;
c_iso  = [176 230 255]/255;
c_nov  = [183 210 205]/255;
cond_names  = {'compared', 'isolated', 'novel'};
cond_colors = {c_comp, c_iso, c_nov};

%% ------------------------------------------------------------------------
%  Parse .asc + merge with behaviour   (kept from S_gaze_reinstat_prep.m)
%  ------------------------------------------------------------------------
[all_fix, all_sac, all_blk, all_sum, all_behav] = deal({});
screen_rect = [];   % filled from the GAZE_COORDS message

for s_idx = 1:length(sub_list)
    sid = sub_list(s_idx); fprintf('subj %d\n', sid);
    s_fldr = sprintf('sub%03d', sid); f_path = fullfile(behav_dir, s_fldr);

    mat_file = fullfile(f_path, sprintf('sub%03d_concat.mat', sid));
    if isfile(mat_file)
        load(mat_file, 'final_data_output');
        r1 = final_data_output.results_1_back_all;
        r1.task = repmat({'1_back'}, height(r1), 1); r1.subj_id = repmat(sid, height(r1), 1); r1.trial_id = (1:height(r1))';
        r2 = final_data_output.results_2_back_all;
        r2.task = repmat({'2_back'}, height(r2), 1); r2.subj_id = repmat(sid, height(r2), 1); r2.trial_id = (1:height(r2))';

        f1 = fieldnames(r1); f2 = fieldnames(r2); all_f = unique([f1; f2]);
        for i = 1:length(all_f)
            if ~ismember(all_f{i}, f1), r1.(all_f{i}) = repmat({''}, height(r1), 1); end
            if ~ismember(all_f{i}, f2), r2.(all_f{i}) = repmat({''}, height(r2), 1); end
        end
        all_behav = [all_behav; {r1}; {r2}];
    end

    acc_tr_t1 = 0; acc_tr_t2 = 0;
    for b = 1:4
        for t_idx = 1:length(tasks)
            tsk = tasks(t_idx); t_lbl = task_lbls{t_idx};
            fn = fullfile(f_path, sprintf('%03d_%01d_%01d.asc', sid, tsk, b));
            if ~isfile(fn), continue; end

            fid = fopen(fn); [f_dat, s_dat, b_dat] = deal(cell(2000,1));
            [n_f, n_s, n_b, blk_max_raw_tr] = deal(0); [curr_tr, st_t, st_id] = deal(0, -1, 'NA');

            while ~feof(fid)
                ln = fgetl(fid);
                if startsWith(ln, 'MSG')
                    if contains(ln, 'GAZE_COORDS') && isempty(screen_rect)
                        gc = sscanf(ln, 'MSG %*d GAZE_COORDS %f %f %f %f');
                        if numel(gc) == 4, screen_rect = gc(:)'; end
                    elseif contains(ln, 'TRIALID')
                        raw_tr = sscanf(ln, 'MSG %*d TRIALID %d');
                        if raw_tr > 0
                            blk_max_raw_tr = max(blk_max_raw_tr, raw_tr);
                            curr_tr = raw_tr + (tsk==1)*acc_tr_t1 + (tsk==2)*acc_tr_t2;
                            all_sum{end+1,1} = {sid, b, t_lbl, curr_tr, 'NA', NaN, 0, 0, 0}; st_t = -1;
                        end
                    elseif contains(ln, 'STIM_ONSET') && curr_tr > 0
                        tmp = textscan(ln, '%*s %f %*s %s'); st_t = tmp{1}; st_id = tmp{2}{1};
                        all_sum{end,1}{5} = st_id; all_sum{end,1}{6} = st_t;
                    end
                elseif curr_tr > 0 && st_t > 0
                    if startsWith(ln, 'EFIX')
                        d = sscanf(ln, 'EFIX %c %d %d %d %f %f %d');
                        if numel(d) == 7
                            n_f = n_f+1; f_dat{n_f} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7)};
                            all_sum{end,1}{7} = all_sum{end,1}{7} + 1;
                        end
                    elseif startsWith(ln, 'ESACC')
                        d = sscanf(strrep(ln, '.', 'NaN'), 'ESACC %c %d %d %d %f %f %f %f');
                        if numel(d) == 8
                            n_s = n_s+1; s_dat{n_s} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7), d(8)};
                            all_sum{end,1}{8} = all_sum{end,1}{8} + 1;
                        end
                    elseif startsWith(ln, 'EBLINK')
                        d = sscanf(ln, 'EBLINK %c %d %d %d');
                        if numel(d) == 4
                            n_b = n_b+1; b_dat{n_b} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4)};
                            all_sum{end,1}{9} = all_sum{end,1}{9} + 1;
                        end
                    end
                end
            end
            fclose(fid);

            if tsk == 1, acc_tr_t1 = acc_tr_t1 + blk_max_raw_tr; else, acc_tr_t2 = acc_tr_t2 + blk_max_raw_tr; end
            if n_f > 0, all_fix = [all_fix; f_dat(1:n_f)]; end
            if n_s > 0, all_sac = [all_sac; s_dat(1:n_s)]; end
            if n_b > 0, all_blk = [all_blk; b_dat(1:n_b)]; end
        end
    end
end

%% ------------------------------------------------------------------------
%  Sample-validity pass (data loss per file)
%  ------------------------------------------------------------------------
% The online detectors in C_run_1_back.m / D_run_2_back.m treat a missing
% sample as "not a break" (out_since = NaN), so heavy data loss shows up as
% gate_timeout while fix_broken stays at 0%. Counting missing samples
% directly is what separates a tracker problem from subject behaviour.
qual = table();
for s_idx = 1:length(sub_list)
    sid = sub_list(s_idx);
    f_path = fullfile(behav_dir, sprintf('sub%03d', sid));
    for b = 1:4
        for t_idx = 1:length(tasks)
            fn = fullfile(f_path, sprintf('%03d_%01d_%01d.asc', sid, tasks(t_idx), b));
            if ~isfile(fn), continue; end
            txt = fileread(fn);
            n_samp = numel(regexp(txt, '\n\d+\t'));          % any sample line
            n_miss = numel(regexp(txt, '\n\d+\t\s*\.\s*\t')); % x field is '.'
            qual = [qual; table(sid, b, string(task_lbls{t_idx}), n_samp, n_miss, ...
                100*n_miss/max(n_samp,1), 'VariableNames', ...
                {'subj_id','block','task','n_samples','n_missing','pct_missing'})];
        end
    end
end

v_fix = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','x','y','pupil'};
v_sac = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','sx','sy','ex','ey'};
v_blk = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur'};
v_sum = {'subj_id','block','task','trial_id','stim_id','onset_time','n_fix','n_sac','n_blink'};

fix     = cell2table(vertcat(all_fix{:}), 'VariableNames', v_fix);
sac     = cell2table(vertcat(all_sac{:}), 'VariableNames', v_sac);
blk     = cell2table(vertcat(all_blk{:}), 'VariableNames', v_blk);
sum_tab = cell2table(vertcat(all_sum{:}), 'VariableNames', v_sum);

all_behav = vertcat(all_behav{:});
M = outerjoin(fix, all_behav, 'Keys', {'subj_id','task','trial_id','block'}, 'MergeKeys', true, 'Type', 'inner');
M = removevars(M, 'stim_id_all_behav'); M = renamevars(M, 'stim_id_fix', 'stim_id');
M.resp_key = cellstr(M.resp_key);
M.resp_key(strcmp(M.resp_key,'NA')) = {'none'};
M.correct = strcmp(cellstr(M.corr_resp), M.resp_key);

% attach stimulus onset time so fixation latency can be recovered
M = outerjoin(M, sum_tab(:, {'subj_id','task','trial_id','block','onset_time'}), ...
    'Keys', {'subj_id','task','trial_id','block'}, 'MergeKeys', true, 'Type', 'left');
M.t_rel = M.onset - M.onset_time;   % fixation start, ms relative to stimulus onset

eye_data = struct('fixations', fix, 'saccades', sac, 'blinks', blk, ...
                  'summary', sum_tab, 'merged', M);
save(fullfile(out_dir, 'group_eye_movement.mat'), '-struct', 'eye_data');

%% ------------------------------------------------------------------------
%  Deviation from screen centre
%  ------------------------------------------------------------------------
if isempty(screen_rect)
    warning('No GAZE_COORDS message found; falling back to 1920x1080.');
    screen_rect = [0 0 1919 1079];
end
cx = (screen_rect(1) + screen_rect(3)) / 2;
cy = (screen_rect(2) + screen_rect(4)) / 2;
fprintf('\nScreen %g x %g, centre (%.1f, %.1f)\n', ...
    screen_rect(3)+1, screen_rect(4)+1, cx, cy);

M.dev = hypot(M.x - cx, M.y - cy);
M = M(isfinite(M.dev), :);
M.condition = string(M.condition);
M.task = string(M.task);

% keep only the three analysed conditions (drops 1-back 'repeat' and 2-back junk)
keep = ismember(M.condition, string(cond_names));
Mc = M(keep, :);

fprintf('fixations parsed: %d | in the 3 conditions: %d\n', height(M), height(Mc));
fprintf('  by condition: ');
for k = 1:numel(cond_names)
    fprintf('%s=%d  ', cond_names{k}, sum(Mc.condition == cond_names{k}));
end
fprintf('\n  by task:      1_back=%d  2_back=%d\n', ...
    sum(Mc.task=="1_back"), sum(Mc.task=="2_back"));

%% ------------------------------------------------------------------------
%  FIGURE 1 -- fixation overlay, 1-back and 2-back side by side
%  ------------------------------------------------------------------------
% Note: 'novel' exists only in the 2-back, so the 1-back panel legitimately
% shows two conditions and the 2-back three.
th = linspace(0, 2*pi, 400);
tasks_plot = ["1_back", "2_back"];

% Axes span the full display (1920 x 1080), so the panels are true-shaped
% screens and the size of the restricted cluster is shown in context.
sx0 = screen_rect(1); sy0 = screen_rect(2);
sx1 = screen_rect(3) + 1; sy1 = screen_rect(4) + 1;

% Confirm the assumed stimulus size against an actual image on disk
stim_files = dir(fullfile(base_dir, 'stimulus', 'stim_final', 'mst_*.png'));
if ~isempty(stim_files)
    info = imfinfo(fullfile(stim_files(1).folder, stim_files(1).name));
    if info.Width ~= IMG_SIZE_PX || info.Height ~= IMG_SIZE_PX
        warning(['Stimulus %s is %dx%d px but IMG_SIZE_PX = %d. ' ...
                 'The drawn image frame will not match the real one.'], ...
                 stim_files(1).name, info.Width, info.Height, IMG_SIZE_PX);
    end
    fprintf('  stimulus size on disk: %d x %d px (frame drawn at %d x %d)\n', ...
        info.Width, info.Height, IMG_SIZE_PX, IMG_SIZE_PX);
end
img_rect = [cx - IMG_SIZE_PX/2, cy - IMG_SIZE_PX/2, IMG_SIZE_PX, IMG_SIZE_PX];

f1 = figure('Visible','off','Position',[100 100 980 1010],'Color','w');
tl = tiledlayout(f1, 2, 1, 'TileSpacing','compact', 'Padding','compact');

hleg = gobjects(1, numel(cond_names)); leg_lbl = cell(1, numel(cond_names));
for t = 1:numel(tasks_plot)
    ax = nexttile(tl); hold(ax,'on'); axis(ax,'equal');
    mt = Mc.task == tasks_plot(t);

    % stimulus extent, drawn first so fixations sit on top of it
    phI = plot(ax, img_rect(1) + [0 IMG_SIZE_PX IMG_SIZE_PX 0 0], ...
                   img_rect(2) + [0 0 IMG_SIZE_PX IMG_SIZE_PX 0], ...
                   '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.5);

    % larger radius first so the smaller sits on top
    ph2 = plot(ax, cx + FIX_TOL_PX*cos(th),      cy + FIX_TOL_PX*sin(th),      'k--', 'LineWidth', 1.5);
    ph1 = plot(ax, cx + FIX_GATE_TOL_PX*cos(th), cy + FIX_GATE_TOL_PX*sin(th), 'k-',  'LineWidth', 1.5);

    for k = 1:numel(cond_names)
        m = mt & Mc.condition == cond_names{k};
        if ~any(m), continue; end
        s = scatter(ax, Mc.x(m), Mc.y(m), 14, cond_colors{k}, 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor','none');
        if t == numel(tasks_plot) || ~isgraphics(hleg(k))
            hleg(k) = s; leg_lbl{k} = cond_names{k};
        end
    end
    plot(ax, cx, cy, 'k+', 'MarkerSize', 12, 'LineWidth', 1.5);

    set(ax, 'YDir','reverse');   % screen coordinates: y grows downward
    xlim(ax, [sx0, sx1]); ylim(ax, [sy0, sy1]);   % full display extent
    ylabel(ax,'y (px)');
    if t == numel(tasks_plot), xlabel(ax,'x (px)'); end
    title(ax, sprintf('%s   (n=%d fixations, %.1f%% within %d px)', ...
        strrep(tasks_plot(t),'_','-'), sum(mt), ...
        100*mean(Mc.dev(mt) <= FIX_TOL_PX), FIX_TOL_PX));
    box(ax,'on');

    if t == numel(tasks_plot)
        keepL = isgraphics(hleg);
        legend(ax, [hleg(keepL), ph1, ph2, phI], ...
            [cellfun(@(c) sprintf('%s', c), leg_lbl(keepL), 'UniformOutput', false), ...
             {sprintf('onset gate (%d px)', FIX_GATE_TOL_PX), ...
              sprintf('break window (%d px)', FIX_TOL_PX), ...
              sprintf('image frame (%dx%d px)', IMG_SIZE_PX, IMG_SIZE_PX)}], ...
            'Location','northeastoutside');
    end
end
title(tl, sprintf('sub%03d: fixations during stimulus viewing (full %g x %g display)', ...
    sub_list(1), sx1-sx0, sy1-sy0), 'FontWeight','bold');
exportgraphics(f1, fullfile(fig_dir, sprintf('sub%03d_fixation_overlay.png', sub_list(1))), 'Resolution', 200);

% per-task x condition counts, for the record
fprintf('\n  fixations by task x condition:\n');
for t = 1:numel(tasks_plot)
    fprintf('    %-7s ', strrep(tasks_plot(t),'_','-'));
    for k = 1:numel(cond_names)
        fprintf('%s=%-5d ', cond_names{k}, sum(Mc.task==tasks_plot(t) & Mc.condition==cond_names{k}));
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------------
%  FIGURE 2 -- deviation by condition
%  ------------------------------------------------------------------------
% trial-level summary: fixations within a trial are not independent, so the
% trial mean is the primary unit; fixation level is reported alongside.
[gt, tkey] = findgroups(Mc(:, {'task','trial_id','block','condition'}));
trial_dev = tkey;
trial_dev.mean_dev = splitapply(@mean, Mc.dev, gt);
trial_dev.max_dev  = splitapply(@max,  Mc.dev, gt);
trial_dev.n_fix    = splitapply(@numel, Mc.dev, gt);

f2 = figure('Visible','off','Position',[100 100 1050 420],'Color','w');

% (a) trial-level distributions
subplot(1,3,1); hold on;
for k = 1:numel(cond_names)
    v = trial_dev.mean_dev(trial_dev.condition == cond_names{k});
    jit = (rand(numel(v),1)-0.5)*0.28;
    scatter(k+jit, v, 16, cond_colors{k}, 'filled', 'MarkerFaceAlpha', 0.5);
    plot(k+[-0.3 0.3], [median(v) median(v)], 'k-', 'LineWidth', 2);
end
yline(FIX_GATE_TOL_PX, 'k-'); yline(FIX_TOL_PX, 'k--');
xticks(1:numel(cond_names)); xticklabels(cond_names); xlim([0.4 numel(cond_names)+0.6]);
ylabel('mean deviation per trial (px)'); title('trial level'); box on;

% (b) fixation-level ECDF
subplot(1,3,2); hold on;
for k = 1:numel(cond_names)
    v = sort(Mc.dev(Mc.condition == cond_names{k}));
    plot(v, (1:numel(v))/numel(v), 'Color', cond_colors{k}, 'LineWidth', 2);
end
xline(FIX_GATE_TOL_PX, 'k-'); xline(FIX_TOL_PX, 'k--');
xlabel('deviation (px)'); ylabel('cumulative proportion');
title('fixation level (ECDF)'); xlim([0 400]); box on;
legend(cond_names, 'Location','southeast');

% (c) proportion outside each radius, with Wilson CIs
subplot(1,3,3); hold on;
for k = 1:numel(cond_names)
    v = Mc.dev(Mc.condition == cond_names{k});
    for r = 1:2
        rad = FIX_GATE_TOL_PX*(r==1) + FIX_TOL_PX*(r==2);
        [pp, lo, hi] = wilson_ci(sum(v > rad), numel(v));
        xpos = k + (r-1.5)*0.22;
        bar(xpos, pp, 0.2, 'FaceColor', cond_colors{k}, ...
            'FaceAlpha', 1 - 0.45*(r-1), 'EdgeColor','k');
        plot([xpos xpos], [lo hi], 'k-', 'LineWidth', 1.2);
    end
end
xticks(1:numel(cond_names)); xticklabels(cond_names); xlim([0.4 numel(cond_names)+0.6]);
ylabel('proportion of fixations outside'); title('left bar: >100 px   right: >150 px'); box on;

exportgraphics(f2, fullfile(fig_dir, sprintf('sub%03d_gaze_deviation.png', sub_list(1))), 'Resolution', 200);

%% ------------------------------------------------------------------------
%  FIGURE 3 -- per-block drift check
%  ------------------------------------------------------------------------
f3 = figure('Visible','off','Position',[100 100 1000 420],'Color','w');
blocks = unique(Mc.block)';

subplot(1,2,1); hold on; axis equal;
plot(cx + FIX_TOL_PX*cos(th),      cy + FIX_TOL_PX*sin(th),      'k--');
plot(cx + FIX_GATE_TOL_PX*cos(th), cy + FIX_GATE_TOL_PX*sin(th), 'k-');
bcol = lines(numel(blocks));
for i = 1:numel(blocks)
    m = Mc.block == blocks(i);
    plot(median(Mc.x(m)), median(Mc.y(m)), 'o', 'MarkerSize', 11, ...
        'MarkerFaceColor', bcol(i,:), 'MarkerEdgeColor','k');
end
plot(cx, cy, 'k+', 'MarkerSize', 12, 'LineWidth', 1.5);
set(gca,'YDir','reverse'); xlim([cx-250 cx+250]); ylim([cy-250 cy+250]);
xlabel('x (px)'); ylabel('y (px)'); title('median fixation position per block'); box on;
legend([{'break window','onset gate'}, arrayfun(@(b) sprintf('block %d', b), blocks, 'UniformOutput', false)], ...
    'Location','eastoutside');

subplot(1,2,2); hold on;
for i = 1:numel(blocks)
    v = Mc.dev(Mc.block == blocks(i));
    jit = (rand(numel(v),1)-0.5)*0.28;
    scatter(blocks(i)+jit, v, 8, bcol(i,:), 'filled', 'MarkerFaceAlpha', 0.25);
    plot(blocks(i)+[-0.3 0.3], [median(v) median(v)], 'k-', 'LineWidth', 2);
end
yline(FIX_GATE_TOL_PX,'k-'); yline(FIX_TOL_PX,'k--');
xticks(blocks); xlabel('block'); ylabel('fixation deviation (px)');
title('deviation by block'); box on;

exportgraphics(f3, fullfile(fig_dir, sprintf('sub%03d_block_drift.png', sub_list(1))), 'Resolution', 200);

%% ------------------------------------------------------------------------
%  STATISTICS
%  ------------------------------------------------------------------------
diary(fullfile(res_dir, sprintf('sub%03d_gaze_check.txt', sub_list(1))));
fprintf('\n========================================================\n');
fprintf(' sub%03d -- restricted-viewing manipulation check\n', sub_list(1));
fprintf('========================================================\n');
fprintf('Screen centre (%.1f, %.1f); radii: gate %d px, break %d px\n\n', ...
    cx, cy, FIX_GATE_TOL_PX, FIX_TOL_PX);

% --- overall containment -------------------------------------------------
fprintf('--- Overall containment (all fixations, both tasks) ---\n');
[pp, lo, hi] = wilson_ci(sum(Mc.dev <= FIX_GATE_TOL_PX), height(Mc));
fprintf('  within %d px : %6.2f%% [%.2f, %.2f]  (%d/%d)\n', FIX_GATE_TOL_PX, ...
    100*pp, 100*lo, 100*hi, sum(Mc.dev <= FIX_GATE_TOL_PX), height(Mc));
[pp, lo, hi] = wilson_ci(sum(Mc.dev <= FIX_TOL_PX), height(Mc));
fprintf('  within %d px : %6.2f%% [%.2f, %.2f]  (%d/%d)\n', FIX_TOL_PX, ...
    100*pp, 100*lo, 100*hi, sum(Mc.dev <= FIX_TOL_PX), height(Mc));
fprintf('  median deviation %.1f px (IQR %.1f-%.1f)\n\n', median(Mc.dev), ...
    prctile(Mc.dev,25), prctile(Mc.dev,75));

% --- by condition --------------------------------------------------------
fprintf('--- Deviation by condition ---\n');
fprintf('%-10s %7s %8s %10s %10s %12s %12s\n', 'condition','nTrial','nFix', ...
    'medTrial','IQR','within100','within150');
for k = 1:numel(cond_names)
    mF = Mc.condition == cond_names{k};
    mT = trial_dev.condition == cond_names{k};
    v  = trial_dev.mean_dev(mT);
    fprintf('%-10s %7d %8d %10.1f %10s %11.1f%% %11.1f%%\n', cond_names{k}, ...
        sum(mT), sum(mF), median(v), ...
        sprintf('%.0f-%.0f', prctile(v,25), prctile(v,75)), ...
        100*mean(Mc.dev(mF) <= FIX_GATE_TOL_PX), ...
        100*mean(Mc.dev(mF) <= FIX_TOL_PX));
end

% Kruskal-Wallis on trial-level means, then pairwise rank-sum (Holm)
grp = trial_dev.condition;
p_kw = kruskalwallis(trial_dev.mean_dev, grp, 'off');
df_kw = numel(unique(grp)) - 1;
fprintf('\n  Kruskal-Wallis across conditions (trial level): chi2(%d) = %.2f, p = %.4f\n', ...
    df_kw, kw_chi2(trial_dev.mean_dev, grp), p_kw);

pairs = nchoosek(1:numel(cond_names), 2);
praw = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    a = trial_dev.mean_dev(trial_dev.condition == cond_names{pairs(i,1)});
    b = trial_dev.mean_dev(trial_dev.condition == cond_names{pairs(i,2)});
    praw(i) = ranksum(a, b);
end
padj = holm(praw);
fprintf('  pairwise rank-sum (Holm-adjusted):\n');
for i = 1:size(pairs,1)
    fprintf('    %-9s vs %-9s : p = %.4f  (adj %.4f)\n', ...
        cond_names{pairs(i,1)}, cond_names{pairs(i,2)}, praw(i), padj(i));
end

% --- by task -------------------------------------------------------------
fprintf('\n--- Deviation by task ---\n');
for t = ["1_back","2_back"]
    m = Mc.task == t;
    fprintf('  %-7s n=%6d  median %.1f px  within100 %.1f%%  within150 %.1f%%\n', ...
        t, sum(m), median(Mc.dev(m)), ...
        100*mean(Mc.dev(m) <= FIX_GATE_TOL_PX), 100*mean(Mc.dev(m) <= FIX_TOL_PX));
end
mt = trial_dev.mean_dev(ismember(trial_dev.task, "1_back"));
mt2 = trial_dev.mean_dev(ismember(trial_dev.task, "2_back"));
fprintf('  rank-sum 1-back vs 2-back (trial level): p = %.4f\n', ranksum(mt, mt2));

% --- per block: drift vs behaviour ---------------------------------------
fprintf('\n--- Per block: offset vs spread ---\n');
fprintf(['  OFFSET = distance of the block''s median fixation from screen centre\n' ...
         '           (a constant shift => calibration drift, not the subject).\n' ...
         '  SPREAD = median distance of fixations from the block''s OWN median\n' ...
         '           (inflated => the subject really was moving their eyes).\n\n']);
fprintf('%7s %8s %9s %9s %9s %9s %9s %11s %11s\n', 'block','nFix','dX','dY', ...
    'OFFSET','SPREAD','medDev','within100','within150');
blk_stats = table();
for i = 1:numel(blocks)
    m = Mc.block == blocks(i);
    mx = median(Mc.x(m)); my = median(Mc.y(m));
    offset = hypot(mx-cx, my-cy);
    spread = median(hypot(Mc.x(m)-mx, Mc.y(m)-my));
    fprintf('%7d %8d %9.1f %9.1f %9.1f %9.1f %9.1f %10.1f%% %10.1f%%\n', ...
        blocks(i), sum(m), mx-cx, my-cy, offset, spread, median(Mc.dev(m)), ...
        100*mean(Mc.dev(m) <= FIX_GATE_TOL_PX), 100*mean(Mc.dev(m) <= FIX_TOL_PX));
    blk_stats = [blk_stats; table(blocks(i), sum(m), mx-cx, my-cy, offset, spread, ...
        median(Mc.dev(m)), 'VariableNames', {'block','n_fix','dx','dy','offset','spread','med_dev'})];
end

% Contrast the suspect block against the rest on each component separately.
oth = ~ismember(blk_stats.block, 3);
if any(blk_stats.block == 3) && any(oth)
    b3 = blk_stats(blk_stats.block == 3, :);
    fprintf('\n  block 3 offset %.1f px vs %.1f px in other blocks (%.1fx)\n', ...
        b3.offset, mean(blk_stats.offset(oth)), b3.offset / mean(blk_stats.offset(oth)));
    fprintf('  block 3 spread %.1f px vs %.1f px in other blocks (%.1fx)\n', ...
        b3.spread, mean(blk_stats.spread(oth)), b3.spread / mean(blk_stats.spread(oth)));
    % is block 3's spread actually different, fixation-level?
    d3  = hypot(Mc.x(Mc.block==3)-median(Mc.x(Mc.block==3)), ...
                Mc.y(Mc.block==3)-median(Mc.y(Mc.block==3)));
    dOt = [];
    for i = 1:numel(blocks)
        if blocks(i) == 3, continue; end
        m = Mc.block == blocks(i);
        dOt = [dOt; hypot(Mc.x(m)-median(Mc.x(m)), Mc.y(m)-median(Mc.y(m)))];
    end
    fprintf('  rank-sum on spread (block 3 vs rest, fixation level): p = %.4g\n', ranksum(d3, dOt));
end

% --- sample validity: the mechanism behind the above --------------------
fprintf('\n--- Data loss per run (missing gaze samples) ---\n');
fprintf('%7s %9s %12s %12s %13s\n', 'block','task','n_samples','n_missing','pct_missing');
for i = 1:height(qual)
    fprintf('%7d %9s %12d %12d %12.2f%%\n', qual.block(i), qual.task(i), ...
        qual.n_samples(i), qual.n_missing(i), qual.pct_missing(i));
end
fprintf(['\n  A run with heavy data loss cannot pass the onset gate (no stable\n' ...
         '  valid fixation to confirm) yet also cannot trigger fix_broken,\n' ...
         '  because invalid samples are explicitly excluded from break\n' ...
         '  detection. High gate_timeout + 0%% fix_broken is therefore the\n' ...
         '  signature of a tracker problem, not of the subject looking away.\n']);

diary off;

res = struct('merged', Mc, 'trial_dev', trial_dev, 'blk_stats', blk_stats, ...
             'quality', qual, 'screen_rect', screen_rect, ...
             'centre', [cx cy], 'radii', [FIX_GATE_TOL_PX FIX_TOL_PX]);
save(fullfile(res_dir, sprintf('sub%03d_gaze_check.mat', sub_list(1))), 'res');
fprintf('\nSaved figures to %s and results to %s\n', fig_dir, res_dir);

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------
function [p, lo, hi] = wilson_ci(k, n)
% Wilson score interval for a binomial proportion (95%).
z = 1.959963984540054;
if n == 0, p = NaN; lo = NaN; hi = NaN; return; end
p = k / n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr - hw); hi = min(1, ctr + hw);
end

function chi2 = kw_chi2(x, g)
% Kruskal-Wallis H statistic, recomputed so it can be printed alongside p.
[~, ~, gi] = unique(g);
n = numel(x); r = tiedrank(x);
k = max(gi); s = 0;
for j = 1:k
    rj = r(gi == j);
    s = s + numel(rj) * (mean(rj) - (n+1)/2)^2;
end
chi2 = 12 / (n*(n+1)) * s;
end

function padj = holm(p)
% Holm-Bonferroni step-down adjustment.
[ps, idx] = sort(p(:));
m = numel(ps);
adj = ps .* (m:-1:1)';
adj = cummax(adj);
padj = nan(size(ps));
padj(idx) = min(adj, 1);
end
