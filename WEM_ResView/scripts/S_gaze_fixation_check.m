%==========================================================================
%   WMEM_ResView -- restricted-viewing manipulation check (eye data)
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
%==========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
%  Config
%  ------------------------------------------------------------------------
sub_list = [601, 602, 603];
tasks = [1, 2]; task_lbls = {'1_back', '2_back'};
FS = struct('tick',16,'lab',18,'ttl',18,'anno',13);   % shared type scale (matches NeurWEM_plt)
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

% Stimulus presentation window. The run scripts hold the image on screen for
% p.timing.image_dur (A_subject_setup.m) after STIM_ONSET, with no offset
% message, so the ITI is bounded by this duration rather than a marker.
IMG_DUR_MS = 1500;       % stimulus presentation duration (image_dur = 1.5 s)

% Stimulus extent. The run scripts call Screen('DrawTexture', ..., [], [], 0)
% with an empty destination rect, so each image is drawn at its native size,
% centred. Verified against a real stimulus further down.
IMG_SIZE_PX = 400;

% Condition colours (c_comp = NYU violet, consistent with NeurWEM_plt group figs)
c_comp = [87 6 140]/255;
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

% Restrict to fixations that began while the image was on screen. The parser
% keeps every fixation from STIM_ONSET to the next TRIALID, so without this
% two kinds of non-viewing fixations leak in: a fixation-cross fixation that
% ends just after onset (t_rel < 0) and any fixation in the post-image tail
% (t_rel >= IMG_DUR_MS). Both belong to the ITI, not to stimulus viewing.
n_before = height(M);
M = M(M.t_rel >= 0 & M.t_rel < IMG_DUR_MS, :);   % NaN t_rel (no onset) also dropped
fprintf('fixations during image presentation: %d of %d kept (%d dropped as pre-stim / ITI)\n', ...
    height(M), n_before, n_before - height(M));

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
%  Per-subject metrics  (subject is the unit for the group check)
%  ------------------------------------------------------------------------
usub = unique(Mc.subj_id)';  nS = numel(usub);
blocks = unique(Mc.block)';  nB = numel(blocks);
tasks_plot = ["1_back", "2_back"];

med_cond  = nan(nS, numel(cond_names));   % median deviation (px) by condition
w100_cond = nan(nS, numel(cond_names));   % fraction within 100 px by condition
w150_cond = nan(nS, numel(cond_names));   % fraction within 150 px by condition
med_task  = nan(nS, numel(tasks_plot));   % median deviation by task
med_block = nan(nS, nB);                  % median deviation by block (drift)
w150_all  = nan(nS, 1);  w100_all = nan(nS, 1);  med_all = nan(nS, 1);
for si = 1:nS
    ms = Mc.subj_id == usub(si);
    med_all(si)  = median(Mc.dev(ms));
    w100_all(si) = mean(Mc.dev(ms) <= FIX_GATE_TOL_PX);
    w150_all(si) = mean(Mc.dev(ms) <= FIX_TOL_PX);
    for k = 1:numel(cond_names)
        mk = ms & Mc.condition == cond_names{k};
        if ~any(mk), continue; end
        med_cond(si,k)  = median(Mc.dev(mk));
        w100_cond(si,k) = mean(Mc.dev(mk) <= FIX_GATE_TOL_PX);
        w150_cond(si,k) = mean(Mc.dev(mk) <= FIX_TOL_PX);
    end
    for t = 1:numel(tasks_plot)
        mt = ms & Mc.task == tasks_plot(t);
        if any(mt), med_task(si,t) = median(Mc.dev(mt)); end
    end
    for b = 1:nB
        mb = ms & Mc.block == blocks(b);
        if any(mb), med_block(si,b) = median(Mc.dev(mb)); end
    end
end
min_n_stat = 3;   % minimum N of subjects before paired tests / brackets are drawn

%% ------------------------------------------------------------------------
%  FIGURE 1 -- fixation overlay, 1-back and 2-back side by side
%  ------------------------------------------------------------------------
% Note: 'novel' exists only in the 2-back, so the 1-back panel legitimately
% shows two conditions and the 2-back three. Fixations are pooled across the
% N subjects here; per-subject summaries carry the group statistics below.
th = linspace(0, 2*pi, 400);

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
        s = scatter(ax, Mc.x(m), Mc.y(m), 12, cond_colors{k}, 'filled', ...
            'MarkerFaceAlpha', 0.30, 'MarkerEdgeColor','none');
        if t == numel(tasks_plot) || ~isgraphics(hleg(k))
            hleg(k) = s; leg_lbl{k} = cond_names{k};
        end
    end
    plot(ax, cx, cy, 'k+', 'MarkerSize', 12, 'LineWidth', 1.5);

    set(ax, 'YDir','reverse', 'FontSize', FS.tick);   % screen coordinates: y grows downward
    xlim(ax, [sx0, sx1]); ylim(ax, [sy0, sy1]);   % full display extent
    ylabel(ax,'y (px)', 'FontSize', FS.lab);
    if t == numel(tasks_plot), xlabel(ax,'x (px)', 'FontSize', FS.lab); end
    title(ax, sprintf('%s   (%d fixations, %.1f%% within %d px)', ...
        strrep(tasks_plot(t),'_','-'), sum(mt), ...
        100*mean(Mc.dev(mt) <= FIX_TOL_PX), FIX_TOL_PX), 'FontSize', FS.ttl);
    box(ax,'on');

    if t == numel(tasks_plot)
        keepL = isgraphics(hleg);
        legend(ax, [hleg(keepL), ph1, ph2, phI], ...
            [cellfun(@(c) sprintf('%s', c), leg_lbl(keepL), 'UniformOutput', false), ...
             {sprintf('onset gate (%d px)', FIX_GATE_TOL_PX), ...
              sprintf('break window (%d px)', FIX_TOL_PX), ...
              sprintf('image frame (%dx%d px)', IMG_SIZE_PX, IMG_SIZE_PX)}], ...
            'Location','northeastoutside', 'FontSize', FS.anno, 'Box','off');
    end
end
title(tl, sprintf('N = %d subjects: fixations during stimulus viewing (full %g x %g display)', ...
    nS, sx1-sx0, sy1-sy0), 'FontWeight','bold', 'FontSize', FS.ttl);
exportgraphics(f1, fullfile(fig_dir, 'group_fixation_overlay.png'), 'Resolution', 200);

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
%  FIGURE 2 -- deviation by condition (group: one dot per subject)
%  ------------------------------------------------------------------------
cpairs = {[1 2],[2 3],[1 3]};
f2 = figure('Visible','off','Position',[80 80 1400 480],'Color','w');
subplot(1,3,1);
paired_plot(med_cond, cond_names, cond_colors, 'median deviation (px)', ...
    'gaze deviation by condition', min_n_stat, cpairs, FS, true);
yline(FIX_GATE_TOL_PX,'k-'); yline(FIX_TOL_PX,'k--'); nice_yticks(4);
subplot(1,3,2);
paired_plot(w100_cond, cond_names, cond_colors, sprintf('fraction within %d px', FIX_GATE_TOL_PX), ...
    'containment (onset gate)', min_n_stat, cpairs, FS, true);
ylim([0 1.05]); nice_yticks(4);
subplot(1,3,3);
paired_plot(w150_cond, cond_names, cond_colors, sprintf('fraction within %d px', FIX_TOL_PX), ...
    'containment (break window)', min_n_stat, cpairs, FS, true);
ylim([0 1.05]); nice_yticks(4);
exportgraphics(f2, fullfile(fig_dir, 'group_gaze_deviation.png'), 'Resolution', 200);

%% ------------------------------------------------------------------------
%  FIGURE 3 -- deviation by task and by block (group: one dot per subject)
%  ------------------------------------------------------------------------
c_1b = [97 125 184]/255; c_2b = c_comp;             % task colours
bc = lines(nB); block_cols = arrayfun(@(i) bc(i,:), 1:nB, 'UniformOutput', false);
block_lbls = arrayfun(@(b) sprintf('block %d', b), blocks, 'UniformOutput', false);

f3 = figure('Visible','off','Position',[80 80 1250 480],'Color','w');
subplot(1,2,1);
paired_plot(med_task, {'1-back','2-back'}, {c_1b, c_2b}, 'median deviation (px)', ...
    'gaze deviation by task', min_n_stat, {[1 2]}, FS, true);
yline(FIX_GATE_TOL_PX,'k-'); yline(FIX_TOL_PX,'k--'); nice_yticks(4);
subplot(1,2,2);
% blocks are ordinal, so the paired_plot trend line is a genuine drift check
paired_plot(med_block, block_lbls, block_cols, 'median deviation (px)', ...
    'gaze deviation by block (drift)', min_n_stat, {}, FS, false);
yline(FIX_GATE_TOL_PX,'k-'); yline(FIX_TOL_PX,'k--'); nice_yticks(4);
exportgraphics(f3, fullfile(fig_dir, 'group_block_drift.png'), 'Resolution', 200);

%% ------------------------------------------------------------------------
%  STATISTICS  (group: the subject is the unit of analysis)
%  ------------------------------------------------------------------------
diary(fullfile(res_dir, 'group_gaze_check.txt'));
fprintf('\n========================================================\n');
fprintf(' GROUP -- restricted-viewing manipulation check (N = %d)\n', nS);
fprintf(' subjects: %s\n', num2str(usub));
fprintf('========================================================\n');
fprintf('Screen centre (%.1f, %.1f); radii: gate %d px, break %d px\n\n', ...
    cx, cy, FIX_GATE_TOL_PX, FIX_TOL_PX);

% --- overall containment (per subject, summarised across subjects) -------
fprintf('--- Overall containment (per subject, both tasks) ---\n');
fprintf('  within %3d px : %5.1f%% +/- %.1f  [range %.1f-%.1f]\n', FIX_GATE_TOL_PX, ...
    100*mean(w100_all), 100*std(w100_all), 100*min(w100_all), 100*max(w100_all));
fprintf('  within %3d px : %5.1f%% +/- %.1f  [range %.1f-%.1f]\n', FIX_TOL_PX, ...
    100*mean(w150_all), 100*std(w150_all), 100*min(w150_all), 100*max(w150_all));
fprintf('  median deviation : %.1f px +/- %.1f  [range %.1f-%.1f]\n', ...
    mean(med_all), std(med_all), min(med_all), max(med_all));
fprintf('  per subject:\n');
for si = 1:nS
    fprintf('    sub%03d  medDev %5.1f px  within100 %5.1f%%  within150 %5.1f%%\n', ...
        usub(si), med_all(si), 100*w100_all(si), 100*w150_all(si));
end

% --- by condition (per-subject means; Friedman + pairwise signed-rank) ---
fprintf('\n--- Deviation by condition (per subject) ---\n');
fprintf('%-10s %16s %12s %12s\n', 'condition','medDev(px)','within100','within150');
for k = 1:numel(cond_names)
    fprintf('%-10s %8.1f +/-%-5.1f %8.1f%% %11.1f%%\n', cond_names{k}, ...
        mean(med_cond(:,k),'omitnan'), std(med_cond(:,k),'omitnan'), ...
        100*mean(w100_cond(:,k),'omitnan'), 100*mean(w150_cond(:,k),'omitnan'));
end
[chi2_c, df_c, p_c] = friedman_safe(med_cond);
fprintf('\n  Friedman across conditions (subject level): chi2(%d) = %.2f, p = %.4f\n', ...
    df_c, chi2_c, p_c);
pairs = nchoosek(1:numel(cond_names), 2);
praw = nan(size(pairs,1),1);
for i = 1:size(pairs,1)
    a = med_cond(:,pairs(i,1)); b = med_cond(:,pairs(i,2));
    ok = ~isnan(a) & ~isnan(b);
    if sum(ok) >= min_n_stat, praw(i) = signrank(a(ok), b(ok)); end
end
padj = holm(praw);
fprintf('  pairwise signed-rank (Holm-adjusted):\n');
for i = 1:size(pairs,1)
    fprintf('    %-9s vs %-9s : p = %.4f  (adj %.4f)\n', ...
        cond_names{pairs(i,1)}, cond_names{pairs(i,2)}, praw(i), padj(i));
end

% --- by task (paired across subjects) ------------------------------------
fprintf('\n--- Deviation by task (per subject) ---\n');
for t = 1:numel(tasks_plot)
    fprintf('  %-7s median %.1f px +/- %.1f\n', strrep(tasks_plot(t),'_','-'), ...
        mean(med_task(:,t),'omitnan'), std(med_task(:,t),'omitnan'));
end
ok = ~isnan(med_task(:,1)) & ~isnan(med_task(:,2));
if sum(ok) >= min_n_stat
    fprintf('  signed-rank 1-back vs 2-back: p = %.4f\n', signrank(med_task(ok,1), med_task(ok,2)));
else
    fprintf('  signed-rank 1-back vs 2-back: n=%d too few\n', sum(ok));
end

% --- by block: drift across the group ------------------------------------
fprintf('\n--- Deviation by block (per subject; drift check) ---\n');
fprintf('%7s %16s %12s %12s\n', 'block','medDev(px)','within100','within150');
for b = 1:nB
    mb = Mc.block == blocks(b);
    fprintf('%7d %8.1f +/-%-5.1f %8.1f%% %11.1f%%\n', blocks(b), ...
        mean(med_block(:,b),'omitnan'), std(med_block(:,b),'omitnan'), ...
        100*mean(Mc.dev(mb) <= FIX_GATE_TOL_PX), 100*mean(Mc.dev(mb) <= FIX_TOL_PX));
end
[chi2_b, df_b, p_b] = friedman_safe(med_block);
fprintf('\n  Friedman across blocks (subject level): chi2(%d) = %.2f, p = %.4f\n', ...
    df_b, chi2_b, p_b);

% --- data loss: the mechanism behind gate timeouts -----------------------
fprintf('\n--- Data loss (missing gaze samples) ---\n');
fprintf('%8s %13s %12s\n', 'subject','mean pct','max pct');
for si = 1:nS
    q = qual(qual.subj_id == usub(si), :);
    fprintf('  sub%03d %11.2f%% %11.2f%%\n', usub(si), mean(q.pct_missing), max(q.pct_missing));
end
fprintf('  overall: %.2f%% of samples missing\n', 100*sum(qual.n_missing)/max(sum(qual.n_samples),1));
fprintf(['\n  A run with heavy data loss cannot pass the onset gate yet also\n' ...
         '  cannot trigger fix_broken (invalid samples are excluded from break\n' ...
         '  detection). High gate_timeout + 0%% fix_broken is the signature of a\n' ...
         '  tracker problem, not of the subject looking away.\n']);

diary off;

res = struct('merged', Mc, 'usub', usub, ...
             'med_cond', med_cond, 'w100_cond', w100_cond, 'w150_cond', w150_cond, ...
             'med_task', med_task, 'med_block', med_block, ...
             'med_all', med_all, 'w100_all', w100_all, 'w150_all', w150_all, ...
             'quality', qual, 'screen_rect', screen_rect, ...
             'centre', [cx cy], 'radii', [FIX_GATE_TOL_PX FIX_TOL_PX]);
save(fullfile(res_dir, 'group_gaze_check.mat'), 'res');
fprintf('\nSaved figures to %s and results to %s\n', fig_dir, res_dir);

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------
function [chi2, df, p] = friedman_safe(X)
% Friedman test across the columns of X (rows = subjects), on complete cases.
% Returns the chi-square statistic, its df, and p so all three can be printed.
    X = X(all(~isnan(X),2), :);
    if size(X,1) < 2 || size(X,2) < 2, chi2 = NaN; df = size(X,2)-1; p = NaN; return; end
    [p, tbl] = friedman(X, 1, 'off');
    chi2 = tbl{2,5}; df = tbl{2,3};
end

function paired_plot(M, lvllbl, cols, ylbl, ttl, min_n, pairs, FS, no_trend)
% Per-subject paired plot (matches NeurWEM_plt/S_group_mstback.m): jittered
% coloured dots per subject, grey lines linking each subject across levels, a
% thin box/whisker per level, signed-rank brackets, and (unless suppressed) a
% dashed trend across ordered levels. M is [nS x nL].
    if nargin < 8 || isempty(FS), FS = struct('tick',16,'lab',18,'ttl',18,'anno',13); end
    if nargin < 9 || isempty(no_trend), no_trend = false; end
    [nS, nL] = size(M);
    if nargin < 7 || isempty(pairs)
        pairs = arrayfun(@(k) [k k+1], 1:nL-1, 'UniformOutput', false);
    end
    hold on; jw = 0.07; bw = 0.34;
    X = (1:nL) + (rand(nS,nL)-0.5)*2*jw;

    for s = 1:nS
        plot(X(s,:), M(s,:), '-', 'Color', [0.55 0.55 0.55 0.45], 'LineWidth', 0.6);
    end
    for L = 1:nL
        d = M(:,L); d = d(~isnan(d));
        if numel(d) >= 2
            q = quantile(d,[.25 .5 .75]); iqrv = q(3)-q(1);
            lo_w = max(min(d), q(1)-1.5*iqrv); hi_w = min(max(d), q(3)+1.5*iqrv);
            plot([L L],[lo_w q(1)],'k-','LineWidth',0.9);
            plot([L L],[q(3) hi_w],'k-','LineWidth',0.9);
            plot(L+[-.06 .06],[lo_w lo_w],'k-','LineWidth',0.9);
            plot(L+[-.06 .06],[hi_w hi_w],'k-','LineWidth',0.9);
            rectangle('Position',[L-bw/2, q(1), bw, max(iqrv,eps)], ...
                'EdgeColor','k','LineWidth',1.1,'FaceColor','none');
            plot(L+[-bw/2 bw/2],[q(2) q(2)],'k-','LineWidth',1.8);
        end
        scatter(X(:,L), M(:,L), 34, cols{L}, 'filled', 'MarkerFaceAlpha',0.85, ...
            'MarkerEdgeColor',[.2 .2 .2], 'LineWidth',0.3);
    end

    all_v = M(~isnan(M));
    yr = range(all_v); if isempty(yr)||yr==0, yr = 1; end

    base = max(all_v) + 0.10*yr; step = 0.10*yr; k = 0;
    if nS >= min_n
        for pp = 1:numel(pairs)
            ij = pairs{pp};
            ok = ~isnan(M(:,ij(1))) & ~isnan(M(:,ij(2)));
            if sum(ok) < min_n, continue; end
            pv = signrank(M(ok,ij(1)), M(ok,ij(2))); s = stars(pv); if isempty(s), s='ns'; end
            y = base + k*step; k = k + 1;
            plot([ij(1) ij(1) ij(2) ij(2)], [y-0.02*yr y y y-0.02*yr], 'k-','LineWidth',0.9);
            text(mean(ij), y, s, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FS.tick);
        end
    end

    if ~no_trend && nL >= 3 && nS >= min_n
        xa = repmat(1:nL, nS, 1); xa = xa(:); ya = M(:);
        ok = ~isnan(ya); xa = xa(ok); ya = ya(ok);
        if numel(unique(xa)) >= 2
            [r, pv] = corr(xa, ya); cf = polyfit(xa, ya, 1);
            xx = [0.7 nL+0.3];
            plot(xx, polyval(cf, xx), 'k--', 'LineWidth', 1.6);
            text(0.7, min(all_v), sprintf('R^2 = %.2f, p = %.2g\ny = %.2f + %.2f x', ...
                r^2, pv, cf(2), cf(1)), 'FontSize', FS.anno, 'VerticalAlignment','bottom');
        end
    end

    set(gca,'XTick',1:nL,'XTickLabel',lvllbl,'FontSize',FS.tick);
    xlim([0.4 nL+0.6]); ylabel(ylbl,'FontSize',FS.lab);
    if ~isempty(ttl), title(ttl,'FontSize',FS.ttl); end
    if nS >= min_n && k > 0
        ylim([min(all_v)-0.06*yr, base + k*step + 0.05*yr]);
    end
    box off; hold off;
end

function s = stars(p), s = repmat('*', 1, (p<0.05)+(p<0.01)+(p<0.001)); end

function nice_yticks(target_n)
% reduce the y-axis to ~target_n nicely-rounded ticks
    yl = ylim; rng = yl(2)-yl(1);
    if ~isfinite(rng) || rng <= 0, return; end
    raw = rng/max(target_n,1);
    mag = 10^floor(log10(raw));
    steps = [1 2 2.5 5 10]*mag;
    step = steps(find(steps >= raw, 1, 'first'));
    if isempty(step), step = steps(end); end
    t = ceil(yl(1)/step)*step : step : floor(yl(2)/step)*step;
    if numel(t) >= 2, set(gca,'YTick',t); end
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
