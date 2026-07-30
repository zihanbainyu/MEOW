clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEM_ResView -- group behavioural results (subject = unit of analysis)
% Aesthetic matches NeurWEM_plt/scripts/S_group_mstback.m (paired_plot dots,
% fixed-size panels, draw_matrix_se confusions). Conditions: compared /
% isolated / novel; post-task test is old/new recognition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---- config ----
base_dir = '..';
data_dir = fullfile(base_dir, 'data');
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');
if ~exist(res_dir,'dir'), mkdir(res_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
min_rt   = 0.150;
min_n_stat = 5;     % minimum N before paired tests / trend are drawn
FS = struct('tick',20,'lab',20,'ttl',20,'anno',20);   % enlarged, consistent type scale

% Fixed panel geometry (pixels): each paired-plot panel is nCond*col wide by h
% tall, so two-group panels are EXACTLY the same size in every figure.
GEO = struct('col',150, 'h',380, 'ml',110, 'mr',40, 'mb',72, 'mt',58, 'hg',105, 'vg',120);

c_same = [97 125 184]/255; c_sim = [255 191 205]/255; c_new = [219 219 219]/255;
c_comp = [87 6 140]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255;   % c_comp = NYU violet
cnd_nm   = {'compared','isolated','novel'};
cnd_cols = {c_comp, c_iso, c_nov};
rec_nm   = {'compared','isolated'};        % only studied conditions are "old"
rec_cols = {c_comp, c_iso};

%% ---- find subjects with an n-back concat ----
dd = dir(fullfile(data_dir, 'sub*'));
subs = [];
for i = 1:numel(dd)
    id = sscanf(dd(i).name, 'sub%d');
    f  = fullfile(data_dir, dd(i).name, sprintf('sub%03d_concat.mat', id));
    if ~isfile(f), continue; end
    L = load(f, 'final_data_output');
    if isfield(L,'final_data_output') && ...
       all(isfield(L.final_data_output, {'results_1_back_all','results_2_back_all'}))
        subs(end+1) = id; %#ok<SAGROW>
    end
end
subs = sort(subs);
nS = numel(subs);
assert(nS >= 1, 'No subjects with an n-back concat found in %s', data_dir);
fprintf('Group WEM_ResView: %d subjects [%s]\n', nS, num2str(subs));

%% ---- per-subject metrics ----
one_acc  = nan(nS,3);   % 1-back accuracy: same / similar / new
one_rt   = nan(nS,2);   % 1-back median RT (correct): same / similar
goal_acc = nan(nS,3);   % 2-back accuracy by goal (pooled): AB / AA / AN
di      = nan(nS,3);    % 2-back LDI (similar discrimination): compared/isolated/novel
dpr      = nan(nS,3);   % 2-back d'  (same detection)       : compared/isolated/novel
rt_lure  = nan(nS,3);   % 2-back median RT (correct AB): by condition
rt_targ  = nan(nS,3);   % 2-back median RT (correct AA): by condition
conf1_all = nan(3,3,nS);                       % 1-back confusion
conf2_all = {nan(3,3,nS), nan(3,3,nS), nan(3,3,nS)};   % 2-back confusion: {comp, iso, nov}
% post-task recognition (old / new by condition)
has_rec  = false(nS,1);
rec_d    = nan(nS,2);   % recognition d': compared / isolated
rec_hit  = nan(nS,2);   % hit rate p(old|old): compared / isolated
rec_far  = nan(nS,1);   % false-alarm rate (pooled foils)
for si = 1:nS
    f = fullfile(data_dir, sprintf('sub%03d', subs(si)), sprintf('sub%03d_concat.mat', subs(si)));
    L = load(f, 'final_data_output'); fdo = L.final_data_output;
    m = nback_metrics(fdo, min_rt);
    one_acc(si,:)  = m.one_acc;
    one_rt(si,:)   = m.one_rt;
    goal_acc(si,:) = m.goal_acc;
    di(si,:)      = m.ldi;
    dpr(si,:)      = m.dpr;
    rt_lure(si,:)  = m.rt_lure;
    rt_targ(si,:)  = m.rt_targ;
    conf1_all(:,:,si)    = m.conf1;
    conf2_all{1}(:,:,si) = m.conf2{1};
    conf2_all{2}(:,:,si) = m.conf2{2};
    conf2_all{3}(:,:,si) = m.conf2{3};
    if isfield(fdo, 'results_recognition')
        rm = rec_metrics(fdo, min_rt, rec_nm); has_rec(si) = true;
        rec_d(si,:) = rm.d; rec_hit(si,:) = rm.hit; rec_far(si) = rm.far;
    end
end
nREC = sum(has_rec);

% grand-mean confusion matrices (+/- SE across subjects)
[C1,  SE1]  = grand_conf(conf1_all);
[C2c, SE2c] = grand_conf(conf2_all{1});
[C2i, SE2i] = grand_conf(conf2_all{2});
[C2n, SE2n] = grand_conf(conf2_all{3});

%% ---- report ----
diary_file = fullfile(res_dir, 'group_behav_report.txt');
if exist(diary_file,'file'), delete(diary_file); end
diary(diary_file);
fprintf('================================================================\n');
fprintf(' GROUP WEM_ResView  --  N = %d   (%s)\n', nS, datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('================================================================\n');
grp_line('1-back same',    one_acc(:,1));
grp_line('1-back similar', one_acc(:,2));
grp_line('1-back new',     one_acc(:,3));
for c = 1:3, grp_line(sprintf('DI %s', cnd_nm{c}), di(:,c)); end
for c = 1:3, grp_line(sprintf('d''  %s', cnd_nm{c}), dpr(:,c)); end
fprintf('\n*** PRIMARY: 2-back similar discrimination (DI), across conditions ***\n');
paired_line('DI compared vs isolated', di(:,1), di(:,2), min_n_stat);
paired_line('DI compared vs novel   ', di(:,1), di(:,3), min_n_stat);
paired_line('DI isolated vs novel   ', di(:,2), di(:,3), min_n_stat);
fprintf('\n-- secondary: same detection (d'') --\n');
paired_line('d''  compared vs isolated', dpr(:,1), dpr(:,2), min_n_stat);
paired_line('d''  compared vs novel   ', dpr(:,1), dpr(:,3), min_n_stat);
paired_line('d''  isolated vs novel   ', dpr(:,2), dpr(:,3), min_n_stat);

fprintf('\n-- POST-TASK RECOGNITION (old/new; N with rec = %d) --\n', nREC);
if nREC >= 1
    grp_line('rec d'' compared',  rec_d(:,1));
    grp_line('rec d'' isolated',  rec_d(:,2));
    grp_line('rec hit compared',  rec_hit(:,1));
    grp_line('rec hit isolated',  rec_hit(:,2));
    grp_line('rec false-alarm',   rec_far);
    paired_line('rec d'' compared vs isolated', rec_d(:,1), rec_d(:,2), min_n_stat);
else
    fprintf('  (no subject has results_recognition)\n');
end
diary off;

%% ---- figures ----
rlbl = {'exp. same','exp. similar','exp. new'};    % confusion rows (presented)
clbl = {'resp. same','resp. similar','resp. new'};  % confusion cols (response)
p3 = {[1 2],[2 3],[1 3]};

% FIGURE 1 -- 1-back: accuracy (3) + RT (2)
[f, L] = panel_grid(GEO, [3 2]);   set(f,'Name','Figure 1: 1-back');
mk_axes(f, L(1,:));
paired_plot(one_acc, {'same','similar','new'}, {c_same,c_sim,c_new}, ...
    'accuracy', '1-back accuracy', min_n_stat, p3, FS);
ylim([0 1.08]); yline(1/3,'k:','chance'); nice_yticks(4);
mk_axes(f, L(2,:));
paired_plot(one_rt, {'same','similar'}, {c_same,c_sim}, ...
    'RT (s)', '1-back RT (correct)', min_n_stat, {[1 2]}, FS);
nice_yticks(4);
save_fig(f, fig_dir, 'group_behav_fig1_1back');

% FIGURE 2 -- 1-back confusion matrix (group mean +/- SE)
f = figure('color','w','Position',[80 80 560 540],'Name','Figure 2: 1-back confusion');
draw_matrix_se(C1, SE1, {c_same,c_sim,c_new}, rlbl, clbl, FS);
title(sprintf('1-back confusions  (N = %d)', nS), 'FontSize', FS.ttl);
save_fig(f, fig_dir, 'group_behav_fig2_1back_confusion');

% FIGURE 3 -- 2-back indices: DI, d', and their RTs (compared/isolated/novel)
[f, L] = panel_grid(GEO, [3 3; 3 3]);   set(f,'Name','Figure 3: 2-back');
mk_axes(f, L(1,:));
paired_plot(di, cnd_nm, cnd_cols, 'DI', ...
    'similar discrimination index (DI)', min_n_stat, p3, FS);
yline(0,'k-'); nice_yticks(4);
mk_axes(f, L(2,:));
paired_plot(dpr, cnd_nm, cnd_cols, 'd''', ...
    'same detection (d'')', min_n_stat, p3, FS);
yline(0,'k-'); nice_yticks(4);
mk_axes(f, L(3,:));
paired_plot(rt_lure, cnd_nm, cnd_cols, 'RT (s)', ...
    'RT: similar discrimination', min_n_stat, p3, FS);
nice_yticks(4);
mk_axes(f, L(4,:));
paired_plot(rt_targ, cnd_nm, cnd_cols, 'RT (s)', ...
    'RT: same detection', min_n_stat, p3, FS);
nice_yticks(4);
save_fig(f, fig_dir, 'group_behav_fig3_2back');

% FIGURE 4 -- 2-back confusion matrices (compared / isolated / novel)
f = figure('color','w','Position',[40 40 1620 560],'Name','Figure 4: 2-back confusion');
Cs = {C2c, C2i, C2n}; SEs = {SE2c, SE2i, SE2n};
for c = 1:3
    subplot(1,3,c);
    draw_matrix_se(Cs{c}, SEs{c}, {c_same,c_sim,c_new}, rlbl, clbl, FS);
    title(sprintf('2-back: %s', cnd_nm{c}), 'FontSize', FS.ttl);
end
save_fig(f, fig_dir, 'group_behav_fig4_2back_confusion');

% FIGURE 5 -- post-task recognition (d' and hit rate by condition)
if nREC >= 1
    [f, L] = panel_grid(GEO, [2 2]);   set(f,'Name','Figure 5: recognition');
    mk_axes(f, L(1,:));
    paired_plot(rec_d, rec_nm, rec_cols, 'd''', ...
        'recognition d''', min_n_stat, {[1 2]}, FS);
    yline(0,'k-'); nice_yticks(4);
    mk_axes(f, L(2,:));
    paired_plot(rec_hit, rec_nm, rec_cols, 'hit rate', ...
        'recognition hit rate', min_n_stat, {[1 2]}, FS);
    ylim([0 1.05]); yline(mean(rec_far,'omitnan'),'r--','FA'); nice_yticks(4);
    save_fig(f, fig_dir, 'group_behav_fig5_recognition');
end

fprintf('\nsaved report + figures to %s , %s\n', res_dir, fig_dir);

%% ==================== local functions ====================
function m = nback_metrics(fdo, min_rt)
    r1 = recode(fdo.results_1_back_all);
    r2 = recode(fdo.results_2_back_all);
    % 1-back accuracy: same / similar / new  (isolated items count as "new")
    i_sam = r1.condition=="repeat"   & strcmp(r1.corr_resp,'j');
    i_sim = r1.condition=="compared" & r1.identity=="B";
    i_new = (r1.condition=="compared" & r1.identity=="A") | ...
             r1.condition=="isolated" | ...
            (r1.condition=="repeat"   & strcmp(r1.corr_resp,'none'));
    m.one_acc = [pmean(r1.correct,i_sam), pmean(r1.correct,i_sim), pmean(r1.correct,i_new)];
    % 1-back RT on correct trials (median): same / similar
    v1 = r1.rt > min_rt;
    m.one_rt = [median(r1.rt(i_sam & r1.correct & v1),'omitnan'), ...
                median(r1.rt(i_sim & r1.correct & v1),'omitnan')];
    m.conf1 = conf_mat(r1.resp_key, {i_sam, i_sim, i_new}, {'j','k','none'});
    % 2-back goals
    real2 = ~contains(r2.goal, "JUNK");
    pan = false(height(r2),1);
    for i = 1:height(r2)-2
        if strcmp(r2.goal(i),'A-N') && r2.block(i+2)==r2.block(i), pan(i+2)=true; end
    end
    aa = real2 & strcmp(r2.goal,'A-A') & strcmp(r2.corr_resp,'j');
    ab = real2 & strcmp(r2.goal,'A-B') & strcmp(r2.corr_resp,'k');
    an = real2 & pan & strcmp(r2.corr_resp,'none');
    m.goal_acc = [pmean(r2.correct,ab), pmean(r2.correct,aa), pmean(r2.correct,an)];
    % LDI, d', RT, and confusion by condition (compared / isolated / novel)
    v2 = r2.rt > min_rt;
    conds = {'compared','isolated','novel'};
    m.conf2 = {nan(3,3), nan(3,3), nan(3,3)};
    for c = 1:3
        cm = strcmp(r2.condition, conds{c});
        m.ldi(c) = pkey(r2,ab&cm,'k') - pkey(r2,an&cm,'k');
        nh = sum(aa&cm); nf = sum(an&cm);
        m.dpr(c) = zc(pkey(r2,aa&cm,'j'),nh) - zc(pkey(r2,an&cm,'j'),nf);
        m.rt_lure(c) = median(r2.rt(ab&cm & r2.correct & v2),'omitnan');   % similar discrimination
        m.rt_targ(c) = median(r2.rt(aa&cm & r2.correct & v2),'omitnan');   % same detection
        m.conf2{c} = conf_mat(r2.resp_key, {aa&cm, ab&cm, an&cm}, {'j','k','none'});
    end
end

function m = rec_metrics(fdo, min_rt, rec_nm)
% Post-task recognition: old (j) vs new (withhold). d' = z(hit) - z(FA), with
% hit rate by study condition and a shared false-alarm rate over fresh foils.
    T = recode(fdo.results_recognition);
    tt = string(T.trial_type);
    is_old = tt=="old"; is_new = ~is_old;
    n_new = sum(is_new);
    m.far = mean(strcmp(T.resp_key(is_new),'j') & T.rt(is_new) > min_rt);
    cc = string(T.condition);
    m.d = nan(1,numel(rec_nm)); m.hit = nan(1,numel(rec_nm));
    for c = 1:numel(rec_nm)
        mk = is_old & cc==rec_nm{c};
        nh = sum(mk); if nh==0, continue; end
        h  = mean(strcmp(T.resp_key(mk),'j') & T.rt(mk) > min_rt);
        m.hit(c) = h;
        m.d(c)   = zc(h,nh) - zc(m.far,n_new);
    end
end

function T = recode(T)
    T.resp_key = cellstr(T.resp_key); T.resp_key(strcmp(T.resp_key,'NA')) = {'none'};
    T.correct  = strcmp(cellstr(T.corr_resp), T.resp_key);
end
function p = pmean(correct, mask), if sum(mask)==0, p=NaN; else, p=mean(correct(mask)); end, end
function p = pkey(T, mask, key), if sum(mask)==0, p=NaN; else, p=mean(strcmp(T.resp_key(mask), key)); end, end
function z = zc(p, n)
    if isnan(p)||n==0, z=NaN; return; end
    z = norminv(max(1/(2*n), min(1-1/(2*n), p)));
end

function grp_line(lbl, v)
    v = v(~isnan(v));
    if isempty(v), fprintf('  %-22s (no data)\n', lbl); return; end
    fprintf('  %-22s mean %.3f  SD %.3f  [range %.3f-%.3f]  n=%d\n', ...
        lbl, mean(v), std(v), min(v), max(v), numel(v));
end

function paired_line(lbl, a, b, min_n)
    ok = ~isnan(a) & ~isnan(b); a = a(ok); b = b(ok); n = numel(a);
    if n < min_n
        fprintf('  %-30s diff %+.3f (n=%d, too few for a test)\n', lbl, mean(b-a), n);
        return;
    end
    p = signrank(a, b);
    fprintf('  %-30s %+.3f vs %+.3f, diff %+.3f, signed-rank p=%.4f %s\n', ...
        lbl, mean(a), mean(b), mean(b-a), p, stars(p));
end

function paired_plot(M, lvllbl, cols, ylbl, ttl, min_n, pairs, FS)
% M: [nS x nL] metric per subject per level. Paired-plot style.
    if nargin < 8 || isempty(FS), FS = struct('tick',16,'lab',18,'ttl',18,'anno',13); end
    [nS, nL] = size(M);
    if nargin < 7 || isempty(pairs)
        pairs = arrayfun(@(k) [k k+1], 1:nL-1, 'UniformOutput', false);
    end
    hold on; jw = 0.07; bw = 0.34;
    dot_area = 90;   % per-subject dot size (scatter marker area, points^2)
    X = (1:nL) + (rand(nS,nL)-0.5)*2*jw;         % jittered x per subject/level

    % gray lines connecting each subject across levels
    for s = 1:nS
        y = M(s,:);
        plot(X(s,:), y, '-', 'Color', [0.55 0.55 0.55 0.45], 'LineWidth', 0.6);
    end
    % box + whiskers per level, then coloured subject dots
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
        scatter(X(:,L), M(:,L), dot_area, cols{L}, 'filled', 'MarkerFaceAlpha',0.85, ...
            'MarkerEdgeColor',[.2 .2 .2], 'LineWidth',0.3);
    end

    all_v = M(~isnan(M));
    yr = range(all_v); if isempty(yr)||yr==0, yr = 1; end

    % significance brackets (paired signed-rank), only with enough subjects
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

    % dashed linear trend (only for >=3 ordered levels, enough N)
    if nL >= 3 && nS >= min_n
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

function M = conf_mat(resp, masks, keys)
% per-subject confusion: row = presented (mask), col = response key; proportions
    M = nan(numel(masks), numel(keys));
    for r = 1:numel(masks)
        n = sum(masks{r}); if n==0, continue; end
        for c = 1:numel(keys)
            M(r,c) = mean(strcmp(resp(masks{r}), keys{c}));
        end
    end
end

function [M, SE] = grand_conf(A)
% grand mean (+/- SE across subjects) of a stack of confusion matrices [r x c x nS]
    M  = mean(A, 3, 'omitnan');
    n  = sum(~isnan(A), 3);
    SE = std(A, 0, 3, 'omitnan') ./ sqrt(max(n,1));
    SE(n < 2) = NaN;
end

function draw_matrix_se(mat, se, cols, ylbl, xlbl, FS)
% confusion matrix: cell shaded by value, mean + (+/-SE) text
    hold on; nR = size(mat,1); nC = size(mat,2);
    for r = 1:nR
        for c = 1:nC
            v = mat(r,c); if isnan(v), v = 0; end
            s = se(r,c);
            t_c = (v > 0.5)*[1 1 1] + (v <= 0.5)*[0.2 0.2 0.2];
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], ...
                cols{r}, 'EdgeColor','none', 'FaceAlpha', v);
            text(c, r-0.10, sprintf('%.2f', v), 'HorizontalAlignment','center', ...
                'Color',t_c, 'FontWeight','bold', 'FontSize',FS.lab);
            if ~isnan(s)
                text(c, r+0.22, sprintf('\\pm%.2f', s), 'HorizontalAlignment','center', ...
                    'Color',t_c, 'FontSize',FS.anno);
            end
            if r == 1
                text(c, 0.32, xlbl{c}, 'HorizontalAlignment','center', ...
                    'Color',cols{c}, 'FontWeight','bold', 'FontSize',FS.tick);
            end
        end
        text(0.42, r, ylbl{r}, 'HorizontalAlignment','right', ...
            'Color',cols{r}, 'FontWeight','bold', 'FontSize',FS.tick);
    end
    axis ij equal; xlim([0.05 3.85]); ylim([0.15 3.55]);
    set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none'); hold off;
end

function nice_yticks(target_n)
% reduce y-axis to ~target_n nicely-rounded ticks
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

function [f, L] = panel_grid(GEO, nLs)
% Build a figure whose panels have a FIXED pixel size: each cell is
% nLs(r,c)*GEO.col wide by GEO.h tall. Two-group panels are therefore exactly
% the same size in every figure. nLs is a [nrow x ncol] matrix of condition
% counts; L is a [nrow*ncol x 4] array of [left bottom w h], row-major.
    [nrow, ncol] = size(nLs);
    col=GEO.col; h=GEO.h; ml=GEO.ml; mr=GEO.mr; mb=GEO.mb; mt=GEO.mt; hg=GEO.hg; vg=GEO.vg;
    figW = max(ml + sum(nLs*col,2) + hg*(ncol-1) + mr);
    figH = mb + nrow*h + vg*(nrow-1) + mt;
    f = figure('color','w','Position',[40 40 figW figH]);
    L = zeros(nrow*ncol, 4); idx = 0;
    for r = 1:nrow
        b = mb + (nrow-r)*(h+vg);        % row 1 sits at the top
        x = ml;
        for c = 1:ncol
            w = nLs(r,c)*col;
            idx = idx+1; L(idx,:) = [x, b, w, h];
            x = x + w + hg;
        end
    end
end

function ax = mk_axes(f, pos)
% Axes at a fixed pixel rectangle; making it current so paired_plot draws here.
    ax = axes('Parent', f, 'Units','pixels', 'Position', pos);
end

function save_fig(f, fig_dir, name)
    set(f,'Renderer','painters','PaperPositionMode','auto');
    pos = get(f,'Position'); set(f,'PaperUnits','points','PaperSize',pos(3:4));
    print(f, fullfile(fig_dir,[name '.png']), '-dpng','-r150');
end
