clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restricted-viewing (WEM_ResView) vs free-viewing (Experiment_1) gaze.
%
% READ-ONLY on Experiment_1: this loads its already-extracted combined fixation
% table and never writes to that project. Both files are fixation-level with the
% same v_fix schema, so no re-extraction is needed.
%
% Between-groups (different participants), so contrasts are rank-sum (Mann-
% Whitney) with Cliff's delta -- legitimate at small restricted N only when the
% groups separate cleanly, which the manipulation should produce.
%
% Aesthetic matches S_group_behav.m (fixed-size panels, dot + box, FS scale).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---- config ----
exp1_file = fullfile('..','..','Experiment_1','results','group_eye_movement_combined.mat');  % free-viewing (READ ONLY)
wem_file  = fullfile('..','data','eye_movement_data','group_eye_movement.mat');               % restricted
res_dir   = fullfile('..','results');
fig_dir   = fullfile('..','figures');
if ~exist(res_dir,'dir'), mkdir(res_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

FS  = struct('tick',20,'lab',20,'ttl',20,'anno',20);
GEO = struct('col',150, 'h',380, 'ml',110, 'mr',40, 'mb',72, 'mt',58, 'hg',105, 'vg',120);
SCREEN = [1920 1080];   % display size (px); fixations outside it are track-loss artifacts
c_free = [140 140 140]/255;      % free viewing (grey)
c_rstr = [87 6 140]/255;         % restricted viewing (NYU violet)
grp_lbl  = {'free','restricted'};
grp_cols = {c_free, c_rstr};

%% ---- load both fixation tables ----
Lf = load(exp1_file, 'Mw');           F = Lf.Mw;         % free-viewing fixations
Lr = load(wem_file,  'merged');       R = Lr.merged;     % restricted fixations (presentation window)
F.task = string(F.task);  R.task = string(R.task);

% identical on-screen filter for both groups: drop track-loss / out-of-monitor
% fixations so the comparison is not distorted by artifacts on one side.
onscr = @(T) T.x>=0 & T.x<=SCREEN(1)-1 & T.y>=0 & T.y<=SCREEN(2)-1 & isfinite(T.x) & isfinite(T.y);
kf = onscr(F);  kr = onscr(R);
fprintf('on-screen filter (%dx%d): free kept %d/%d (%.1f%% dropped) | restricted kept %d/%d (%.1f%% dropped)\n', ...
    SCREEN(1), SCREEN(2), sum(kf), numel(kf), 100*mean(~kf), sum(kr), numel(kr), 100*mean(~kr));
F = F(kf, :);  R = R(kr, :);
fprintf('free (Exp1): %d fixations, %d subjects | restricted (ResView): %d fixations, %d subjects\n', ...
    height(F), numel(unique(F.subj_id)), height(R), numel(unique(R.subj_id)));

%% ---- common center = restricted group's median fixation (the fixation cross) ----
% The restricted group was held at the central fixation target, so their median
% fixation is the screen/stimulus center. Applied to BOTH groups. Group medians
% are printed so a monitor mismatch would be visible.
cx = median(R.x, 'omitnan'); cy = median(R.y, 'omitnan');
fprintf('common center (restricted median): (%.1f, %.1f)\n', cx, cy);
fprintf('free-group median fixation       : (%.1f, %.1f)   [should match if same monitor]\n', ...
    median(F.x,'omitnan'), median(F.y,'omitnan'));
fprintf('free x-range [%.0f %.0f]  y-range [%.0f %.0f]\n', ...
    min(F.x), max(F.x), min(F.y), max(F.y));

%% ---- per-subject metrics for each group ----
[Pf, subs_f] = persubj(F, cx, cy);
[Pr, subs_r] = persubj(R, cx, cy);

M = struct();
M.nfix = {Pf.nfix, Pr.nfix};   % mean fixations per trial
M.rate = {Pf.rate, Pr.rate};   % median fixations per second (duration-robust)
M.sacc = {Pf.sacc, Pr.sacc};   % saccade amplitude (px) -- exploration
M.scan = {Pf.scan, Pr.scan};   % scanpath length (px)  -- exploration
M.area = {Pf.area, Pr.area};   % explored area (px^2)   -- exploration
M.disp = {Pf.disp, Pr.disp};   % within-trial dispersion (RMS px from centroid)
M.dev  = {Pf.dev,  Pr.dev};    % deviation from center (px)

%% ---- report ----
diary_file = fullfile(res_dir, 'compare_gaze_freeviewing_report.txt');
if exist(diary_file,'file'), delete(diary_file); end
diary(diary_file);
fprintf('================================================================\n');
fprintf(' GAZE: restricted (N=%d) vs free viewing (N=%d)   (%s)\n', ...
    numel(subs_r), numel(subs_f), datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('================================================================\n');
fprintf(['NOTE: restricted fixations are the presentation-window set; this\n' ...
         'assumes the free-viewing Mw is likewise viewing-period fixations. If\n' ...
         'Mw spans the whole trial epoch, free counts are inflated (the contrast\n' ...
         'direction still holds; magnitude does not).\n']);
fprintf('\n-- sampling --');
report_metric('fixations per trial',        M.nfix);
report_metric('fixations per second',        M.rate);
fprintf('\n-- exploration --');
report_metric('saccade amplitude (px)',      M.sacc);
report_metric('scanpath length (px)',        M.scan);
report_metric('explored area (px^2)',        M.area);
fprintf('\n-- spread --');
report_metric('within-trial dispersion (px)', M.disp);
report_metric('deviation from center (px)',  M.dev);
diary off;

%% ---- figure: exploration (top row) + sampling/spread (bottom row) ----
[f, L] = panel_grid(GEO, [2 2 2; 2 2 2]);   set(f,'Name','free vs restricted gaze');
mk_axes(f, L(1,:));
group_compare_plot(M.sacc, grp_lbl, grp_cols, 'saccade amp (px)', 'saccade amplitude', FS);
nice_yticks(4);
mk_axes(f, L(2,:));
group_compare_plot(M.scan, grp_lbl, grp_cols, 'scanpath (px)', 'scanpath length', FS);
nice_yticks(4);
mk_axes(f, L(3,:));
group_compare_plot(M.area, grp_lbl, grp_cols, 'area (px^2)', 'explored area', FS);
nice_yticks(4);
mk_axes(f, L(4,:));
group_compare_plot(M.rate, grp_lbl, grp_cols, 'fixations / s', 'sampling rate', FS);
nice_yticks(4);
mk_axes(f, L(5,:));
group_compare_plot(M.disp, grp_lbl, grp_cols, 'dispersion (px)', 'within-trial spread', FS);
nice_yticks(4);
mk_axes(f, L(6,:));
group_compare_plot(M.dev, grp_lbl, grp_cols, 'deviation (px)', 'distance from center', FS);
nice_yticks(4);
save_fig(f, fig_dir, 'compare_gaze_freeviewing');

%% ---- figures: fixation-density overlays (free vs restricted) ----
% "Where did gaze land?" is a density question: with ~10^5 overlapping fixations,
% colouring by condition or subject overplots into mud, so a single-hue density
% map carries the spatial pattern. Reference elements are drawn to true on-screen
% scale at the stimulus centre. The fixation TARGET and the gaze windows apply
% only to the restricted (gaze-contingent) task; free viewing has neither.
scx = (SCREEN(1)-1)/2; scy = (SCREEN(2)-1)/2;   % true stimulus centre
cfg = struct('binpx',20, 'IMG',400, 'GATE',100, 'BREAK',150, ...   % px; GATE/BREAK are the enforced radii (main.m)
             'SCREEN',SCREEN, 'FS',FS, 'fig_dir',fig_dir);
draw_density(F, scx, scy, 'free', cfg, ...
    sprintf('free viewing: fixation density  (N = %d, %d fixations)', numel(subs_f), height(F)), ...
    'free_viewing_fixation_density');
draw_density(R, scx, scy, 'restricted', cfg, ...
    sprintf('restricted viewing: fixation density  (N = %d, %d fixations)', numel(subs_r), height(R)), ...
    'resview_fixation_density');

fprintf('\nsaved report + figures to %s , %s\n', res_dir, fig_dir);

%% ==================== local functions ====================
function [P, subs] = persubj(T, cx, cy)
% Per-subject gaze summaries from a fixation table, aggregated over trials.
% Sampling  : nfix (mean fixations/trial), rate (median fixations/s, duration-
%             robust via the observed viewing span).
% Exploration: sacc (median saccade amplitude = inter-fixation step, px),
%             scan (median scanpath length = total gaze travel/trial, px),
%             area (median explored area = fixation convex hull, px^2).
% Spread    : disp (RMS px from trial centroid), dev (deviation from center).
    [g, gk] = findgroups(T(:, {'subj_id','task','trial_id'}));
    nfix_t = splitapply(@numel, T.x, g);
    span_t = splitapply(@(on,du) (max(on+du)-min(on))/1000, T.onset, T.dur, g);   % viewing span (s)
    rate_t = nfix_t ./ max(span_t, eps);
    disp_t = splitapply(@(x,y) sqrt(mean((x-mean(x)).^2 + (y-mean(y)).^2)), T.x, T.y, g);
    scan_t = splitapply(@(x,y,on) scanpath_len(x,y,on), T.x, T.y, T.onset, g);
    sacc_t = splitapply(@(x,y,on) sacc_amp(x,y,on),     T.x, T.y, T.onset, g);
    area_t = splitapply(@(x,y) hull_area(x,y),          T.x, T.y, g);
    dev    = hypot(T.x - cx, T.y - cy);
    subs = unique(T.subj_id)';
    z = nan(numel(subs),1);
    P = struct('nfix',z,'rate',z,'sacc',z,'scan',z,'area',z,'disp',z,'dev',z);
    for i = 1:numel(subs)
        tr = gk.subj_id == subs(i);
        P.nfix(i) = mean(nfix_t(tr), 'omitnan');
        P.rate(i) = median(rate_t(tr), 'omitnan');
        P.sacc(i) = median(sacc_t(tr), 'omitnan');
        P.scan(i) = median(scan_t(tr), 'omitnan');
        P.area(i) = median(area_t(tr), 'omitnan');
        P.disp(i) = median(disp_t(tr), 'omitnan');
        P.dev(i)  = median(dev(T.subj_id == subs(i)), 'omitnan');
    end
end

function L = scanpath_len(x,y,on)
% total gaze travel within a trial (fixations ordered by onset)
    [~,o] = sort(on); x = x(o); y = y(o);
    if numel(x) < 2, L = 0; return; end
    L = sum(hypot(diff(x), diff(y)));
end

function a = sacc_amp(x,y,on)
% mean saccade amplitude = mean inter-fixation step; NaN if no saccade made
    [~,o] = sort(on); x = x(o); y = y(o);
    if numel(x) < 2, a = NaN; return; end
    a = mean(hypot(diff(x), diff(y)));
end

function A = hull_area(x,y)
% area of the fixation convex hull (spatial extent explored); 0 if < 3 points
    if numel(x) < 3, A = 0; return; end
    pts = unique([x(:) y(:)], 'rows');
    if size(pts,1) < 3, A = 0; return; end
    try, [~, A] = convhull(pts(:,1), pts(:,2)); catch, A = 0; end   % collinear -> 0
end

function report_metric(lbl, vals)
    a = vals{1}(~isnan(vals{1})); b = vals{2}(~isnan(vals{2}));
    fprintf('\n%s\n', lbl);
    fprintf('  free       median %.2f  [IQR %.2f-%.2f]  mean %.2f  n=%d\n', ...
        median(a), quantile(a,.25), quantile(a,.75), mean(a), numel(a));
    fprintf('  restricted median %.2f  [IQR %.2f-%.2f]  mean %.2f  n=%d\n', ...
        median(b), quantile(b,.25), quantile(b,.75), mean(b), numel(b));
    if numel(a) >= 2 && numel(b) >= 2
        p = ranksum(a, b); d = cliffs_delta(a, b);
        fprintf('  rank-sum: p = %.4f %s | Cliff''s delta = %+.2f (free - restricted)\n', ...
            p, stars(p), d);
    else
        fprintf('  (too few for a rank-sum test)\n');
    end
end

function d = cliffs_delta(a, b)
% Cliff's delta: P(a>b) - P(a<b). Robust nonparametric effect size in [-1, 1].
    a = a(:); b = b(:); s = 0;
    for i = 1:numel(a), s = s + sum(sign(a(i) - b)); end
    d = s / (numel(a)*numel(b));
end

function group_compare_plot(vals, labels, cols, ylbl, ttl, FS)
% Independent-groups dot + box plot (no connecting lines), matched to the
% paired_plot cosmetics elsewhere. Two groups get a rank-sum bracket.
    if nargin < 6 || isempty(FS), FS = struct('tick',20,'lab',20,'ttl',20,'anno',20); end
    hold on; nG = numel(vals); dot_area = 90; jw = 0.09; bw = 0.34; all_v = [];
    for gi = 1:nG
        d = vals{gi}(:); d = d(~isnan(d)); all_v = [all_v; d]; %#ok<AGROW>
        if numel(d) >= 2
            q = quantile(d,[.25 .5 .75]); iqrv = q(3)-q(1);
            lo_w = max(min(d), q(1)-1.5*iqrv); hi_w = min(max(d), q(3)+1.5*iqrv);
            plot([gi gi],[lo_w q(1)],'k-','LineWidth',0.9);
            plot([gi gi],[q(3) hi_w],'k-','LineWidth',0.9);
            plot(gi+[-.06 .06],[lo_w lo_w],'k-','LineWidth',0.9);
            plot(gi+[-.06 .06],[hi_w hi_w],'k-','LineWidth',0.9);
            rectangle('Position',[gi-bw/2, q(1), bw, max(iqrv,eps)], ...
                'EdgeColor','k','LineWidth',1.1,'FaceColor','none');
            plot(gi+[-bw/2 bw/2],[q(2) q(2)],'k-','LineWidth',1.8);
        end
        x = gi + (rand(numel(d),1)-0.5)*2*jw;
        scatter(x, d, dot_area, cols{gi}, 'filled', 'MarkerFaceAlpha',0.85, ...
            'MarkerEdgeColor',[.2 .2 .2], 'LineWidth',0.3);
    end
    if isempty(all_v), all_v = 0; end
    yr = range(all_v); if yr==0, yr = 1; end
    if nG == 2
        a = vals{1}(~isnan(vals{1})); b = vals{2}(~isnan(vals{2}));
        if numel(a) >= 2 && numel(b) >= 2
            p = ranksum(a,b); dd = cliffs_delta(a,b); s = stars(p); if isempty(s), s='ns'; end
            y = max(all_v) + 0.10*yr;
            plot([1 1 2 2],[y-0.02*yr y y y-0.02*yr],'k-','LineWidth',0.9);
            text(1.5, y, sprintf('%s  \\delta=%.2f', s, dd), 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize',FS.anno);
            ylim([min(all_v)-0.06*yr, y + 0.16*yr]);
        end
    end
    set(gca,'XTick',1:nG,'XTickLabel',labels,'FontSize',FS.tick);
    xlim([0.4 nG+0.6]); ylabel(ylbl,'FontSize',FS.lab);
    if ~isempty(ttl), title(ttl,'FontSize',FS.ttl); end
    box off; hold off;
end

function s = stars(p), s = repmat('*', 1, (p<0.05)+(p<0.01)+(p<0.001)); end

function nice_yticks(target_n)
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
    [nrow, ncol] = size(nLs);
    col=GEO.col; h=GEO.h; ml=GEO.ml; mr=GEO.mr; mb=GEO.mb; mt=GEO.mt; hg=GEO.hg; vg=GEO.vg;
    figW = max(ml + sum(nLs*col,2) + hg*(ncol-1) + mr);
    figH = mb + nrow*h + vg*(nrow-1) + mt;
    f = figure('color','w','Position',[40 40 figW figH]);
    L = zeros(nrow*ncol, 4); idx = 0;
    for r = 1:nrow
        b = mb + (nrow-r)*(h+vg); x = ml;
        for c = 1:ncol
            w = nLs(r,c)*col;
            idx = idx+1; L(idx,:) = [x, b, w, h];
            x = x + w + hg;
        end
    end
end

function ax = mk_axes(f, pos)
    ax = axes('Parent', f, 'Units','pixels', 'Position', pos);
end

function draw_density(T, scx, scy, mode, cfg, ttl, fname)
% Fixation-density heatmap over the full screen. Reference geometry (image
% frame, and for restricted the gaze windows + fixation target) is drawn to
% true on-screen scale in data coordinates, centred at (scx,scy).
    xe = 0:cfg.binpx:cfg.SCREEN(1); ye = 0:cfg.binpx:cfg.SCREEN(2);
    D  = gsmooth(histcounts2(T.x, T.y, xe, ye), 1.0);
    f = figure('color','w','Position',[80 80 960 640],'Name',fname);
    ax = axes(f); hold(ax,'on');
    imagesc(ax, xe(1:end-1)+cfg.binpx/2, ye(1:end-1)+cfg.binpx/2, D');   % transpose: image rows = y
    colormap(ax, jelly_ramp());
    cmax = quantile(D(D>0), 0.99); if isempty(cmax) || cmax <= 0, cmax = max(D(:)); end
    if cmax > 0, caxis(ax, [0 cmax]); end                % cap so the peak doesn't wash it out
    set(ax,'YDir','reverse'); axis(ax,'image');
    xlim(ax,[0 cfg.SCREEN(1)]); ylim(ax,[0 cfg.SCREEN(2)]);
    % stimulus frame -- on screen during viewing in both tasks
    rectangle(ax,'Position',[scx-cfg.IMG/2, scy-cfg.IMG/2, cfg.IMG, cfg.IMG], ...
        'EdgeColor',[.15 .15 .15],'LineWidth',1.5);
    % 100 px bound -- drawn on both: the enforced gate radius in the restricted
    % task, and a reference on the free map (gaze visibly spills beyond it).
    th = linspace(0, 2*pi, 200);
    plot(ax, scx+cfg.GATE*cos(th), scy+cfg.GATE*sin(th), 'k-', 'LineWidth',1.2);
    if strcmp(mode,'restricted')
        % real fixation target to scale: outer disc 36 px, central dot 12 px (main.m)
        rectangle(ax,'Position',[scx-18, scy-18, 36, 36],'Curvature',[1 1], ...
            'EdgeColor','k','LineWidth',1.5,'FaceColor','none');
        rectangle(ax,'Position',[scx-6,  scy-6,  12, 12],'Curvature',[1 1], ...
            'EdgeColor','none','FaceColor','k');
    end
    cb = colorbar(ax); cb.Label.String = 'fixation density';
    cb.FontSize = cfg.FS.tick; cb.Label.FontSize = cfg.FS.lab;
    set(ax,'FontSize',cfg.FS.tick);
    xlabel(ax,'x (px)','FontSize',cfg.FS.lab); ylabel(ax,'y (px)','FontSize',cfg.FS.lab);
    title(ax, ttl, 'FontSize', cfg.FS.ttl); box(ax,'on');
    save_fig(f, cfg.fig_dir, fname);
end

function cmap = jelly_ramp()
% white -> jelly blue -> deep blue, for the density heatmap (clean single-hue)
    stops = [1.00 1.00 1.00;    % white  (low density)
             0.83 0.91 0.99;    % pale blue
             0.53 0.75 0.96;    % sky blue
             0.22 0.49 0.90;    % jelly blue
             0.06 0.22 0.52];   % deep blue (high density)
    x  = linspace(0, 1, size(stops,1));
    xi = linspace(0, 1, 256)';
    cmap = [interp1(x,stops(:,1),xi), interp1(x,stops(:,2),xi), interp1(x,stops(:,3),xi)];
end

function S = gsmooth(D, sigma)
% small Gaussian blur via conv2 (no Image Processing Toolbox dependency)
    r = max(1, ceil(3*sigma)); [xx,yy] = meshgrid(-r:r, -r:r);
    k = exp(-(xx.^2 + yy.^2)/(2*sigma^2)); k = k/sum(k(:));
    S = conv2(D, k, 'same');
end

function save_fig(f, fig_dir, name)
    set(f,'Renderer','painters','PaperPositionMode','auto');
    pos = get(f,'Position'); set(f,'PaperUnits','points','PaperSize',pos(3:4));
    print(f, fullfile(fig_dir,[name '.png']), '-dpng','-r150');
end
