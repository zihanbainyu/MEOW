clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NeurWEM_plt -- GROUP analysis, MST-Back (1-back + 2-back) only
%
% Auto-detects every subject with an n-back concat, computes per-subject
% behavioural metrics, and draws them in the paired-plot style: per-level box
% + jittered subject dots + gray lines connecting each subject across levels,
% significance brackets, and a dashed linear trend with an R^2 / p / equation
% ("stat_cor") annotation for the >=3-level panels.
%
% Stats scale with N: paired tests and the trend are only drawn once enough
% subjects are present (they are skipped, not faked, at small N).
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

c_same = [97 125 184]/255; c_sim = [255 191 205]/255; c_new = [219 219 219]/255;
c_comp = [180 174 211]/255; c_nov = [183 210 205]/255;
% ordered goal colours (lure -> target -> foil = hard -> easy)
c_ab = [214 96 77]/255; c_aa = [103 169 207]/255; c_an = [90 180 172]/255;

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
        subs(end+1) = id; %#ok<AGROW>
    end
end
subs = sort(subs);
nS = numel(subs);
assert(nS >= 1, 'No subjects with an n-back concat found in %s', data_dir);
fprintf('Group MST-Back: %d subjects [%s]\n', nS, num2str(subs));

%% ---- per-subject metrics ----
one_acc  = nan(nS,3);   % 1-back accuracy: same / similar / new
goal_acc = nan(nS,3);   % 2-back accuracy by goal: AB(lure) / AA(target) / AN(foil)
ldi      = nan(nS,2);   % 2-back LDI: compared / novel
dpr      = nan(nS,2);   % 2-back d' : compared / novel
for si = 1:nS
    f = fullfile(data_dir, sprintf('sub%03d', subs(si)), sprintf('sub%03d_concat.mat', subs(si)));
    L = load(f, 'final_data_output'); fdo = L.final_data_output;
    m = nback_metrics(fdo, min_rt);
    one_acc(si,:)  = m.one_acc;
    goal_acc(si,:) = m.goal_acc;
    ldi(si,:)      = m.ldi;
    dpr(si,:)      = m.dpr;
end

%% ---- report ----
diary_file = fullfile(res_dir, 'group_mstback_report.txt');
if exist(diary_file,'file'), delete(diary_file); end
diary(diary_file);
fprintf('================================================================\n');
fprintf(' GROUP MST-BACK  --  N = %d   (%s)\n', nS, datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('================================================================\n');
grp_line('1-back same',    one_acc(:,1));
grp_line('1-back similar', one_acc(:,2));
grp_line('1-back new',     one_acc(:,3));
grp_line('2-back AB acc',  goal_acc(:,1));
grp_line('2-back AA acc',  goal_acc(:,2));
grp_line('2-back AN acc',  goal_acc(:,3));
grp_line('LDI compared',   ldi(:,1));    grp_line('LDI novel', ldi(:,2));
grp_line('d''  compared',  dpr(:,1));    grp_line('d''  novel', dpr(:,2));
fprintf('\n*** PRIMARY OUTCOME: 2-back lure discrimination (LDI), compared vs novel ***\n');
paired_line('LDI compared vs novel', ldi(:,1), ldi(:,2), min_n_stat);
fprintf('    (per-subject compared - novel: %s)\n', ...
    strjoin(compose('%+.3f', ldi(:,1)-ldi(:,2))', ', '));
fprintf('\n-- secondary --\n');
paired_line('d''  compared vs novel', dpr(:,1), dpr(:,2), min_n_stat);
diary off;

%% ---- figures (reference paired-plot style) ----
% PRIMARY: 2-back lure discrimination (LDI), compared vs novel -- its own panel
f = figure('color','w','Position',[100 100 560 560],'Name','PRIMARY: LDI compared vs novel');
paired_plot(ldi, {'compared','novel'}, {c_comp,c_nov}, 'LDI  ( p(sim|lure) - p(sim|foil) )', ...
    '2-back lure discrimination', min_n_stat, {[1 2]});
yline(0,'k-');
save_fig(f, fig_dir, 'group_mstback_LDI');

f = figure('color','w','Position',[60 60 1300 430],'Name','MST-Back accuracy');
subplot(1,2,1);
paired_plot(one_acc, {'same','similar','new'}, {c_same,c_sim,c_new}, ...
    'accuracy', '1-back accuracy', min_n_stat, {[1 2],[2 3],[1 3]});
ylim([0 1.05]); yline(1/3,'k:','chance');
subplot(1,2,2);
paired_plot(goal_acc, {'AB (lure)','AA (target)','AN (foil)'}, {c_ab,c_aa,c_an}, ...
    'accuracy', '2-back accuracy by goal', min_n_stat, {[1 2],[2 3],[1 3]});
ylim([0 1.05]); yline(1/3,'k:','chance');
save_fig(f, fig_dir, 'group_mstback_accuracy');

f = figure('color','w','Position',[60 60 900 430],'Name','MST-Back 2-back indices');
subplot(1,2,1);
paired_plot(ldi, {'compared','novel'}, {c_comp,c_nov}, 'LDI', ...
    '2-back lure discrimination', min_n_stat, {[1 2]});
yline(0,'k-');
subplot(1,2,2);
paired_plot(dpr, {'compared','novel'}, {c_comp,c_nov}, 'd''', ...
    '2-back target detection', min_n_stat, {[1 2]});
yline(0,'k-');
save_fig(f, fig_dir, 'group_mstback_indices');

fprintf('\nsaved report + figures to %s , %s\n', res_dir, fig_dir);

%% ==================== local functions ====================
function m = nback_metrics(fdo, min_rt)
    r1 = recode(fdo.results_1_back_all);
    r2 = recode(fdo.results_2_back_all);
    % 1-back accuracy: same / similar / new
    i_sam = r1.condition=="repeat"   & strcmp(r1.corr_resp,'j');
    i_sim = r1.condition=="compared" & r1.identity=="B";
    i_new = (r1.condition=="compared" & r1.identity=="A") | ...
            (r1.condition=="repeat"   & strcmp(r1.corr_resp,'none')) | ...
             r1.condition=="filler";
    m.one_acc = [pmean(r1.correct,i_sam), pmean(r1.correct,i_sim), pmean(r1.correct,i_new)];
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
    % LDI and d' by condition
    conds = {'compared','novel'};
    for c = 1:2
        cm = strcmp(r2.condition, conds{c});
        m.ldi(c) = pkey(r2,ab&cm,'k') - pkey(r2,an&cm,'k');
        nh = sum(aa&cm); nf = sum(an&cm);
        m.dpr(c) = zc(pkey(r2,aa&cm,'j'),nh) - zc(pkey(r2,an&cm,'j'),nf);
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

function paired_plot(M, lvllbl, cols, ylbl, ttl, min_n, pairs)
% M: [nS x nL] metric per subject per level. Paired-plot style.
    [nS, nL] = size(M);
    if nargin < 7 || isempty(pairs)
        pairs = arrayfun(@(k) [k k+1], 1:nL-1, 'UniformOutput', false);
    end
    hold on; jw = 0.07; bw = 0.34;
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
        scatter(X(:,L), M(:,L), 34, cols{L}, 'filled', 'MarkerFaceAlpha',0.85, ...
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
            text(mean(ij), y, s, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11);
        end
    end

    % dashed linear trend + stat_cor (only for >=3 ordered levels, enough N)
    if nL >= 3 && nS >= min_n
        xa = repmat(1:nL, nS, 1); xa = xa(:); ya = M(:);
        ok = ~isnan(ya); xa = xa(ok); ya = ya(ok);
        if numel(unique(xa)) >= 2
            [r, pv] = corr(xa, ya); cf = polyfit(xa, ya, 1);
            xx = [0.7 nL+0.3];
            plot(xx, polyval(cf, xx), 'k--', 'LineWidth', 1.6);
            xt = 0.7; yt = min(all_v);
            text(xt, yt, sprintf('R^2 = %.2f, p = %.2g\ny = %.2f + %.2f x', ...
                r^2, pv, cf(2), cf(1)), 'FontSize', 10, 'VerticalAlignment','bottom');
        end
    end

    set(gca,'XTick',1:nL,'XTickLabel',lvllbl,'FontSize',12);
    xlim([0.4 nL+0.6]); ylabel(ylbl,'FontSize',13);
    if ~isempty(ttl), title(ttl,'FontSize',13); end
    if nS >= min_n && k > 0
        ylim([min(all_v)-0.06*yr, base + k*step + 0.05*yr]);
    end
    box off; hold off;
end

function s = stars(p), s = repmat('*', 1, (p<0.05)+(p<0.01)+(p<0.001)); end

function save_fig(f, fig_dir, name)
    set(f,'Renderer','painters','PaperPositionMode','auto');
    pos = get(f,'Position'); set(f,'PaperUnits','points','PaperSize',pos(3:4));
    print(f, fullfile(fig_dir,[name '.png']), '-dpng','-r150');
end
