% %% A1-B2 cross-pair reinstatement: A1 encoding gaze vs B2 retrieval gaze
% % Correlates fixation maps from A1 (1-back encoding) with B2 (2-back retrieval),
% % matched via base_id. Baseline-corrected by other B1 encoding maps.
% clear; clc; close all;
% 
% base_dir = '..';
% res_dir = fullfile(base_dir, 'results');
% 
% load(fullfile(res_dir, 'group_eye_movement_combined.mat'), 'Mw');
% 
% spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);
% 
% comp_idx = strcmp(Mw.condition, 'compared');
% iso_idx  = strcmp(Mw.condition, 'isolated');
% ab_idx   = strcmp(Mw.goal, 'A-B');
% a_idx    = strcmp(Mw.identity, 'A');
% b_idx    = strcmp(Mw.identity, 'B');
% task_1b_idx = strcmp(Mw.task, '1_back');
% task_2b_idx = strcmp(Mw.task, '2_back');
% 
% %% Build trial tables
% tr_1b_a_comp = unique(Mw(task_1b_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
% tr_1b_a_iso  = unique(Mw(task_1b_idx & a_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
% tr_1b_b_comp = unique(Mw(task_1b_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
% tr_1b_b_iso  = unique(Mw(task_1b_idx & b_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
% tr_2b_b_comp = unique(Mw(task_2b_idx & ab_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
% tr_2b_b_iso  = unique(Mw(task_2b_idx & ab_idx & b_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
% 
% %% Match A1-B2 pairs via base_id
% tr_1b_a_comp.base_id = regexprep(tr_1b_a_comp.stim_id, '_A_', '_BASE_');
% tr_1b_a_iso.base_id  = regexprep(tr_1b_a_iso.stim_id,  '_A_', '_BASE_');
% tr_2b_b_comp.base_id = regexprep(tr_2b_b_comp.stim_id, '_B_', '_BASE_');
% tr_2b_b_iso.base_id  = regexprep(tr_2b_b_iso.stim_id,  '_B_', '_BASE_');
% 
% pairs_a1b2_comp = innerjoin(tr_1b_a_comp, tr_2b_b_comp, 'Keys', {'subj_id','base_id'}, ...
%     'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
% pairs_a1b2_comp.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id_a','base_id','tr_2b_b','stim_id_b'};
% 
% pairs_a1b2_iso = innerjoin(tr_1b_a_iso, tr_2b_b_iso, 'Keys', {'subj_id','base_id'}, ...
%     'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
% pairs_a1b2_iso.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id_a','base_id','tr_2b_b','stim_id_b'};
% 
% fprintf('A1-B2 pairs: comp=%d, iso=%d\n', height(pairs_a1b2_comp), height(pairs_a1b2_iso));
% 
% %% Add B2 correctness
% pairs_a1b2_comp.correct = nan(height(pairs_a1b2_comp), 1);
% for i = 1:height(pairs_a1b2_comp)
%     sid = pairs_a1b2_comp.subj_id(i);
%     tr  = pairs_a1b2_comp.tr_2b_b(i);
%     match_row = Mw.subj_id == sid & Mw.trial_id == tr & task_2b_idx & b_idx;
%     if any(match_row)
%         pairs_a1b2_comp.correct(i) = Mw.correct(find(match_row, 1));
%     end
% end
% 
% pairs_a1b2_iso.correct = nan(height(pairs_a1b2_iso), 1);
% for i = 1:height(pairs_a1b2_iso)
%     sid = pairs_a1b2_iso.subj_id(i);
%     tr  = pairs_a1b2_iso.tr_2b_b(i);
%     match_row = Mw.subj_id == sid & Mw.trial_id == tr & task_2b_idx & b_idx;
%     if any(match_row)
%         pairs_a1b2_iso.correct(i) = Mw.correct(find(match_row, 1));
%     end
% end
% 
% fprintf('B2 correctness matched: comp=%d/%d, iso=%d/%d\n', ...
%     sum(~isnan(pairs_a1b2_comp.correct)), height(pairs_a1b2_comp), ...
%     sum(~isnan(pairs_a1b2_iso.correct)), height(pairs_a1b2_iso));
% 
% %% Build baseline pools: other B1 encoding trials per subject per condition
% unique_subjs = unique([pairs_a1b2_comp.subj_id; pairs_a1b2_iso.subj_id]);
% subj_baseline_b_comp = struct();
% subj_baseline_b_iso  = struct();
% for s_idx = 1:length(unique_subjs)
%     sid = unique_subjs(s_idx);
%     skey = sprintf('s%d', sid);
%     b_comp = tr_1b_b_comp(tr_1b_b_comp.subj_id == sid, :);
%     b_iso  = tr_1b_b_iso(tr_1b_b_iso.subj_id == sid, :);
%     if height(b_comp) > 0, subj_baseline_b_comp.(skey) = b_comp; end
%     if height(b_iso)  > 0, subj_baseline_b_iso.(skey)  = b_iso;  end
% end
% 
% %% Compute A1-B2 compared
% % Match score: corr(A1_map, B2_map)
% % Baseline: mean(corr(otherB1_map, B2_map))
% fprintf('\nComputing A1-B2 compared (%d pairs)...\n', height(pairs_a1b2_comp));
% results_a1b2_comp = zeros(height(pairs_a1b2_comp), 7);
% n_a1b2_comp = 0;
% 
% for i = 1:height(pairs_a1b2_comp)
%     if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_a1b2_comp)); end
%     pair = pairs_a1b2_comp(i,:);
%     skey = sprintf('s%d', pair.subj_id);
%     if ~isfield(subj_baseline_b_comp, skey), continue; end
% 
%     fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
%     fix_b2 = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_b, :);
%     if height(fix_a) < 2 || height(fix_b2) < 2, continue; end
% 
%     map_a  = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
%     map_b2 = create_fixation_map(fix_b2.x, fix_b2.y, fix_b2.dur, spatial_params);
%     if sum(map_a(:)) == 0 || sum(map_b2(:)) == 0, continue; end
% 
%     match_score = corr(map_a(:), map_b2(:));
% 
%     % Baseline: correlate B2 map with other B1 encoding maps
%     pool = subj_baseline_b_comp.(skey);
%     % Exclude the partner B from baseline
%     partner_b_stim = regexprep(pair.stim_id_a, '_A_', '_B_');
%     baseline_pool = pool(~strcmp(pool.stim_id, partner_b_stim), :);
%     if height(baseline_pool) < 1, continue; end
% 
%     mismatch_scores = nan(height(baseline_pool), 1);
%     for j = 1:height(baseline_pool)
%         fix_other = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
%         if height(fix_other) < 2, continue; end
%         map_other = create_fixation_map(fix_other.x, fix_other.y, fix_other.dur, spatial_params);
%         if sum(map_other(:)) == 0, continue; end
%         mismatch_scores(j) = corr(map_other(:), map_b2(:));
%     end
%     baseline_score = mean(mismatch_scores, 'omitnan');
% 
%     n_a1b2_comp = n_a1b2_comp + 1;
%     results_a1b2_comp(n_a1b2_comp,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_b, match_score, baseline_score, match_score - baseline_score, pair.correct];
% end
% results_a1b2_comp = results_a1b2_comp(1:n_a1b2_comp,:);
% fprintf('A1-B2 compared complete: %d valid pairs\n\n', n_a1b2_comp);
% 
% %% Compute A1-B2 isolated
% fprintf('Computing A1-B2 isolated (%d pairs)...\n', height(pairs_a1b2_iso));
% results_a1b2_iso = zeros(height(pairs_a1b2_iso), 7);
% n_a1b2_iso = 0;
% 
% for i = 1:height(pairs_a1b2_iso)
%     if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_a1b2_iso)); end
%     pair = pairs_a1b2_iso(i,:);
%     skey = sprintf('s%d', pair.subj_id);
%     if ~isfield(subj_baseline_b_iso, skey), continue; end
% 
%     fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
%     fix_b2 = Mw(task_2b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_2b_b, :);
%     if height(fix_a) < 2 || height(fix_b2) < 2, continue; end
% 
%     map_a  = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
%     map_b2 = create_fixation_map(fix_b2.x, fix_b2.y, fix_b2.dur, spatial_params);
%     if sum(map_a(:)) == 0 || sum(map_b2(:)) == 0, continue; end
% 
%     match_score = corr(map_a(:), map_b2(:));
% 
%     pool = subj_baseline_b_iso.(skey);
%     partner_b_stim = regexprep(pair.stim_id_a, '_A_', '_B_');
%     baseline_pool = pool(~strcmp(pool.stim_id, partner_b_stim), :);
%     if height(baseline_pool) < 1, continue; end
% 
%     mismatch_scores = nan(height(baseline_pool), 1);
%     for j = 1:height(baseline_pool)
%         fix_other = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(j), :);
%         if height(fix_other) < 2, continue; end
%         map_other = create_fixation_map(fix_other.x, fix_other.y, fix_other.dur, spatial_params);
%         if sum(map_other(:)) == 0, continue; end
%         mismatch_scores(j) = corr(map_other(:), map_b2(:));
%     end
%     baseline_score = mean(mismatch_scores, 'omitnan');
% 
%     n_a1b2_iso = n_a1b2_iso + 1;
%     results_a1b2_iso(n_a1b2_iso,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_2b_b, match_score, baseline_score, match_score - baseline_score, pair.correct];
% end
% results_a1b2_iso = results_a1b2_iso(1:n_a1b2_iso,:);
% fprintf('A1-B2 isolated complete: %d valid pairs\n\n', n_a1b2_iso);
% 
% %% Convert to tables and save
% results_a1b2_comp = array2table(results_a1b2_comp, 'VariableNames', {'subj_id','tr_1b_a','tr_2b_b','match_score','baseline_score','reinst_index','correct'});
% results_a1b2_iso  = array2table(results_a1b2_iso,  'VariableNames', {'subj_id','tr_1b_a','tr_2b_b','match_score','baseline_score','reinst_index','correct'});
% 
% reinstat_res_a1b2 = struct();
% reinstat_res_a1b2.a1b2_compared = results_a1b2_comp;
% reinstat_res_a1b2.a1b2_isolated = results_a1b2_iso;
% reinstat_res_a1b2.spatial_params = spatial_params;
% 
% save(fullfile(res_dir, 'gaze_reinstat_res_a1b2.mat'), 'reinstat_res_a1b2');
% fprintf('Saved gaze_reinstat_res_a1b2.mat\n');
% fprintf('Done! A1-B2: comp=%d, iso=%d\n', n_a1b2_comp, n_a1b2_iso);



% 
% %% Helper
% function map = create_fixation_map(x, y, dur, params)
%     in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & y >= params.roi_y(1) & y <= params.roi_y(2);
%     x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
%     if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end
%     x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
%     y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
%     map = zeros(params.grid_y, params.grid_x);
%     for i = 1:length(x)
%         x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
%         y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
%         if ~isempty(x_idx) && ~isempty(y_idx), map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i); end
%     end
%     if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
%     if sum(map(:)) > 0, map = map / sum(map(:)); end
% end


base_dir = '..';
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');

load(fullfile(res_dir, 'gaze_reinstat_res_a1b2.mat'), 'reinstat_res_a1b2');

res_comp = reinstat_res_a1b2.a1b2_compared;
res_iso  = reinstat_res_a1b2.a1b2_isolated;

% exc = [609, 606, 608, 618];
exc = [];
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255;

%% 2x2 ANOVA + figure
sids = unique([res_comp.subj_id; res_iso.subj_id]);
ns   = length(sids);
sm   = nan(ns, 4);
for s = 1:ns, sid = sids(s);
    for j = 1:4, [tbl,acc] = deal({res_comp,res_comp,res_iso,res_iso},{1,0,1,0});
        d = tbl{j}.reinst_index(tbl{j}.subj_id == sid & tbl{j}.correct == acc{j});
        if ~isempty(d), sm(s,j) = mean(d,'omitnan'); end, end, end

for e = exc
    idx = find(sids == e);
    if ~isempty(idx), sm(idx,:) = []; sids(idx) = []; end
end

v = all(~isnan(sm), 2); sm = sm(v,:);
fprintf('A1-B2 SME n=%d\n', sum(v));

fprintf('\n=== A1-B2 by B2 accuracy ===\n');
t = table(sm(:,1),sm(:,2),sm(:,3),sm(:,4),'VariableNames',{'cc','ci','ic','ii'});
within = table({'comp';'comp';'iso';'iso'},{'corr';'incorr';'corr';'incorr'},'VariableNames',{'cond','correctness'});
rm  = fitrm(t,'cc-ii~1','WithinDesign',within);
tbl = ranova(rm,'WithinModel','cond*correctness');
tn  = string(tbl.Properties.RowNames);
rows   = {'(Intercept):cond','(Intercept):correctness','(Intercept):cond:correctness'};
labels = {'condition','correctness','interaction'};
for i = 1:3, idx = find(contains(tn,rows{i}),1); if ~isempty(idx)
    eta = tbl.SumSq(idx)/(tbl.SumSq(idx)+tbl.SumSq(idx+1));
    fprintf('  %s: F(%d,%d)=%.2f, p=%.4f, eta2=%.3f %s\n',labels{i},tbl.DF(idx),tbl.DF(idx+1),tbl.F(idx),tbl.pValue(idx),eta,...
        repmat('*',1,(tbl.pValue(idx)<.05)+(tbl.pValue(idx)<.01)+(tbl.pValue(idx)<.001))); end, end

nv = size(sm,1);
[~,p1,~,s1] = ttest(sm(:,1),sm(:,2)); [~,p2,~,s2] = ttest(sm(:,3),sm(:,4));
[~,pc,~,sc] = ttest(mean(sm(:,1:2),2),mean(sm(:,3:4),2));
fprintf('  comp corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s1.df,s1.tstat,p1,s1.tstat/sqrt(nv),repmat('*',1,(p1<.05)+(p1<.01)+(p1<.001)));
fprintf('  iso  corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s2.df,s2.tstat,p2,s2.tstat/sqrt(nv),repmat('*',1,(p2<.05)+(p2<.01)+(p2<.001)));
fprintf('  comp>iso:         t(%d)=%.2f, p=%.4f, d=%.2f %s\n',sc.df,sc.tstat,pc,sc.tstat/sqrt(nv),repmat('*',1,(pc<.05)+(pc<.01)+(pc<.001)));
fprintf('  comp corr:  M=%.4f SD=%.4f\n  comp incorr: M=%.4f SD=%.4f\n  iso  corr:  M=%.4f SD=%.4f\n  iso  incorr: M=%.4f SD=%.4f\n',...
    mean(sm(:,1)),std(sm(:,1)),mean(sm(:,2)),std(sm(:,2)),mean(sm(:,3)),std(sm(:,3)),mean(sm(:,4)),std(sm(:,4)));

cols = {c_comp, c_comp*.5+.5, c_iso, c_iso*.5+.5};
xp = [1 2 3.5 4.5]; yl = [-0.2 0.5];

figure('Color','w','Position',[100 100 500 500]); hold on;
jt = -0.15 - rand(size(sm)) * .2;
for s = 1:size(sm,1), for g = [1 3]
    plot(xp(g:g+1)+jt(s,g:g+1), sm(s,g:g+1), '-', 'Color', [.7 .7 .7 .4], 'LineWidth', .5);
end, end
for i = 1:4, d = sm(:,i); d_c = d(~isnan(d)); if isempty(d_c), continue; end
    x = xp(i); [f,xi] = ksdensity(d_c); f = f/max(f)*.55;
    patch([x+f, x*ones(1,length(f))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',.5);
    scatter(x+jt(~isnan(d),i), d_c, 35, cols{i}, 'filled', 'MarkerFaceAlpha', .7);
    q = quantile(d_c,[.25 .5 .75]);
    rectangle('Position',[x-.1,q(1),.2,q(3)-q(1)],'FaceColor',cols{i},'EdgeColor','k','LineWidth',1.2);
    plot([x-.1,x+.1],[q(2) q(2)],'k-','LineWidth',2.5); end
yline(0,'r--','LineWidth',2); ylim(yl); set(gca,'YTick',[-.2 0 .2 .4]);
st = diff(yl)*.06; lv = 0;
for g = [1 3], [~,p] = ttest(sm(:,g),sm(:,g+1)); lv = lv+1; yp = yl(2)-st*(2-lv);
    if p<.001,tx='***';elseif p<.01,tx='**';elseif p<.05,tx='*';else,tx='n.s.';end
    plot([xp(g) xp(g) xp(g+1) xp(g+1)],[yp-st*.3,yp,yp,yp-st*.3],'k-','LineWidth',1.2);
    text(mean(xp(g:g+1)),yp+st*.15,tx,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
set(gca,'XTick',xp,'XTickLabel',{'correct','incorrect','correct','incorrect'},'FontSize',20);
ylabel('Gaze Similarity (r)','FontSize',20); xlim([.3 5.2]); box off;
print(gcf, fullfile(fig_dir,'sme_a1b2.pdf'),'-dpdf','-vector');


function raincloud(mat, cols, xlbls, ylbl, ttl, ylims)
    [n_rows, n_grps] = size(mat); hold on;
    d_v = mat(:); d_v = d_v(~isnan(d_v)); 
    mn = min(d_v); mx = max(d_v); rng = mx-mn; if rng==0, rng=1; end
    auto_lim = [mn-(rng*0.15), mx+(rng*0.15)];
    jit = -0.15 - (rand(size(mat)) * 0.20); x_c = repmat(1:n_grps, n_rows, 1) + jit;
    plot(x_c', mat', '-', 'Color', [0.7 0.7 0.7, 0.4], 'LineWidth', 0.5);
    for i = 1:n_grps
        d = mat(:,i); d = d(~isnan(d)); if isempty(d), continue; end
        [f, xi] = ksdensity(d); f = f/max(f)*0.4;
        patch([i+f, i*ones(1,length(f))], [xi, fliplr(xi)], cols{i}, 'EdgeColor','none','FaceAlpha',0.5);
        scatter(x_c(:,i), mat(:,i), 20, cols{i}, 'filled', 'MarkerFaceAlpha',0.6);
        q = quantile(d, [0.25, 0.5, 0.75]);
        rectangle('Position', [i-0.06, q(1), 0.12, q(3)-q(1)], 'FaceColor', cols{i}, 'EdgeColor','k', 'LineWidth',1.2);
        plot([i-0.06, i+0.06], [q(2), q(2)], 'k-', 'LineWidth', 2);
    end
    set(gca, 'XTick', 1:n_grps, 'XTickLabel', xlbls, 'FontSize', 20);
    if ~isempty(ttl), title(ttl, 'FontSize', 20); end
    ylabel(ylbl,'FontSize',20); xlim([0.2, n_grps+0.8]);
    if nargin>5 && ~isempty(ylims), ylim(ylims); else, ylim(auto_lim); end
    grid off; set(gca,'GridAlpha',0.1); box off; hold off;
end

function add_sig(data, pairs)
    [~, n_grps] = size(data);
    cl = ylim; y_top = cl(2); rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng * 0.08; line_lvl = 0; hold on;
    real_mx = max(data(:)); if isnan(real_mx), real_mx=y_top; end
    base = max(y_top, real_mx + step*0.5);
    for i = 1:size(pairs,1)
        c1 = pairs(i,1); c2 = pairs(i,2);
        if c2==0, [~, p]=ttest(data(:,c1)); paired=false; else, [~, p]=ttest(data(:,c1),data(:,c2)); paired=true; end
        if p < 0.05
            line_lvl=line_lvl+1; y_p = base + (line_lvl*step);
            if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
            if paired, plot([c1,c1,c2,c2],[y_p-step*0.3,y_p,y_p,y_p-step*0.3],'k-','LineWidth',1.2);
               text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',20,'FontWeight','bold');
            else, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
        end
    end
    if line_lvl>0, ylim([cl(1), base+(line_lvl*step)+step]); end; hold off;
end

function draw_matrix(mat, se_mat, cols, ylbl, xlbl)
    hold on;
    for r = 1:3
        for c = 1:3
            v = mat(r,c); s = se_mat(r,c);
            t_c = (v > 0.5)*[1 1 1] + (v <= 0.5)*[0.2 0.2 0.2];
            patch([c-0.48 c+0.48 c+0.48 c-0.48], [r-0.48 r-0.48 r+0.48 r+0.48], cols{r}, 'EdgeColor','none', 'FaceAlpha',v);
            text(c, r-0.08, sprintf('%.2f', v), 'HorizontalAlignment','center', 'Color',t_c, 'FontWeight','bold', 'FontSize',10);
            text(c, r+0.18, sprintf('±%.2f', s), 'HorizontalAlignment','center', 'Color',t_c, 'FontSize',8);
            if r == 1, text(c, 0.4, xlbl{c}, 'HorizontalAlignment','center', 'Color',cols{c}, 'FontWeight','bold', 'FontSize',10); end
        end
        text(0.4, r, ylbl{r}, 'HorizontalAlignment','right', 'Color',cols{r}, 'FontWeight','bold', 'FontSize',10);
    end
    axis ij equal; xlim([0.2 3.8]); ylim([0.3 3.5]); set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); hold off;
end

function [p_val, obs_diff] = run_permutation(cond_A, cond_B, n_perms)
    valid = ~isnan(cond_A) & ~isnan(cond_B);
    diffs = cond_A(valid) - cond_B(valid);
    obs_diff = mean(diffs); null_dist = zeros(n_perms, 1);
    for i = 1:n_perms
        signs = sign(rand(size(diffs)) - 0.5);
        null_dist(i) = mean(diffs .* signs);
    end
    p_val = mean(abs(null_dist) >= abs(obs_diff));
end

function add_sig_perm(data, pairs, pvals)
    cl = ylim; rng = cl(2)-cl(1); if rng==0, rng=1; end
    step = rng * 0.08; line_lvl = 0; hold on;
    base = max(cl(2), max(data(:)) + step*0.5);
    for i = 1:size(pairs,1)
        p = pvals(i); c1 = pairs(i,1); c2 = pairs(i,2);
        if p >= 0.05, continue; end
        line_lvl = line_lvl+1; y_p = base + (line_lvl*step);
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; else, txt = '*'; end
        if c2 == 0, text(c1, y_p, txt, 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
        else, plot([c1 c1 c2 c2], [y_p-step*0.3 y_p y_p y_p-step*0.3], 'k-', 'LineWidth', 1.2);
              text(mean([c1 c2]), y_p+step*0.1, txt, 'HorizontalAlignment','center','FontSize',14,'FontWeight','bold'); end
    end
    if line_lvl>0, ylim([cl(1), base+(line_lvl*step)+step]); end; hold off;
end

function [p_val, obs_diff, null_dist] = run_permutation_with_dist(cond_A, cond_B, n_perms)
    valid = ~isnan(cond_A) & ~isnan(cond_B);
    diffs = cond_A(valid) - cond_B(valid);
    obs_diff = mean(diffs); null_dist = zeros(n_perms, 1);
    for i = 1:n_perms
        signs = sign(rand(size(diffs)) - 0.5);
        null_dist(i) = mean(diffs .* signs);
    end
    p_val = mean(abs(null_dist) >= abs(obs_diff));
end

function [res, posthoc] = run_rm_anova(lbl, d1, d2, d3, within)
    t = table(d1, d2, d3, 'VariableNames', {'compared','isolated','novel'});
    rm = fitrm(t, 'compared-novel~1', 'WithinDesign', within);
    tbl = ranova(rm); eps = epsilon(rm); m = mauchly(rm);
    res.F = tbl.F(1); res.df1 = tbl.DF(1); res.df2_gg = tbl.DF(2)*eps.GreenhouseGeisser;
    res.p_gg = tbl.pValueGG(1); res.eps = eps.GreenhouseGeisser;
    res.eta2_p = tbl.SumSq(1) / (tbl.SumSq(1) + tbl.SumSq(2));
    fprintf('\n%s: F(%d,%.1f) = %.2f, p = %.4f, eta2_p = %.3f, eps = %.3f\n', ...
        lbl, res.df1, res.df2_gg, res.F, res.p_gg, res.eta2_p, res.eps);
    n = length(d1); subj = (1:n)';
    bf_tbl = table(subj, d1, d2, d3, 'VariableNames', {'subj','c1','c2','c3'});
    bf_long = stack(bf_tbl, {'c1','c2','c3'}, 'NewDataVariableName','y', 'IndexVariableName','cond');
    bf_long.subj = categorical(bf_long.subj); bf_long.cond = categorical(bf_long.cond);
    res.bf10 = bf.anova(bf_long, 'y~cond', 'treatAsRandom', {'subj'}, 'verbose', false);
    if res.bf10 < 1, fprintf('  BF01 = %.2f\n', 1/res.bf10); else, fprintf('  BF10 = %.2f\n', res.bf10); end
    cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    [~,p1,~,s1] = ttest(d1,d2); [~,p2,~,s2] = ttest(d2,d3); [~,p3,~,s3] = ttest(d1,d3);
    adj = holm_bonf([p1 p2 p3]); df = s1.df;
    posthoc.comp_iso = struct('t',s1.tstat,'df',df,'d',cd(d1,d2),'p',p1,'p_adj',adj(1));
    posthoc.iso_nov  = struct('t',s2.tstat,'df',df,'d',cd(d2,d3),'p',p2,'p_adj',adj(2));
    posthoc.comp_nov = struct('t',s3.tstat,'df',df,'d',cd(d1,d3),'p',p3,'p_adj',adj(3));
    fprintf('  post-hoc:\n');
    fprintf('    comp-iso: t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s1.tstat, cd(d1,d2), adj(1), repmat('*',1,(adj(1)<.05)+(adj(1)<.01)+(adj(1)<.001)));
    fprintf('    iso-nov:  t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s2.tstat, cd(d2,d3), adj(2), repmat('*',1,(adj(2)<.05)+(adj(2)<.01)+(adj(2)<.001)));
    fprintf('    comp-nov: t(%d) = %.2f, d = %.2f, p_adj = %.4f %s\n', df, s3.tstat, cd(d1,d3), adj(3), repmat('*',1,(adj(3)<.05)+(adj(3)<.01)+(adj(3)<.001)));
end

function adj = holm_bonf(p)
    n = length(p); [ps, idx] = sort(p); as = min(1, ps .* (n:-1:1));
    for i = 2:n, as(i) = max(as(i), as(i-1)); end
    adj(idx) = as;
end