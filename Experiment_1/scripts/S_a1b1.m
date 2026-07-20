% Subsequent memory: A1-B1 encoding alignment predicting B2 accuracy
clear; clc; close all;

base_dir = '..';
res_dir = fullfile(base_dir, 'results');
fig_dir = fullfile(base_dir, 'figures');

% load(fullfile(res_dir, 'group_eye_movement_combined.mat'), 'Mw');
% 
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255; c_nov = [183 210 205]/255; 
c_sim  = [255 191 205]/255; c_same = [97 125 184]/255; c_new = [219 219 219]/255;
cohend = @(x,y) mean(x-y,'omitnan') / std(x-y,'omitnan');
% 
% comp_idx = strcmp(Mw.condition, 'compared');
% iso_idx = strcmp(Mw.condition, 'isolated');
% ab_idx = strcmp(Mw.goal, 'A-B');
% a_idx = strcmp(Mw.identity, 'A');
% b_idx = strcmp(Mw.identity, 'B');
% task_1b_idx = strcmp(Mw.task, '1_back');
% task_2b_idx = strcmp(Mw.task, '2_back');

% spatial_params = struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);
% 
% % Get 1-back A and B trials
% tr_one_b_a_comp = unique(Mw(task_1b_idx & a_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
% tr_one_b_a_iso  = unique(Mw(task_1b_idx & a_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
% tr_one_b_b_comp = unique(Mw(task_1b_idx & b_idx & comp_idx, {'subj_id','trial_id','stim_id'}));
% tr_one_b_b_iso  = unique(Mw(task_1b_idx & b_idx & iso_idx,  {'subj_id','trial_id','stim_id'}));
% 
% % Match A1-B1 pairs via base_id
% tr_one_b_a_comp.base_id = regexprep(tr_one_b_a_comp.stim_id, '_A_', '_BASE_');
% tr_one_b_a_iso.base_id  = regexprep(tr_one_b_a_iso.stim_id,  '_A_', '_BASE_');
% tr_one_b_b_comp.base_id = regexprep(tr_one_b_b_comp.stim_id, '_B_', '_BASE_');
% tr_one_b_b_iso.base_id  = regexprep(tr_one_b_b_iso.stim_id,  '_B_', '_BASE_');
% 
% pairs_ab_comp = innerjoin(tr_one_b_a_comp, tr_one_b_b_comp, 'Keys', {'subj_id','base_id'}, ...
%     'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
% pairs_ab_comp.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id_a','base_id','tr_1b_b','stim_id_b'};
% 
% pairs_ab_iso = innerjoin(tr_one_b_a_iso, tr_one_b_b_iso, 'Keys', {'subj_id','base_id'}, ...
%     'LeftVariables', {'subj_id','trial_id','stim_id','base_id'}, 'RightVariables', {'trial_id','stim_id'});
% pairs_ab_iso.Properties.VariableNames = {'subj_id','tr_1b_a','stim_id_a','base_id','tr_1b_b','stim_id_b'};
% 
% fprintf('A1-B1 encoding pairs: comp=%d, iso=%d\n', height(pairs_ab_comp), height(pairs_ab_iso));
% 
% %% Add B2 lure discrimination correctness
% pairs_ab_comp.correct = nan(height(pairs_ab_comp), 1);
% for i = 1:height(pairs_ab_comp)
%     sid = pairs_ab_comp.subj_id(i);
%     b_rows = Mw(Mw.subj_id == sid & task_2b_idx & b_idx & comp_idx & ab_idx, :);
%     b_base = regexprep(b_rows.stim_id, '_B_', '_BASE_');
%     match = strcmp(b_base, pairs_ab_comp.base_id{i});
%     if any(match)
%         pairs_ab_comp.correct(i) = b_rows.correct(find(match, 1));
%     end
% end
% 
% pairs_ab_iso.correct = nan(height(pairs_ab_iso), 1);
% for i = 1:height(pairs_ab_iso)
%     sid = pairs_ab_iso.subj_id(i);
%     b_rows = Mw(Mw.subj_id == sid & task_2b_idx & b_idx & iso_idx & ab_idx, :);
%     b_base = regexprep(b_rows.stim_id, '_B_', '_BASE_');
%     match = strcmp(b_base, pairs_ab_iso.base_id{i});
%     if any(match)
%         pairs_ab_iso.correct(i) = b_rows.correct(find(match, 1));
%     end
% end
% 
% fprintf('B2 correctness matched: comp=%d/%d, iso=%d/%d\n', ...
%     sum(~isnan(pairs_ab_comp.correct)), height(pairs_ab_comp), ...
%     sum(~isnan(pairs_ab_iso.correct)), height(pairs_ab_iso));
% fprintf('B2 correct rate: comp=%.1f%%, iso=%.1f%%\n', ...
%     100*mean(pairs_ab_comp.correct, 'omitnan'), ...
%     100*mean(pairs_ab_iso.correct, 'omitnan'));
% 
% %% Filter to pairs with B2 accuracy matched only
% pairs_ab_comp = pairs_ab_comp(~isnan(pairs_ab_comp.correct), :);
% pairs_ab_iso  = pairs_ab_iso(~isnan(pairs_ab_iso.correct), :);
% fprintf('Pairs with B2 accuracy: comp=%d, iso=%d\n', height(pairs_ab_comp), height(pairs_ab_iso));
% 
% %% Build baseline pools: other B1 encoding trials per subject per condition
% unique_subjs = unique([pairs_ab_comp.subj_id; pairs_ab_iso.subj_id]);
% subj_baseline_b_comp = struct();
% subj_baseline_b_iso  = struct();
% for s_idx = 1:length(unique_subjs)
%     sid = unique_subjs(s_idx);
%     skey = sprintf('s%d', sid);
%     b_comp = tr_one_b_b_comp(tr_one_b_b_comp.subj_id == sid, :);
%     b_iso  = tr_one_b_b_iso(tr_one_b_b_iso.subj_id == sid, :);
%     if height(b_comp) > 0, subj_baseline_b_comp.(skey) = b_comp; end
%     if height(b_iso)  > 0, subj_baseline_b_iso.(skey)  = b_iso;  end
% end
% 
% n_base_max = 40; % subsample baseline to equate across subjects
% rng(42);
% 
% %% Compute A1-B1 gaze similarity with other-B baseline correction
% fprintf('\nComputing A1-B1 gaze similarity (baseline-corrected, n_base=%d)...\n', n_base_max);
% 
% results_comp = nan(height(pairs_ab_comp), 7); % subj_id, tr_1b_a, tr_1b_b, match, baseline, reinst_index, correct
% n_comp = 0;
% for i = 1:height(pairs_ab_comp)
%     if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_ab_comp)); end
%     pair = pairs_ab_comp(i,:);
%     skey = sprintf('s%d', pair.subj_id);
%     if ~isfield(subj_baseline_b_comp, skey), continue; end
% 
%     fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
%     fix_b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_b, :);
%     if height(fix_a) < 2 || height(fix_b) < 2, continue; end
%     map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
%     map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);
%     if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end
% 
%     match_score = corr(map_a(:), map_b(:));
% 
%     % Baseline: correlate B1 map with subsampled other B1 encoding maps
%     pool = subj_baseline_b_comp.(skey);
%     baseline_pool = pool(pool.trial_id ~= pair.tr_1b_b, :);
%     if height(baseline_pool) < 1, continue; end
%     n_base = min(n_base_max, height(baseline_pool));
%     base_idx = randsample(height(baseline_pool), n_base);
%     mismatch_scores = nan(n_base, 1);
%     for j = 1:n_base
%         fix_other = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(base_idx(j)), :);
%         if height(fix_other) < 2, continue; end
%         map_other = create_fixation_map(fix_other.x, fix_other.y, fix_other.dur, spatial_params);
%         if sum(map_other(:)) == 0, continue; end
%         mismatch_scores(j) = corr(map_other(:), map_b(:));
%     end
%     baseline_score = mean(mismatch_scores, 'omitnan');
% 
%     n_comp = n_comp + 1;
%     results_comp(n_comp,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_1b_b, match_score, baseline_score, match_score - baseline_score, pair.correct];
% end
% results_comp = results_comp(1:n_comp,:);
% 
% results_iso = nan(height(pairs_ab_iso), 7);
% n_iso = 0;
% for i = 1:height(pairs_ab_iso)
%     if mod(i, 50) == 0, fprintf('  %d/%d\n', i, height(pairs_ab_iso)); end
%     pair = pairs_ab_iso(i,:);
%     skey = sprintf('s%d', pair.subj_id);
%     if ~isfield(subj_baseline_b_iso, skey), continue; end
% 
%     fix_a = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_a, :);
%     fix_b = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == pair.tr_1b_b, :);
%     if height(fix_a) < 2 || height(fix_b) < 2, continue; end
%     map_a = create_fixation_map(fix_a.x, fix_a.y, fix_a.dur, spatial_params);
%     map_b = create_fixation_map(fix_b.x, fix_b.y, fix_b.dur, spatial_params);
%     if sum(map_a(:)) == 0 || sum(map_b(:)) == 0, continue; end
% 
%     match_score = corr(map_a(:), map_b(:));
% 
%     pool = subj_baseline_b_iso.(skey);
%     baseline_pool = pool(pool.trial_id ~= pair.tr_1b_b, :);
%     if height(baseline_pool) < 1, continue; end
%     n_base = min(n_base_max, height(baseline_pool));
%     base_idx = randsample(height(baseline_pool), n_base);
%     mismatch_scores = nan(n_base, 1);
%     for j = 1:n_base
%         fix_other = Mw(task_1b_idx & Mw.subj_id == pair.subj_id & Mw.trial_id == baseline_pool.trial_id(base_idx(j)), :);
%         if height(fix_other) < 2, continue; end
%         map_other = create_fixation_map(fix_other.x, fix_other.y, fix_other.dur, spatial_params);
%         if sum(map_other(:)) == 0, continue; end
%         mismatch_scores(j) = corr(map_other(:), map_b(:));
%     end
%     baseline_score = mean(mismatch_scores, 'omitnan');
% 
%     n_iso = n_iso + 1;
%     results_iso(n_iso,:) = [pair.subj_id, pair.tr_1b_a, pair.tr_1b_b, match_score, baseline_score, match_score - baseline_score, pair.correct];
% end
% results_iso = results_iso(1:n_iso,:);
% 
% fprintf('Valid pairs: comp=%d, iso=%d\n', n_comp, n_iso);
% 
% %% Convert to tables and save for Subtitle 1 figure
% results_comp = array2table(results_comp, 'VariableNames', {'subj_id','tr_1b_a','tr_1b_b','match_score','baseline_score','reinst_index','correct'});
% results_iso  = array2table(results_iso,  'VariableNames', {'subj_id','tr_1b_a','tr_1b_b','match_score','baseline_score','reinst_index','correct'});
% 
% reinstat_res_ab = struct();
% reinstat_res_ab.ab_compared = results_comp;
% reinstat_res_ab.ab_isolated = results_iso;
% reinstat_res_ab.spatial_params = spatial_params;
% save(fullfile(res_dir, 'gaze_reinstat_res_ab.mat'), 'reinstat_res_ab');
% fprintf('Saved gaze_reinstat_res_ab.mat\n');
% 
% 
% %% Subject-level means by condition x accuracy (baseline-corrected)
% subj_ids = unique([results_comp.subj_id; results_iso.subj_id]);
% n_subj = length(subj_ids);
% sim_comp_corr = nan(n_subj, 1); sim_comp_incr = nan(n_subj, 1);
% sim_iso_corr  = nan(n_subj, 1); sim_iso_incr  = nan(n_subj, 1);
% sim_comp_all  = nan(n_subj, 1); sim_iso_all   = nan(n_subj, 1);
% 
% for s = 1:n_subj
%     sid = subj_ids(s);
%     cc = results_comp.reinst_index(results_comp.subj_id == sid & results_comp.correct == 1);
%     ci = results_comp.reinst_index(results_comp.subj_id == sid & results_comp.correct == 0);
%     ic = results_iso.reinst_index(results_iso.subj_id == sid & results_iso.correct == 1);
%     ii_t = results_iso.reinst_index(results_iso.subj_id == sid & results_iso.correct == 0);
%     ca = results_comp.reinst_index(results_comp.subj_id == sid);
%     ia = results_iso.reinst_index(results_iso.subj_id == sid);
%     if ~isempty(cc), sim_comp_corr(s) = mean(cc, 'omitnan'); end
%     if ~isempty(ci), sim_comp_incr(s) = mean(ci, 'omitnan'); end
%     if ~isempty(ic), sim_iso_corr(s)  = mean(ic, 'omitnan'); end
%     if ~isempty(ii_t), sim_iso_incr(s) = mean(ii_t, 'omitnan'); end
%     if ~isempty(ca), sim_comp_all(s) = mean(ca, 'omitnan'); end
%     if ~isempty(ia), sim_iso_all(s)  = mean(ia, 'omitnan'); end
% end
% 
% %% 2x2 ANOVA
% v = ~isnan(sim_comp_corr) & ~isnan(sim_comp_incr) & ~isnan(sim_iso_corr) & ~isnan(sim_iso_incr);
% fprintf('\nSubjects with all 4 cells: %d\n', sum(v));
% 
% [stat_results.sme_ab_2x2, stat_results.sme_ab_2x2_ph] = run_2x2_anova_correctness( ...
%     'Subsequent memory: A1-B1 encoding alignment (baseline-corrected) by B2 accuracy', ...
%     sim_comp_corr(v), sim_comp_incr(v), sim_iso_corr(v), sim_iso_incr(v), ...
%     sim_comp_all(v), sim_iso_all(v));
% 
% %% Figure: A1-B1 SME
% clear; clc; close all;
% base_dir = '..'; res_dir = fullfile(base_dir,'results'); fig_dir = fullfile(base_dir,'figures');
% c_comp = [180 174 211]/255; c_iso = [176 230 255]/255;
% load(fullfile(res_dir,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab');
% rc = reinstat_res_ab.ab_compared; ri = reinstat_res_ab.ab_isolated;
% sids = unique([rc.subj_id; ri.subj_id]); ns = length(sids);
% sm = nan(ns,4);
% for s=1:ns, sid=sids(s);
%     for j=1:4, [tbl,acc]= deal({rc,rc,ri,ri},{1,0,1,0});
%     d=tbl{j}.reinst_index(tbl{j}.subj_id==sid & tbl{j}.correct==acc{j});
%     if ~isempty(d), sm(s,j)=mean(d,'omitnan'); end, end, end
% v=all(~isnan(sm),2); sm=sm(v,:);
% fprintf('A1-B1 SME n=%d\n', sum(v));
% sim_comp_all = mean(sm(:,1:2), 2, 'omitnan');
% sim_iso_all  = mean(sm(:,3:4), 2, 'omitnan');
% [stat_ab.sme_2x2, stat_ab.sme_2x2_ph] = run_2x2_anova_correctness('A1-B1', sm(:,1), sm(:,2), sm(:,3), sm(:,4), sim_comp_all, sim_iso_all);
% cols={c_comp,c_comp*.5+.5,c_iso,c_iso*.5+.5}; xp=[1 2 3 4];
% figure('color','w','position',[100 100 450 400]);  hold on;
% jt=-0.15-rand(size(sm))*.2;
% for s=1:size(sm,1), for g=[1 3], plot(xp(g:g+1)+jt(s,g:g+1),sm(s,g:g+1),'-','Color',[.7 .7 .7 .4],'LineWidth',.5); end, end
% for i=1:4, d=sm(:,i); x=xp(i); [f,xi]=ksdensity(d); f=f/max(f)*.55;
%     patch([x+f,x*ones(1,length(f))],[xi,fliplr(xi)],cols{i},'EdgeColor','none','FaceAlpha',.5);
%     scatter(x+jt(:,i),d,35,cols{i},'filled','MarkerFaceAlpha',.7);
%     q=quantile(d,[.25 .5 .75]); rectangle('Position',[x-.1,q(1),.2,q(3)-q(1)],'FaceColor',cols{i},'EdgeColor','k','LineWidth',1.2);
%     plot([x-.1,x+.1],[q(2) q(2)],'k-','LineWidth',2.5); end
% yline(0,'r--','LineWidth',2);
% cl=ylim; st=(cl(2)-cl(1))*.08; bs=max(sm(:))+st*.5; lv=0;
% for g=[1 3], [~,p]=ttest(sm(:,g),sm(:,g+1),'Tail','left'); lv=lv+1; yp=bs+lv*st;
%     if p<.001,tx='***';elseif p<.01,tx='**';elseif p<.05,tx='*';else,tx='n.s.';end
%     plot([xp(g) xp(g) xp(g+1) xp(g+1)],[yp-st*.3,yp,yp,yp-st*.3],'k-','LineWidth',1.2);
%     text(mean(xp(g:g+1)),yp+st*.1,tx,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
% ylim([cl(1),bs+lv*st+st]);
% set(gca,'XTick',xp,'XTickLabel',{'Correct','Incorrect','Correct','Incorrect'},'FontSize',20);
% ylabel('Gaze Similarity(r)','FontSize',20); xlim([.3 4.7]); box off;
% print(gcf,fullfile(fig_dir,'sme_a1b1_v2.pdf'),'-dpdf','-vector');

%% Figure: A1-A2 SME
load(fullfile(res_dir,'gaze_reinstat_res_aa.mat'),'reinstat_res_aa');
rc=reinstat_res_aa.aa_compared; ri=reinstat_res_aa.aa_isolated;
sids=unique([rc.subj_id;ri.subj_id]); ns=length(sids); sm=nan(ns,4);
for s=1:ns, sid=sids(s);
    for j=1:4, [tbl,acc]=deal({rc,rc,ri,ri},{1,0,1,0});
    d=tbl{j}.reinst_index(tbl{j}.subj_id==sid & tbl{j}.correct==acc{j});
    if ~isempty(d), sm(s,j)=mean(d,'omitnan'); end, end, end
v=all(~isnan(sm),2); sm=sm(v,:);
fprintf('A1-A2 SME n=%d\n',sum(v));

cols={c_comp,c_comp*.5+.5,c_iso,c_iso*.5+.5}; xp=[1 2 3.5 4.5];
yl=[-0.2 0.5];

figure('color','w','position',[100 100 500 500]); hold on;
jt=-0.15-rand(size(sm))*.2;
for s=1:size(sm,1), for g=[1 3], plot(xp(g:g+1)+jt(s,g:g+1),sm(s,g:g+1),'-','Color',[.7 .7 .7 .4],'LineWidth',.5); end, end
for i=1:4, d=sm(:,i); x=xp(i); [f,xi]=ksdensity(d); f=f/max(f)*.55;
    patch([x+f,x*ones(1,length(f))],[xi,fliplr(xi)],cols{i},'EdgeColor','none','FaceAlpha',.5);
    scatter(x+jt(:,i),d,35,cols{i},'filled','MarkerFaceAlpha',.7);
    q=quantile(d,[.25 .5 .75]); rectangle('Position',[x-.1,q(1),.2,q(3)-q(1)],'FaceColor',cols{i},'EdgeColor','k','LineWidth',1.2);
    plot([x-.1,x+.1],[q(2) q(2)],'k-','LineWidth',2.5); end
yline(0,'r--','LineWidth',2); ylim(yl); set(gca,'YTick',[-0.2 0 0.2 0.4]);
st=diff(yl)*.06; lv=0;
for g=[1 3], [~,p]=ttest(sm(:,g),sm(:,g+1)); lv=lv+1; yp=yl(2)-st*(2-lv);
    if p<.001,tx='***';elseif p<.01,tx='**';elseif p<.05,tx='*';else,tx='n.s.';end
    plot([xp(g) xp(g) xp(g+1) xp(g+1)],[yp-st*.3,yp,yp,yp-st*.3],'k-','LineWidth',1.2);
    text(mean(xp(g:g+1)),yp+st*.15,tx,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
set(gca,'XTick',xp,'XTickLabel',{'Correct','Incorrect','Correct','Incorrect'},'FontSize',20);
ylabel('Gaze Similarity (r)','FontSize',20); xlim([.3 5.2]); box off;
print(gcf,fullfile(fig_dir,'sme_a1a2.pdf'),'-dpdf','-vector');

t=table(sm(:,1),sm(:,2),sm(:,3),sm(:,4),'VariableNames',{'cc','ci','ic','ii'});
within=table({'comp';'comp';'iso';'iso'},{'corr';'incorr';'corr';'incorr'},'VariableNames',{'cond','correctness'});
rm=fitrm(t,'cc-ii~1','WithinDesign',within); tbl=ranova(rm,'WithinModel','cond*correctness');
tn=string(tbl.Properties.RowNames); rows={'(Intercept):cond','(Intercept):correctness','(Intercept):cond:correctness'};
labels={'condition','correctness','interaction'};
for i=1:3, idx=find(contains(tn,rows{i}),1); if ~isempty(idx)
    eta=tbl.SumSq(idx)/(tbl.SumSq(idx)+tbl.SumSq(idx+1));
    fprintf('  %s: F(%d,%d)=%.2f, p=%.4f, eta2=%.3f %s\n',labels{i},tbl.DF(idx),tbl.DF(idx+1),tbl.F(idx),tbl.pValue(idx),eta,...
        repmat('*',1,(tbl.pValue(idx)<.05)+(tbl.pValue(idx)<.01)+(tbl.pValue(idx)<.001))); end, end
nv=size(sm,1);
[~,p1,~,s1]=ttest(sm(:,1),sm(:,2)); [~,p2,~,s2]=ttest(sm(:,3),sm(:,4));
fprintf('  comp corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f\n',s1.df,s1.tstat,p1,s1.tstat/sqrt(nv));
fprintf('  iso  corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f\n',s2.df,s2.tstat,p2,s2.tstat/sqrt(nv));
fprintf('  comp corr:  M=%.4f SD=%.4f\n  comp incorr: M=%.4f SD=%.4f\n  iso  corr:  M=%.4f SD=%.4f\n  iso  incorr: M=%.4f SD=%.4f\n',...
    mean(sm(:,1)),std(sm(:,1)),mean(sm(:,2)),std(sm(:,2)),mean(sm(:,3)),std(sm(:,3)),mean(sm(:,4)),std(sm(:,4)));
% 
% %% 2x2 ANOVA + post-hocs
% fprintf('\n=== A1-A2 reinstatement by B2 accuracy ===\n');
% t=table(sm(:,1),sm(:,2),sm(:,3),sm(:,4),'VariableNames',{'cc','ci','ic','ii'});
% within=table({'comp';'comp';'iso';'iso'},{'corr';'incorr';'corr';'incorr'},'VariableNames',{'cond','correctness'});
% rm=fitrm(t,'cc-ii~1','WithinDesign',within);
% tbl=ranova(rm,'WithinModel','cond*correctness');
% rows={'(Intercept):cond','(Intercept):correctness','(Intercept):cond:correctness'};
% labels={'condition','correctness','interaction'};
% tn=string(tbl.Properties.RowNames);
% for i=1:3, idx=find(contains(tn,rows{i}),1);
%     if ~isempty(idx), f=tbl.F(idx); pv=tbl.pValue(idx); df1=tbl.DF(idx); df2=tbl.DF(idx+1);
%     eta=tbl.SumSq(idx)/(tbl.SumSq(idx)+tbl.SumSq(idx+1));
%     fprintf('  %s: F(%d,%d)=%.2f, p=%.4f, eta2=%.3f %s\n',labels{i},df1,df2,f,pv,eta,...
%         repmat('*',1,(pv<.05)+(pv<.01)+(pv<.001))); end, end
% nv=size(sm,1);
% [~,p1,~,s1]=ttest(sm(:,1),sm(:,2)); [~,p2,~,s2]=ttest(sm(:,3),sm(:,4));
% [~,pc,~,sc]=ttest(mean(sm(:,1:2),2),mean(sm(:,3:4),2));
% [~,p3,~,s3]=ttest(sm(:,1)-sm(:,2),sm(:,3)-sm(:,4));
% fprintf('  comp corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s1.df,s1.tstat,p1,s1.tstat/sqrt(nv),repmat('*',1,(p1<.05)+(p1<.01)+(p1<.001)));
% fprintf('  iso corr>incorr:  t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s2.df,s2.tstat,p2,s2.tstat/sqrt(nv),repmat('*',1,(p2<.05)+(p2<.01)+(p2<.001)));
% fprintf('  comp>iso (all):   t(%d)=%.2f, p=%.4f, d=%.2f %s\n',sc.df,sc.tstat,pc,sc.tstat/sqrt(nv),repmat('*',1,(pc<.05)+(pc<.01)+(pc<.001)));
% fprintf('  interaction:      t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s3.df,s3.tstat,p3,s3.tstat/sqrt(nv),repmat('*',1,(p3<.05)+(p3<.01)+(p3<.001)));
% fprintf('  comp corr: M=%.4f SD=%.4f\n',mean(sm(:,1)),std(sm(:,1)));
% fprintf('  comp incorr: M=%.4f SD=%.4f\n',mean(sm(:,2)),std(sm(:,2)));
% fprintf('  iso corr: M=%.4f SD=%.4f\n',mean(sm(:,3)),std(sm(:,3)));
% fprintf('  iso incorr: M=%.4f SD=%.4f\n',mean(sm(:,4)),std(sm(:,4)));

%% Helper function
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

function [res, posthoc] = run_2x2_anova_correctness(lbl, cc, ci, ic, ii, cond_comp, cond_iso)
    fprintf('\n%s:\n', lbl);
    t = table(cc, ci, ic, ii, 'VariableNames', {'cc','ci','ic','ii'});
    within = table({'comp';'comp';'iso';'iso'}, {'corr';'incorr';'corr';'incorr'}, 'VariableNames', {'cond','correctness'});
    rm = fitrm(t, 'cc-ii~1', 'WithinDesign', within);
    tbl = ranova(rm, 'WithinModel', 'cond*correctness'); eps = epsilon(rm); m = mauchly(rm);
    rows = {'(Intercept):cond', '(Intercept):correctness', '(Intercept):cond:correctness'};
    labels = {'condition', 'correctness', 'interaction'};
    term_names = string(tbl.Properties.RowNames);
    for i = 1:3
        idx = find(contains(term_names, rows{i}), 1);
        if ~isempty(idx)
            f = tbl.F(idx); p_gg = tbl.pValueGG(idx); df1 = tbl.DF(idx); df2_gg = tbl.DF(idx+1) * eps.GreenhouseGeisser;
            mauch_p = m.pValue(1); ep = eps.GreenhouseGeisser; s_ef = tbl.SumSq(idx); s_er = tbl.SumSq(idx+1); eta = s_ef/(s_ef+s_er);
            fprintf('  %s: F(%d,%.1f)=%.2f, p=%.4f (GG), mauchly p=%.3f, eps=%.3f, eta=%.2f %s\n', ...
                labels{i}, df1, df2_gg, f, p_gg, mauch_p, ep, eta, repmat('*',1,(p_gg<0.05)+(p_gg<0.01)+(p_gg<0.001)));
            if i==1, res.cond.F=f; res.cond.p=p_gg; res.cond.eta=eta; res.cond.mauchly_p=mauch_p; res.cond.eps=ep; end
            if i==2, res.corr.F=f; res.corr.p=p_gg; res.corr.eta=eta; res.corr.mauchly_p=mauch_p; res.corr.eps=ep; end
            if i==3, res.int.F=f; res.int.p=p_gg; res.int.eta=eta; res.int.mauchly_p=mauch_p; res.int.eps=ep; end
        end
    end
    fprintf('  condition post-hoc:\n');
    [~,pc,~,sc] = ttest(cond_comp, cond_iso);
    posthoc.cond_comp_iso.t=sc.tstat; posthoc.cond_comp_iso.p=pc; posthoc.cond_comp_iso.d=sc.tstat/sqrt(sc.df+1);
    fprintf('    comp>iso: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', sc.df, sc.tstat, posthoc.cond_comp_iso.d, pc, repmat('*',1,(pc<0.05)+(pc<0.01)+(pc<0.001)));
    fprintf('  correctness post-hoc:\n');
    [~,p1,~,s1] = ttest(cc, ci, 'Tail','left'); [~,p2,~,s2] = ttest(ic, ii, 'Tail','left'); [~,p3,~,s3] = ttest(cc-ci, ic-ii);
    df = s1.df;
    posthoc.comp_corr_incorr.t=s1.tstat; posthoc.comp_corr_incorr.p=p1; posthoc.comp_corr_incorr.d=s1.tstat/sqrt(df+1);
    posthoc.iso_corr_incorr.t=s2.tstat; posthoc.iso_corr_incorr.p=p2; posthoc.iso_corr_incorr.d=s2.tstat/sqrt(df+1);
    posthoc.comp_vs_iso.t=s3.tstat; posthoc.comp_vs_iso.p=p3; posthoc.comp_vs_iso.d=s3.tstat/sqrt(df+1);
    fprintf('    comp corr>incorr: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s1.tstat, posthoc.comp_corr_incorr.d, p1, repmat('*',1,(p1<0.05)+(p1<0.01)+(p1<0.001)));
    fprintf('    iso corr>incorr: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s2.tstat, posthoc.iso_corr_incorr.d, p2, repmat('*',1,(p2<0.05)+(p2<0.01)+(p2<0.01)+(p2<0.001)));
    fprintf('    (comp-iso) diff: t(%d)=%.2f, d=%.2f, p=%.4f %s\n', df, s3.tstat, posthoc.comp_vs_iso.d, p3, repmat('*',1,(p3<0.05)+(p3<0.01)+(p3<0.001)));
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