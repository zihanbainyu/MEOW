% %% A2-B2 retrieval gaze similarity with other-B1 baseline correction
clear; clc; close all;
base_dir='..'; res_dir=fullfile(base_dir,'results'); fig_dir=fullfile(base_dir,'figures');
c_comp=[180 174 211]/255; c_iso=[176 230 255]/255;
% 
% load(fullfile(res_dir,'group_eye_movement_combined.mat'),'Mw');
% spatial_params=struct('xres',1920,'yres',1080,'roi_x',[760,1160],'roi_y',[340,740],'grid_x',20,'grid_y',20,'sigma',2);
% 
% comp_idx=strcmp(Mw.condition,'compared'); iso_idx=strcmp(Mw.condition,'isolated');
% ab_idx=strcmp(Mw.goal,'A-B'); a_idx=strcmp(Mw.identity,'A'); b_idx=strcmp(Mw.identity,'B');
% task_1b_idx=strcmp(Mw.task,'1_back'); task_2b_idx=strcmp(Mw.task,'2_back');
% 
% tr_2b_a_comp=unique(Mw(task_2b_idx&ab_idx&a_idx&comp_idx,{'subj_id','trial_id','stim_id'}));
% tr_2b_a_iso=unique(Mw(task_2b_idx&ab_idx&a_idx&iso_idx,{'subj_id','trial_id','stim_id'}));
% tr_2b_b_comp=unique(Mw(task_2b_idx&ab_idx&b_idx&comp_idx,{'subj_id','trial_id','stim_id'}));
% tr_2b_b_iso=unique(Mw(task_2b_idx&ab_idx&b_idx&iso_idx,{'subj_id','trial_id','stim_id'}));
% tr_1b_b_comp=unique(Mw(task_1b_idx&b_idx&comp_idx,{'subj_id','trial_id','stim_id'}));
% tr_1b_b_iso=unique(Mw(task_1b_idx&b_idx&iso_idx,{'subj_id','trial_id','stim_id'}));
% 
% tr_2b_a_comp.base_id=regexprep(tr_2b_a_comp.stim_id,'_A_','_BASE_');
% tr_2b_a_iso.base_id=regexprep(tr_2b_a_iso.stim_id,'_A_','_BASE_');
% tr_2b_b_comp.base_id=regexprep(tr_2b_b_comp.stim_id,'_B_','_BASE_');
% tr_2b_b_iso.base_id=regexprep(tr_2b_b_iso.stim_id,'_B_','_BASE_');
% 
% pairs_comp=innerjoin(tr_2b_a_comp,tr_2b_b_comp,'Keys',{'subj_id','base_id'},...
%     'LeftVariables',{'subj_id','trial_id','stim_id','base_id'},'RightVariables',{'trial_id','stim_id'});
% pairs_comp.Properties.VariableNames={'subj_id','tr_2b_a','stim_id_a','base_id','tr_2b_b','stim_id_b'};
% pairs_iso=innerjoin(tr_2b_a_iso,tr_2b_b_iso,'Keys',{'subj_id','base_id'},...
%     'LeftVariables',{'subj_id','trial_id','stim_id','base_id'},'RightVariables',{'trial_id','stim_id'});
% pairs_iso.Properties.VariableNames={'subj_id','tr_2b_a','stim_id_a','base_id','tr_2b_b','stim_id_b'};
% fprintf('A2-B2 pairs: comp=%d, iso=%d\n',height(pairs_comp),height(pairs_iso));
% 
% %% B2 correctness
% pairs_comp.correct=nan(height(pairs_comp),1);
% for i=1:height(pairs_comp)
%     row=Mw(Mw.subj_id==pairs_comp.subj_id(i)&Mw.trial_id==pairs_comp.tr_2b_b(i)&task_2b_idx&b_idx,:);
%     if height(row)>=1, pairs_comp.correct(i)=row.correct(1); end
% end
% pairs_iso.correct=nan(height(pairs_iso),1);
% for i=1:height(pairs_iso)
%     row=Mw(Mw.subj_id==pairs_iso.subj_id(i)&Mw.trial_id==pairs_iso.tr_2b_b(i)&task_2b_idx&b_idx,:);
%     if height(row)>=1, pairs_iso.correct(i)=row.correct(1); end
% end
% fprintf('B2 correctness: comp=%d/%d, iso=%d/%d\n',sum(~isnan(pairs_comp.correct)),height(pairs_comp),sum(~isnan(pairs_iso.correct)),height(pairs_iso));
% 
% %% Baseline pools: other B1 encoding
% unique_subjs=unique([pairs_comp.subj_id;pairs_iso.subj_id]);
% bl_comp=struct(); bl_iso=struct();
% for s=1:length(unique_subjs), sid=unique_subjs(s); sk=sprintf('s%d',sid);
%     bc=tr_1b_b_comp(tr_1b_b_comp.subj_id==sid,:); bi=tr_1b_b_iso(tr_1b_b_iso.subj_id==sid,:);
%     if height(bc)>0, bl_comp.(sk)=bc; end; if height(bi)>0, bl_iso.(sk)=bi; end
% end
% 
% %% Compute A2-B2 compared
% fprintf('\nComputing A2-B2 compared (%d pairs)...\n',height(pairs_comp));
% res_comp=nan(height(pairs_comp),7); nc=0;
% for i=1:height(pairs_comp)
%     if mod(i,50)==0, fprintf('  %d/%d\n',i,height(pairs_comp)); end
%     p=pairs_comp(i,:); sk=sprintf('s%d',p.subj_id);
%     if ~isfield(bl_comp,sk), continue; end
%     fa=Mw(task_2b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==p.tr_2b_a,:);
%     fb=Mw(task_2b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==p.tr_2b_b,:);
%     if height(fa)<2||height(fb)<2, continue; end
%     ma=create_fixation_map(fa.x,fa.y,fa.dur,spatial_params);
%     mb=create_fixation_map(fb.x,fb.y,fb.dur,spatial_params);
%     if sum(ma(:))==0||sum(mb(:))==0, continue; end
%     match=corr(ma(:),mb(:));
%     pool=bl_comp.(sk); mm=nan(height(pool),1);
%     for j=1:height(pool)
%         fo=Mw(task_1b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==pool.trial_id(j),:);
%         if height(fo)<2, continue; end
%         mo=create_fixation_map(fo.x,fo.y,fo.dur,spatial_params);
%         if sum(mo(:))>0, mm(j)=corr(mo(:),mb(:)); end
%     end
%     nc=nc+1; res_comp(nc,:)=[p.subj_id,p.tr_2b_a,p.tr_2b_b,match,mean(mm,'omitnan'),match-mean(mm,'omitnan'),p.correct];
% end
% res_comp=res_comp(1:nc,:); fprintf('A2-B2 compared: %d valid\n',nc);
% 
% %% Compute A2-B2 isolated
% fprintf('Computing A2-B2 isolated (%d pairs)...\n',height(pairs_iso));
% res_iso=nan(height(pairs_iso),7); ni=0;
% for i=1:height(pairs_iso)
%     if mod(i,50)==0, fprintf('  %d/%d\n',i,height(pairs_iso)); end
%     p=pairs_iso(i,:); sk=sprintf('s%d',p.subj_id);
%     if ~isfield(bl_iso,sk), continue; end
%     fa=Mw(task_2b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==p.tr_2b_a,:);
%     fb=Mw(task_2b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==p.tr_2b_b,:);
%     if height(fa)<2||height(fb)<2, continue; end
%     ma=create_fixation_map(fa.x,fa.y,fa.dur,spatial_params);
%     mb=create_fixation_map(fb.x,fb.y,fb.dur,spatial_params);
%     if sum(ma(:))==0||sum(mb(:))==0, continue; end
%     match=corr(ma(:),mb(:));
%     pool=bl_iso.(sk); mm=nan(height(pool),1);
%     for j=1:height(pool)
%         fo=Mw(task_1b_idx&Mw.subj_id==p.subj_id&Mw.trial_id==pool.trial_id(j),:);
%         if height(fo)<2, continue; end
%         mo=create_fixation_map(fo.x,fo.y,fo.dur,spatial_params);
%         if sum(mo(:))>0, mm(j)=corr(mo(:),mb(:)); end
%     end
%     ni=ni+1; res_iso(ni,:)=[p.subj_id,p.tr_2b_a,p.tr_2b_b,match,mean(mm,'omitnan'),match-mean(mm,'omitnan'),p.correct];
% end
% res_iso=res_iso(1:ni,:); fprintf('A2-B2 isolated: %d valid\n',ni);
% 
% %% Save
% res_comp=array2table(res_comp,'VariableNames',{'subj_id','tr_2b_a','tr_2b_b','match_score','baseline_score','reinst_index','correct'});
% res_iso=array2table(res_iso,'VariableNames',{'subj_id','tr_2b_a','tr_2b_b','match_score','baseline_score','reinst_index','correct'});
% reinstat_res_a2b2=struct('a2b2_compared',res_comp,'a2b2_isolated',res_iso,'spatial_params',spatial_params);
% save(fullfile(res_dir,'gaze_reinstat_res_a2b2.mat'),'reinstat_res_a2b2');
% fprintf('Saved gaze_reinstat_res_a2b2.mat\n');

base_dir = '..';
res_dir  = fullfile(base_dir, 'results');
fig_dir  = fullfile(base_dir, 'figures');

load(fullfile(res_dir, 'gaze_reinstat_res_a2b2.mat'), 'reinstat_res_a2b2');

res_comp = reinstat_res_a2b2.a2b2_compared;
res_iso  = reinstat_res_a2b2.a2b2_isolated;

% exc = [609, 606, 608, 618];
exc = [];
c_comp = [180 174 211]/255; c_iso = [176 230 255]/255;

%% SME: 2x2 ANOVA + figure
sids=unique([res_comp.subj_id;res_iso.subj_id]); ns=length(sids); sm=nan(ns,4);
for s=1:ns, sid=sids(s);
    for j=1:4, [tbl,acc]=deal({res_comp,res_comp,res_iso,res_iso},{1,0,1,0});
    d=tbl{j}.reinst_index(tbl{j}.subj_id==sid&tbl{j}.correct==acc{j});
    if ~isempty(d), sm(s,j)=mean(d,'omitnan'); end, end, end
v=all(~isnan(sm),2); sm=sm(v,:); fprintf('A2-B2 SME n=%d\n',sum(v));

fprintf('\n=== A2-B2 by B2 accuracy ===\n');
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
[~,pc,~,sc]=ttest(mean(sm(:,1:2),2),mean(sm(:,3:4),2));
fprintf('  comp corr>incorr: t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s1.df,s1.tstat,p1,s1.tstat/sqrt(nv),repmat('*',1,(p1<.05)+(p1<.01)+(p1<.001)));
fprintf('  iso corr>incorr:  t(%d)=%.2f, p=%.4f, d=%.2f %s\n',s2.df,s2.tstat,p2,s2.tstat/sqrt(nv),repmat('*',1,(p2<.05)+(p2<.01)+(p2<.001)));
fprintf('  comp>iso:         t(%d)=%.2f, p=%.4f, d=%.2f %s\n',sc.df,sc.tstat,pc,sc.tstat/sqrt(nv),repmat('*',1,(pc<.05)+(pc<.01)+(pc<.001)));
fprintf('  comp corr: M=%.4f SD=%.4f\n  comp incorr: M=%.4f SD=%.4f\n  iso corr: M=%.4f SD=%.4f\n  iso incorr: M=%.4f SD=%.4f\n',...
    mean(sm(:,1)),std(sm(:,1)),mean(sm(:,2)),std(sm(:,2)),mean(sm(:,3)),std(sm(:,3)),mean(sm(:,4)),std(sm(:,4)));

cols={c_comp,c_comp*.5+.5,c_iso,c_iso*.5+.5}; xp=[1 2 3.5 4.5]; yl=[-0.2 0.5];
figure('color','w','position',[100 100 500 500]); hold on;
jt=-0.15-rand(size(sm))*.2;
for s=1:size(sm,1), for g=[1 3], plot(xp(g:g+1)+jt(s,g:g+1),sm(s,g:g+1),'-','Color',[.7 .7 .7 .4],'LineWidth',.5); end, end
for i=1:4, d=sm(:,i); d_c=d(~isnan(d)); if isempty(d_c), continue; end
    x=xp(i); [f,xi]=ksdensity(d_c); f=f/max(f)*.55;
    patch([x+f,x*ones(1,length(f))],[xi,fliplr(xi)],cols{i},'EdgeColor','none','FaceAlpha',.5);
    scatter(x+jt(~isnan(d),i),d_c,35,cols{i},'filled','MarkerFaceAlpha',.7);
    q=quantile(d_c,[.25 .5 .75]); rectangle('Position',[x-.1,q(1),.2,q(3)-q(1)],'FaceColor',cols{i},'EdgeColor','k','LineWidth',1.2);
    plot([x-.1,x+.1],[q(2) q(2)],'k-','LineWidth',2.5); end
yline(0,'r--','LineWidth',2); ylim(yl); set(gca,'YTick',[-.2 0 .2 .4]);
st=diff(yl)*.06; lv=0;
for g=[1 3], [~,p]=ttest(sm(:,g),sm(:,g+1)); lv=lv+1; yp=yl(2)-st*(2-lv);
    if p<.001,tx='***';elseif p<.01,tx='**';elseif p<.05,tx='*';else,tx='n.s.';end
    plot([xp(g) xp(g) xp(g+1) xp(g+1)],[yp-st*.3,yp,yp,yp-st*.3],'k-','LineWidth',1.2);
    text(mean(xp(g:g+1)),yp+st*.15,tx,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold'); end
set(gca,'XTick',xp,'XTickLabel',{'Correct','Incorrect','Correct','Incorrect'},'FontSize',20);
ylabel('Gaze Similarity (r)','FontSize',20); xlim([.3 5.2]); box off;
print(gcf,fullfile(fig_dir,'sme_a2b2.pdf'),'-dpdf','-vector');

function map=create_fixation_map(x,y,dur,params)
    in=x>=params.roi_x(1)&x<=params.roi_x(2)&y>=params.roi_y(1)&y<=params.roi_y(2);
    x=x(in);y=y(in);dur=dur(in);
    if isempty(x), map=zeros(params.grid_y,params.grid_x); return; end
    xb=linspace(params.roi_x(1),params.roi_x(2),params.grid_x+1);
    yb=linspace(params.roi_y(1),params.roi_y(2),params.grid_y+1);
    map=zeros(params.grid_y,params.grid_x);
    for i=1:length(x), xi=find(x(i)>=xb(1:end-1)&x(i)<xb(2:end),1); yi=find(y(i)>=yb(1:end-1)&y(i)<yb(2:end),1);
        if ~isempty(xi)&&~isempty(yi), map(yi,xi)=map(yi,xi)+dur(i); end, end
    if params.sigma>0, map=imgaussfilt(map,params.sigma); end
    if sum(map(:))>0, map=map/sum(map(:)); end
end