%% Figure: Same-item vs cross-item gaze similarity
clear; clc; close all;
base_dir='..'; res_dir=fullfile(base_dir,'results'); fig_dir=fullfile(base_dir,'figures');

load(fullfile(res_dir,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab');
load(fullfile(res_dir,'gaze_reinstat_res_a1b2.mat'),'reinstat_res_a1b2');
load(fullfile(res_dir,'gaze_reinstat_res_a2b2.mat'),'reinstat_res_a2b2');
load(fullfile(res_dir,'gaze_reinstat_res_aa.mat'),'reinstat_res_aa');
load(fullfile(res_dir,'gaze_reinstat_res_full.mat'),'reinstat_res');
load(fullfile(res_dir,'gaze_same_item_within.mat'),'same_item_res');

ab_all=[reinstat_res_ab.ab_compared;reinstat_res_ab.ab_isolated];
a1b2_all=[reinstat_res_a1b2.a1b2_compared;reinstat_res_a1b2.a1b2_isolated];
a2b2_all=[reinstat_res_a2b2.a2b2_compared;reinstat_res_a2b2.a2b2_isolated];
ba_all=[reinstat_res.ba_compared;reinstat_res.ba_isolated];
aa_all=[reinstat_res_aa.aa_compared;reinstat_res_aa.aa_isolated];
bb_all=[reinstat_res.bb_compared;reinstat_res.bb_isolated];
same_within_all=[same_item_res.a1a1;same_item_res.a2a2];

sids=unique([ab_all.subj_id;a1b2_all.subj_id;a2b2_all.subj_id;ba_all.subj_id;aa_all.subj_id;bb_all.subj_id;same_within_all.subj_id]);
ns=length(sids);

% Per-subject means for each pair type
sm=nan(ns,7);
all_data={ab_all,a1b2_all,a2b2_all,ba_all,aa_all,bb_all,same_within_all};
for s=1:ns, sid=sids(s); for j=1:7
    d=all_data{j}.reinst_index(all_data{j}.subj_id==sid);
    if ~isempty(d), sm(s,j)=mean(d,'omitnan'); end
end, end

v=all(~isnan(sm),2); sm=sm(v,:); n=sum(v);
fprintf('N=%d\n',n);

% Same-item: A1A2 (col5), B1B2 (col6), plus same-response A1A1/A2A2 trials (col7)
% Similar-item: A1B1 (col1), A1B2 (col2), A2B2 (col3), B1A2 (col4)
subj_same_item=mean(sm(:,[5 6 7]),2);
subj_sim_item=mean(sm(:,[1 2 3 4]),2);

[~,p_same,~,s_same]=ttest(subj_same_item);
[~,p_sim,~,s_sim]=ttest(subj_sim_item);
fprintf('\nSame-item:    M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_same_item),std(subj_same_item),s_same.df,s_same.tstat,p_same,mean(subj_same_item)/std(subj_same_item));
fprintf('Similar-item: M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_sim_item),std(subj_sim_item),s_sim.df,s_sim.tstat,p_sim,mean(subj_sim_item)/std(subj_sim_item));
[~,p_same_sim,~,s_same_sim]=ttest(subj_same_item,subj_sim_item);
fprintf('Same-item vs similar-item: t(%d)=%.2f, p=%.4f\n',s_same_sim.df,s_same_sim.tstat,p_same_sim);

n_perm=1000; rng(42);
subj_pooled=mean([subj_same_item,subj_sim_item],2);
null_means=nan(n_perm,1);
for p=1:n_perm, null_means(p)=mean(subj_pooled.*((rand(n,1)>.5)*2-1)); end

pp=nan(2,1);
for i=1:2, d={subj_same_item,subj_sim_item}; om=mean(d{i});
    pd=nan(n_perm,1); for p=1:n_perm, pd(p)=mean(d{i}.*((rand(n,1)>.5)*2-1)); end
    pp(i)=mean(abs(pd)>=abs(om)); end
fprintf('Perm: same-item=%.4f, similar-item=%.4f\n',pp(1),pp(2));

data_plot={subj_same_item,subj_sim_item};
plot_pp=pp;
labels={'Same-item','Similar-item','Null'};
c_obs=[.03 .03 .03]; c_box=[.30 .30 .30]; c_null=[.82 .82 .82];
xp=[1 1.82]; x_null=2.64;

figure('color','w','position',[100 100 620 560]); hold on;
for i=1:2
    d=data_plot{i}; x=xp(i);
    [f,xi]=ksdensity(d); f=f/max(f)*.24;
    patch([x+f,x-fliplr(f)],[xi,fliplr(xi)],c_obs,'EdgeColor','none','FaceAlpha',.68);
    jit=(rand(size(d))-.5)*.10;
    scatter(x+jit,d,34,c_obs,'filled','MarkerFaceAlpha',.82);
    q=quantile(d,[.25 .5 .75]); rectangle('Position',[x-.055,q(1),.11,q(3)-q(1)],'FaceColor',c_box,'EdgeColor','k','LineWidth',1.7);
    plot([x-.055,x+.055],[q(2) q(2)],'w-','LineWidth',2.6);
    if plot_pp(i)<.001,tx='***';elseif plot_pp(i)<.01,tx='**';elseif plot_pp(i)<.05,tx='*';else,tx='n.s.';end
    text(x,max(d)+range([min(d);max(d)])*.15,tx,'HorizontalAlignment','center','FontSize',18,'FontWeight','bold');
end

[f,xi]=ksdensity(null_means); f=f/max(f)*.22;
patch([x_null+f,x_null-fliplr(f)],[xi,fliplr(xi)],c_null,'EdgeColor','none','FaceAlpha',.6);
q=quantile(null_means,[.25 .5 .75]); rectangle('Position',[x_null-.05,q(1),.10,q(3)-q(1)],'FaceColor',c_null,'EdgeColor','k','LineWidth',1.2);
plot([x_null-.05,x_null+.05],[q(2) q(2)],'k-','LineWidth',2);

yline(0,'k--','LineWidth',1.5);
set(gca,'XTick',[xp,x_null],'XTickLabel',labels,'FontSize',16);
yl=ylim; set(gca,'YTick',linspace(yl(1),yl(2),4));
ylabel('Gaze similarity (baseline-corrected r)','FontSize',18);
xlim([.55 3.1]); box off;
print(gcf,fullfile(fig_dir,'fig_gaze_tracks_similarity_3cond.pdf'),'-dpdf','-vector');