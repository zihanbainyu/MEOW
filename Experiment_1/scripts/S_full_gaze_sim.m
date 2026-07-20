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

% Same-item: mean of A1A2 (col5) and B1B2 (col6)
% Cross-item: mean of A1B1 (col1), A1B2 (col2), A2B2 (col3), B1A2 (col4)
% Same-item across tasks: A1A2 (col5), B1B2 (col6)
% Similar-item across tasks: A1B2 (col2), B1A2 (col4)
% Similar-item within tasks: A1B1 (col1), A2B2 (col3)
% Same-item within tasks: A1A1 and A2A2 same-response trials (col7)
subj_same_across=mean(sm(:,[5 6]),2);
subj_sim_across=mean(sm(:,[2 4]),2);
subj_sim_within=mean(sm(:,[1 3]),2);
subj_same_within=sm(:,7);

[~,p1,~,s1]=ttest(subj_same_across); [~,p2,~,s2]=ttest(subj_sim_across); [~,p3,~,s3]=ttest(subj_sim_within); [~,p4,~,s4]=ttest(subj_same_within);
fprintf('\nSame-item across:    M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_same_across),std(subj_same_across),s1.df,s1.tstat,p1,mean(subj_same_across)/std(subj_same_across));
fprintf('Similar-item across: M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_sim_across),std(subj_sim_across),s2.df,s2.tstat,p2,mean(subj_sim_across)/std(subj_sim_across));
fprintf('Similar-item within: M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_sim_within),std(subj_sim_within),s3.df,s3.tstat,p3,mean(subj_sim_within)/std(subj_sim_within));
fprintf('Same-item within:    M=%.4f SD=%.4f, t(%d)=%.2f, p=%.4f, d=%.2f\n',mean(subj_same_within),std(subj_same_within),s4.df,s4.tstat,p4,mean(subj_same_within)/std(subj_same_within));
[~,p12,~,s12]=ttest(subj_same_across,subj_sim_across); [~,p13,~,s13]=ttest(subj_same_across,subj_sim_within); [~,p14,~,s14]=ttest(subj_same_across,subj_same_within); [~,p23,~,s23]=ttest(subj_sim_across,subj_sim_within); [~,p34,~,s34]=ttest(subj_sim_within,subj_same_within);
fprintf('Same vs sim-across:   t(%d)=%.2f, p=%.4f\n',s12.df,s12.tstat,p12);
fprintf('Same vs sim-within:   t(%d)=%.2f, p=%.4f\n',s13.df,s13.tstat,p13);
fprintf('Same across vs within: t(%d)=%.2f, p=%.4f\n',s14.df,s14.tstat,p14);
fprintf('Sim-across vs within: t(%d)=%.2f, p=%.4f\n',s23.df,s23.tstat,p23);
fprintf('Sim-within vs same-within: t(%d)=%.2f, p=%.4f\n',s34.df,s34.tstat,p34);

n_perm=1000; rng(42);
subj_pooled=mean([subj_same_across,subj_sim_across,subj_sim_within,subj_same_within],2);
null_means=nan(n_perm,1);
for p=1:n_perm, null_means(p)=mean(subj_pooled.*((rand(n,1)>.5)*2-1)); end

pp=nan(4,1);
for i=1:4, d={subj_same_across,subj_sim_across,subj_sim_within,subj_same_within}; om=mean(d{i});
    pd=nan(n_perm,1); for p=1:n_perm, pd(p)=mean(d{i}.*((rand(n,1)>.5)*2-1)); end
    pp(i)=mean(abs(pd)>=abs(om)); end
fprintf('Perm: same=%.4f, sim-across=%.4f, sim-within=%.4f, same-within=%.4f\n',pp(1),pp(2),pp(3),pp(4));

data_plot={subj_same_within,subj_sim_within,subj_same_across,subj_sim_across};
plot_pp=pp([4 3 1 2]);
labels={'Same-item within','Similar-item within','Same-item across','Similar-item across','Null'};
c_obs=[.03 .03 .03]; c_box=[.30 .30 .30]; c_null=[.82 .82 .82];
xp=[1 1.82 2.64 3.46]; x_null=4.28;

figure('color','w','position',[100 100 620 560]); hold on;
for i=1:4
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
ylabel('Gaze similarity (baseline-corrected r)','FontSize',18);
xlim([.55 4.75]); box off;
print(gcf,fullfile(fig_dir,'fig_gaze_tracks_similarity.pdf'),'-dpdf','-vector');