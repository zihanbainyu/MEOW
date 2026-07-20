clear; clc; close all;
base_dir='..'; res_dir=fullfile(base_dir,'results'); fig_dir=fullfile(base_dir,'figures');
c_comp=[180 174 211]/255; c_iso=[176 230 255]/255; exc=609;

%% A1A2
load(fullfile(res_dir,'cumu_reinstat_aa_cor.mat'));
nf_aa=n_fix_to_plot; xc=(1:nf_aa)'-mean((1:nf_aa)');
sl_aa_c=trial_slopes(cumulative_results_comp,5,nf_aa,xc);
sl_aa_i=trial_slopes(cumulative_results_iso,5,nf_aa,xc);
[aa_sc,aa_si,~]=subj_slopes(cumulative_results_comp.subj_id,sl_aa_c,cumulative_results_iso.subj_id,sl_aa_i,exc);
print_stats('A1A2',aa_sc,aa_si);
valid_aa=~any(isnan(subj_traj_comp),2)&~any(isnan(subj_traj_iso),2); valid_aa(common_subjs==exc)=false;
tc_aa=subj_traj_comp(valid_aa,:); ti_aa=subj_traj_iso(valid_aa,:);

%% B1A2
load(fullfile(res_dir,'cumu_reinstat_ba.mat'));
nf_ba=n_fix_to_plot;
sl_ba_c=trial_slopes(cumulative_results_comp,4,nf_ba,xc);
sl_ba_i=trial_slopes(cumulative_results_iso,4,nf_ba,xc);
[ba_sc,ba_si,~]=subj_slopes(cumulative_results_comp.subj_id,sl_ba_c,cumulative_results_iso.subj_id,sl_ba_i,exc);
print_stats('B1A2',ba_sc,ba_si);
valid_ba=~any(isnan(subj_traj_comp),2)&~any(isnan(subj_traj_iso),2);
tc_ba=subj_traj_comp(valid_ba,:); ti_ba=subj_traj_iso(valid_ba,:);

%% 2x2 ANOVA
run_anova('Cumulative Slope ANOVA',aa_sc,aa_si,ba_sc,ba_si);

%% Figure: A1A2 + B1A2 side by side + interaction
figure('color','w','position',[50 50 1200 400]);

yl=[-.02 .09];
yt=[-.01 0 .02 .04 .06 .08];

subplot('Position',[.06 .15 .27 .75]); hold on;
plot_traj_ax(tc_ba,ti_ba,ba_sc,ba_si,nf_ba,c_comp,c_iso,'B1-A2',yl,yt);

subplot('Position',[.38 .15 .27 .75]); hold on;
plot_traj_ax(tc_aa,ti_aa,aa_sc,aa_si,nf_aa,c_comp,c_iso,'A1-A2',yl,yt);

subplot('Position',[.74 .15 .22 .75]); hold on;
ns=length(ba_sc); jx=0.12;
rng(42); jt_c=0.06*(rand(ns,1)-0.5); jt_i=0.06*(rand(ns,1)-0.5);
% Per-subject connecting lines: comp B1A2→A1A2, iso B1A2→A1A2
for s=1:ns
    p=plot([1+jx+jt_c(s) 2-jx+jt_c(s)],[ba_sc(s) aa_sc(s)],'-','Color',[c_comp 0.18],'LineWidth',0.8); p.Annotation.LegendInformation.IconDisplayStyle='off';
    p=plot([1-jx+jt_i(s) 2+jx+jt_i(s)],[ba_si(s) aa_si(s)],'-','Color',[c_iso 0.18],'LineWidth',0.8); p.Annotation.LegendInformation.IconDisplayStyle='off';
end
sc1=scatter(1+jx+jt_c,ba_sc,35,c_comp,'filled','MarkerFaceAlpha',.5); sc1.Annotation.LegendInformation.IconDisplayStyle='off';
sc2=scatter(1-jx+jt_i,ba_si,35,c_iso,'filled','MarkerFaceAlpha',.5); sc2.Annotation.LegendInformation.IconDisplayStyle='off';
sc3=scatter(2-jx+jt_c,aa_sc,35,c_comp,'filled','MarkerFaceAlpha',.5); sc3.Annotation.LegendInformation.IconDisplayStyle='off';
sc4=scatter(2+jx+jt_i,aa_si,35,c_iso,'filled','MarkerFaceAlpha',.5); sc4.Annotation.LegendInformation.IconDisplayStyle='off';

% Mean ± SE
ms=[mean(ba_sc),mean(aa_sc); mean(ba_si),mean(aa_si)];
se=[std(ba_sc)/sqrt(ns),std(aa_sc)/sqrt(ns); std(ba_si)/sqrt(ns),std(aa_si)/sqrt(ns)];
errorbar([1+jx 2-jx],ms(1,:),se(1,:),'o-','Color',c_comp,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',c_comp,'CapSize',8,'DisplayName','Compared');
errorbar([1-jx 2+jx],ms(2,:),se(2,:),'s-','Color',c_iso,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',c_iso,'CapSize',8,'DisplayName','Isolated');
yline(0,'k--','LineWidth',1.5);
set(gca,'XTick',[1 2],'XTickLabel',{'B1-A2','A1-A2'},'FontSize',18,'LineWidth',1.5);
ylabel('Cumulative slope','FontSize',18); xlim([.5 2.5]); box off;
legend({'Compared','Isolated'},'FontSize',14,'Box','off','Location','best');
title('Interaction','FontSize',20,'FontWeight','normal');

print(gcf,fullfile(fig_dir,'gaze_cumu_aa_ba.pdf'),'-dpdf','-vector');

%% functions
function sl=trial_slopes(tbl,col0,nf,xc)
    sl=nan(height(tbl),1);
    for i=1:height(tbl), y=tbl{i,col0:col0+nf-1}'; v=~isnan(y); if sum(v)>=2, sl(i)=xc(v)\y(v); end, end
end

function [sc,si,common]=subj_slopes(sid_c,sl_c,sid_i,sl_i,exc)
    agg=@(s,sl) grpstats(table(s(~isnan(sl)),sl(~isnan(sl)),'VariableNames',{'subj_id','slope'}),'subj_id','mean','DataVars','slope');
    gc=agg(sid_c,sl_c); gc.Properties.VariableNames{'mean_slope'}='slope';
    gi=agg(sid_i,sl_i); gi.Properties.VariableNames{'mean_slope'}='slope';
    common=intersect(gc.subj_id,gi.subj_id); common(common==exc)=[];
    [~,ic]=ismember(common,gc.subj_id); [~,ii]=ismember(common,gi.subj_id);
    sc=gc.slope(ic); si=gi.slope(ii);
end

function print_stats(lbl,sc,si)
    fprintf('\ncumulative %s\n',lbl);
    [~,pc,~,tc]=ttest(sc,0,'Tail','right'); [~,pi,~,ti]=ttest(si,0,'Tail','right');
    fprintf('Comp: M=%.4f SD=%.4f, t(%d)=%.3f, p=%.4f, d=%.3f\n',mean(sc),std(sc),tc.df,tc.tstat,pc,mean(sc)/std(sc));
    fprintf('Iso:  M=%.4f SD=%.4f, t(%d)=%.3f, p=%.4f, d=%.3f\n',mean(si),std(si),ti.df,ti.tstat,pi,mean(si)/std(si));
    [~,pd,cd,td]=ttest(sc,si,'Tail','right'); fprintf('Diff: t(%d)=%.3f, p=%.4f, d=%.3f, 95%%CI=[%.4f,%.4f]\n',td.df,td.tstat,pd,mean(sc-si)/std(sc-si),cd(1),cd(2));
    [pw,~,ws]=signrank(sc,si); fprintf('Wilcoxon: z=%.3f, p=%.4f\n',ws.zval,pw);
end

function plot_traj_ax(tc,ti,sc,si,nf,cc,ci,ttl,yl,yt)
    mc=mean(tc,1,'omitnan'); se_c=std(tc,0,1,'omitnan')/sqrt(size(tc,1));
    mi=mean(ti,1,'omitnan'); se_i=std(ti,0,1,'omitnan')/sqrt(size(ti,1));
    fn=1:nf; xl=linspace(1,nf,50); xlc=xl-mean(1:nf);
    % SE bands (subtle context)
    p=fill([fn,fliplr(fn)],[mc+se_c,fliplr(mc-se_c)],cc,'FaceAlpha',.15,'EdgeColor','none'); p.Annotation.LegendInformation.IconDisplayStyle='off';
    p=fill([fn,fliplr(fn)],[mi+se_i,fliplr(mi-se_i)],ci,'FaceAlpha',.15,'EdgeColor','none'); p.Annotation.LegendInformation.IconDisplayStyle='off';
    % Mean trajectory (thin dashed, raw context)
    p=plot(fn,mc,'--','Color',cc,'LineWidth',1.5); p.Annotation.LegendInformation.IconDisplayStyle='off';
    p=plot(fn,mi,'--','Color',ci,'LineWidth',1.5); p.Annotation.LegendInformation.IconDisplayStyle='off';
    % Mean regression line (thick — primary visual claim)
    plot(xl,mean(sc)*xlc+mean(mc),'-','Color',cc,'LineWidth',3);
    plot(xl,mean(si)*xlc+mean(mi),'-','Color',ci,'LineWidth',3);
    yline(0,'k--','LineWidth',1.5);
    xlabel('Cumulative fixation number','FontSize',18); ylabel('Gaze similarity (r)','FontSize',18);
    title(ttl,'FontSize',20,'FontWeight','normal');
    legend({'Compared','Isolated'},'Location','best','FontSize',16,'Box','off');
    xlim([0 nf+.5]); ylim(yl); set(gca,'XTick',1:nf,'YTick',yt,'FontSize',18,'LineWidth',1.5); box off;
end

function run_anova(lbl,ac,ai,bc,bi)
    within=table({'A1A2';'A1A2';'B1A2';'B1A2'},{'Comp';'Iso';'Comp';'Iso'},'VariableNames',{'PairType','TaskCond'});
    t=table(ac,ai,bc,bi,'VariableNames',{'ac','ai','bc','bi'});
    rm=fitrm(t,'ac-bi~1','WithinDesign',within); tbl=ranova(rm,'WithinModel','PairType*TaskCond');
    df1=tbl.DF(3); df2=tbl.DF(4);
    fprintf('\n%s:\n',lbl);
    for k=[3 5 7], eta=tbl.SumSq(k)/(tbl.SumSq(k)+tbl.SumSq(k+1));
        fprintf('  %s: F(%d,%d)=%.2f, p=%.4f, eta2=%.3f %s\n',tbl.Properties.RowNames{k},df1,df2,tbl.F(k),tbl.pValue(k),eta,...
            repmat('*',1,(tbl.pValue(k)<.05)+(tbl.pValue(k)<.01)+(tbl.pValue(k)<.001))); end
    cd=@(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
    [~,p1,~,s1]=ttest(ac,ai,'Tail','right'); [~,p2,~,s2]=ttest(bc,bi,'Tail','right'); [~,p3,~,s3]=ttest(ac,bc); [~,p4,~,s4]=ttest(ai,bi);
    fprintf('  A1A2 comp-iso:  t(%d)=%.2f, d=%.2f, p=%.4f %s\n',s1.df,s1.tstat,cd(ac,ai),p1,repmat('*',1,(p1<.05)+(p1<.01)+(p1<.001)));
    fprintf('  B1A2 comp-iso:  t(%d)=%.2f, d=%.2f, p=%.4f %s\n',s2.df,s2.tstat,cd(bc,bi),p2,repmat('*',1,(p2<.05)+(p2<.01)+(p2<.001)));
    fprintf('  comp A1A2-B1A2: t(%d)=%.2f, d=%.2f, p=%.4f %s\n',s3.df,s3.tstat,cd(ac,bc),p3,repmat('*',1,(p3<.05)+(p3<.01)+(p3<.001)));
    fprintf('  iso  A1A2-B1A2: t(%d)=%.2f, d=%.2f, p=%.4f %s\n',s4.df,s4.tstat,cd(ai,bi),p4,repmat('*',1,(p4<.05)+(p4<.01)+(p4<.001)));
end