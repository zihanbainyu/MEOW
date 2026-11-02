function make_roi_2x2()
% MAKE_ROI_2X2  One house-style vector PDF per ROI from the condition-level 2x2 GLM.
% GLM: MST condition x accuracy 2x2 (old+similar [EM] vs novel  x  correct vs incorrect),
%      condition-level first-level betas (NOT LS-S single-trial), N=7.
% Each figure: region name > anatomical icon > correct vs incorrect bars, per-subject
% dots + mean crossbar; bracket = correct > incorrect (one-tailed paired t).
% Source: m3_roi.mat (fields <roi>_em_c, <roi>_em_i per subject).
here=fileparts(mfilename('fullpath')); S=load(fullfile(here,'m3_roi.mat'));

% ROI: posterior hippocampus
roi_panel(S,here,'HC_post','Posterior hippocampus','roi_HCpost');

% ROI: angular gyrus
roi_panel(S,here,'Angular','Angular gyrus','roi_Angular');

% ROI: precuneus
roi_panel(S,here,'Precuneus','Precuneus','roi_Precuneus');

% ROI: visual cortex
roi_panel(S,here,'Visual','Visual cortex','roi_Visual');
end

function roi_panel(S,here,roi,roilab,icon)
N=numel(S.subs); INK=[0.149 0.149 0.184]; CORR=[0.768 0.306 0.321]; INC=[0.357 0.561 0.788]; GREY=[0.72 0.72 0.76];
ec=S.([roi '_em_c'])(:); ei=S.([roi '_em_i'])(:); dat={ec,ei}; x=[1 2]; col=[CORR;INC];
f=figure('Color','w','Units','inches','Position',[1 1 4.4 4.8]);
ip=fullfile(here,'icons',[icon '.png']);
if isfile(ip), ia=axes(f,'Position',[0.30 0.70 0.40 0.20]); imshow(imread(ip),'Parent',ia); axis(ia,'image','off'); end
ax=axes(f,'Position',[0.17 0.14 0.78 0.50]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
for k=1:2
    v=dat{k}; rng(round(x(k)*7)); jit=(rand(N,1)-0.5)*0.22;
    scatter(ax,x(k)+jit,v,28,col(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
    plot(ax,[x(k)-0.28 x(k)+0.28],[mean(v) mean(v)],'-','Color',INK,'LineWidth',2.0);
end
yline(ax,0,':','Color',GREY,'LineWidth',0.8);
d=ec-ei; t=mean(d)/(std(d)/sqrt(N)); p=1-tcdf(t,N-1);
allv=[ec;ei]; dmin=min(allv); dmax=max(allv); yr=dmax-dmin; yy=dmax+0.11*yr;
plot(ax,[1 1 2 2],[yy-0.025*yr yy yy yy-0.025*yr],'-','Color',INK,'LineWidth',0.8);
if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
text(ax,1.5,yy+0.01*yr,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
tk=niceticks(dmin,dmax,5);
ax.XTick=x; ax.XTickLabel={'correct','incorrect'}; ax.YTick=tk;
xlim(ax,[0.4 2.6]); ylim(ax,[min(dmin,tk(1))-0.06*yr max(dmax+0.40*yr, tk(end)+0.01*yr)]);
title(ax,roilab,'FontName','Helvetica','FontWeight','normal','FontSize',13);
ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',11);
box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
exportgraphics(f,fullfile(here,sprintf('fig_%s_2x2.pdf',lower(roi))),'ContentType','vector');
close(f);
end

function tk=niceticks(dmin,dmax,n)
raw=(dmax-dmin)/(n-1); mag=10^floor(log10(raw)); cand=[1 2 2.5 5 10]*mag;
step=cand(find(cand>=raw,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1);
if tk(end)<dmax, step=cand(find(cand>=(dmax-dmin)/(n-1)*1.5,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1); end
end
