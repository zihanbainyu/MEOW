function make_mst_idx()
% MST behaviour, GLM-panel style, single figure, shared y. Two measures side by side:
% Recognition (old|old vs old|foil) and Lure discrimination (sim|lure vs sim|foil).
% Per-subject dots + mean crossbar; bracket over each pair = one-sample t of the
% difference (index) vs chance. y = proportion. Helvetica, vector PDF.
here=fileparts(mfilename('fullpath')); S=load(fullfile(here,'mst_idx.mat'));
N=numel(S.subs); INK=[0.149 0.149 0.184]; TARG=[0.768 0.306 0.321]; FOIL=[0.357 0.561 0.788]; GREY=[0.72 0.72 0.76];
dat={S.old_old(:),S.foil_old(:),S.lure_sim(:),S.foil_sim(:)};
x=[0 0.8 2.0 2.8]; col=[TARG;FOIL;TARG;FOIL];
xlab={'old | old','old | foil','sim | lure','sim | foil'};
f=figure('Color','w','Units','inches','Position',[1 1 4.8 4.8]);
ax=axes(f,'Position',[0.15 0.17 0.80 0.73]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
for k=1:4
    v=dat{k}; rng(round((x(k)+1)*7)); jit=(rand(N,1)-0.5)*0.20;
    scatter(ax,x(k)+jit,v,28,col(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
    plot(ax,[x(k)-0.26 x(k)+0.26],[mean(v) mean(v)],'-','Color',INK,'LineWidth',2.0);
end
    function bracket(i,j,d)
        [~,p]=ttest(d,0); yy=0.90;                 % one-sample t of index vs chance
        plot(ax,[x(i) x(i) x(j) x(j)],[yy-0.02 yy yy yy-0.02],'-','Color',INK,'LineWidth',0.8);
        if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
        text(ax,(x(i)+x(j))/2,yy+0.015,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
    end
bracket(1,2,S.REC(:)); bracket(3,4,S.LDI(:));
% cross comparison: Recognition index vs Lure-discrimination index (paired)
[~,pc]=ttest(S.REC(:),S.LDI(:)); yc=1.05;
plot(ax,[0.4 0.4 2.4 2.4],[yc-0.02 yc yc yc-0.02],'-','Color',INK,'LineWidth',0.8);
if pc<0.001, pct='p < .001'; else, pct=sprintf('p = %.3f',pc); end
text(ax,1.4,yc+0.015,pct,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
ax.XTick=x; ax.XTickLabel=xlab; ax.XTickLabelRotation=0;
xlim(ax,[-0.5 3.3]); ylim(ax,[0 1.2]); ax.YTick=0:0.25:1;
title(ax,'MST performance','FontName','Helvetica','FontWeight','normal','FontSize',14);
text(ax,0.40,-0.20,'Recognition','Units','data','HorizontalAlignment','center','FontName','Helvetica','FontSize',12,'Color',INK);
text(ax,2.40,-0.20,'Lure discrimination','Units','data','HorizontalAlignment','center','FontName','Helvetica','FontSize',12,'Color',INK);
ylabel(ax,'proportion of responses','FontName','Helvetica','FontSize',12);
box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
exportgraphics(f,fullfile(here,'fig_mst_idx.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(here,'fig_mst_idx.png'),'Resolution',300);
disp('wrote fig_mst_idx.pdf / .png');
end
