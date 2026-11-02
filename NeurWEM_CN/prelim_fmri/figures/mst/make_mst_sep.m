function make_mst_sep()
% MST canonical separation signal, GLM-panel style (one figure, ROI panels).
% Trial-type tuning by behavioural response: repeat (old-hit), lure-correct,
% novel (foil-CR). Canonical contrast bracketed = lure-correct vs repeat.
here=fileparts(mfilename('fullpath')); S=load(fullfile(here,'sep_roi.mat'));
N=numel(S.subs); INK=[0.149 0.149 0.184]; SIG=[0.768 0.306 0.321]; CMP=[0.357 0.561 0.788]; GREY=[0.72 0.72 0.76];
rois={'HC_post','HC_ant','Parahipp','Visual'};
disp_={'Posterior hippocampus','Anterior hippocampus','Parahippocampal','Visual cortex'};
xpos=[1 2 3]; col=[CMP;SIG;GREY];   % repeat / lure-correct / novel
f=figure('Color','w','Units','inches','Position',[1 1 11.0 3.6]);
tl=tiledlayout(f,1,4,'TileSpacing','compact','Padding','compact');
for r=1:4
    rn=rois{r};
    oh=S.([rn '_old_hit'])(:); lc=S.([rn '_lure_cor'])(:); fc=S.([rn '_foil_cr'])(:);
    dat={oh,lc,fc}; ax=nexttile(tl); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',11);
    for k=1:3
        v=dat{k}; rng(round(xpos(k)*7)); jit=(rand(N,1)-0.5)*0.22;
        scatter(ax,xpos(k)+jit,v,26,col(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
        plot(ax,[xpos(k)-0.26 xpos(k)+0.26],[mean(v) mean(v)],'-','Color',INK,'LineWidth',1.8);
    end
    yline(ax,0,':','Color',GREY,'LineWidth',0.8);
    p=ttest_rel_greater(lc,oh);                    % canonical: lure-correct > repeat
    allv=[oh;lc;fc]; ytop=max(allv); yr=max(allv)-min(allv)+1e-6; yy=ytop+0.10*yr;
    plot(ax,[1 1 2 2],[yy-0.03*yr yy yy yy-0.03*yr],'-','Color',INK,'LineWidth',0.8);
    if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
    text(ax,1.5,yy+0.015*yr,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9,'Color',INK);
    ax.XTick=xpos; ax.XTickLabel={'repeat','lure-correct','novel'}; ax.XTickLabelRotation=22;
    xlim(ax,[0.5 3.5]); ylim(ax,[min(allv)-0.08*yr ytop+0.42*yr]);
    ax.YAxis.Exponent=0; ax.YTickLabelRotation=0;
    title(ax,disp_{r},'FontName','Helvetica','FontWeight','normal','FontSize',12,'HorizontalAlignment','left');
    set(get(ax,'Title'),'Units','normalized','Position',[0 1.03 0]);
    box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
    if r==1, ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',12); end
end
exportgraphics(f,fullfile(here,'fig_mst_separation.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(here,'fig_mst_separation.png'),'Resolution',300);
disp('wrote fig_mst_separation.pdf / .png');
end
function p=ttest_rel_greater(a,b)
d=a-b; n=numel(d); t=mean(d)/(std(d)/sqrt(n));
p=1-tcdf(t,n-1);
end
