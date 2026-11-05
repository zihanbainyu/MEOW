function make_mst_m3()
% MST memory-success regions (old + similar trials), GLM style with small ROI icons.
% Per region: name (text) > small anatomical icon > correct vs incorrect panel;
% bracket = correct > incorrect (one-tailed paired t). Each panel own y-range, 4 ticks.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..'); S=load(fullfile(root,'data','mst','m3_roi.mat'));
N=numel(S.subs); INK=[0.149 0.149 0.184]; CORR=[0.768 0.306 0.321]; INC=[0.357 0.561 0.788]; GREY=[0.72 0.72 0.76];
rois={'HC_post','Angular','Precuneus','Visual'}; icons={'roi_HCpost','roi_Angular','roi_Precuneus','roi_Visual'};
disp_={'Posterior hippocampus','Angular gyrus','Precuneus','Visual cortex'};
nR=numel(rois); x=[1 2]; col=[CORR;INC];
f=figure('Color','w','Units','inches','Position',[1 1 10.2 4.4]);
Lm=0.085; Rr=0.985; gap=0.022; colw=((Rr-Lm)-(nR-1)*gap)/nR; iw=0.55*colw;
for r=1:nR
    x0=Lm+(r-1)*(colw+gap);
    annotation(f,'textbox',[x0 0.845 colw 0.055],'String',disp_{r},'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',11,'Color',INK);
    img=imread(fullfile(root,'icons',[icons{r} '.png']));
    ia=axes(f,'Position',[x0+(colw-iw)/2 0.645 iw 0.20]); imshow(img,'Parent',ia); axis(ia,'image','off');
    rn=rois{r}; ec=S.([rn '_em_c'])(:); ei=S.([rn '_em_i'])(:); dat={ec,ei};
    ax=axes(f,'Position',[x0 0.14 colw 0.50]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
    for k=1:2
        v=dat{k}; rng(round(x(k)*7)); jit=(rand(N,1)-0.5)*0.22;
        scatter(ax,x(k)+jit,v,28,col(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
        plot(ax,[x(k)-0.28 x(k)+0.28],[mean(v) mean(v)],'-','Color',INK,'LineWidth',2.0);
    end
    yline(ax,0,':','Color',GREY,'LineWidth',0.8);
    d=ec-ei; n=numel(d); t=mean(d)/(std(d)/sqrt(n)); p=1-tcdf(t,n-1);
    allv=[ec;ei]; dmin=min(allv); dmax=max(allv); yr=dmax-dmin; yy=dmax+0.11*yr;
    plot(ax,[1 1 2 2],[yy-0.025*yr yy yy yy-0.025*yr],'-','Color',INK,'LineWidth',0.8);
    if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
    text(ax,1.5,yy+0.01*yr,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
    % own range, same number of ticks across panels
    tk=niceticks(dmin,dmax,5);
    ax.XTick=x; ax.XTickLabel={'correct','incorrect'}; ax.YTick=tk;
    xlim(ax,[0.4 2.6]); ylim(ax,[min(dmin,tk(1))-0.06*yr max(dmax+0.40*yr, tk(end)+0.01*yr)]);
    box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
    if r==1, ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',11); end
end
exportgraphics(f,fullfile(root,'figures','mst','fig_mst_memory.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(root,'figures','mst','fig_mst_memory.png'),'Resolution',300);
disp('wrote fig_mst_memory.pdf / .png');
end
function tk=niceticks(dmin,dmax,n)
% n nice ticks spanning the data (same count across panels)
raw=(dmax-dmin)/(n-1); mag=10^floor(log10(raw)); cand=[1 2 2.5 5 10]*mag;
step=cand(find(cand>=raw,1));
lo=floor(dmin/step)*step; tk=lo+step*(0:n-1);
if tk(end)<dmax, step=cand(find(cand>=(dmax-dmin)/(n-1)*1.5,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1); end
end
