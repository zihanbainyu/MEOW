function make_mst_combined()
% MST summary panel: (A) behaviour (recognition + lure discrimination) and
% (B) memory-success ROIs (correct vs incorrect, old+similar), side by side.
here=fileparts(mfilename('fullpath'));
I=load(fullfile(here,'mst_idx.mat')); M=load(fullfile(here,'m3_roi.mat'));
INK=[0.149 0.149 0.184]; CORR=[0.768 0.306 0.321]; INC=[0.357 0.561 0.788]; GREY=[0.72 0.72 0.76];
N=numel(I.subs);
f=figure('Color','w','Units','inches','Position',[1 1 15 4.7]);

% ===================== Panel A : behaviour =====================
axA=axes(f,'Position',[0.055 0.17 0.185 0.63]); hold(axA,'on'); set(axA,'FontName','Helvetica','FontSize',10);
datA={I.old_old(:),I.foil_old(:),I.lure_sim(:),I.foil_sim(:)}; xA=[0 1.0 2.4 3.4]; colA=[CORR;INC;CORR;INC];
xlabA={'old | old','old | foil','sim | lure','sim | foil'};
for k=1:4
    v=datA{k}; rng(round((xA(k)+1)*7)); jit=(rand(N,1)-0.5)*0.20;
    scatter(axA,xA(k)+jit,v,26,colA(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
    plot(axA,[xA(k)-0.26 xA(k)+0.26],[mean(v) mean(v)],'-','Color',INK,'LineWidth',2.0);
end
    function pbracket(ax,xi,xj,d,yy)
        [~,p]=ttest(d,0);
        plot(ax,[xi xi xj xj],[yy-0.02 yy yy yy-0.02],'-','Color',INK,'LineWidth',0.8);
        if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
        text(ax,(xi+xj)/2,yy+0.012,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9,'Color',INK);
    end
pbracket(axA,xA(1),xA(2),I.REC(:),0.90); pbracket(axA,xA(3),xA(4),I.LDI(:),0.90);
cL=(xA(1)+xA(2))/2; cR=(xA(3)+xA(4))/2;
[~,pc]=ttest(I.REC(:),I.LDI(:)); yc=1.05;
plot(axA,[cL cL cR cR],[yc-0.02 yc yc yc-0.02],'-','Color',INK,'LineWidth',0.8);
if pc<0.001, pct='p < .001'; else, pct=sprintf('p = %.3f',pc); end
text(axA,(cL+cR)/2,yc+0.012,pct,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9,'Color',INK);
axA.XTick=xA; axA.XTickLabel=xlabA; axA.XTickLabelRotation=0;
xlim(axA,[-0.6 4.0]); ylim(axA,[0 1.2]); axA.YTick=0:0.25:1;
ylabel(axA,'proportion of responses','FontName','Helvetica','FontSize',11);
text(axA,cL,-0.18,'Recognition','Units','data','HorizontalAlignment','center','FontName','Helvetica','FontSize',10.5,'Color',INK);
text(axA,cR,-0.18,'Lure discrimination','Units','data','HorizontalAlignment','center','FontName','Helvetica','FontSize',10.5,'Color',INK);
title(axA,'MST performance','FontName','Helvetica','FontWeight','normal','FontSize',12);
box(axA,'off'); axA.TickLength=[0 0]; axA.XColor=INK; axA.YColor=INK; axA.LineWidth=0.8;

% ===================== Panel B : memory-success ROIs =====================
rois={'HC_post','Angular','Precuneus','Visual'}; icons={'roi_HCpost','roi_Angular','roi_Precuneus','roi_Visual'};
disp_={'Posterior hippocampus','Angular gyrus','Precuneus','Visual cortex'};
nR=4; xB=[1 2]; colB=[CORR;INC];
Lm=0.345; Rr=0.985; gp=0.018; colw=((Rr-Lm)-(nR-1)*gp)/nR; iw=0.55*colw;
for r=1:nR
    x0=Lm+(r-1)*(colw+gp);
    annotation(f,'textbox',[x0 0.855 colw 0.05],'String',disp_{r},'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',10.5,'Color',INK);
    img=imread(fullfile(here,'icons',[icons{r} '.png']));
    ia=axes(f,'Position',[x0+(colw-iw)/2 0.66 iw 0.19]); imshow(img,'Parent',ia); axis(ia,'image','off');
    rn=rois{r}; ec=M.([rn '_em_c'])(:); ei=M.([rn '_em_i'])(:); dat={ec,ei};
    ax=axes(f,'Position',[x0 0.17 colw 0.47]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10);
    for k=1:2
        v=dat{k}; rng(round(xB(k)*7)); jit=(rand(N,1)-0.5)*0.22;
        scatter(ax,xB(k)+jit,v,26,colB(k,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85);
        plot(ax,[xB(k)-0.28 xB(k)+0.28],[mean(v) mean(v)],'-','Color',INK,'LineWidth',2.0);
    end
    yline(ax,0,':','Color',GREY,'LineWidth',0.8);
    d=ec-ei; nn=numel(d); t=mean(d)/(std(d)/sqrt(nn)); p=1-tcdf(t,nn-1);
    allv=[ec;ei]; dmin=min(allv); dmax=max(allv); yr=dmax-dmin; yy=dmax+0.11*yr;
    plot(ax,[1 1 2 2],[yy-0.025*yr yy yy yy-0.025*yr],'-','Color',INK,'LineWidth',0.8);
    if p<0.001, pt='p < .001'; else, pt=sprintf('p = %.3f',p); end
    text(ax,1.5,yy+0.01*yr,pt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9,'Color',INK);
    tk=niceticks(dmin,dmax,5);
    ax.XTick=xB; ax.XTickLabel={'correct','incorrect'}; ax.YTick=tk;
    xlim(ax,[0.4 2.6]); ylim(ax,[min(dmin,tk(1))-0.06*yr max(dmax+0.40*yr, tk(end)+0.01*yr)]);
    box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
    if r==1, ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',11); end
end
% panel letters
annotation(f,'textbox',[0.005 0.92 0.04 0.06],'String','A','EdgeColor','none','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',INK);
annotation(f,'textbox',[0.30 0.92 0.04 0.06],'String','B','EdgeColor','none','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',INK);
exportgraphics(f,fullfile(here,'fig_mst_panel.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(here,'fig_mst_panel.png'),'Resolution',300);
disp('wrote fig_mst_panel.pdf / .png');
end
function tk=niceticks(dmin,dmax,n)
raw=(dmax-dmin)/(n-1); mag=10^floor(log10(raw)); cand=[1 2 2.5 5 10]*mag;
step=cand(find(cand>=raw,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1);
if tk(end)<dmax, step=cand(find(cand>=(dmax-dmin)/(n-1)*1.5,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1); end
end
