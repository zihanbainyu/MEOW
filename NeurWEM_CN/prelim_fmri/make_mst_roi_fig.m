function make_mst_roi_fig()
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data','mst_roi.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
rois={'Hippocampus','Parahipp','Fusiform','Visual','Angular','Precuneus','DLPFC','IPS'};
roilab={'Hipp','PHC','Fusi','Vis','Ang','Prec','DLPFC','IPS'};
cons={'recognition','lure_vs_new','discrimination','error'};
conlab={'Recognition (old - new)','Similar-discrim (lure - new)','Resolve competitor (lure - old)','Error (incorrect - correct)'};
subs=unique(T.subj); ns=numel(subs); nr=numel(rois);
col=[0.20 0.40 0.70];
f=figure('Color','w','Position',[20 20 1250 640]);
for k=1:4
    subplot(2,2,k); hold on;
    M=nan(ns,nr);
    for r=1:nr, for s=1:ns
        v=T.value(T.subj==subs(s)&strcmp(T.contrast,cons{k})&strcmp(T.roi,rois{r}));
        if ~isempty(v), M(s,r)=v(1); end
    end, end
    for r=1:nr
        mu=nanmean(M(:,r));
        c=col; if mu<0, c=[0.75 0.30 0.25]; end
        bar(r,mu,0.6,'FaceColor',c,'FaceAlpha',0.25,'EdgeColor',c,'LineWidth',1);
        jit=(rand(1,ns)-0.5)*0.25; scatter(r+jit,M(:,r),30,c,'filled','MarkerFaceAlpha',0.9,'MarkerEdgeColor','w');
    end
    yline(0,'k-','LineWidth',0.75);
    set(gca,'XTick',1:nr,'XTickLabel',roilab,'FontSize',12,'LineWidth',1,'TickDir','out');
    xlim([0.3 nr+0.7]); box off;
    title(conlab{k},'FontSize',14,'FontWeight','normal');
    if mod(k,2)==1, ylabel('% signal change','FontSize',12); end
end
sgtitle('Final MST (pure EM): regional contrast values  (N=3, points = subjects)','FontSize',15,'FontWeight','normal');
savepdf(f,fullfile(base,'figures','fig_mst_roi.pdf')); close(f);
fprintf('contrast        ROI        mean(N=3)\n');
for k=1:4, for r=1:nr
    M=[]; for s=1:ns, v=T.value(T.subj==subs(s)&strcmp(T.contrast,cons{k})&strcmp(T.roi,rois{r})); if ~isempty(v), M(end+1)=v(1); end, end %#ok<AGROW>
    fprintf('%-15s %-10s %+.3f\n',cons{k},rois{r},mean(M));
end, end
fprintf('wrote fig_mst_roi.pdf\n');
end
function savepdf(f,path)
set(f,'Units','points'); pos=get(f,'Position');
set(f,'PaperUnits','points','PaperPosition',[0 0 pos(3) pos(4)],'PaperSize',[pos(3) pos(4)]);
print(f,path,'-dpdf','-vector');
end
