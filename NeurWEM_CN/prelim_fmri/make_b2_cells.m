function make_b2_cells()
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data','b2_cells_z.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22];
rois={'HC','HC_ant','HC_post','Angular','DLPFC','IPS'};
roilab={'Hippocampus','Hippocampus (ant)','Hippocampus (post)','Angular','DLPFC','IPS'};
cells={'compared','correct',1,cBlue; 'compared','incorrect',2,cRed; ...
       'novel','correct',3,cBlue; 'novel','incorrect',4,cRed};
f=figure('Color','w','Position',[20 20 1200 720]);
for r=1:numel(rois)
    subplot(2,3,r); hold on;
    allp=[];
    for k=1:size(cells,1)
        d=T.z(strcmp(T.roi,rois{r})&strcmp(T.condition,cells{k,1})&strcmp(T.acc,cells{k,2}));
        house_box(cells{k,3},d,cells{k,4}); allp=[allp;d];
    end
    lo=prctile(allp,1); hi=prctile(allp,99); pad=0.12*(hi-lo+eps); ylim([lo-pad hi+pad]);
    yline(0,'k:','LineWidth',0.75);
    plot([2.5 2.5],ylim,'-','Color',[0.82 0.82 0.82],'LineWidth',0.75);
    set(gca,'XTick',[1.5 3.5],'XTickLabel',{'compared','novel'},'FontSize',13,'LineWidth',1,'TickDir','out');
    xlim([0.4 4.6]); box off;
    title(roilab{r},'FontSize',15,'FontWeight','normal');
    if r==1||r==4, ylabel('B2 single-trial \beta (within-subj z)','FontSize',13); end
    % per-cell mean markers (thick, for the effect)
    for k=1:size(cells,1)
        d=T.z(strcmp(T.roi,rois{r})&strcmp(T.condition,cells{k,1})&strcmp(T.acc,cells{k,2}));
        if ~isempty(d), plot(cells{k,3}+[-0.33 0.33],[mean(d) mean(d)],'-','Color',cells{k,4}*0.65,'LineWidth',3); end
    end
end
annotation('textbox',[0.30 0.955 0.4 0.04],'String',...
    '\color[rgb]{0.16,0.44,0.74}\bfcorrect\rm    \color[rgb]{0.80,0.25,0.22}\bfincorrect\rm\color{black}   (thick line = cell mean)',...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',12,'Interpreter','tex');
sgtitle('B2 trials by condition \times accuracy (within-subject z, pooled across 5 subjects)','FontSize',16,'FontWeight','normal');
savepdf(f,fullfile(base,'figures','fig_b2_cells_byROI.pdf'));
fprintf('wrote fig_b2_cells_byROI.pdf\n');
for r=1:numel(rois)
    row='';
    for k=1:size(cells,1)
        d=T.z(strcmp(T.roi,rois{r})&strcmp(T.condition,cells{k,1})&strcmp(T.acc,cells{k,2}));
        row=[row sprintf(' %s-%s=%.2f',cells{k,1}(1),cells{k,2}(1),mean(d))];
    end
    fprintf('%-18s%s\n',roilab{r},row);
end
close(f);
end

function savepdf(f,path)
set(f,'Units','points'); pos=get(f,'Position');
set(f,'PaperUnits','points','PaperPosition',[0 0 pos(3) pos(4)],'PaperSize',[pos(3) pos(4)]);
print(f,path,'-dpdf','-vector');
end

function house_box(x,d,c)
d=d(~isnan(d)); if isempty(d), return; end
w=0.30; q=quantile(d,[0.25 0.5 0.75]);
patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',0.14,'EdgeColor',c,'LineWidth',1.1);
plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2);
jit=(rand(size(d))-0.5)*0.26;
scatter(x+jit,d,14,c,'filled','MarkerFaceAlpha',0.40,'MarkerEdgeColor','none');
end
