function make_sep3_fig()
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data','b2_separation.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22]; cGray=[0.55 0.55 0.55];
nets={'Hippocampus','Visual','WM_network'}; netlab={'Hippocampus','Visual areas','WM network'};
subs=unique(T.subj); ns=numel(subs);
cells={'compared',1,cBlue; 'compared',0,cRed; 'novel',1,cBlue; 'novel',0,cRed};
f=figure('Color','w','Position',[30 30 1200 470]);
for k=1:3
    subplot(1,3,k); hold on;
    M=nan(ns,4);
    for s=1:ns
        for c=1:4
            M(s,c)=mean(T.dist(T.subj==subs(s)&strcmp(T.roi,nets{k})&strcmp(T.condition,cells{c,1})&T.correct==cells{c,2}),'omitnan');
        end
    end
    for c=1:4, house_box(c,M(:,c),cells{c,3}); end
    for s=1:ns
        plot([1 2],M(s,1:2),'-','Color',[cGray 0.4],'LineWidth',0.6);
        plot([3 4],M(s,3:4),'-','Color',[cGray 0.4],'LineWidth',0.6);
    end
    plot([2.5 2.5],ylim,'-','Color',[0.85 0.85 0.85],'LineWidth',0.75);
    set(gca,'XTick',[1.5 3.5],'XTickLabel',{'compared','novel'},'FontSize',13,'LineWidth',1,'TickDir','out');
    xlim([0.4 4.6]); box off;
    title(netlab{k},'FontSize',16,'FontWeight','normal');
    if k==1, ylabel('A2\leftrightarrowB2 neural distance (1 - r)','FontSize',14); end
    pc=signrank_safe(M(:,1),M(:,2));
    yy=max(M(:))*1.01; if ~isnan(pc), text(1.5,yy,sprintf('d=%+.2f',nanmean(M(:,1)-M(:,2))/nanstd(M(:,1)-M(:,2))),'HorizontalAlignment','center','FontSize',12,'Color',[0.3 0.3 0.3]); end
    fprintf('%-12s compared c=%.3f i=%.3f (d=%+.2f)   novel c=%.3f i=%.3f (d=%+.2f)\n',netlab{k},...
        nanmean(M(:,1)),nanmean(M(:,2)),nanmean(M(:,1)-M(:,2))/nanstd(M(:,1)-M(:,2)),...
        nanmean(M(:,3)),nanmean(M(:,4)),nanmean(M(:,3)-M(:,4))/nanstd(M(:,3)-M(:,4)));
end
text(0.5,0.985,'\color[rgb]{0.16,0.44,0.74}correct   \color[rgb]{0.80,0.25,0.22}incorrect','Units','normalized','Parent',gca,'HorizontalAlignment','center','FontSize',13);
sgtitle('Pairmate neural separation (A2 vs B2) at the probe, by outcome and condition','FontSize',16,'FontWeight','normal');
savepdf(f,fullfile(base,'figures','fig_pairmate_sep_3net.pdf')); close(f);
fprintf('wrote fig_pairmate_sep_3net.pdf\n');
end
function savepdf(f,path)
set(f,'Units','points'); pos=get(f,'Position');
set(f,'PaperUnits','points','PaperPosition',[0 0 pos(3) pos(4)],'PaperSize',[pos(3) pos(4)]);
print(f,path,'-dpdf','-vector');
end
function house_box(x,d,c)
d=d(~isnan(d)); if isempty(d), return; end
w=0.30; q=quantile(d,[0.25 0.5 0.75]);
patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',0.15,'EdgeColor',c,'LineWidth',1.2);
plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2.5);
jit=(rand(size(d))-0.5)*0.22; scatter(x+jit,d,44,c,'filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','w','LineWidth',0.5);
end
function p=signrank_safe(a,b)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m); if numel(a)<2, p=NaN; return; end
try, p=signrank(a,b); catch, [~,p]=ttest(a,b); end
end
