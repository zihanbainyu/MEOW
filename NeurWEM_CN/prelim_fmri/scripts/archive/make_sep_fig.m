function make_sep_fig()
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data','b2_separation.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22]; cGray=[0.55 0.55 0.55];
rois={'EarlyVis','LOC','Fusiform','HC','HC_post','Angular','Precuneus','DLPFC','IPS'};
roilab={'Early vis','LOC','Fusiform','Hipp','post Hipp','Angular','Precun','DLPFC','IPS'};
subs=unique(T.subj); ns=numel(subs); nr=numel(rois);
% subject-level mean A2-B2 distance per ROI x correctness
C=nan(ns,nr); I=nan(ns,nr);
for r=1:nr
    for s=1:ns
        C(s,r)=mean(T.dist(T.subj==subs(s)&strcmp(T.roi,rois{r})&T.correct==1),'omitnan');
        I(s,r)=mean(T.dist(T.subj==subs(s)&strcmp(T.roi,rois{r})&T.correct==0),'omitnan');
    end
end
f=figure('Color','w','Position',[30 30 1150 560]); hold on;
xt=[]; xl={};
for r=1:nr
    x0=(r-1)*3;
    house_box(x0+1,C(:,r),cBlue); house_box(x0+2,I(:,r),cRed);
    for s=1:ns, if ~isnan(C(s,r))&&~isnan(I(s,r)), plot([x0+1 x0+2],[C(s,r) I(s,r)],'-','Color',[cGray 0.4],'LineWidth',0.6); end, end
    xt=[xt x0+1.5]; xl{end+1}=roilab{r}; %#ok<AGROW>
    p=signrank_safe(C(:,r),I(:,r));
    yy=max([C(:,r);I(:,r)])+0.02; if ~isnan(p)&&p<0.05, text(x0+1.5,yy,star(p),'HorizontalAlignment','center','FontSize',15); end
end
set(gca,'XTick',xt,'XTickLabel',xl,'FontSize',13,'LineWidth',1,'TickDir','out');
ylabel('A2\leftrightarrowB2 neural distance (1 - r)','FontSize',15);
title('Pairmate neural separation at the 2-back probe, by discrimination outcome','FontSize',17,'FontWeight','normal');
xlim([0.2 nr*3-0.2]); box off;
text(0.012,0.985,'\color[rgb]{0.16,0.44,0.74}correct   \color[rgb]{0.80,0.25,0.22}incorrect','Units','normalized','VerticalAlignment','top','FontSize',15);
savepdf(f,fullfile(base,'figures','fig_pairmate_separation.pdf')); close(f);
fprintf('ROI          corr   inc    diff(c-i)  d     p\n');
for r=1:nr
    df=C(:,r)-I(:,r);
    fprintf('%-11s %5.3f %5.3f   %+5.3f   %+5.2f  %.3f\n',roilab{r},nanmean(C(:,r)),nanmean(I(:,r)),nanmean(df),nanmean(df)/nanstd(df),signrank_safe(C(:,r),I(:,r)));
end
fprintf('wrote fig_pairmate_separation.pdf\n');
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
jit=(rand(size(d))-0.5)*0.22; scatter(x+jit,d,46,c,'filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','w','LineWidth',0.5);
end
function p=signrank_safe(a,b)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m); if numel(a)<2, p=NaN; return; end
try, p=signrank(a,b); catch, [~,p]=ttest(a,b); end
end
function s=star(p)
if isnan(p), s=''; elseif p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s=''; end
end
