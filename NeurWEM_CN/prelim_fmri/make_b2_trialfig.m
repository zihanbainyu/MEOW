function make_b2_trialfig(roi)
if nargin<1, roi='HC'; end
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data',sprintf('b2_trialbetas_%s.csv',roi)));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22];
subs=unique(T.subj); ns=numel(subs);
f=figure('Color','w','Position',[30 30 260*ns 460]);
for s=1:ns
    subplot(1,ns,s); hold on;
    C=T.beta(T.subj==subs(s)&strcmp(T.acc,'correct'));
    I=T.beta(T.subj==subs(s)&strcmp(T.acc,'incorrect'));
    house_box(1,C,cBlue); house_box(2,I,cRed);
    allv=[C;I]; lo=prctile(allv,2); hi=prctile(allv,98); pad=0.15*(hi-lo+eps);
    ylim([lo-pad hi+pad]);
    set(gca,'XTick',[1 2],'XTickLabel',{'correct','incorrect'},'FontSize',14,'LineWidth',1,'TickDir','out');
    xlim([0.4 2.6]); box off; yline(0,'k:','LineWidth',0.75);
    title(sprintf('sub-%d',subs(s)),'FontSize',16,'FontWeight','normal');
    if s==1, ylabel(sprintf('%s single-trial \\beta (B2)',roi),'FontSize',15); end
end
sgtitle('B2 discrimination trials: single-trial betas by accuracy','FontSize',17,'FontWeight','normal');
savepdf(f,fullfile(base,'figures',sprintf('fig_b2_trialbetas_%s.pdf',roi)));
fprintf('wrote fig_b2_trialbetas_%s.pdf\n',roi);
for s=1:ns
    C=T.beta(T.subj==subs(s)&strcmp(T.acc,'correct')); I=T.beta(T.subj==subs(s)&strcmp(T.acc,'incorrect'));
    [~,p]=ttest2(C,I);
    fprintf('sub-%d  correct=%.2f (n=%d)  incorrect=%.2f (n=%d)  diff=%.2f  p=%.3f\n',subs(s),mean(C),numel(C),mean(I),numel(I),mean(C)-mean(I),p);
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
w=0.30;
q=quantile(d,[0.25 0.5 0.75]);
patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',0.15,'EdgeColor',c,'LineWidth',1.2);
plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2.5);
jit=(rand(size(d))-0.5)*0.26;
scatter(x+jit,d,26,c,'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
end
