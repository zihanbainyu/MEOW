function make_roi_paired()
% Pair each ROI locator with its RSA-mechanism result (correct vs incorrect).
% Separation = A1-B1 (1-back pairmate similarity), Completion = B2-B1, Predictive recall = A2-B1.
base='/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/viewer';
T=readtable(fullfile(base,'data','fig3_rsa.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cB=[0.16 0.44 0.74]; cR=[0.80 0.25 0.22]; cG=[0.55 0.55 0.55];
% roi_png | csv_roi | mechanism | title
rows={'HCpost','HC_post','separation','Separation  (A1\midB1, 1-back)';
      'HCant','HC_ant','completion','Completion  (B2\midB1)';
      'Angular','Angular','predrecall','Predictive recall  (A2\midB1)'};
subs=unique(T.subj);
f=figure('Color','w','Position',[30 30 820 820]);
for k=1:3
    % --- left: ROI locator image ---
    subplot(3,3,(k-1)*3+1);
    img=imread(fullfile(base,'roi',['roi_' rows{k,1} '.png']));
    imshow(img); axis image off;
    % --- right: correct vs incorrect box ---
    subplot(3,3,[(k-1)*3+2 (k-1)*3+3]); hold on;
    C=nan(numel(subs),1); I=nan(numel(subs),1);
    for s=1:numel(subs)
        C(s)=mean(T.pair_r(T.subj==subs(s)&strcmp(T.roi,rows{k,2})&strcmp(T.mechanism,rows{k,3})&T.correct==1),'omitnan');
        I(s)=mean(T.pair_r(T.subj==subs(s)&strcmp(T.roi,rows{k,2})&strcmp(T.mechanism,rows{k,3})&T.correct==0),'omitnan');
    end
    hb(1,C,cB); hb(2,I,cR);
    for s=1:numel(subs), if ~isnan(C(s))&&~isnan(I(s)), plot([1 2],[C(s) I(s)],'-','Color',[cG 0.4],'LineWidth',0.6); end, end
    yline(0,'k:','LineWidth',0.75);
    set(gca,'XTick',[1 2],'XTickLabel',{'correct','incorrect'},'FontSize',13,'LineWidth',1,'TickDir','out');
    ylabel('pattern similarity (r)','FontSize',12); xlim([0.5 2.5]); box off;
    d=nanmean(C-I)/nanstd(C-I); p=sr(C,I);
    title(sprintf('%s     d=%+.2f, p=%.2f',rows{k,4},d,p),'FontSize',13,'FontWeight','normal');
end
sgtitle('RSA mechanisms: ROI + effect  (N=5, points = subjects)','FontSize',15,'FontWeight','normal');
savepdf(f,fullfile(base,'figures','fig_roi_paired_rsa.pdf')); close(f);
fprintf('wrote fig_roi_paired_rsa.pdf\n');
end
function savepdf(f,p)
set(f,'Units','points'); q=get(f,'Position');
set(f,'PaperUnits','points','PaperPosition',[0 0 q(3) q(4)],'PaperSize',[q(3) q(4)]); print(f,p,'-dpdf','-vector');
end
function hb(x,d,c)
d=d(~isnan(d)); if isempty(d), return; end
q=quantile(d,[.25 .5 .75]); w=0.3;
patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',.15,'EdgeColor',c,'LineWidth',1.2); hold on;
plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2.5);
scatter(x+(rand(size(d))-.5)*.2,d,44,c,'filled','MarkerFaceAlpha',.85,'MarkerEdgeColor','w');
end
function p=sr(a,b)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m); if numel(a)<2, p=NaN; return; end
try, p=signrank(a,b); catch, [~,p]=ttest(a,b); end
end
