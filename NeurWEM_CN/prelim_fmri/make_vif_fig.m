function make_vif_fig()
base='/Users/bai/Downloads/NeurWEM_prelim';
T=readtable(fullfile(base,'data','lss_vif_n6.csv'));
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22]; cGray=[0.55 0.55 0.55];
subs=unique(T.subj); ns=numel(subs);
cap=@(v) min(max(v,1),1000);                      % clip pathological 1e8 singulars for display
f=figure('Color','w','Position',[40 40 1050 560]); hold on;
xt=[]; xl={};
for s=1:ns
    x0=(s-1)*3;
    L=log10(cap(T.vif_lss(T.subj==subs(s))));
    A=log10(cap(T.vif_lsa(T.subj==subs(s))));
    house_box(x0+1,L,cBlue); house_box(x0+2,A,cRed);
    xt=[xt x0+1.5]; xl{end+1}=sprintf('sub-%d',subs(s)); %#ok<AGROW>
end
yline(log10(5),'k--','LineWidth',1);
text(0.3,log10(5),'  VIF = 5','VerticalAlignment','bottom','FontSize',12,'Color',[0.3 0.3 0.3]);
set(gca,'XTick',xt,'XTickLabel',xl,'FontSize',14,'LineWidth',1,'TickDir','out');
yt=[1 2 5 10 30 100 1000]; set(gca,'YTick',log10(yt),'YTickLabel',arrayfun(@num2str,yt,'uni',0));
ylabel('VIF of single-trial regressor','FontSize',15);
title('Single-trial design collinearity: LS-S vs LS-A (same designs)','FontSize',18,'FontWeight','normal');
xlim([0.2 ns*3-0.2]); box off;
text(0.015,0.97,'\color[rgb]{0.16,0.44,0.74}LS-S (used)   \color[rgb]{0.80,0.25,0.22}LS-A','Units','normalized','VerticalAlignment','top','FontSize',15);
savepdf(f,fullfile(base,'figures','fig_lss_vif_n6.pdf')); close(f);
% report
fprintf('subj  medSOA  LSS[med 90th]   LSA[med 90th]   %%LSS>5  %%LSA>5\n');
for s=1:ns
    d=T(T.subj==subs(s),:);
    fprintf('%4d   %5.2f   %4.2f %4.2f      %5.1f %6.1f     %4.1f   %4.1f\n',subs(s),median(d.soa),...
        median(d.vif_lss),quantile(d.vif_lss,.9),median(d.vif_lsa),quantile(d.vif_lsa,.9),...
        100*mean(d.vif_lss>5),100*mean(d.vif_lsa>5));
end
fprintf('wrote fig_lss_vif_n6.pdf\n');
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
jit=(rand(size(d))-0.5)*0.26;
scatter(x+jit,d,8,c,'filled','MarkerFaceAlpha',0.18,'MarkerEdgeColor','none');
end
