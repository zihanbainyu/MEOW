function make_prelim_figs()
base='/Users/bai/Downloads/NeurWEM_prelim';
D=fullfile(base,'data'); F=fullfile(base,'figures');
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
cBlue=[0.16 0.44 0.74]; cRed=[0.80 0.25 0.22]; cGray=[0.55 0.55 0.55];

%% ---------- Figure 1b: graded similarity in early visual (same>similar>null) ----------
T=readtable(fullfile(D,'fig1b_graded.csv'));
roi='EarlyVis'; cats={'same','similar','null'};
M=nan(0,3); subs=unique(T.subj);
for s=1:numel(subs)
    row=nan(1,3);
    for c=1:3
        v=T.mean_r(T.subj==subs(s)&strcmp(T.roi,roi)&strcmp(T.category,cats{c}));
        if ~isempty(v), row(c)=v(1); end
    end
    M(end+1,:)=row; %#ok<AGROW>
end
f=figure('Color','w','Position',[50 50 560 620]); hold on;
cols=[0.15 0.15 0.15; 0.45 0.45 0.45; 0.78 0.78 0.78];
for c=1:3, house_box(c, M(:,c), cols(c,:)); end
for s=1:size(M,1), plot(1:3, M(s,:),'-','Color',[cGray 0.5],'LineWidth',0.75); end
set(gca,'XTick',1:3,'XTickLabel',{'Same','Similar','Null'},'FontSize',18,'LineWidth',1,'TickDir','out');
ylabel('Pattern similarity (r)','FontSize',18);
title('Early visual cortex','FontSize',20,'FontWeight','normal');
xlim([0.4 3.6]); box off; yline(0,'k--','LineWidth',1);
p=signrank_safe(M(:,1),M(:,2)); p2=signrank_safe(M(:,2),M(:,3));
annotate_bracket(1,2,max(M(:))*1.05,star(p)); annotate_bracket(2,3,max(M(:))*1.12,star(p2));
hold off; savepdf(f,fullfile(F,'fig1b_graded_similarity.pdf')); close(f);
fprintf('fig1b: same=%.3f similar=%.3f null=%.3f (n=%d)\n',nanmean(M),size(M,1));

%% ---------- Figure 2b: univariate discrimination effect (correct-incorrect) x condition ----------
T=readtable(fullfile(D,'fig2_glm.csv'));
rois={'HC','Angular','DLPFC','IPS'}; roilab={'Hippocampus','Angular','DLPFC','IPS'};
subs=unique(T.subj);
eff=struct();
for r=1:numel(rois)
    for cc=1:2
        condn={'compared','novel'}; co=condn{cc};
        d=nan(numel(subs),1);
        for s=1:numel(subs)
            v=T.effect_pscb(T.subj==subs(s)&strcmp(T.roi,rois{r})&strcmp(T.condition,co));
            if ~isempty(v), d(s)=v(1); end
        end
        eff.(sprintf('r%d_c%d',r,cc))=d;
    end
end
f=figure('Color','w','Position',[50 50 1050 620]); hold on;
xt=[]; xl={};
for r=1:numel(rois)
    x0=(r-1)*3;
    dc=eff.(sprintf('r%d_c%d',r,1)); dn=eff.(sprintf('r%d_c%d',r,2));
    house_box(x0+1, dc, cBlue); house_box(x0+2, dn, cRed);
    for s=1:numel(dc), if ~isnan(dc(s))&&~isnan(dn(s)), plot([x0+1 x0+2],[dc(s) dn(s)],'-','Color',[cGray 0.4],'LineWidth',0.6); end, end
    xt=[xt x0+1.5]; xl{end+1}=roilab{r}; %#ok<AGROW>
end
yline(0,'k--','LineWidth',1);
set(gca,'XTick',xt,'XTickLabel',xl,'FontSize',17,'LineWidth',1,'TickDir','out');
ylabel('Discrimination effect (correct - incorrect, % signal)','FontSize',17);
title('Univariate GLM: accuracy modulation by condition','FontSize',20,'FontWeight','normal');
xlim([0.2 numel(rois)*3-0.2]); box off;
text(0.02,0.98,'\color[rgb]{0.16,0.44,0.74}compared   \color[rgb]{0.80,0.25,0.22}novel',...
    'Units','normalized','VerticalAlignment','top','FontSize',16);
hold off; savepdf(f,fullfile(F,'fig2_interaction.pdf')); close(f);
for r=1:numel(rois)
    dc=eff.(sprintf('r%d_c%d',r,1)); dn=eff.(sprintf('r%d_c%d',r,2));
    fprintf('fig2 %-10s compared d+=%.3f  novel d+=%.3f  interaction(c-n)=%.3f\n',rois{r},nanmean(dc),nanmean(dn),nanmean(dc-dn));
end

%% ---------- Figure 3: multivariate mechanisms ----------
T=readtable(fullfile(D,'fig3_rsa.csv'));
mechs={'separation','completion','predrecall'};
mroi={'HC_post','HC_ant','Angular'};  % a-priori region per mechanism
mlab={'Separation (A1\midB1)','Completion (B2\midB1)','Predictive recall (A2\midB1)'};
mrlab={'post. hippocampus','ant. hippocampus','angular'};
subs=unique(T.subj);
f=figure('Color','w','Position',[40 40 1250 560]);
for m=1:3
    subplot(1,3,m); hold on;
    C=nan(numel(subs),1); I=nan(numel(subs),1);
    for s=1:numel(subs)
        mc=T.pair_r(T.subj==subs(s)&strcmp(T.roi,mroi{m})&strcmp(T.mechanism,mechs{m})&T.correct==1);
        mi=T.pair_r(T.subj==subs(s)&strcmp(T.roi,mroi{m})&strcmp(T.mechanism,mechs{m})&T.correct==0);
        if ~isempty(mc), C(s)=mean(mc); end
        if ~isempty(mi), I(s)=mean(mi); end
    end
    house_box(1,C,cBlue); house_box(2,I,cRed);
    for s=1:numel(C), if ~isnan(C(s))&&~isnan(I(s)), plot([1 2],[C(s) I(s)],'-','Color',[cGray 0.4],'LineWidth',0.6); end, end
    yline(0,'k--','LineWidth',1);
    set(gca,'XTick',[1 2],'XTickLabel',{'correct','incorrect'},'FontSize',15,'LineWidth',1,'TickDir','out');
    if m==1, ylabel('Pattern similarity (r)','FontSize',16); end
    title({mlab{m};['\rm\fontsize{13}' mrlab{m}]},'FontSize',16,'FontWeight','normal');
    xlim([0.4 2.6]); box off;
    p=signrank_safe(C,I);
    yy=max([C;I])*1.05; annotate_bracket(1,2,yy,star(p));
    fprintf('fig3 %-11s %-9s correct=%.3f incorrect=%.3f d=%.2f p=%.3f\n',mechs{m},mroi{m},nanmean(C),nanmean(I),cohen_d(C,I),p);
end
savepdf(f,fullfile(F,'fig3_rsa_mechanisms.pdf')); close(f);
fprintf('ALL FIGURES WRITTEN to %s\n',F);
end

function savepdf(f,path)
set(f,'Units','points'); pos=get(f,'Position');
set(f,'PaperUnits','points','PaperPosition',[0 0 pos(3) pos(4)],'PaperSize',[pos(3) pos(4)]);
print(f,path,'-dpdf','-vector');
end

function house_box(x,d,c)
d=d(~isnan(d)); if isempty(d), return; end
w=0.28;
if numel(d)>=2
    q=quantile(d,[0.25 0.5 0.75]);
    patch([x-w x+w x+w x-w],[q(1) q(1) q(3) q(3)],c,'FaceAlpha',0.18,'EdgeColor',c,'LineWidth',1.2);
    plot([x-w x+w],[q(2) q(2)],'-','Color',c,'LineWidth',2.5);
end
jit=(rand(size(d))-0.5)*0.22;
scatter(x+jit,d,46,c,'filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','w','LineWidth',0.5);
end

function p=signrank_safe(a,b)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m);
if numel(a)<2, p=NaN; return; end
try, p=signrank(a,b); catch, [~,p]=ttest(a,b); end
end

function d=cohen_d(a,b)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m); df=a-b;
if numel(df)<2||std(df)==0, d=NaN; else, d=mean(df)/std(df); end
end

function s=star(p)
if isnan(p), s=''; elseif p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s='n.s.'; end
end

function annotate_bracket(x1,x2,y,txt)
if isempty(txt)||isnan(y), return; end
plot([x1 x2],[y y],'k-','LineWidth',1);
text((x1+x2)/2,y,txt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16);
end
