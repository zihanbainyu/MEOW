function nm_emwm_behav(outdir)
% Additional EM x WM behavioral analyses (reviewer-driven), subjects 101-109.
% Produces four house-style figures + a short console report.
% Depends on nm_load_long.m, nm_metrics.m on the path.
if nargin<1 || isempty(outdir), outdir=pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end
set(0,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
subs=101:109; nS=numel(subs);

% palette (matches S_group_mstback)
c_comp=[87 6 140]/255; c_nov=[183 210 205]/255;
c_same=[97 125 184]/255; c_sim=[255 191 205]/255; c_new=[219 219 219]/255;
c_rem=[60 120 70]/255; c_forg=[200 110 60]/255;

% ---- gather per-subject metrics ----
dSim=nan(nS,2); cSim=nan(nS,2); di=nan(nS,2); na=nan(nS,2); rt=nan(nS,2); miss=nan(nS,2);
goalAcc=nan(nS,3,2);                 % subj x [AA AB AN] x [compared novel]
confLV=nan(nS,3);                    % 2-back A-B p("same" confusion) by lure l2 l3 l4
recByLV=nan(nS,3);                   % recognition accuracy by lure l2 l3 foil
mst_d=nan(nS,1); rec_d=nan(nS,1);
has_rec=false(nS,1); has_mst=false(nS,1);
% WM->EM coupling (recognition): mean 2-back B acc/RT for later-remembered vs forgotten
wmAcc_rf=nan(nS,2); wmRT_rf=nan(nS,2);   % cols: [remembered forgotten]
D=cell(nS,1);
for i=1:nS
    Dc=nm_load_long(subs(i)); D{i}=Dc; if isempty(Dc), continue; end
    M=nm_metrics(Dc);
    dSim(i,:)=M.dSim; cSim(i,:)=M.cSim; di(i,:)=M.di; na(i,:)=M.na; rt(i,:)=M.rt; miss(i,:)=M.miss;
    [goalAcc(i,:,:),confLV(i,:)]=goal_and_conf(Dc);
    if ~isempty(Dc.mst),  has_mst(i)=true; mst_d(i)=mst_dprime(Dc.mst); end
    if ~isempty(Dc.rec),  has_rec(i)=true;
        [rec_d(i),recByLV(i,:)]=rec_sdt(Dc.rec);
        [wmAcc_rf(i,:),wmRT_rf(i,:)]=wm_to_em(Dc);
    end
end

% ================= REPORT =================
fprintf('\n===== EM x WM additional behavior  (N=%d: %s) =====\n',nS,num2str(subs));
pr('2-back d-prime (similar)  compared vs novel', dSim);
pr('2-back criterion c         compared vs novel', cSim);
pr('2-back DI (p[sim|lure]-p[sim|foil])          ', di);
pr('no-response rate            compared vs novel', na);
pr('median correct RT (s)       compared vs novel', rt);
fprintf('  -> sensitivity vs bias: if the advantage is real it sits in d-prime with c matched.\n');
fprintf('\nWM->EM coupling (recognition subjects, N=%d): 2-back B-acc remembered vs forgotten\n',sum(has_rec));
prc('  2-back acc  rem vs forg', wmAcc_rf(has_rec,:));
prc('  2-back RT   rem vs forg', wmRT_rf(has_rec,:));

% ================= FIGURE A: sensitivity vs bias =================
fA=figure('Color','w','Position',[40 60 1180 300],'Name','FigA sensitivity vs bias');
tl=tiledlayout(fA,1,4,'Padding','compact','TileSpacing','compact');
nexttile; pp(dSim,{'compared','novel'},{c_comp,c_nov},'d'' (similar)','Sensitivity',subs,{[1 2]});
nexttile; pp(cSim,{'compared','novel'},{c_comp,c_nov},'criterion c','Response bias',subs,{[1 2]});
nexttile; pp(na,{'compared','novel'},{c_comp,c_nov},'no-response rate','Missingness (matched?)',subs,{[1 2]});
nexttile; pp(rt,{'compared','novel'},{c_comp,c_nov},'median correct RT (s)','Speed (no tradeoff?)',subs,{[1 2]});
title(tl,'2-back similar discrimination: is compared>novel sensitivity or bias?','FontWeight','normal','FontSize',14);
savef(fA,outdir,'emwm_figA_sensitivity_vs_bias');

% ================= FIGURE B: WM -> EM coupling =================
fB=figure('Color','w','Position',[40 60 640 320],'Name','FigB WM->EM coupling');
tl=tiledlayout(fB,1,2,'Padding','compact','TileSpacing','compact');
nexttile; pp(wmAcc_rf(has_rec,:),{'remembered','forgotten'},{c_rem,c_forg},'2-back accuracy','WM accuracy by later memory',subs(has_rec),{[1 2]});
nexttile; pp(wmRT_rf(has_rec,:),{'remembered','forgotten'},{c_rem,c_forg},'2-back RT (s)','WM RT by later memory',subs(has_rec),{[1 2]});
title(tl,sprintf('Compared items: WM processing predicts later recognition  (N=%d)',sum(has_rec)),'FontWeight','normal','FontSize',13);
savef(fB,outdir,'emwm_figB_wm_to_em_coupling');

% ================= FIGURE C: similarity gradient =================
fC=figure('Color','w','Position',[40 60 680 320],'Name','FigC similarity gradient');
tl=tiledlayout(fC,1,2,'Padding','compact','TileSpacing','compact');
nexttile; grad(confLV,{'l2','l3','l4'},'p("same" confusion)','2-back: A-B confusability',subs,c_sim);
nexttile; grad(recByLV,{'l2','l3','foil'},'recognition accuracy','later memory by similarity',subs(has_rec),c_same,has_rec);
title(tl,'Mnemonic-similarity gradient (harder l2 \rightarrow easier)','FontWeight','normal','FontSize',13);
savef(fC,outdir,'emwm_figC_similarity_gradient');

% ================= FIGURE D: pairmate specificity + individual diffs =================
fD=figure('Color','w','Position',[40 60 760 320],'Name','FigD specificity + indiv diff');
tl=tiledlayout(fD,1,2,'Padding','compact','TileSpacing','compact');
nexttile; goalfig(goalAcc,{'A-A','A-B','A-N'},{c_comp,c_nov},subs);
nexttile; scatterID(di(:,1)-di(:,2), mst_d, has_mst, subs, c_comp);
title(tl,'Pairmate specificity (left) and WM benefit \times EM ability (right)','FontWeight','normal','FontSize',13);
savef(fD,outdir,'emwm_figD_specificity_indivdiff');

fprintf('\nwrote 4 figures to %s\n',outdir);
end

% ============================ helpers ============================
function [ga,conf]=goal_and_conf(Dc)
two=Dc.two; g=string(two.goal); cond=string(two.condition); cr=string(two.corr_resp); rk=string(two.resp_key);
blk=two.block; pan=false(height(two),1);
for i=1:height(two)-2, if g(i)=="A-N" && blk(i+2)==blk(i), pan(i+2)=true; end, end
aa=g=="A-A"&cr=="1"; ab=g=="A-B"&cr=="2"; an=pan&cr=="none";
conds={'compared','novel'}; ga=nan(3,2);
for c=1:2
    cm=cond==conds{c};
    ga(1,c)=accm(two.correct,aa&cm); ga(2,c)=accm(two.correct,ab&cm); ga(3,c)=accm(two.correct,an&cm);
end
% A-B confusion = calling the similar pairmate "same" (resp '1'), by lure level
lv=string(two.lure); conf=nan(1,3); lvls={'l2','l3','l4'};
for k=1:3
    m=ab & lv==lvls{k};
    if sum(m)>0, conf(k)=mean(rk(m)=="1"); end
end
end

function [d,byLV]=rec_sdt(T)
tt=string(T.trial_type); rk=string(T.resp_key); cond=string(T.condition); lv=string(T.lure);
old=tt=="old"; new=tt=="new";
H=mean(rk(old)=="1"); FA=mean(rk(new)=="1");
nH=sum(old); nF=sum(new);
d=zc(H,nH)-zc(FA,nF);
byLV=nan(1,3); lvls={'l2','l3','foil'};
for k=1:3
    m=old & lv==lvls{k};
    if sum(m)>0, byLV(k)=mean(T.correct(m)); end
end
end

function d=mst_dprime(T)
tt=string(T.trial_type); rk=string(T.resp_key);
old=tt=="old"; new=tt=="new";
H=mean(rk(old)=="1"); FA=mean(rk(new)=="1");
if (H==0 && FA==0) || (H==1 && FA==1), d=NaN; return; end   % degenerate: no "old" endorsements (withholding) -> undefined
d=zc(H,sum(old))-zc(FA,sum(new));
end

function [acc_rf,rt_rf]=wm_to_em(Dc)
% Link each compared item's 2-back B-trial (A-B, identity B) processing to its
% later recognition fate (tested as "old": correct = hit, else forgotten).
two=Dc.two; g=string(two.goal); id=string(two.identity); cond=string(two.condition);
isB = cond=="compared" & g=="A-B" & id=="B";
wmItem=two.item(isB); wmAcc=double(two.correct(isB)); wmRT=two.rt(isB);
rec=Dc.rec; rtt=string(rec.trial_type); rcond=string(rec.condition);
isOld = rcond=="compared" & rtt=="old";
emItem=rec.item(isOld); emHit=double(rec.correct(isOld));
accR=[]; accF=[]; rtR=[]; rtF=[];
for j=1:numel(wmItem)
    it=wmItem(j); k=find(emItem==it,1);
    if isempty(k), continue; end
    if emHit(k)==1, accR(end+1)=wmAcc(j); rtR(end+1)=wmRT(j); %#ok<AGROW>
    else,           accF(end+1)=wmAcc(j); rtF(end+1)=wmRT(j); %#ok<AGROW>
    end
end
acc_rf=[mn(accR) mn(accF)];
rt_rf =[median(rtR(rtR>0.15),'omitnan') median(rtF(rtF>0.15),'omitnan')];
end

function pp(M,labels,cols,ylab,ttl,subs,pairs)
[nS,nL]=size(M); hold on; bw=0.34; jw=0.06;
scol=lines(max(nS,1)); mks={'o','s','d','^','v','p','h','>','<'};
X=(1:nL)+(rand(nS,nL)-0.5)*2*jw;
for s=1:nS
    if all(~isnan(M(s,:))), plot(X(s,:),M(s,:),'-','Color',[scol(s,:) 0.45],'LineWidth',1.0); end
end
for L=1:nL
    d=M(:,L); d=d(~isnan(d));
    if numel(d)>=2
        q=quantile(d,[.25 .5 .75]); iqrv=q(3)-q(1);
        lo=max(min(d),q(1)-1.5*iqrv); hi=min(max(d),q(3)+1.5*iqrv);
        plot([L L],[lo q(1)],'k-','LineWidth',0.9); plot([L L],[q(3) hi],'k-','LineWidth',0.9);
        rectangle('Position',[L-bw/2 q(1) bw max(iqrv,eps)],'EdgeColor','k','LineWidth',1.0);
        plot(L+[-bw/2 bw/2],[q(2) q(2)],'-','Color',cols{L},'LineWidth',2.6);
    end
    for s=1:nS
        if ~isnan(M(s,L))
            scatter(X(s,L),M(s,L),58,cols{L},mks{mod(s-1,9)+1},'filled',...
                'MarkerFaceAlpha',0.85,'MarkerEdgeColor',[.25 .25 .25],'LineWidth',0.3);
        end
    end
end
set(gca,'XTick',1:nL,'XTickLabel',labels,'TickDir','out','FontSize',12,'LineWidth',1);
xlim([0.4 nL+0.6]); ylabel(ylab,'FontSize',12); box off;
if ~isempty(ttl), title(ttl,'FontSize',13,'FontWeight','normal'); end
allv=M(~isnan(M)); if isempty(allv), return; end
yr=range(allv); if yr==0, yr=1; end
if nargin>=7 && ~isempty(pairs) && nS>=5
    base=max(allv)+0.10*yr;
    for pp_i=1:numel(pairs)
        ij=pairs{pp_i}; ok=~isnan(M(:,ij(1)))&~isnan(M(:,ij(2)));
        if sum(ok)<5, continue; end
        pv=signrank(M(ok,ij(1)),M(ok,ij(2)));
        plot([ij(1) ij(1) ij(2) ij(2)],[base-0.02*yr base base base-0.02*yr],'k-','LineWidth',0.9);
        text(mean(ij),base,sprintf('p=%.3f',pv),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11);
    end
    ylim([min(allv)-0.08*yr base+0.12*yr]);
end
hold off;
end

function grad(M,labels,ylab,ttl,subs,col,mask)
if nargin>=7 && ~isempty(mask), M=M(mask,:); end
keep=any(~isnan(M),1); M=M(:,keep); labels=labels(keep);   % drop empty lure levels
[nS,nL]=size(M); hold on;
for s=1:nS
    if any(~isnan(M(s,:))), plot(1:nL,M(s,:),'-','Color',[0.6 0.6 0.6 0.5],'LineWidth',0.9,'Marker','o','MarkerSize',4,'MarkerFaceColor','w'); end
end
mu=mean(M,1,'omitnan'); se=std(M,0,1,'omitnan')./sqrt(sum(~isnan(M),1));
errorbar(1:nL,mu,se,'-','Color',col,'LineWidth',2.4,'Marker','o','MarkerSize',8,'MarkerFaceColor',col,'CapSize',0);
set(gca,'XTick',1:nL,'XTickLabel',labels,'TickDir','out','FontSize',12,'LineWidth',1);
xlim([0.6 nL+0.4]); ylabel(ylab,'FontSize',12); box off;
title(ttl,'FontSize',13,'FontWeight','normal');
% slope test across ordered levels
x=repmat(1:nL,nS,1); y=M; ok=~isnan(y);
if sum(ok(:))>=4 && numel(unique(x(ok)))>=2
    [r,pv]=corr(x(ok),y(ok)); yl=ylim;
    text(0.65,yl(1)+0.04*range(yl),sprintf('slope r=%.2f, p=%.3f',r,pv),'FontSize',10,'VerticalAlignment','bottom');
end
hold off;
end

function goalfig(GA,labels,cols,subs)
% GA: subj x 3(goal) x 2(cond). Plot compared vs novel accuracy per goal.
nS=size(GA,1); hold on; offs=[-0.16 0.16]; cc={cols{1},cols{2}};
for c=1:2
    for L=1:3
        xx=L+offs(c); d=squeeze(GA(:,L,c)); d=d(~isnan(d));
        if numel(d)>=2
            q=quantile(d,[.25 .5 .75]);
            rectangle('Position',[xx-0.12 q(1) 0.24 max(q(3)-q(1),eps)],'EdgeColor',cc{c},'LineWidth',1.1);
            plot(xx+[-0.12 0.12],[q(2) q(2)],'-','Color',cc{c},'LineWidth',2.2);
        end
        scatter(xx+(rand(nS,1)-0.5)*0.12,squeeze(GA(:,L,c)),40,cc{c},'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor',[.3 .3 .3],'LineWidth',0.3);
    end
end
set(gca,'XTick',1:3,'XTickLabel',labels,'TickDir','out','FontSize',12,'LineWidth',1);
xlim([0.5 3.5]); ylabel('2-back accuracy','FontSize',12); box off;
title('Benefit by goal (purple=compared, sage=novel)','FontSize',13,'FontWeight','normal');
hold off;
end

function scatterID(wmBenefit, emAbility, mask, subs, col)
x=wmBenefit(mask); y=emAbility(mask); sid=subs(mask);
hold on; scatter(x,y,70,col,'filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor',[.25 .25 .25]);
for i=1:numel(x), text(x(i),y(i),sprintf('  %d',sid(i)),'FontSize',9,'Color',[.3 .3 .3]); end
ok=~isnan(x)&~isnan(y);
if sum(ok)>=3
    cf=polyfit(x(ok),y(ok),1); xx=[min(x(ok)) max(x(ok))];
    plot(xx,polyval(cf,xx),'k--','LineWidth',1.4);
    [r,pv]=corr(x(ok),y(ok)); yl=ylim; xl=xlim;
    text(xl(1)+0.04*range(xl),yl(2)-0.06*range(yl),sprintf('r=%.2f, p=%.2g (N=%d, exploratory)',r,pv,sum(ok)),'FontSize',10,'VerticalAlignment','top');
end
set(gca,'TickDir','out','FontSize',12,'LineWidth',1); box off;
xlabel('WM benefit: DI compared - novel','FontSize',12); ylabel('EM ability: MST d'' (old)','FontSize',12);
title('WM benefit vs EM ability','FontSize',13,'FontWeight','normal'); hold off;
end

function savef(f,outdir,name)
set(f,'Renderer','painters','PaperPositionMode','auto');
pos=get(f,'Position'); set(f,'PaperUnits','points','PaperSize',pos(3:4));
print(f,fullfile(outdir,[name '.png']),'-dpng','-r170');
print(f,fullfile(outdir,[name '.pdf']),'-dpdf','-vector');
end

function a=accm(c,m), if sum(m)==0, a=NaN; else, a=mean(c(m)); end, end
function v=mn(x), if isempty(x), v=NaN; else, v=mean(x); end, end
function z=zc(p,n), if isnan(p)||n==0, z=NaN; return; end, z=norminv(max(1/(2*n),min(1-1/(2*n),p))); end
function pr(lbl,M), fprintf('  %-46s %.3f / %.3f  (diff %+.3f)\n',lbl,mean(M(:,1),'omitnan'),mean(M(:,2),'omitnan'),mean(M(:,1)-M(:,2),'omitnan')); end
function prc(lbl,M)
ok=~isnan(M(:,1))&~isnan(M(:,2));
if sum(ok)>=5, pv=signrank(M(ok,1),M(ok,2)); ps=sprintf('signed-rank p=%.3f',pv); else, ps='(n<5)'; end
fprintf('  %-26s %.3f / %.3f  (diff %+.3f, %s)\n',lbl,mean(M(:,1),'omitnan'),mean(M(:,2),'omitnan'),mean(M(:,1)-M(:,2),'omitnan'),ps);
end
