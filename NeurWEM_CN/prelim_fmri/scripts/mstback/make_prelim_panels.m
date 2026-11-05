function make_prelim_panels()
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
dat=fullfile(root,'data','mstback'); datm=fullfile(root,'data','mst'); out=fullfile(root,'figures','mstback'); outm=fullfile(root,'figures','mst'); addpath(fullfile(root,'..','scripts'));
W4=5; W2=3.9; H=3.8;
T=readtable(fullfile(dat,'roi_2x2_cells.csv')); M=load(fullfile(datm,'m3_roi.mat'));
HT=readtable(fullfile(dat,'rsa_hyp_cells.csv')); ST=readtable(fullfile(dat,'rsa_hyp_stats.csv'));
RT=readtable(fullfile(dat,'rsa_reinstate_hcsplit.csv'));

%%%%%%%%%%%%% MST-back
% two-back behaviour: similar-discrimination index, compared vs novel
behav_combined(here,[4.4 4.4],'di','di_one','Similar-discrimination index','Similar discrimination','fig_behav_2back');

% two-back 2x2  ROI: anterior hippocampus
roi_2x2(T,here,[W4 H],'HC_ant','Anterior hippocampus',[1 0],'roi_HCant');

% two-back 2x2  ROI: angular gyrus
roi_2x2(T,here,[W4 H],'Angular','Angular gyrus',[1 0],'roi_Angular');

% two-back 2x2  ROI: DLPFC
roi_2x2(T,here,[W4 H],'DLPFC','DLPFC',[0 1],'roi_DLPFC');

% two-back 2x2  ROI: superior parietal (IPS)
roi_2x2(T,here,[W4 H],'IPS','Superior parietal (IPS)',[0 1],'roi_IPS');

% two-back 2x2  ROI: visual cortex
roi_2x2(T,here,[W4 H],'Visual','Visual cortex',[1 1],'roi_Visual');

% two-back 2x2  ROI: ventral temporal
roi_2x2(T,here,[W4 H],'VentralTemporal','Ventral temporal',[1 1],'roi_VentralTemporal');

%%%%%%%%%%%%% MST-back control: same-detection (goal A-A, identical repeats)
% two-back behaviour: same-detection d', compared vs novel
behav_combined(here,[4.4 4.4],'dSame','dSame_one','same-detection d''','Same detection','fig_behav_2back_same');

% repeat-detection 2x2  ROI: anterior hippocampus
roi_2x2(T,here,[W4 H],'HC_ant','Anterior hippocampus (repeat detection)',[0 0],'roi_HCant','same');

% repeat-detection 2x2  ROI: angular gyrus
roi_2x2(T,here,[W4 H],'Angular','Angular gyrus (repeat detection)',[0 0],'roi_Angular','same');

% repeat-detection 2x2  ROI: DLPFC
roi_2x2(T,here,[W4 H],'DLPFC','DLPFC (repeat detection)',[0 0],'roi_DLPFC','same');

% repeat-detection 2x2  ROI: superior parietal (IPS)
roi_2x2(T,here,[W4 H],'IPS','Superior parietal (repeat detection)',[0 0],'roi_IPS','same');

% repeat-detection 2x2  ROI: visual cortex
roi_2x2(T,here,[W4 H],'Visual','Visual cortex (repeat detection)',[0 0],'roi_Visual','same');

% repeat-detection 2x2  ROI: ventral temporal
roi_2x2(T,here,[W4 H],'VentralTemporal','Ventral temporal (repeat detection)',[0 0],'roi_VentralTemporal','same');

%%%%%%%%%%%%% MST
% MST behaviour: recognition & lure discrimination index
behav_mst(datm,here,[W4 H]);

% memory success  ROI: posterior hippocampus
mem_roi(M,here,[W2 H],'HC_post','Posterior hippocampus','roi_HCpost');

% memory success  ROI: angular gyrus
mem_roi(M,here,[W2 H],'Angular','Angular gyrus','roi_Angular');

% memory success  ROI: precuneus
mem_roi(M,here,[W2 H],'Precuneus','Precuneus','roi_Precuneus');

% memory success  ROI: visual cortex
mem_roi(M,here,[W2 H],'Visual','Visual cortex','roi_Visual');

%%%%%%%%%%%%% RSA reinstatement indices (item-specific, Fisher z; one figure per similarity, ROIs share a y-axis)
RROI={'HC_ant','HC_post','Angular','DLPFC','IPS','Visual','VentralTemporal'};
RLAB={'Ant HC','Post HC','Angular','DLPFC','IPS','Visual','Vent temp'};

% B2-B1 similarity (compared)
rsahyp_strip(HT,ST,here,RROI,RLAB,'completion','compared','b2b1','B2-B1 neural similarity');

% B1-A1 similarity (compared)
rsahyp_strip(HT,ST,here,RROI,RLAB,'separation','compared','b1a1','A1-B1 neural similarity');

% A2-B1 similarity (compared)
rsahyp_strip(HT,ST,here,RROI,RLAB,'predictive','compared','a2b1','A2-B1 neural similarity');

% A2-B2 similarity (compared)
rsahyp_strip(HT,ST,here,RROI,RLAB,'wm_discrim','compared','a2b2_compared','A2-B2 neural similarity (compared)');

% A2-B2 similarity (novel)
rsahyp_strip(HT,ST,here,RROI,RLAB,'wm_discrim','novel','a2b2_novel','A2-B2 neural similarity (novel)');

%%%%%%%%%%%%% Anterior-HC - cortex reinstatement coupling (B2-B1, trial-wise within subject)
% all trials, coupling vs zero (above chance)
rsa_couple0(RT,here,'HC_ant',{'Angular','DLPFC','VentralTemporal','Visual','vmPFC'},{'Angular','DLPFC','Vent temp','Visual','vmPFC'},'coupling_hcant_vs0','Anterior HC - cortex coupling (B2-B1)');
% split by accuracy, correct vs incorrect
rsa_couple(RT,here,'HC_ant',{'Angular','DLPFC','VentralTemporal','Visual','vmPFC'},{'Angular','DLPFC','Vent temp','Visual','vmPFC'},'coupling_hcant','Anterior HC - cortex coupling (B2-B1)');
end

function roi_2x2(T,here,SZ,roi,roilab,tails,icon,pre)
if nargin<8, pre='sim'; end
INK=[0.149 0.149 0.184]; PUR=[0.58 0.52 0.78]; GRY=[0.60 0.60 0.63];
s=unique(T.subj); g=@(c)arrayfun(@(x)getv(T,x,roi,c),s);
cc=g([pre '_cmp_hit']); ci=g([pre '_cmp_miss']); nc=g([pre '_nov_hit']); ni=g([pre '_nov_miss']);
x=[1 2 3.3 4.3]; D={cc,ci,nc,ni}; C=[PUR;PUR;GRY;GRY];
[f,ax]=canvas(here,SZ,icon,roilab);
for k=1:4, boxdot(ax,x(k),D{k},C(k,:),INK); end
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
allv=cell2mat(D(:)); ymax=max(allv,[],'omitnan'); ymin=min(allv,[],'omitnan'); yr=ymax-ymin;
yb=ymax+0.16*yr;
brklab(ax,x(1),x(2),yb,yr,ttp(cc,ci,tails(1)),'compared',INK);
brklab(ax,x(3),x(4),yb,yr,ttp(nc,ni,tails(2)),'novel',INK);
[~,pa]=ttest(((cc+nc)-(ci+ni))/2); [~,pint]=ttest((cc-ci)-(nc-ni));
annotation(f,'textbox',[0.12 0.785 0.85 0.05],'String',sprintf('accuracy  %s        condition \\times accuracy  %s',pp(pa),pp(pint)),'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',9,'Color',INK);
ax.XTick=x; ax.XTickLabel={'correct','incorrect','correct','incorrect'}; ax.XTickLabelRotation=0;
xlim(ax,[0.4 4.9]); ylim(ax,[ymin-0.10*yr ymax+0.48*yr]);
ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',11);
fn='fig_2x2'; if ~strcmp(pre,'sim'), fn=['fig_2x2' pre]; end
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',sprintf('%s_%s.pdf',fn,lower(roi))),'ContentType','vector'); close(f);
end

function mem_roi(M,here,SZ,roi,roilab,icon)
INK=[0.149 0.149 0.184]; PUR=[0.58 0.52 0.78];
ec=M.([roi '_em_c'])(:); ei=M.([roi '_em_i'])(:);
[f,ax]=canvas(here,SZ,icon,roilab);
boxdot(ax,1,ec,PUR,INK); boxdot(ax,2,ei,PUR,INK);
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
allv=[ec;ei]; ymax=max(allv); ymin=min(allv); yr=ymax-ymin;
brklab(ax,1,2,ymax+0.12*yr,yr,ttp(ec,ei,1),'',INK);
tk=niceticks(ymin,ymax,5);
ax.XTick=[1 2]; ax.XTickLabel={'correct','incorrect'}; ax.YTick=tk; ax.XTickLabelRotation=0;
xlim(ax,[0.4 2.6]); ylim(ax,[min(ymin,tk(1))-0.06*yr max(ymax+0.26*yr,tk(end)+0.01*yr)]);
ylabel(ax,'parameter estimate (\beta)','FontName','Helvetica','FontSize',11);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mst',sprintf('fig_mem_%s.pdf',lower(roi))),'ContentType','vector'); close(f);
end

function rsahyp_strip(HT,ST,here,ROI,LAB,contr,cond,tag,ylab)
INK=[0.149 0.149 0.184]; TARG=[0.768 0.306 0.321]; FOIL=[0.357 0.561 0.788];
nR=numel(ROI); W=1.4+0.92*nR; Hh=3.6;
f=figure('Color','w','Units','inches','Position',[1 1 W Hh]);
annotation(f,'textbox',[0.085 0.885 0.62 0.09],'String',ylab,'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Color',INK);
annotation(f,'textbox',[0.70 0.885 0.13 0.09],'String','correct','EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',TARG);
annotation(f,'textbox',[0.83 0.885 0.15 0.09],'String','incorrect','EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',FOIL);
ax=axes(f,'Position',[0.115 0.145 0.86 0.66]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
allv=[];
for k=1:nR
    cc=getH(HT,ROI{k},contr,cond,'correct'); ci=getH(HT,ROI{k},contr,cond,'incorrect');
    boxdot(ax,k-0.19,cc,TARG,INK); boxdot(ax,k+0.19,ci,FOIL,INK); allv=[allv;cc;ci];
end
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
ymax=max(allv,[],'omitnan'); ymin=min(allv,[],'omitnan'); yr=ymax-ymin;
for k=1:nR
    p=getp(ST,ROI{k},contr,cond);
    if ~isnan(p), text(ax,k,ymax+0.09*yr,pp(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',8,'Color',INK); end
end
ax.XTick=1:nR; ax.XTickLabel=LAB; ax.XTickLabelRotation=0;
xlim(ax,[0.4 nR+0.6]); ylim(ax,[ymin-0.10*yr ymax+0.32*yr]);
ylabel(ax,'pattern similarity (item-specific, Fisher z)','FontName','Helvetica','FontSize',10.5);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',sprintf('fig_rsa_%s.pdf',tag)),'ContentType','vector'); close(f);
end

function rsa_couple(RT,here,seed,CROI,CLAB,tag,ttl)
INK=[0.149 0.149 0.184]; TARG=[0.768 0.306 0.321]; FOIL=[0.357 0.561 0.788];
subs=unique(RT.subj); nR=numel(CROI); ZC=nan(numel(subs),nR); ZI=nan(numel(subs),nR);
for si=1:numel(subs)
    hm=RT.subj==subs(si)&strcmp(RT.roi,seed)&strcmp(RT.cell,'B2-B1_PC'); hv=RT.sim(hm); hf=RT.correct(hm);
    for k=1:nR
        cm=RT.subj==subs(si)&strcmp(RT.roi,CROI{k})&strcmp(RT.cell,'B2-B1_PC'); cv=RT.sim(cm);
        if numel(cv)~=numel(hv), continue; end
        for a=[1 0]
            sel=hf==a&isfinite(hv)&isfinite(cv);
            if sum(sel)>=5 && std(hv(sel))>0 && std(cv(sel))>0
                z=atanh(min(max(corr(hv(sel),cv(sel)),-0.9999),0.9999));
                if a==1, ZC(si,k)=z; else, ZI(si,k)=z; end
            end
        end
    end
end
W=1.4+0.92*nR; Hh=3.6; f=figure('Color','w','Units','inches','Position',[1 1 W Hh]);
annotation(f,'textbox',[0.085 0.885 0.60 0.09],'String',ttl,'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Color',INK);
annotation(f,'textbox',[0.70 0.885 0.13 0.09],'String','correct','EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',TARG);
annotation(f,'textbox',[0.83 0.885 0.15 0.09],'String','incorrect','EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',FOIL);
ax=axes(f,'Position',[0.115 0.145 0.86 0.66]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
for k=1:nR, boxdot(ax,k-0.19,ZC(:,k),TARG,INK); boxdot(ax,k+0.19,ZI(:,k),FOIL,INK); end
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
allv=[ZC(:);ZI(:)]; ymax=max(allv,[],'omitnan'); ymin=min(allv,[],'omitnan'); yr=ymax-ymin;
for k=1:nR
    p=ttp(ZC(:,k),ZI(:,k),1);
    text(ax,k,ymax+0.09*yr,pp(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',8,'Color',INK);
end
ax.XTick=1:nR; ax.XTickLabel=CLAB; ax.XTickLabelRotation=0;
xlim(ax,[0.4 nR+0.6]); ylim(ax,[ymin-0.10*yr ymax+0.32*yr]);
ylabel(ax,'HC_{ant}-cortex coupling (Fisher z)','FontName','Helvetica','FontSize',10.5);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',sprintf('fig_rsa_%s.pdf',tag)),'ContentType','vector'); close(f);
end

function rsa_couple0(RT,here,seed,CROI,CLAB,tag,ttl)
INK=[0.149 0.149 0.184]; PUR=[0.58 0.52 0.78];
subs=unique(RT.subj); nR=numel(CROI); ZA=nan(numel(subs),nR);
for si=1:numel(subs)
    hm=RT.subj==subs(si)&strcmp(RT.roi,seed)&strcmp(RT.cell,'B2-B1_PC'); hv=RT.sim(hm);
    for k=1:nR
        cm=RT.subj==subs(si)&strcmp(RT.roi,CROI{k})&strcmp(RT.cell,'B2-B1_PC'); cv=RT.sim(cm);
        if numel(cv)~=numel(hv), continue; end
        sel=isfinite(hv)&isfinite(cv);
        if sum(sel)>=5 && std(hv(sel))>0 && std(cv(sel))>0
            ZA(si,k)=atanh(min(max(corr(hv(sel),cv(sel)),-0.9999),0.9999));
        end
    end
end
W=1.4+0.92*nR; Hh=3.6; f=figure('Color','w','Units','inches','Position',[1 1 W Hh]);
annotation(f,'textbox',[0.085 0.885 0.80 0.09],'String',[ttl '  (all trials, vs 0)'],'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Color',INK);
ax=axes(f,'Position',[0.115 0.145 0.86 0.66]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
for k=1:nR, boxdot(ax,k,ZA(:,k),PUR,INK); end
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
allv=ZA(:); ymax=max(allv,[],'omitnan'); ymin=min(allv,[],'omitnan'); yr=ymax-ymin;
for k=1:nR
    v=ZA(:,k); v=v(isfinite(v)); if numel(v)<2, continue; end
    [~,p]=ttest(v,0,'Tail','right');
    text(ax,k,ymax+0.09*yr,pp(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',8,'Color',INK);
end
ax.XTick=1:nR; ax.XTickLabel=CLAB; ax.XTickLabelRotation=0;
xlim(ax,[0.4 nR+0.6]); ylim(ax,[ymin-0.10*yr ymax+0.32*yr]);
ylabel(ax,'HC_{ant}-cortex coupling (Fisher z)','FontName','Helvetica','FontSize',10.5);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',sprintf('fig_rsa_%s.pdf',tag)),'ContentType','vector'); close(f);
end

function v=getH(HT,roi,contr,cond,acc)
m=strcmp(HT.roi,roi)&strcmp(HT.contrast,contr)&strcmp(HT.condition,cond)&strcmp(HT.accuracy,acc);
v=HT.index(m);
end

function p=getp(ST,roi,contr,cond)
m=strcmp(ST.roi,roi)&strcmp(ST.contrast,contr)&strcmp(ST.condition,cond);
if any(m), p=ST.p_pred1tail(find(m,1)); else, p=NaN; end
end

function behav_combined(here,SZ,fld2,fld1,ylab,ttl,oname)
INK=[0.149 0.149 0.184]; PUR=[0.58 0.52 0.78]; GRY=[0.60 0.60 0.63]; ONE=[0.42 0.56 0.66];
subs=101:109; V=nan(numel(subs),3);
for i=1:numel(subs)
    Dc=nm_load_long(subs(i)); if isempty(Dc), continue; end
    m=nm_metrics(Dc); d2=m.(fld2); V(i,:)=[m.(fld1) d2(1) d2(2)];
end
% order: 1-back, novel, compared   (V columns are [1-back, compared, novel])
x=[1 1.8 2.6]; one=V(:,1); nov=V(:,3); cmp=V(:,2);
f=figure('Color','w','Units','inches','Position',[1 1 SZ]);
annotation(f,'textbox',[0.08 0.875 0.86 0.08],'String',ttl,'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Color',INK);
ax=axes(f,'Position',[0.165 0.135 0.72 0.72]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
boxdot(ax,x(1),one,ONE,INK); boxdot(ax,x(2),nov,GRY,INK); boxdot(ax,x(3),cmp,PUR,INK);
allv=V(~isnan(V)); ymax=max(allv); ymin=min(allv); yr=ymax-ymin;
if contains(fld2,'Same'), lo=ymin-0.14*yr; else, lo=min(ymin,0)-0.08*yr; yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8); end
brklab(ax,x(2),x(3),ymax+0.08*yr,yr,ttp(nov,cmp,0),'',INK);
brklab(ax,x(1),x(2),ymax+0.22*yr,yr,ttp(one,nov,0),'',INK);
ax.XTick=x; ax.XTickLabel={'1-back','2-back novel','2-back compared'}; ax.XTickLabelRotation=0;
xlim(ax,[0.6 3.0]); ylim(ax,[lo ymax+0.38*yr]);
ylabel(ax,ylab,'FontName','Helvetica','FontSize',11);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',[oname '.pdf']),'ContentType','vector'); close(f);
end

function behav_2back(here,SZ,fld,ylab,ttl,oname)
if nargin<3, fld='di'; ylab='similar-discrimination index'; ttl='Two-back discrimination'; oname='fig_behav_2back'; end
INK=[0.149 0.149 0.184]; PUR=[0.58 0.52 0.78]; GRY=[0.60 0.60 0.63];
subs=101:109; DI=nan(numel(subs),2);
for i=1:numel(subs), Dc=nm_load_long(subs(i)); if isempty(Dc), continue; end; m=nm_metrics(Dc); DI(i,:)=m.(fld); end
[f,ax]=canvas(here,SZ,'',ttl);
boxdot(ax,1,DI(:,1),PUR,INK); boxdot(ax,2,DI(:,2),GRY,INK);
yline(ax,0,':','Color',[0.72 0.72 0.75],'LineWidth',0.8);
allv=DI(~isnan(DI)); ymax=max(allv); ymin=min(allv); yr=ymax-ymin;
brklab(ax,1,2,ymax+0.12*yr,yr,ttp(DI(:,1),DI(:,2),0),'',INK);
ax.XTick=[1 2]; ax.XTickLabel={'compared','novel'}; ax.XTickLabelRotation=0;
xlim(ax,[0.4 2.6]); ylim(ax,[ymin-0.10*yr ymax+0.26*yr]);
ylabel(ax,ylab,'FontName','Helvetica','FontSize',11);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mstback',[oname '.pdf']),'ContentType','vector'); close(f);
end

function behav_mst(dat,here,SZ)
INK=[0.149 0.149 0.184]; TARG=[0.768 0.306 0.321]; FOIL=[0.357 0.561 0.788]; S=load(fullfile(dat,'mst_idx.mat'));
D={S.old_old(:),S.foil_old(:),S.lure_sim(:),S.foil_sim(:)}; x=[1 2 3.3 4.3]; C=[TARG;FOIL;TARG;FOIL];
[f,ax]=canvas(here,SZ,'','MST performance');
for k=1:4, boxdot(ax,x(k),D{k},C(k,:),INK); end
idxbr(ax,x(1),x(2),S.REC(:),0.80,'Recognition',INK);
idxbr(ax,x(3),x(4),S.LDI(:),0.80,'Lure discrimination',INK);
[~,pc]=ttest(S.REC(:),S.LDI(:));
plot(ax,[1.5 1.5 3.8 3.8],[1.10 1.12 1.12 1.10],'-','Color',INK,'LineWidth',0.8);
text(ax,2.65,1.13,pp(pc),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
ax.XTick=x; ax.XTickLabel={'old | old','old | foil','sim | lure','sim | foil'}; ax.XTickLabelRotation=0;
xlim(ax,[0.4 4.9]); ylim(ax,[0 1.28]); ax.YTick=0:0.25:1;
ylabel(ax,'proportion of responses','FontName','Helvetica','FontSize',11);
finish(ax,INK); exportgraphics(f,fullfile(here,'..','..','figures','mst','fig_behav_mst.pdf'),'ContentType','vector'); close(f);
end

function [f,ax]=canvas(here,SZ,icon,ttl)
INK=[0.149 0.149 0.184]; f=figure('Color','w','Units','inches','Position',[1 1 SZ]);
ip=fullfile(here,'..','..','icons',[icon '.png']);
if ~isempty(icon) && isfile(ip)
    ia=axes(f,'Position',[0.06 0.75 0.21 0.23]); imshow(imread(ip),'Parent',ia); axis(ia,'image','off');
end
annotation(f,'textbox',[0.1 0.885 0.8 0.085],'String',ttl,'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontName','Helvetica','FontSize',13,'Color',INK);
ax=axes(f,'Position',[0.17 0.15 0.80 0.58]); hold(ax,'on'); set(ax,'FontName','Helvetica','FontSize',10.5);
end

function boxdot(ax,x,v,c,INK)
v=v(~isnan(v)); if isempty(v), return; end
w=0.15; q=quantile(v,[.25 .5 .75]); iqv=q(3)-q(1);
lo=max(min(v),q(1)-1.5*iqv); hi=min(max(v),q(3)+1.5*iqv);
plot(ax,[x x],[lo q(1)],'-','Color',INK,'LineWidth',0.8); plot(ax,[x x],[q(3) hi],'-','Color',INK,'LineWidth',0.8);
patch(ax,'XData',x+[-w -w w w],'YData',[q(1) q(3) q(3) q(1)],'FaceColor',c,'FaceAlpha',0.22,'EdgeColor',INK,'LineWidth',1.0);
plot(ax,x+[-w w],[q(2) q(2)],'-','Color',INK,'LineWidth',1.8);
rng(round(x*37)); jit=(rand(numel(v),1)-0.5)*0.22;
scatter(ax,x+jit,v,22,c,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.82);
end

function brklab(ax,x1,x2,yb,yr,p,grp,INK)
plot(ax,[x1 x1 x2 x2],[yb-0.03*yr yb yb yb-0.03*yr],'-','Color',INK,'LineWidth',0.8);
text(ax,(x1+x2)/2,yb+0.03*yr,pp(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
if ~isempty(grp), text(ax,(x1+x2)/2,yb+0.24*yr,grp,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',11.5,'Color',INK); end
end

function idxbr(ax,x1,x2,idx,yb,grp,INK)
[~,p]=ttest(idx,0);
plot(ax,[x1 x1 x2 x2],[yb-0.02 yb yb yb-0.02],'-','Color',INK,'LineWidth',0.8);
text(ax,(x1+x2)/2,yb+0.02,pp(p),'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',9.5,'Color',INK);
text(ax,(x1+x2)/2,yb+0.15,grp,'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Helvetica','FontSize',11.5,'Color',INK);
end

function finish(ax,INK)
box(ax,'off'); ax.TickLength=[0 0]; ax.XColor=INK; ax.YColor=INK; ax.LineWidth=0.8;
end

function v=getv(T,s,roi,c)
v=T.beta(T.subj==s & strcmp(T.roi,roi) & strcmp(T.cell,c)); if isempty(v), v=NaN; else, v=v(1); end
end

function p=ttp(a,b,tail)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m);
if tail==1, [~,p]=ttest(a,b,'Tail','right'); else, [~,p]=ttest(a,b); end
end

function s=pp(p)
if p<0.001, s='p < .001'; else, s=sprintf('p = .%03d',round(p*1000)); end
end

function tk=niceticks(dmin,dmax,n)
raw=(dmax-dmin)/(n-1); mag=10^floor(log10(raw)); cand=[1 2 2.5 5 10]*mag;
step=cand(find(cand>=raw,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1);
if tk(end)<dmax, step=cand(find(cand>=(dmax-dmin)/(n-1)*1.5,1)); lo=floor(dmin/step)*step; tk=lo+step*(0:n-1); end
end
