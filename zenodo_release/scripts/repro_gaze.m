function repro_gaze()
% gaze-pattern stats: figure 4 (reinstatement), figure 5 (cumulative),
% supplementary entropy, and the self-selection controls.


lib = fullfile(fileparts(mfilename('fullpath')),'lib');
addpath(lib, genpath(fullfile(lib,'bayesFactor-master')));
DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','data'); end
load(fullfile(DATA,'gaze.mat'),'gaze');

% reinstatement 2x2 (figure 4d-f)
fprintf('\ngaze reinstatement (figure 4)\n');
R = gaze.reinst_ab;
sm = cellmeans(R.ab_compared, R.ab_isolated, 'correct');
report('A1-B1 (pattern separation)', sm, 'left');     % correct < incorrect

F = gaze.reinst_full;
sm = cellmeans(drop(F.bb_compared,[609 606 608]), drop(F.bb_isolated,[609 606 608]), 'correct');
report('B2-B1 (pattern completion)', sm, 'right');     % correct > incorrect

bc = attach(F.ba_compared, F.bb_compared);
bi = attach(F.ba_isolated, F.bb_isolated);
sm = cellmeans(drop(bc,[609 606 608 618]), drop(bi,[609 606 608 618]), 'b2_correct');
report('A2-B1 (predictive recall)', sm, 'right');      % correct > incorrect

% cumulative slopes (figure 5)
fprintf('\ncumulative reinstatement (figure 5)\n');
[aa_c, aa_i] = slopes_subj(gaze.cumu_aa, 5);
[ba_c, ba_i] = slopes_subj(gaze.cumu_ba, 4);
vs_zero('A2-A1', aa_c, aa_i);
vs_zero('A2-B1', ba_c, ba_i);
fprintf('\n2x2 anova: pair type x encoding condition\n');
rm_2x2('cumulative slopes', ba_c, ba_i, aa_c, aa_i, {'pair type','encoding cond','interaction'});
[~,pd,~,sd]=ttest(ba_c,ba_i,'Tail','right');
fprintf('   A2-B1 compared vs isolated slope: t(%d) = %.2f, d = %.2f, %s\n', sd.df, sd.tstat, abs(cohend(ba_c,ba_i)), pstr(pd));

% spatial entropy (supplementary)
entropy(gaze.entropy);

% self-selection control
controls(gaze);

% bin-split of the A1-B1 effect
binsplit(gaze);
end

% reinstatement
function report(label, sm, tail)
rm_2x2(label, sm(:,1), sm(:,2), sm(:,3), sm(:,4));
[~,p1,~,s1]=ttest(sm(:,1),sm(:,2),'Tail',tail);
[~,p2,~,s2]=ttest(sm(:,3),sm(:,4),'Tail',tail);
fprintf('   simple effects (one-tailed)\n');
fprintf('     compared correct vs incorrect: t(%d) = %.2f, %s\n', s1.df, s1.tstat, pstr(p1));
fprintf('     isolated correct vs incorrect: t(%d) = %.2f, %s\n', s2.df, s2.tstat, pstr(p2));
end
function sm = cellmeans(Tc, Ti, av)
sids=unique([Tc.subj_id;Ti.subj_id]); n=numel(sids); sm=nan(n,4);
for s=1:n
    sm(s,1)=cm(Tc,sids(s),av,1); sm(s,2)=cm(Tc,sids(s),av,0);
    sm(s,3)=cm(Ti,sids(s),av,1); sm(s,4)=cm(Ti,sids(s),av,0);
end
end
function v = cm(T,sid,av,val)
d=T.reinst_index(T.subj_id==sid & T.(av)==val); if isempty(d), v=NaN; else, v=mean(d,'omitnan'); end
end
function T = drop(T,exc), T(ismember(T.subj_id,exc),:) = []; end
function ba = attach(ba, bb)
ba.b2_correct=nan(height(ba),1);
for i=1:height(ba)
    m=bb(bb.subj_id==ba.subj_id(i) & bb.tr_1b_b==ba.tr_1b_b(i),:);
    if height(m)==1, ba.b2_correct(i)=m.correct; end
end
end

% cumulative
function [sc, si] = slopes_subj(C, col0)
% per-trial slope over cumulative fixations, then per-subject mean; exclude 609
nf=C.n_fix_to_plot; xc=(1:nf)'-mean((1:nf)');
sl_c=trial_slope(C.cumulative_results_comp,col0,nf,xc);
sl_i=trial_slope(C.cumulative_results_iso, col0,nf,xc);
agg=@(s,sl) grpstats(table(s(~isnan(sl)),sl(~isnan(sl)),'VariableNames',{'subj_id','slope'}),'subj_id','mean','DataVars','slope');
gc=agg(C.cumulative_results_comp.subj_id,sl_c); gc.Properties.VariableNames{'mean_slope'}='slope';
gi=agg(C.cumulative_results_iso.subj_id, sl_i); gi.Properties.VariableNames{'mean_slope'}='slope';
common=intersect(gc.subj_id,gi.subj_id); common(common==609)=[];
[~,ic]=ismember(common,gc.subj_id); [~,ii]=ismember(common,gi.subj_id);
sc=gc.slope(ic); si=gi.slope(ii);
end
function sl = trial_slope(tbl,col0,nf,xc)
sl=nan(height(tbl),1);
for i=1:height(tbl), y=tbl{i,col0:col0+nf-1}'; v=~isnan(y); if sum(v)>=2, sl(i)=xc(v)\y(v); end, end
end
function vs_zero(lbl, sc, si)
[~,pc,~,tc]=ttest(sc,0,'Tail','right'); [~,pii,~,ti]=ttest(si,0,'Tail','right');
fprintf('\n%s slope vs zero (one-tailed)\n', lbl);
fprintf('   compared: M = %.4f, SD = %.4f, t(%d) = %.2f, %s\n', mean(sc), std(sc), tc.df, tc.tstat, pstr(pc));
fprintf('   isolated: M = %.4f, SD = %.4f, t(%d) = %.2f, %s\n', mean(si), std(si), ti.df, ti.tstat, pstr(pii));
end

% entropy
function entropy(E)
[~,ia,ib]=intersect(E.common_subjs_a,E.common_subjs_b);
ac=E.entropy_a_comp_subj(ia); ai=E.entropy_a_iso_subj(ia);
bc=E.entropy_b_comp_subj(ib); bi=E.entropy_b_iso_subj(ib);
fprintf('\nspatial entropy (supplementary)\n  n = %d\n', numel(ia));
t=table(ac,ai,bc,bi,'VariableNames',{'A_comp','A_iso','B_comp','B_iso'});
w=table(categorical({'A';'A';'B';'B'}),categorical({'comp';'iso';'comp';'iso'}),'VariableNames',{'Item','Condition'});
tbl=ranova(fitrm(t,'A_comp-B_iso ~ 1','WithinDesign',w),'WithinModel','Item*Condition');
ep=@(k) tbl.SumSq(k)/(tbl.SumSq(k)+tbl.SumSq(k+1));
fprintf('  item:           F(%d,%d) = %.2f, %s, eta_p2 = %.3f\n', tbl.DF(3),tbl.DF(4),tbl.F(3),pstr(tbl.pValue(3)),ep(3));
fprintf('  condition:      F(%d,%d) = %.2f, %s, eta_p2 = %.3f\n', tbl.DF(5),tbl.DF(6),tbl.F(5),pstr(tbl.pValue(5)),ep(5));
fprintf('  item x cond:    F(%d,%d) = %.2f, %s, eta_p2 = %.3f\n', tbl.DF(7),tbl.DF(8),tbl.F(7),pstr(tbl.pValue(7)),ep(7));
[~,p1,~,s1]=ttest(ac,ai); [~,p2,~,s2]=ttest(bc,bi); [~,p3,~,s3]=ttest(ac,bc); [~,p4,~,s4]=ttest(ai,bi);
q=bh_fdr([p1 p2 p3 p4]);
fprintf('  follow-ups (bh-fdr)\n');
fprintf('    A comp vs iso: t(%d) = %.2f, d = %.2f, %s, p_adj = %.3f\n', s1.df, s1.tstat, abs(cohend(ac,ai)), pstr(p1), q(1));
fprintf('    B comp vs iso: t(%d) = %.2f, d = %.2f, %s, p_adj = %.3f\n', s2.df, s2.tstat, abs(cohend(bc,bi)), pstr(p2), q(2));
fprintf('    comp A vs B:   t(%d) = %.2f, d = %.2f, %s, p_adj = %.3f\n', s3.df, s3.tstat, abs(cohend(ac,bc)), pstr(p3), q(3));
fprintf('    iso  A vs B:   t(%d) = %.2f, d = %.2f, %s, p_adj = %.3f\n', s4.df, s4.tstat, abs(cohend(ai,bi)), pstr(p4), q(4));
n=numel(ac);
y=[ac;ai;bc;bi]; itm=categorical([repmat({'A'},2*n,1);repmat({'B'},2*n,1)]);
cnd=categorical([repmat({'comp'},n,1);repmat({'iso'},n,1);repmat({'comp'},n,1);repmat({'iso'},n,1)]);
B=table(categorical(repmat((1:n)',4,1)),itm,cnd,y,'VariableNames',{'subj','Item','Condition','y'});
bf_full=bf.anova(B,'y ~ Item*Condition','treatAsRandom',{'subj'},'options',bfopts);
bf_add =bf.anova(B,'y ~ Item+Condition','treatAsRandom',{'subj'},'options',bfopts);
fprintf('  bayes: full BF10 = %.2f, interaction (full/additive) = %.2f\n', bf_full, bf_full/bf_add);
end

% self-selection control
function controls(gaze)
pk=string(gaze.img_sim.base_id)+"_"+string(gaze.img_sim.bin);
imgM=containers.Map(cellstr(pk),gaze.img_sim.pix_corr);
binM=containers.Map(cellstr(pk),double(string(gaze.img_sim.bin)=="l1"));
B=gaze.control; B.pair_key=string(B.pair_key);
B.img=mapv(imgM,B.pair_key); B.l1=mapv(binM,B.pair_key);
fprintf('\nself-selection control\n  matched %d / %d b2 trials\n', sum(~isnan(B.img)), height(B));
fprintf('\npixel similarity: 2x2 (condition x b2 accuracy)\n'); cell2x2(B,'img');
fprintf('\nlure bin (proportion l1/hard): 2x2 (condition x b2 accuracy)\n'); cell2x2(B,'l1');
end
function cell2x2(B, var)
subs=unique(B.subj_id); cc=nan(numel(subs),1); ci=cc; ic=cc; ii=cc;
for s=1:numel(subs)
    d=B(B.subj_id==subs(s) & ~isnan(B.(var)),:);
    cc(s)=mean(d.(var)(strcmp(d.condition,'compared') & d.correct==1),'omitnan');
    ci(s)=mean(d.(var)(strcmp(d.condition,'compared') & d.correct==0),'omitnan');
    ic(s)=mean(d.(var)(strcmp(d.condition,'isolated') & d.correct==1),'omitnan');
    ii(s)=mean(d.(var)(strcmp(d.condition,'isolated') & d.correct==0),'omitnan');
end
rm_2x2(var, cc, ci, ic, ii);
[~,p1,~,s1]=ttest(cc,ci,'Tail','left'); [~,p2,~,s2]=ttest(ic,ii,'Tail','left');
fprintf('   compared correct vs incorrect: t(%d) = %.2f, %s\n', s1.df, s1.tstat, pstr(p1));
fprintf('   isolated correct vs incorrect: t(%d) = %.2f, %s\n', s2.df, s2.tstat, pstr(p2));
ca=mean([cc ci],2,'omitnan'); ia=mean([ic ii],2,'omitnan'); v=~isnan(ca)&~isnan(ia);
[~,pc,~,sc]=ttest(ca(v),ia(v));
fprintf('   compared vs isolated (overall): t(%d) = %.2f, %s  [comp %.3f, iso %.3f]\n', sc.df, sc.tstat, pstr(pc), mean(ca(v)), mean(ia(v)));
end

% bin-split of A1-B1
function binsplit(gaze)
kmap=containers.Map(compose('%d_%d',gaze.bin_lookup.subj_id,gaze.bin_lookup.trial_id), gaze.bin_lookup.isl1);
fprintf('\nA1-B1 effect split by lure bin\n');
one_split('compared', gaze.reinst_ab.ab_compared, kmap);
one_split('isolated', gaze.reinst_ab.ab_isolated, kmap);
end
function one_split(name, C, kmap)
C.isl1=nan(height(C),1);
for i=1:height(C), k=sprintf('%d_%d',C.subj_id(i),C.tr_1b_b(i)); if isKey(kmap,k), C.isl1(i)=kmap(k); end, end
subs=unique(C.subj_id); h1=nan(numel(subs),1); m1=h1; h2=h1; m2=h1;  % l1-hit l1-miss l2-hit l2-miss
for s=1:numel(subs)
    d=C(C.subj_id==subs(s) & ~isnan(C.isl1),:);
    h1(s)=mean(d.reinst_index(d.isl1==1 & d.correct==1),'omitnan');
    m1(s)=mean(d.reinst_index(d.isl1==1 & d.correct==0),'omitnan');
    h2(s)=mean(d.reinst_index(d.isl1==0 & d.correct==1),'omitnan');
    m2(s)=mean(d.reinst_index(d.isl1==0 & d.correct==0),'omitnan');
end
fprintf('\n%s: 2x2 (lure bin x b2 accuracy)\n', name);
rm_2x2('bin x accuracy', h1, m1, h2, m2, {'lure bin','accuracy','interaction'});
[~,~,~,st]=ttest(m1-h1,m2-h2);
fprintf('   miss-hit: l1 = %.4f, l2 = %.4f (l1 vs l2: t(%d) = %.2f)\n', mean(m1-h1,'omitnan'), mean(m2-h2,'omitnan'), st.df, st.tstat);
bf_int(h1,m1,h2,m2);
end
function bf_int(h1,m1,h2,m2)
M=[h1 m1 h2 m2]; M=M(all(~isnan(M),2),:); n=size(M,1);
try
    y=[M(:,1);M(:,2);M(:,3);M(:,4)];
    bin=categorical([repmat({'l1'},2*n,1);repmat({'l2'},2*n,1)]);
    acc=categorical([repmat({'hit'},n,1);repmat({'miss'},n,1);repmat({'hit'},n,1);repmat({'miss'},n,1)]);
    T=table(categorical(repmat((1:n)',4,1)),bin,acc,y,'VariableNames',{'subj','bin','acc','y'});
    bfi=bf.anova(T,'y ~ bin*acc','treatAsRandom',{'subj'},'options',bfopts) / ...
        bf.anova(T,'y ~ bin+acc','treatAsRandom',{'subj'},'options',bfopts);
    fprintf('   interaction bf01 = %.2f (evidence for no bin-dependence)\n', 1/bfi);
catch
    fprintf('   interaction bf skipped: bayesFactor toolbox not on path\n');
end
end

function v = mapv(M, keys)
v=nan(numel(keys),1);
for i=1:numel(keys), k=char(keys(i)); if isKey(M,k), v(i)=M(k); end, end
end
