T = readtable('/Users/bai/Documents/GitHub/MEOW/NeurWEM_CN/prelim_fmri/roi_2x2_persubj.csv');
subs = unique(T.subj);
ROIS = {'mPFC','vmPFC','Visual'};
fprintf('%-10s | interaction F(1,8),p,eta2p | compared p | novel p\n','ROI');
fprintf('%s\n', repmat('-',1,72));
for ri = 1:numel(ROIS)
  r = ROIS{ri};
  cc=nan(numel(subs),1); ci=cc; ic=cc; ii=cc;
  for i=1:numel(subs)
    g=@(cond,a) T.z(T.subj==subs(i)&strcmp(T.roi,r)&strcmp(T.condition,cond)&T.acc==a);
    v=g('compared',1); if ~isempty(v), cc(i)=v; end
    v=g('compared',0); if ~isempty(v), ci(i)=v; end
    v=g('novel',1);    if ~isempty(v), ic(i)=v; end
    v=g('novel',0);    if ~isempty(v), ii(i)=v; end
  end
  t = table(cc, ci, ic, ii, 'VariableNames', {'cc','ci','ic','ii'});
  within = table({'comp';'comp';'nov';'nov'}, {'corr';'incorr';'corr';'incorr'}, 'VariableNames', {'cond','correctness'});
  rm = fitrm(t, 'cc-ii~1', 'WithinDesign', within);
  tbl = ranova(rm, 'WithinModel', 'cond*correctness');
  tn = string(tbl.Properties.RowNames);
  idx = find(contains(tn,'(Intercept):cond:correctness'),1);
  F=tbl.F(idx); pI=tbl.pValueGG(idx); eta=tbl.SumSq(idx)/(tbl.SumSq(idx)+tbl.SumSq(idx+1));
  [~,pc]=ttest(cc,ci); [~,pn]=ttest(ic,ii);
  star = repmat('*',1,(pI<0.05)+(pI<0.01)+(pI<0.001));
  fprintf('%-10s | F=%.2f, p=%.3f, eta=%.2f %-3s | p=%.3f | p=%.3f\n', r, F, pI, eta, star, pc, pn);
end
