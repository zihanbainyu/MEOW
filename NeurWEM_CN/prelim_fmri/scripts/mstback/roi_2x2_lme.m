function roi_2x2_lme()
% Condition-level LME on the 2x2 ROI cell betas (cabs): beta ~ condition*accuracy + (1|subject).
% Effect-coded (compared/novel = +/-.5, correct/incorrect = +/-.5) so anova(lme) gives the
% accuracy main effect and condition x accuracy interaction, the mixed-model analogue of the rmANOVA.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
T=readtable(fullfile(root,'data','mstback','roi_2x2_cells.csv'));
rois={'HC_ant','HC_post','Angular','DLPFC','IPS','Visual','VentralTemporal'};
pres={'sim','same'}; prl={'similar-discrimination','same-detection'};
for pp=1:2
  pre=pres{pp};
  fprintf('\n================ %s (2x2 LME) ================\n',prl{pp});
  fprintf('%-16s | accuracy F,p | cond F,p | cond x acc F,p\n','ROI');
  for r=1:numel(rois)
    rn=rois{r}; s=unique(T.subj); rows=[];
    for i=1:numel(s)
      g=@(c) T.beta(T.subj==s(i)&strcmp(T.roi,rn)&strcmp(T.cell,c));
      ch=g([pre '_cmp_hit']); ci=g([pre '_cmp_miss']); nh=g([pre '_nov_hit']); ni=g([pre '_nov_miss']);
      if isempty(ch)||isempty(ci)||isempty(nh)||isempty(ni), continue; end
      rows=[rows; s(i) +.5 +.5 ch(1); s(i) +.5 -.5 ci(1); s(i) -.5 +.5 nh(1); s(i) -.5 -.5 ni(1)];
    end
    if isempty(rows), continue; end
    tbl=table(categorical(rows(:,1)),rows(:,2),rows(:,3),rows(:,4),'VariableNames',{'subj','cond','acc','beta'});
    lme=fitlme(tbl,'beta ~ cond*acc + (1|subj)');
    a=anova(lme); nm=a.Term; F=a.FStat; P=a.pValue;
    gp=@(t) deal(F(strcmp(nm,t)),P(strcmp(nm,t)));
    [Fa,Pa]=gp('acc'); [Fc,Pc]=gp('cond'); [Fi,Pi]=gp('cond:acc');
    fprintf('%-16s | F=%5.2f p=%.3f | F=%5.2f p=%.3f | F=%5.2f p=%.3f\n',rn,Fa,Pa,Fc,Pc,Fi,Pi);
  end
end
end
