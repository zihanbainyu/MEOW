function roi_2x2_lme()
% Condition-level linear mixed-effects model on the 2x2 ROI cell betas (cabs):
%   beta ~ condition * accuracy + (1 | subject)
% Effect-coded (compared/novel = +/-.5, correct/incorrect = +/-.5) so the fixed-effect
% coefficients are the main effects and the condition x accuracy interaction.
% Prints journal-format stats: b [95% CI], t(df), p  (residual df; fitlme default).
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
T=readtable(fullfile(root,'data','mstback','roi_2x2_cells.csv'));
rois={'HC_ant','HC_post','Angular','DLPFC','IPS','Visual','VentralTemporal'};
pres={'sim','same'}; prl={'similar-discrimination','same-detection'};
for pp=1:2
  pre=pres{pp};
  fprintf('\n================ %s :  beta ~ condition*accuracy + (1|subject) ================\n',prl{pp});
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
    lme=fitlme(tbl,'beta ~ cond*acc + (1|subj)'); C=lme.Coefficients; N=numel(unique(rows(:,1)));
    fprintf('%-16s (N=%d)\n',rn,N);
    report(C,'acc',      '  accuracy (correct-incorrect)');
    report(C,'cond',     '  condition (compared-novel)  ');
    report(C,'cond:acc', '  condition x accuracy        ');
    % post-hoc simple effects of accuracy within each condition (coef order: Intercept, cond, acc, cond:acc)
    report_ph(lme,[0 0 1 +0.5],'    post-hoc  accuracy | compared');
    report_ph(lme,[0 0 1 -0.5],'    post-hoc  accuracy | novel   ');
  end
end
end

function report_ph(lme,H,lab)
b=lme.Coefficients.Estimate; V=lme.CoefficientCovariance;
est=H*b; se=sqrt(H*V*H'); t=est/se;
[p2,~,~,df]=coefTest(lme,H);
if est>0, p=p2/2; else, p=1-p2/2; end           % one-tailed: correct > incorrect
tc=tinv(0.975,df);
fprintf('%s:  b = %+.3f [%+.3f, %+.3f], t(%d) = %+.2f, %s (1-tailed)\n',lab,est,est-tc*se,est+tc*se,df,t,pstr(p));
end

function report(C,name,lab)
i=find(strcmp(C.Name,name)); if isempty(i), return; end
fprintf('%s:  b = %+.3f [%+.3f, %+.3f], t(%d) = %+.2f, %s\n',...
        lab, C.Estimate(i), C.Lower(i), C.Upper(i), C.DF(i), C.tStat(i), pstr(C.pValue(i)));
end

function s=pstr(p)
if p<0.001, s='p < .001'; else, s=sprintf('p = %.3f',p); end
end
