function behav_lme()
% Journal-format LME/paired stats for the behavioural comparisons, matching the neural
% reporting (b [95% CI], t(df), p). Paired comparisons: y ~ grp + (1|subject).
% One-sample (vs chance) kept as one-sample t (one value per subject). One-tailed where directional.
here=fileparts(mfilename('fullpath')); root=fullfile(here,'..','..');
addpath(fullfile(root,'..','scripts'));
subs=101:109; di=nan(numel(subs),3); ds=nan(numel(subs),3);   % cols: 1-back, compared, novel
for i=1:numel(subs)
  Dc=nm_load_long(subs(i)); if isempty(Dc), continue; end; m=nm_metrics(Dc);
  di(i,:)=[m.di_one m.di(1) m.di(2)]; ds(i,:)=[m.dSame_one m.dSame(1) m.dSame(2)];
end
fprintf('\n================ MST-back behaviour (LME: y ~ condition + (1|subject)) ================\n');
fprintf('-- similar-discrimination index --\n');
reportpair(di(:,1),di(:,3),1,'  1-back > novel');
reportpair(di(:,2),di(:,3),1,'  compared > novel');
fprintf('-- same-detection d'' --\n');
reportpair(ds(:,1),ds(:,3),1,'  1-back > novel');
reportpair(ds(:,2),ds(:,3),1,'  compared > novel');

S=load(fullfile(root,'data','mst','mst_idx.mat'));
fprintf('\n================ MST behaviour ================\n');
report1(S.REC(:),'  recognition index > chance');
report1(S.LDI(:),'  lure-discrimination index > chance');
reportpair(S.REC(:),S.LDI(:),0,'  recognition vs lure-discrimination');
end

function reportpair(a,b,tail,lab)
m=~isnan(a)&~isnan(b); a=a(m); b=b(m); n=numel(a);
tb=table(categorical([(1:n)';(1:n)']),[ones(n,1);zeros(n,1)],[a;b],'VariableNames',{'subj','grp','y'});
C=fitlme(tb,'y ~ grp + (1|subj)').Coefficients; i=strcmp(C.Name,'grp');
est=C.Estimate(i); df=C.DF(i); t=C.tStat(i); p2=C.pValue(i); lo=C.Lower(i); hi=C.Upper(i);
if tail==1, if est>0, p=p2/2; else, p=1-p2/2; end, ts=' (1-tailed)'; else, p=p2; ts=''; end
fprintf('%-38s: b = %+.3f [%+.3f, %+.3f], t(%d) = %+.2f, %s%s\n',lab,est,lo,hi,df,t,pstr(p),ts);
end

function report1(x,lab)   % one-sample vs chance (0), one-tailed x>0
x=x(~isnan(x)); n=numel(x); mu=mean(x); se=std(x)/sqrt(n); t=mu/se; df=n-1; p=tcdf(-t,df);
ci=mu+[-1 1]*tinv(0.975,df)*se;
fprintf('%-38s: M = %+.3f [%+.3f, %+.3f], t(%d) = %+.2f, %s (1-tailed, one-sample)\n',lab,mu,ci(1),ci(2),df,t,pstr(p));
end

function s=pstr(p)
if p<0.001, s='p < .001'; else, s=sprintf('p = %.3f',p); end
end
