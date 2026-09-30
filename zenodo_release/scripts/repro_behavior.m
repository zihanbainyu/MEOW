function repro_behavior()
% behavioral stats: figure 2 and supplementary figure 1.

lib = fullfile(fileparts(mfilename('fullpath')),'lib');
addpath(lib, genpath(fullfile(lib,'bayesFactor-master')));
DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','data'); end
load(fullfile(DATA,'behavior.mat'),'behavior');
gv = @(f1,f2) arrayfun(@(x) x.stats.(f1).(f2), behavior)';

ldi = {gv('two','ldi_comp'),    gv('two','ldi_iso'),    gv('two','ldi_nov')};
dpr = {gv('two','dprime_comp'), gv('two','dprime_iso'), gv('two','dprime_nov')};
rtL = {gv('two','rt_AB_comp'),  gv('two','rt_AB_iso'),  gv('two','rt_AB_nov')};
rtT = {gv('two','rt_AA_comp'),  gv('two','rt_AA_iso'),  gv('two','rt_AA_nov')};

fprintf('\nbehavioral (figure 2)\n');
fprintf('\ntwo-back one-way rm-anovas (condition: compared/isolated/novel)\n');
rm_oneway('discrimination idx', ldi{:}); posthoc('  discrimination idx', ldi{:}, +1);
rm_oneway('same-item d''',       dpr{:}); % null control
rm_oneway('RT lure',            rtL{:}); posthoc('  RT lure', rtL{:}, -1);
rm_oneway('RT target',          rtT{:}); % null control
fprintf('  bayes factors\n');
bf_report('discrimination idx', ldi{:});
bf_report('same-item d''',       dpr{:});
bf_report('RT lure',            rtL{:});
bf_report('RT target',          rtT{:});

fprintf('\nwm-em correlations (one-tailed)\n');
wm  = (dpr{1}+dpr{2}+dpr{3})/3;
emb = (ldi{1}+ldi{2})/2 - ldi{3};
em  = (gv('rec','d_comp') + gv('rec','d_iso'))/2;
corr1('em benefit x wm-d''', emb, wm);
corr1('em benefit x em-d''', emb, em);
corr1('wm-d'' x em-d''',      wm,  em);

fprintf('\none-back (supplementary figure 1)\n');
as=gv('one','acc_same'); asi=gv('one','acc_sim'); an=gv('one','acc_new');
rs=gv('one','rt_same');  rsi=gv('one','rt_sim');
[~,p1,~,s1]=ttest(as,asi); [~,p2,~,s2]=ttest(asi,an); [~,p3,~,s3]=ttest(as,an); [~,p4,~,s4]=ttest(rs,rsi);
q=bh_fdr([p1 p2 p3 p4]);
prow('acc same vs similar', s1, cohend(as,asi), p1, q(1));
prow('acc similar vs new ', s2, cohend(asi,an), p2, q(2));
prow('acc same vs new    ', s3, cohend(as,an),  p3, q(3));
prow('RT  same vs similar', s4, cohend(rs,rsi), p4, q(4));

fprintf('\nrecognition d''\n');
rc=gv('rec','d_comp'); ri=gv('rec','d_iso');
[~,pp1,~,ss1]=ttest(rc); [~,pp2,~,ss2]=ttest(ri); [~,pp3,~,ss3]=ttest(rc,ri);
qq=bh_fdr([pp1 pp2 pp3]);
prow('compared > 0   ', ss1, cohend(rc), pp1, qq(1));
prow('isolated > 0   ', ss2, cohend(ri), pp2, qq(2));
prow('compared vs iso', ss3, cohend(rc,ri), pp3, qq(3));
end

function posthoc(lbl, xc, xi, xn, sgn)
tail='right'; if sgn<0, tail='left'; end
[~,p1,~,s1]=ttest(xc,xi,'Tail',tail); [~,p2,~,s2]=ttest(xi,xn,'Tail',tail); [~,p3,~,s3]=ttest(xc,xn,'Tail',tail);
q=bh_fdr([p1 p2 p3]);
fprintf('%s post-hocs (one-tailed, bh-fdr)\n', lbl);
prow('  compared vs isolated', s1, cohend(xc,xi), p1, q(1));
prow('  isolated vs novel   ', s2, cohend(xi,xn), p2, q(2));
prow('  compared vs novel   ', s3, cohend(xc,xn), p3, q(3));
end

function corr1(lbl, x, y)
v=~isnan(x)&~isnan(y); [r,p]=corr(x(v),y(v),'Tail','right');
fprintf('  %-22s r(%d) = %.2f, p = %.3f, n = %d\n', lbl, sum(v)-2, r, p, sum(v));
end

function prow(lbl, s, d, p, q)
fprintf('  %-22s t(%d) = %.2f, d = %.2f, %s, p_adj = %.3f\n', lbl, s.df, s.tstat, abs(d), pstr(p), q);
end

function bf_report(lbl, a, b, c)
% one-way bayes factor (subject random); deterministic quadrature
n=numel(a); T=table((1:n)',a,b,c,'VariableNames',{'subj','x1','x2','x3'});
L=stack(T,{'x1','x2','x3'},'NewDataVariableName','y','IndexVariableName','cond');
L.subj=categorical(L.subj); L.cond=categorical(L.cond);
b10=bf.anova(L,'y~cond','treatAsRandom',{'subj'},'options',bfopts);
if b10>=1, fprintf('    %-20s BF10 = %.2f\n', lbl, b10);
else,      fprintf('    %-20s BF01 = %.2f\n', lbl, 1/b10); end
end
