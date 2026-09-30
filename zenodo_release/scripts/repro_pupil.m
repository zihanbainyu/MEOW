function repro_pupil()
% pupil stats: figure 3 and supplementary figure 3.
% the cluster-based permutation (which selects the 1.0-1.5 s window) is
% precomputed in pupil.mat; here we report it and run the window-mean anovas.
lib = fullfile(fileparts(mfilename('fullpath')),'lib');
addpath(lib, genpath(fullfile(lib,'bayesFactor-master')));
DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','data'); end
load(fullfile(DATA,'pupil.mat'),'pupil');
t = pupil.t; win = t >= pupil.win(1) & t <= pupil.win(2);
P = pupil.pup; C = pupil.corr; I = pupil.incorr;

fprintf('\npupil (figure 3)\n');

% cluster-based permutation, lure trials (precomputed during data generation)
fprintf('\ncluster-based permutation, lure trials\n');
report_clusters('omnibus (condition)', pupil.clusters.omnibus);
report_clusters('compared > isolated', pupil.clusters.comp_iso);
report_clusters('compared > novel',    pupil.clusters.comp_nov);
report_clusters('isolated > novel',    pupil.clusters.iso_nov);

% window means
lc = [wm(C.ab_com,win) wm(C.ab_iso,win) wm(C.ab_nov,win)];
li = [wm(I.ab_com,win) wm(I.ab_iso,win) wm(I.ab_nov,win)];
cd = [wm(P.ab_com,win) wm(P.ab_iso,win) wm(P.ab_nov,win)];
tc = [wm(P.aa_com,win) wm(P.aa_iso,win) wm(P.aa_nov,win)];

fprintf('\nlure trials: 3x2 (encoding condition x accuracy)\n');
rm_3x2('pupil lure 3x2', lc(:,1),li(:,1), lc(:,2),li(:,2), lc(:,3),li(:,3));
fprintf('   condition post-hocs (one-tailed)\n');
oh('compared vs isolated', cd(:,1), cd(:,2));
oh('isolated vs novel   ', cd(:,2), cd(:,3));
oh('compared vs novel   ', cd(:,1), cd(:,3));
fprintf('   correct vs incorrect within condition (one-tailed)\n');
oh('compared', lc(:,1), li(:,1));
oh('isolated', lc(:,2), li(:,2));
oh('novel   ', lc(:,3), li(:,3));
fprintf('   bayes factors\n');
fprintf('     condition                      BF10 = %.2f\n', bfow(cd(:,1),cd(:,2),cd(:,3)));
fprintf('     accuracy                       BF10 = %.2f\n', bftwo(mean(lc,2), mean(li,2)));
fprintf('     compared correct vs incorrect  BF10 = %.2f\n', bft(lc(:,1),li(:,1)));
fprintf('     isolated correct vs incorrect  BF01 = %.2f\n', 1/bft(lc(:,2),li(:,2)));
fprintf('     novel correct vs incorrect     BF10 = %.2f\n', bft(lc(:,3),li(:,3)));

fprintf('\nsame-detection trials: one-way anova (supplementary figure 3)\n');
rm_oneway('pupil same-detect', tc(:,1), tc(:,2), tc(:,3));
fprintf('   same-detection condition       BF01 = %.2f\n', 1/bfow(tc(:,1),tc(:,2),tc(:,3)));
end

function b = bfow(c1,c2,c3)
v=~isnan(c1)&~isnan(c2)&~isnan(c3); n=sum(v);
T=table((1:n)',c1(v),c2(v),c3(v),'VariableNames',{'subj','x1','x2','x3'});
L=stack(T,{'x1','x2','x3'},'NewDataVariableName','y','IndexVariableName','cond');
L.subj=categorical(L.subj); L.cond=categorical(L.cond);
b=bf.anova(L,'y~cond','treatAsRandom',{'subj'},'options',bfopts);
end

function b = bftwo(a1,a2)
v=~isnan(a1)&~isnan(a2); n=sum(v);
T=table((1:n)',a1(v),a2(v),'VariableNames',{'subj','x1','x2'});
L=stack(T,{'x1','x2'},'NewDataVariableName','y','IndexVariableName','acc');
L.subj=categorical(L.subj); L.acc=categorical(L.acc);
b=bf.anova(L,'y~acc','treatAsRandom',{'subj'},'options',bfopts);
end

function b = bft(a, c)
d=a-c; d=d(~isnan(d)); b=bf.ttest(d);
end

function m = wm(X, win), m = mean(X(:,win),2,'omitnan'); end

function oh(lbl, a, b)
[~,p,~,s] = ttest(a, b, 'Tail','right');
fprintf('     %-22s t(%d) = %.2f, d = %.2f, %s\n', lbl, s.df, s.tstat, abs(cohend(a,b)), pstr(p));
end

function report_clusters(lbl, M)
% M rows: [t0 t1 mass p]
if isempty(M), fprintf('   %-22s no clusters above threshold\n', lbl); return; end
for i = 1:size(M,1)
    fprintf('   %-22s %.3f-%.3f s, mass = %.1f, %s\n', lbl, M(i,1), M(i,2), M(i,3), pstr(M(i,4)));
end
end
