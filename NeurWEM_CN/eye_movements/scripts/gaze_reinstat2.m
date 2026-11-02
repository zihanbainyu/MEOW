function gaze_reinstat2()
% Faithful port of paper-1 gaze pipeline. A1-B1 / B2-B1 / A2-B1 matched by base_id
% (item+level: _A_/_B_ -> _BASE_), condition = compared/novel, baseline = other
% same-condition B1 maps. Each row = one matched similar pair, tagged with the B2
% probe's discrimination accuracy. Then: do the indices predict B2 success?
here=fileparts(mfilename('fullpath')); root=fileparts(here);
L=load(fullfile(root,'results','gaze_maps.mat'),'T'); T=L.T;
T.base=regexprep(string(T.stim),'_[AB]_','_BASE_');          % item+level key
subs=unique(T.subj)';
rows=[];                                                     % subj base compared b2corr a1b1 b2b1 a2b1
for s=subs
    Ts=T(T.subj==s,:);
    A1=Ts(Ts.task=="1_back"&Ts.identity=="A",:); B1=Ts(Ts.task=="1_back"&Ts.identity=="B",:);
    A2=Ts(Ts.task=="2_back"&Ts.identity=="A",:); B2=Ts(Ts.task=="2_back"&Ts.identity=="B",:);
    for ib=1:height(B1)                                      % anchor on each encoded B1
        b=B1.base(ib); mB1=B1.map{ib};
        p2=B2(B2.base==b,:); if height(p2)<1, continue; end  % must be probed in 2-back
        cond=string(p2.condition(1));
        rk=regexprep(string(p2.respkey(1)),'\.0',''); rs=regexprep(string(p2.resp(1)),'\.0','');
        if ~ismember(rk,["1" "2"]), b2c=NaN; else, b2c=double(rk==rs); end
        poolB1=B1.map(B1.base~=b & B1.condition==p2.condition(1));   % baseline: other same-cond B1
        pa1=A1(A1.base==b,:); pa2=A2(A2.base==b,:);
        a1b1=NaN; if height(pa1)>=1, a1b1=ri(pa1.map{1},mB1,poolB1); end
        b2b1=ri(p2.map{1},mB1,poolB1);
        a2b1=NaN; if height(pa2)>=1, a2b1=ri(pa2.map{1},mB1,poolB1); end
        rows=[rows; s, double(cond=="compared"), b2c, a1b1, b2b1, a2b1]; %#ok<AGROW>
    end
end
R=array2table(rows,'VariableNames',{'subj','compared','b2correct','a1b1','b2b1','a2b1'});
writetable(R,fullfile(root,'results','gaze_pairs.csv'));
% counts per subject
fprintf('pairs per subject (compared, valid a1b1):\n');
for s=subs, fprintf('  sub-%d: total=%d  compared=%d  a1b1valid=%d\n',s,sum(R.subj==s),sum(R.subj==s&R.compared==1),sum(R.subj==s&~isnan(R.a1b1))); end

idx={'a1b1','b2b1','a2b1'}; nm={'A1-B1 separation','B2-B1 completion','A2-B1 predictive'};
fprintf('\n=== Does gaze similarity predict B2 discrimination success? (compared trials) ===\n');
for j=1:3
    v=idx{j}; c=nan(numel(subs),1); i=nan(numel(subs),1);
    for k=1:numel(subs)
        s=subs(k); Rs=R(R.subj==s & R.compared==1 & ~isnan(R.(v)) & ~isnan(R.b2correct),:);
        c(k)=mean(Rs.(v)(Rs.b2correct==1),'omitnan'); i(k)=mean(Rs.(v)(Rs.b2correct==0),'omitnan');
    end
    m=~isnan(c)&~isnan(i); [~,p]=ttest(c(m),i(m));
    fprintf('%-18s correct=%+.4f incorrect=%+.4f diff=%+.4f  t-test p=%.4f (n=%d)\n',nm{j},nanmean(c),nanmean(i),nanmean(c-i),p,sum(m));
end
% trial-level mixed logistic (compared)
fprintf('\n--- trial-level mixed logistic (compared) P(correct) ~ z(index) + (1|subj) ---\n');
for j=1:3
    v=idx{j}; D=R(R.compared==1 & ~isnan(R.(v)) & ~isnan(R.b2correct),:);
    D.z=(D.(v)-mean(D.(v)))/std(D.(v)); D.subj=categorical(D.subj); D.y=D.b2correct;
    try
        mdl=fitglme(D,'y ~ z + (1|subj)','Distribution','Binomial'); r=find(strcmp(mdl.CoefficientNames,'z'));
        fprintf('%-18s trials=%d  OR=%.2f  p=%.4f\n',nm{j},height(D),exp(mdl.Coefficients.Estimate(r)),mdl.Coefficients.pValue(r));
    catch e, fprintf('%-18s ERR %s\n',nm{j},e.message); end
end
disp('done');
end
function r=ri(mx,my,pool)
if isempty(mx)||isempty(my), r=NaN; return; end
match=corr(mx(:),my(:)); b=nan(numel(pool),1);
for i=1:numel(pool), if ~isempty(pool{i}), b(i)=corr(my(:),pool{i}(:)); end, end
r=match-nanmean(b);
end
