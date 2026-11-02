function gaze_accuracy()
% Do the gaze reinstatement indices predict B2 lure-discrimination success?
% Per subject, each index (A1-B1 separation, B2-B1 completion, A2-B1 predictive)
% split by B2 accuracy (correct/incorrect), overall and within condition; paired t + 2x2.
here=fileparts(mfilename('fullpath')); root=fileparts(here);
R=readtable(fullfile(root,'results','gaze_reinstat_trials.csv'));
subs=unique(R.subj)'; idx={'a1b1','b2b1','a2b1'}; nm={'A1-B1 separation','B2-B1 completion','A2-B1 predictive'};
fprintf('N=%d subjects\n\n',numel(subs));
for j=1:3
    v=idx{j};
    [c,i]=deal(nan(numel(subs),1)); [cc,ci,nc,ni]=deal(nan(numel(subs),1));
    for k=1:numel(subs)
        s=subs(k); Rs=R(R.subj==s & ~isnan(R.(v)) & ~isnan(R.b2correct),:);
        c(k)=mean(Rs.(v)(Rs.b2correct==1),'omitnan'); i(k)=mean(Rs.(v)(Rs.b2correct==0),'omitnan');
        cc(k)=mean(Rs.(v)(Rs.b2correct==1 & Rs.compared==1),'omitnan');
        ci(k)=mean(Rs.(v)(Rs.b2correct==0 & Rs.compared==1),'omitnan');
        nc(k)=mean(Rs.(v)(Rs.b2correct==1 & Rs.compared==0),'omitnan');
        ni(k)=mean(Rs.(v)(Rs.b2correct==0 & Rs.compared==0),'omitnan');
    end
    [~,p]=ttest(c,i);
    fprintf('=== %s ===\n',nm{j});
    fprintf('  overall   correct=%+.4f  incorrect=%+.4f  diff=%+.4f  p=%.4f  (n=%d)\n',nanmean(c),nanmean(i),nanmean(c-i),p,sum(~isnan(c-i)));
    pr=pairp(cc,ci); fprintf('  compared  correct=%+.4f  incorrect=%+.4f  diff=%+.4f  p=%.4f\n',nanmean(cc),nanmean(ci),nanmean(cc-ci),pr);
    pr=pairp(nc,ni); fprintf('  novel     correct=%+.4f  incorrect=%+.4f  diff=%+.4f  p=%.4f\n',nanmean(nc),nanmean(ni),nanmean(nc-ni),pr);
    inter=(cc-ci)-(nc-ni); pi=pairp(inter,zeros(size(inter)));
    fprintf('  condition x accuracy interaction  p=%.4f\n\n',pi);
end
end
function p=pairp(a,b)
m=~isnan(a)&~isnan(b); if sum(m)<3, p=NaN; return; end
[~,p]=ttest(a(m),b(m));
end
