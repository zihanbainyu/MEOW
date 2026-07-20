clear; clc;
addpath(genpath('/Users/bai/Documents/GitHub/MEOW/toolbox/bayesFactor-master'));
load('../results/spatial_entropy_results.mat');

% align subjects across A and B (use intersection)
[common_ab, ia, ib] = intersect(common_subjs_a, common_subjs_b);
ea_c = entropy_a_comp_subj(ia); ea_i = entropy_a_iso_subj(ia);
eb_c = entropy_b_comp_subj(ib); eb_i = entropy_b_iso_subj(ib);
n = numel(common_ab); fprintf('n = %d common subjects\n\n', n);

%% 2x2 RM-ANOVA: ItemPosition (A,B) x Condition (compared,isolated)
t = table(ea_c, ea_i, eb_c, eb_i, 'VariableNames', {'A_comp','A_iso','B_comp','B_iso'});
within = table(categorical({'A';'A';'B';'B'}), categorical({'comp';'iso';'comp';'iso'}), ...
    'VariableNames', {'Item','Condition'});
rm = fitrm(t, 'A_comp-B_iso ~ 1', 'WithinDesign', within);
tbl = ranova(rm, 'WithinModel', 'Item*Condition');
disp(tbl);

%% partial eta-squared
get_eta2p = @(SSeff, SSerr) SSeff / (SSeff + SSerr);
SS = tbl.SumSq;
eta2p_item = get_eta2p(SS(3), SS(4));
eta2p_cond = get_eta2p(SS(5), SS(6));
eta2p_int  = get_eta2p(SS(7), SS(8));
fprintf('Item:           F(%d,%d)=%.2f, p=%.4f, eta2_p=%.3f\n', tbl.DF(3), tbl.DF(4), tbl.F(3), tbl.pValue(3), eta2p_item);
fprintf('Condition:      F(%d,%d)=%.2f, p=%.4f, eta2_p=%.3f\n', tbl.DF(5), tbl.DF(6), tbl.F(5), tbl.pValue(5), eta2p_cond);
fprintf('Item*Condition: F(%d,%d)=%.2f, p=%.4f, eta2_p=%.3f\n', tbl.DF(7), tbl.DF(8), tbl.F(7), tbl.pValue(7), eta2p_int);

%% Bayes factor (full model vs null)
subj = (1:n)';
y = [ea_c; ea_i; eb_c; eb_i];
sj = repmat(subj, 4, 1);
itm = categorical([repmat({'A'},2*n,1); repmat({'B'},2*n,1)]);
cnd = categorical([repmat({'comp'},n,1); repmat({'iso'},n,1); repmat({'comp'},n,1); repmat({'iso'},n,1)]);
bf_tbl = table(categorical(sj), itm, cnd, y, 'VariableNames', {'subj','Item','Condition','y'});

bf_full = bf.anova(bf_tbl, 'y ~ Item*Condition', 'treatAsRandom', {'subj'}, 'verbose', false);
bf_item = bf.anova(bf_tbl, 'y ~ Item',            'treatAsRandom', {'subj'}, 'verbose', false);
bf_cond = bf.anova(bf_tbl, 'y ~ Condition',       'treatAsRandom', {'subj'}, 'verbose', false);
bf_add  = bf.anova(bf_tbl, 'y ~ Item+Condition',  'treatAsRandom', {'subj'}, 'verbose', false);

fprintf('\nBayes factors (vs null = subj-only):\n');
fprintf('  Item only:      BF10 = %.2f  (BF01 = %.2f)\n', bf_item, 1/bf_item);
fprintf('  Condition only: BF10 = %.2f  (BF01 = %.2f)\n', bf_cond, 1/bf_cond);
fprintf('  Item+Cond:      BF10 = %.2f  (BF01 = %.2f)\n', bf_add,  1/bf_add);
fprintf('  Full (Item*Cond): BF10 = %.2f  (BF01 = %.2f)\n', bf_full, 1/bf_full);
fprintf('  Interaction (full/additive): BF = %.2f\n', bf_full / bf_add);

%% follow-up paired t-tests (BH-FDR over 4)
cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
[~,p1,~,s1] = ttest(ea_c, ea_i);
[~,p2,~,s2] = ttest(eb_c, eb_i);
[~,p3,~,s3] = ttest(ea_c, eb_c);
[~,p4,~,s4] = ttest(ea_i, eb_i);
q = bh_fdr([p1 p2 p3 p4]);
fprintf('\nFollow-up paired t-tests (BH-FDR over 4):\n');
fprintf('  A: comp vs iso:  t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s1.df, s1.tstat, cd(ea_c,ea_i), p1, q(1));
fprintf('  B: comp vs iso:  t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s2.df, s2.tstat, cd(eb_c,eb_i), p2, q(2));
fprintf('  comp: A vs B:    t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s3.df, s3.tstat, cd(ea_c,eb_c), p3, q(3));
fprintf('  iso:  A vs B:    t(%d)=%.2f, d=%.2f, p=%.4f, q=%.4f\n', s4.df, s4.tstat, cd(ea_i,eb_i), p4, q(4));

fprintf('\nMeans (SD):\n');
fprintf('  A compared: %.3f (%.3f)\n', mean(ea_c), std(ea_c));
fprintf('  A isolated: %.3f (%.3f)\n', mean(ea_i), std(ea_i));
fprintf('  B compared: %.3f (%.3f)\n', mean(eb_c), std(eb_c));
fprintf('  B isolated: %.3f (%.3f)\n', mean(eb_i), std(eb_i));

function q = bh_fdr(p)
    n=numel(p); [ps,idx]=sort(p); qs=ps.*n./(1:n);
    for i=n-1:-1:1, qs(i)=min(qs(i),qs(i+1)); end
    q=zeros(size(p)); q(idx)=min(qs,1);
end
