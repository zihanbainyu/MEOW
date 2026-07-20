clear; clc;
addpath(genpath('/Users/bai/Documents/GitHub/MEOW/toolbox/bayesFactor-master'));
load('../results/spatial_entropy_results.mat');

[~, ia, ib] = intersect(common_subjs_a, common_subjs_b);
ent_comp = (entropy_a_comp_subj(ia) + entropy_b_comp_subj(ib)) / 2;
ent_iso  = (entropy_a_iso_subj(ia)  + entropy_b_iso_subj(ib))  / 2;
n = numel(ent_comp); fprintf('n = %d\n\n', n);

cd = @(x,y) mean(x-y,'omitnan')/std(x-y,'omitnan');
[~,p,~,s] = ttest(ent_comp, ent_iso);
bf10 = bf.ttest(ent_comp, ent_iso);
fprintf('Encoding entropy, aggregated across A and B pair members:\n');
fprintf('  compared: M=%.3f, SD=%.3f\n', mean(ent_comp), std(ent_comp));
fprintf('  isolated: M=%.3f, SD=%.3f\n', mean(ent_iso),  std(ent_iso));
fprintf('  paired t(%d) = %.2f, d = %.2f, p = %.4f\n', s.df, s.tstat, cd(ent_comp,ent_iso), p);
fprintf('  BF10 = %.3f  (BF01 = %.2f)\n', bf10, 1/bf10);
