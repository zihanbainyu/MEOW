% S_gaze_bin_split.m
% Sebastian's check: does the A1-B1 gaze pattern-separation effect (hit vs miss)
% differ by lure bin (l1 = harder/more similar, l2 = easier)?
% 2x2 (lure bin x subsequent B2 accuracy) on A1-B1 gaze similarity, run SEPARATELY
% for the compared and isolated conditions, WITH a Bayes factor on the interaction:
%   BF01 > 1  =  positive evidence that the gaze effect does NOT depend on bin
%              (i.e. is independent of image similarity) -- the key claim, which a
%              non-significant p-value alone cannot support.
% Requires rm_2x2.m, gg_eps.m, pstr.m and the bayesFactor toolbox on the path.
% Run from Experiment_1/scripts.

clear; clc;
base_dir='..'; res_dir=fullfile(base_dir,'results');
bf_path = fullfile(base_dir,'..','toolbox','bayesFactor-master');
if exist(bf_path,'dir'), addpath(genpath(bf_path)); end

S = load(fullfile(res_dir,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab');
load(fullfile(base_dir,'data','eye_movement_data','group_eye_movement_combined.mat'),'Mw');

% bin lookup: (subj_id, 1-back B trial_id) -> is_l1
oneB = unique(Mw(strcmp(Mw.task,'1_back') & strcmp(Mw.identity,'B'), {'subj_id','trial_id','stim_id'}));
oneB.isl1 = double(contains(string(oneB.stim_id),'_l1'));
kmap = containers.Map(compose('%d_%d', oneB.subj_id, oneB.trial_id), oneB.isl1);

bin_split('compared', S.reinstat_res_ab.ab_compared, kmap);
bin_split('isolated', S.reinstat_res_ab.ab_isolated, kmap);

function bin_split(name, C, kmap)
    C.isl1 = nan(height(C),1);
    for i=1:height(C)
        k = sprintf('%d_%d', C.subj_id(i), C.tr_1b_b(i));
        if isKey(kmap,k), C.isl1(i)=kmap(k); end
    end
    subs = unique(C.subj_id);
    h1=nan(numel(subs),1); m1=h1; h2=h1; m2=h1;   % l1-hit l1-miss l2-hit l2-miss
    for s=1:numel(subs)
        d = C(C.subj_id==subs(s) & ~isnan(C.isl1), :);
        h1(s)=mean(d.reinst_index(d.isl1==1 & d.correct==1),'omitnan');
        m1(s)=mean(d.reinst_index(d.isl1==1 & d.correct==0),'omitnan');
        h2(s)=mean(d.reinst_index(d.isl1==0 & d.correct==1),'omitnan');
        m2(s)=mean(d.reinst_index(d.isl1==0 & d.correct==0),'omitnan');
    end
    fprintf('\n=== A1-B1 gaze (%s): 2x2 (lure bin x B2 accuracy) ===\n', upper(name));
    fprintf('   [factor "condition" = lure bin l1/l2 ; "accuracy" = hit/miss]\n');
    rm_2x2('bin x accuracy', h1, m1, h2, m2);
    d_l1 = m1 - h1; d_l2 = m2 - h2;            % pattern-separation effect (miss - hit)
    v = ~isnan(d_l1) & ~isnan(d_l2);
    [~,p,~,st] = ttest(d_l1(v), d_l2(v));
    fprintf('   gaze effect (miss - hit):  L1 = %.4f, L2 = %.4f\n', mean(d_l1(v)), mean(d_l2(v)));
    fprintf('   L1 vs L2 effect difference: t(%d) = %.2f, %s, n = %d\n', st.df, st.tstat, pstr(p), sum(v));

    % ---- Bayes factor for the bin x accuracy INTERACTION ----
    M = [h1 m1 h2 m2]; keep = all(~isnan(M),2); M = M(keep,:); n = size(M,1);
    try
        subj = (1:n)';
        y   = [M(:,1); M(:,2); M(:,3); M(:,4)];
        bin = categorical([repmat({'l1'},2*n,1);  repmat({'l2'},2*n,1)]);
        acc = categorical([repmat({'hit'},n,1); repmat({'miss'},n,1); repmat({'hit'},n,1); repmat({'miss'},n,1)]);
        T = table(categorical(repmat(subj,4,1)), bin, acc, y, 'VariableNames', {'subj','bin','acc','y'});
        bf_full = bf.anova(T, 'y ~ bin*acc', 'treatAsRandom', {'subj'}, 'verbose', false);
        bf_add  = bf.anova(T, 'y ~ bin+acc', 'treatAsRandom', {'subj'}, 'verbose', false);
        bf_int  = bf_full / bf_add;   % BF10 for the interaction term
        fprintf('   interaction BF10 = %.3f  ->  BF01 = %.2f  (evidence for NO bin-dependence)\n', bf_int, 1/bf_int);
    catch ME
        fprintf('   [BF skipped: %s]\n', ME.message);
    end
end