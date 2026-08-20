function repro_controls()
% REPRO_CONTROLS  Self-selection control for the A1-B1 pattern-separation effect.
% Tests whether objective stimulus similarity between the paired exemplars
% differs by encoding condition or subsequent B2 accuracy. Two measures:
%   (1) pixel correlation between the grayscale A and B images (graded)
%   (2) MST lure-difficulty bin  (l1 = harder, l2 = easier)
% Both are analysed with the same 2x2 (condition x B2 accuracy) used for gaze.
% A stimulus-driven accuracy effect that is equal across conditions cannot
% explain the compared-specific gaze effect.
%
% Source data:
%   results/image_similarity_AB.csv                 (built by make_image_similarity)
%   data/eye_movement_data/group_eye_movement_combined.mat  (Mw)

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
STIM = getenv('MEOW_STIM'); if isempty(STIM), STIM = fullfile('..','stimulus','stim_final'); end
csv  = fullfile(DATA,'image_similarity_AB.csv');
if ~isfile(csv), make_image_similarity(STIM, csv); end     % self-contained

T  = readtable(csv);
pk = string(T.base_id) + "_" + string(T.bin);
simMap = containers.Map(cellstr(pk), T.pix_corr);
binMap = containers.Map(cellstr(pk), double(string(T.bin)=="l1"));

MW = getenv('MEOW_MW');
if isempty(MW), MW = fullfile('..','data','eye_movement_data','group_eye_movement_combined.mat'); end
load(MW, 'Mw');
isB2 = strcmp(Mw.task,'2_back') & strcmp(Mw.identity,'B') & strcmp(Mw.goal,'A-B') & ...
       (strcmp(Mw.condition,'compared') | strcmp(Mw.condition,'isolated'));
B2 = unique(Mw(isB2, {'subj_id','trial_id','stim_id','condition','correct'}));
key = regexprep(regexprep(string(B2.stim_id),'_B_','_'), '\.png$','');
B2.img = getmap(simMap, key); B2.l1 = getmap(binMap, key);

fprintf('\n================ SELF-SELECTION CONTROL (image similarity) ================\n');
fprintf('matched %d / %d B2 trials\n', sum(~isnan(B2.img)), height(B2));

fprintf('\n-- objective pixel similarity: 2x2 (condition x B2 accuracy) --\n');
run_cell(B2, 'img');
fprintf('\n-- MST lure bin (proportion l1/hard): 2x2 (condition x B2 accuracy) --\n');
run_cell(B2, 'l1');
end

% ================= local helpers =================
function run_cell(B2, var)
subs = unique(B2.subj_id); cc=nan(numel(subs),1); ci=cc; ic=cc; ii=cc;
for s = 1:numel(subs)
    d = B2(B2.subj_id==subs(s) & ~isnan(B2.(var)), :);
    cc(s)=mean(d.(var)(strcmp(d.condition,'compared') & d.correct==1),'omitnan');
    ci(s)=mean(d.(var)(strcmp(d.condition,'compared') & d.correct==0),'omitnan');
    ic(s)=mean(d.(var)(strcmp(d.condition,'isolated') & d.correct==1),'omitnan');
    ii(s)=mean(d.(var)(strcmp(d.condition,'isolated') & d.correct==0),'omitnan');
end
rm_2x2(var, cc, ci, ic, ii);
[~,p1,~,s1]=ttest(cc,ci,'Tail','left'); [~,p2,~,s2]=ttest(ic,ii,'Tail','left');
fprintf('   compared correct vs incorrect: t(%d) = %.2f, %s\n', s1.df, s1.tstat, pstr(p1));
fprintf('   isolated correct vs incorrect: t(%d) = %.2f, %s\n', s2.df, s2.tstat, pstr(p2));
comp_all=mean([cc ci],2,'omitnan'); iso_all=mean([ic ii],2,'omitnan'); v=~isnan(comp_all)&~isnan(iso_all);
[~,pc,~,sc]=ttest(comp_all(v),iso_all(v));
fprintf('   compared vs isolated (overall): t(%d) = %.2f, %s  [M_comp = %.3f, M_iso = %.3f]\n', ...
        sc.df, sc.tstat, pstr(pc), mean(comp_all(v)), mean(iso_all(v)));
end

function v = getmap(M, keys)
v = nan(numel(keys),1);
for i=1:numel(keys), k=char(keys(i)); if isKey(M,k), v(i)=M(k); end, end
end
