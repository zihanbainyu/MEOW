function make_group_data()
% one-time generator. reads the raw individual-level files once and writes the
% consolidated group data used by everything else. after this, the repro scripts
% never touch the raw files.
%   raw (not deposited): eye table Mw, pupil trial table.
%   out: data/behavior.mat, data/gaze.mat, data/pupil.mat

raw = fullfile(getenv('HOME'),'Downloads','backup_etdata');   % individual-level raw
src = fullfile('..','..','Experiment_1','results');           % existing group intermediates
out = fullfile('..','data');
if ~exist(out,'dir'), mkdir(out); end

% behavior: per-subject summaries (one-back, two-back, recognition)
S = load(fullfile(src,'all_subjs_stats.mat'),'all_subjs');
behavior = S.all_subjs; %#ok<NASGU>
save(fullfile(out,'behavior.mat'),'behavior');
fprintf('behavior.mat: %d subjects\n', numel(S.all_subjs));

% gaze: reinstatement, cumulative slopes, entropy, image similarity, control maps
gaze = struct();
t = load(fullfile(src,'gaze_reinstat_res_ab.mat'),'reinstat_res_ab'); gaze.reinst_ab   = t.reinstat_res_ab;
t = load(fullfile(src,'gaze_reinstat_res_full.mat'),'reinstat_res');  gaze.reinst_full = t.reinstat_res;
gaze.cumu_aa = load(fullfile(src,'cumu_reinstat_aa_cor.mat'));
gaze.cumu_ba = load(fullfile(src,'cumu_reinstat_ba.mat'));
gaze.entropy = load(fullfile(src,'spatial_entropy_results.mat'));
gaze.img_sim = imgsim(fullfile('..','..','Experiment_1','stimulus','stim_final'));

% control + lure-bin maps derived from the raw eye table (used here only)
M = load(fullfile(raw,'group_eye_movement_combined.mat'),'Mw'); Mw = M.Mw; clear M
isB2 = strcmp(Mw.task,'2_back') & strcmp(Mw.identity,'B') & strcmp(Mw.goal,'A-B') & ...
       (strcmp(Mw.condition,'compared') | strcmp(Mw.condition,'isolated'));
b2  = unique(Mw(isB2, {'subj_id','stim_id','condition','correct'}));
key = regexprep(regexprep(string(b2.stim_id),'_B_','_'),'\.png$','');   % mst_###_B_l# -> mst_###_l#
gaze.control = table(b2.subj_id, key, string(b2.condition), double(b2.correct), ...
                     'VariableNames',{'subj_id','pair_key','condition','correct'});
ob = unique(Mw(strcmp(Mw.task,'1_back') & strcmp(Mw.identity,'B'), {'subj_id','trial_id','stim_id'}));
gaze.bin_lookup = table(ob.subj_id, ob.trial_id, double(contains(string(ob.stim_id),'_l1')), ...
                        'VariableNames',{'subj_id','trial_id','isl1'});
save(fullfile(out,'gaze.mat'),'gaze');
fprintf('gaze.mat: %d control rows, %d bin rows\n', height(gaze.control), height(gaze.bin_lookup));
clear Mw

% pupil: per-subject baseline-corrected time-courses by condition and correctness
L = load(fullfile(raw,'all_trials_pupil.mat'),'all_preprocessed');
pupil = extract_pupil(L.all_preprocessed);
pupil.clusters = pupil_clusters(pupil);   % cluster-based permutation, computed once here
save(fullfile(out,'pupil.mat'),'pupil');
fprintf('pupil.mat: %d subjects, %d samples\n', numel(pupil.subj), numel(pupil.t));
end

function P = extract_pupil(ap)
% per-subject condition time-courses (a-b lure and a-a target), plus a-b lure
% split by correctness. window means are computed downstream.
v = ap(ap.preprocess_success, :);
maxs = max(cellfun(@length, v.pupil_preprocessed));
sr   = v.sample_rate(1);
subj = unique(v.subj_id); n = numel(subj);
conds = {'compared','isolated','novel'}; goals = {'A-B','A-A'}; gt = {'ab','aa'};
pup = struct(); pc = struct(); pin = struct();
for g=1:2, for k=1:3, pup.(sprintf('%s_%s',gt{g},conds{k}(1:3))) = nan(n,maxs); end, end
for k=1:3, pc.(sprintf('ab_%s',conds{k}(1:3)))=nan(n,maxs); pin.(sprintf('ab_%s',conds{k}(1:3)))=nan(n,maxs); end
for s=1:n
    sv = v(v.subj_id==subj(s),:);
    for g=1:2, for k=1:3
        pup.(sprintf('%s_%s',gt{g},conds{k}(1:3)))(s,:) = tracemean(sv(strcmp(sv.condition,conds{k}) & strcmp(sv.goal,goals{g}),:), maxs);
    end, end
    for k=1:3
        lure = sv(strcmp(sv.condition,conds{k}) & strcmp(sv.goal,'A-B') & strcmp(cellstr(sv.corr_resp),'k'),:);
        pc.(sprintf('ab_%s',conds{k}(1:3)))(s,:)  = tracemean(lure(lure.correct==1,:), maxs);
        pin.(sprintf('ab_%s',conds{k}(1:3)))(s,:) = tracemean(lure(lure.correct==0,:), maxs);
    end
end
P = struct('pup',pup,'corr',pc,'incorr',pin,'t',(0:maxs-1)/sr,'subj',subj,'win',[1.0 1.5]);
end

function a = tracemean(tr, maxs)
% baseline-correct each trial (subtract mean baseline) and average across trials
if height(tr)==0, a = nan(1,maxs); return; end
m = nan(height(tr),maxs);
for i=1:height(tr)
    p = tr.pupil_preprocessed{i} - mean(tr.baseline_pupil_preprocessed{i},'omitnan');
    m(i,1:numel(p)) = p(:)';
end
a = mean(m,1,'omitnan');
end

function T = imgsim(stim_dir)
% objective a-b similarity: grayscale pixel correlation per pair (bins l1/l2)
A = dir(fullfile(stim_dir,'mst_*_A_l*.png'));
rows = cell(0,3);
for i=1:numel(A)
    b = fullfile(stim_dir, strrep(A(i).name,'_A_l','_B_l'));
    if ~isfile(b), continue; end
    tok = regexp(A(i).name,'(mst_\d+)_A_(l\d)','tokens','once');
    if isempty(tok), continue; end
    va = double(togray(imread(fullfile(stim_dir,A(i).name))));
    vb = double(togray(imread(b)));
    if ~isequal(size(va),size(vb)), vb = double(imresize(togray(imread(b)),size(va))); end
    rows(end+1,:) = {tok{1}, tok{2}, corr(va(:),vb(:))}; %#ok<AGROW>
end
T = cell2table(rows,'VariableNames',{'base_id','bin','pix_corr'});
end

function g = togray(im)
if ndims(im)==3, g = rgb2gray(im); else, g = im; end
end
