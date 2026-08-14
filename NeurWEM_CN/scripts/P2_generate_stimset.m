%==========================================================================
%                  Generate stimulus set
%==========================================================================
% author: Zihan Bai, Michelmann lab at NYU
% email: zihan.bai@nyu.edu
%
% Selection (fMRI design):
%   condition pairs : 120 from lure bin 2 (-> l2) + 120 from bin 3 (-> l3)
%                     = 240 pairs. In A_subject_setup (P3) each level's 120
%                     splits into 60 compared + 60 novel, giving 120 compared +
%                     120 novel overall; both A and B of every pair are shown.
%   repeat pairs    : 80 from lure bin 4 (-> l4), real A-B pairs whose A is
%                     shown twice (A-A) in the 1-back and whose B is the MST lure.
%   foils           : N_FOIL single images (A only) drawn from the remaining
%                     pool -- new/filler items for the 1-back, recognition and
%                     MST. No condition pairs are recycled as foils.
% All selections are balanced by set. Same per-set machinery as the old l1/l2.
clear;
clc;
rng('shuffle');

%% SETUP
%--------------------------------------------------------------------------
% lure bin -> level label -> number of pairs
sel_bins   = [2      3      4    ];
sel_labels = {'l2', 'l3', 'l4'};
sel_counts = [120    120    80   ];   % l2 + l3 = 240 condition pairs, l4 = 80 repeat
% Each level's 120 splits into 60 compared + 60 novel in A_subject_setup (P3),
% giving 120 compared + 120 novel overall. Foils are separate single images.
n_experimental_total = 240;

base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/stim_norm/');                 % normalized images (P1 output; no .txt here)
bin_dir = fullfile(base_dir, 'stimulus/stim_Stark_etal');   % lure-bin labels stay with the raw Stark set
out_dir = fullfile(base_dir, 'stimulus/stim_pool');        % lean pool: only the stimuli actually used
% Foils consumed per subject across all three tasks: 1-back fillers (30 x 4 =
% 120) + recognition-new (80) + MST-new (40) = 240, ALL drawn as single A-images
% from the remaining pool (no condition pairs recycled). P2 supplies exactly 240
% dedicated foils -- no buffer. (The 2-back is built to spend zero padding foils;
% if an unlucky block spends any, setup generation errors on an empty foil pool
% and you re-run it.)
N_FOIL = 240;
mst_set_folders = {'Set 1','Set 2','Set 3','Set 4','Set 5','Set 6'};
lure_bin_files  = {'Set1 bins.txt','Set2 bins.txt','Set3 bins.txt',...
    'Set4 bins.txt','Set5 bins.txt','Set6 bins.txt'};
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% LOAD ALL STIMULUS PAIRS AND LURE-BIN LABELS
%--------------------------------------------------------------------------
fprintf('Loading all available pairs and lure-bin information...\n');
master_pair_list = [];
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    current_dir = fullfile(in_dir, set_name);                     % P1-normalized PNGs
    bin_table = readmatrix(fullfile(bin_dir, lure_bin_files{s})); % bins from the raw Stark set
    if size(bin_table,2) < 2
        error('Lure-bin file for %s must have 2 columns.', set_name);
    end
    for j = 1:192
        path_a = fullfile(current_dir, sprintf('%03da.png', j));
        path_b = fullfile(current_dir, sprintf('%03db.png', j));
        if exist(path_a, 'file') && exist(path_b, 'file')
            lure_bin = bin_table(j, 2);
            master_pair_list = [master_pair_list; {path_a, path_b, lure_bin, set_name}];
        end
    end
end
fprintf('Found %d total valid pairs across all sets.\n', size(master_pair_list,1));
all_bins = cell2mat(master_pair_list(:,3));

%% SELECT PAIRS PER BIN (balanced by set), using the same machinery for each
%--------------------------------------------------------------------------
selected = cell(1, numel(sel_bins));   % selected{k} = pairs for level sel_labels{k}
nSets = numel(mst_set_folders);
for k = 1:numel(sel_bins)
    this_bin  = sel_bins(k);
    n_total   = sel_counts(k);

    % Distribute n_total across the sets as evenly as possible. When it does
    % not divide evenly (e.g. l4 = 80 over 6 sets -> 14,14,13,13,13,13), a
    % random subset of sets takes one extra pair, so per-set counts differ by
    % at most one and the total is EXACT -- the pool holds no buffer, and every
    % subject draws from the same fixed pool.
    base    = floor(n_total / nSets);
    n_extra = n_total - base * nSets;
    per_set_counts = base * ones(1, nSets);
    per_set_counts(randperm(nSets, n_extra)) = base + 1;

    fprintf('\nSelecting %d bin-%d pairs (%s), per-set counts [%s]...\n', ...
        n_total, this_bin, sel_labels{k}, num2str(per_set_counts));
    this_bin_pairs = master_pair_list(all_bins==this_bin, :);
    picked = [];
    for s = 1:nSets
        set_name = mst_set_folders{s};
        set_pairs = this_bin_pairs(strcmp(this_bin_pairs(:,4), set_name), :);
        if size(set_pairs,1) < per_set_counts(s)
            error('%s has insufficient bin-%d pairs (need %d, found %d).', ...
                set_name, this_bin, per_set_counts(s), size(set_pairs,1));
        end
        idx = randperm(size(set_pairs,1), per_set_counts(s));
        picked = [picked; set_pairs(idx,:)];
    end
    selected{k} = picked;
    fprintf('  -> %d %s pairs selected.\n', size(picked,1), sel_labels{k});
end

%% SELECT INSTRUCTION PAIR (a spare experimental-difficulty pair)
%--------------------------------------------------------------------------
% one leftover bin-2 pair, not used for the experiment
used_l2 = [string(selected{1}(:,1)); string(selected{1}(:,2))];
bin2_all = master_pair_list(all_bins==2, :);
is_l2_avail = ~ismember(string(bin2_all(:,1)), used_l2) & ...
              ~ismember(string(bin2_all(:,2)), used_l2);
remaining_l2_pool = bin2_all(is_l2_avail, :);
assert(height(remaining_l2_pool) >= 1, 'No remaining bin-2 pairs for instruction!');
instr_pair_exp = remaining_l2_pool(randi(size(remaining_l2_pool,1)), :);

%% DEFINE POOL OF ALL REMAINING PAIRS (for foils / practice)
%--------------------------------------------------------------------------
used_all = string([]);
for k = 1:numel(selected)
    used_all = [used_all; string(selected{k}(:,1)); string(selected{k}(:,2))];
end
used_all = unique([used_all; string(instr_pair_exp(:,1)); string(instr_pair_exp(:,2))]);
is_pair_available = ~ismember(string(master_pair_list(:,1)), used_all) & ...
                    ~ismember(string(master_pair_list(:,2)), used_all);
available_pool = master_pair_list(is_pair_available, :);
fprintf('\nFound %d remaining pairs available for foils/practice.\n', size(available_pool,1));

%% SET ASIDE INSTRUCTION/PRACTICE FOILS
%--------------------------------------------------------------------------
n_instr_foil = 1;
n_prac_pairs = 18;                          % practice pool: 1-back uses objects 1-9, 2-back uses 10-18 (disjoint, all used)
n_total_setaside = n_instr_foil + n_prac_pairs;
assert(size(available_pool,1) >= n_total_setaside, ...
    'Not enough remaining foils for practice! Need %d, have %d', ...
    n_total_setaside, size(available_pool,1));
idx_setaside = randperm(size(available_pool,1), n_total_setaside);
setaside_pairs = available_pool(idx_setaside, :);
available_pool(idx_setaside, :) = [];
instr_pair_foil = setaside_pairs(1:n_instr_foil, :);
prac_pairs      = setaside_pairs(n_instr_foil+1:end, :);

%% MAIN FOILS = everything left
%--------------------------------------------------------------------------
assert(size(available_pool,1) >= N_FOIL, ...
    'Need %d foils, only %d pairs remain.', N_FOIL, size(available_pool,1));
foil_idx = randperm(size(available_pool,1), N_FOIL);
main_foils = available_pool(foil_idx, :);   % exactly N_FOIL, chosen at random
n_foils = N_FOIL;
fprintf('Assigned %d foils (A-image only).\n', n_foils);

%% COPY EXPERIMENTAL / REPEAT PAIRS (l2, l3, l4)
%--------------------------------------------------------------------------
fprintf('\nCopying selected pairs...\n');
file_idx = 1;
for k = 1:numel(selected)
    lab = sel_labels{k};
    for i = 1:size(selected{k},1)
        copyfile(selected{k}{i,1}, fullfile(out_dir, sprintf('mst_%03d_A_%s.png', file_idx, lab)));
        copyfile(selected{k}{i,2}, fullfile(out_dir, sprintf('mst_%03d_B_%s.png', file_idx, lab)));
        file_idx = file_idx + 1;
    end
end

%% COPY MAIN FOILS
%--------------------------------------------------------------------------
fprintf('Copying %d foils (A-image only)...\n', n_foils);
for i = 1:n_foils
    copyfile(main_foils{i,1}, fullfile(out_dir, sprintf('mst_%03d_A_foil.png', i)));
end

%% COPY INSTRUCTION / PRACTICE STIMULI
%--------------------------------------------------------------------------
fprintf('\nCopying instruction and practice files...\n');
copyfile(instr_pair_exp{1,1},  fullfile(out_dir, 'instr_A.png'));
copyfile(instr_pair_exp{1,2},  fullfile(out_dir, 'instr_B.png'));
copyfile(instr_pair_foil{1,1}, fullfile(out_dir, 'instr_A_foil.png'));
copyfile(instr_pair_foil{1,2}, fullfile(out_dir, 'instr_B_foil.png'));
for i = 1:size(prac_pairs,1)
    copyfile(prac_pairs{i,1}, fullfile(out_dir, sprintf('prac_%02d_A.png', i)));
    copyfile(prac_pairs{i,2}, fullfile(out_dir, sprintf('prac_%02d_B.png', i)));
end

%% SUMMARY
%--------------------------------------------------------------------------
fprintf('\n------------------------------------------------------------\n');
fprintf('Final output folder: %s\n', out_dir);
fprintf('  -> l2 (bin 2): %d | l3 (bin 3): %d  [experimental = %d]\n', ...
    sel_counts(1), sel_counts(2), n_experimental_total);
fprintf('  -> l4 (bin 4): %d  [repeat / MST]\n', sel_counts(3));
fprintf('  -> Foils: %d (A-image only, exact)\n', n_foils);
fprintf('  -> Instruction: 2 pairs | Practice: %d pairs\n', size(prac_pairs,1));
fprintf('------------------------------------------------------------\n');
