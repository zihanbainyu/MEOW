%==========================================================================
%                  Generate stimulus set
%==========================================================================
% author: Zihan Bai, Michelmann lab at NYU
% email: zihan.bai@nyu.edu
clear;
clc;
rng('shuffle');

%% SETUP
%--------------------------------------------------------------------------
n_experimental_pairs_total = 360; % 180 from Bin 1, 180 from Bin 2
n_bin1_pairs = 180;
n_bin2_pairs = 180;
base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/stim_norm/');
out_dir = fullfile(base_dir, 'stimulus/stim_final');
mst_set_folders = {'Set 1','Set 2','Set 3','Set 4','Set 5','Set 6'};
lure_bin_files  = {'Set1_bins.txt','Set2_bins.txt','Set3_bins.txt',...
    'Set4_bins.txt','Set5_bins.txt','Set6_bins.txt'};
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
%% LOAD ALL STIMULUS PAIRS AND LURE-BIN LABELS
%--------------------------------------------------------------------------
fprintf('Loading all available pairs and lure-bin information...\n');
master_pair_list = [];
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    current_dir = fullfile(in_dir, set_name);
    bin_table = readmatrix(fullfile(current_dir, lure_bin_files{s}));
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
%% SELECT EXPERIMENTAL PAIRS (BIN=1 and BIN=2)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d total experimental pairs (%d from Bin 1, %d from Bin 2)...\n', ...
    n_experimental_pairs_total, n_bin1_pairs, n_bin2_pairs);
% balance by set within each bin
n_bin1_per_set = floor(n_bin1_pairs / numel(mst_set_folders)); % 180/6 = 30
n_bin2_per_set = floor(n_bin2_pairs / numel(mst_set_folders)); % 180/6 = 30
experimental_pairs_l1 = [];
experimental_pairs_l2 = [];
% get Bin 1 Pairs, balanced by set
fprintf('Selecting %d Bin 1 pairs (%d per set)...\n', n_bin1_pairs, n_bin1_per_set);
all_bin1_pairs = master_pair_list(all_bins==1, :);
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(all_bin1_pairs(:,4), set_name);
    set_bin1_pairs = all_bin1_pairs(in_set, :);
    if size(set_bin1_pairs,1) < n_bin1_per_set
        error('%s has insufficient bin-1 pairs (need %d, found %d).', ...
            set_name, n_bin1_per_set, size(set_bin1_pairs,1));
    end
    idx = randperm(size(set_bin1_pairs,1), n_bin1_per_set);
    experimental_pairs_l1 = [experimental_pairs_l1; set_bin1_pairs(idx,:)];
end
% get Bin 2 Pairs, balanced by set
fprintf('Selecting %d Bin 2 pairs (%d per set)...\n', n_bin2_pairs, n_bin2_per_set);
all_bin2_pairs = master_pair_list(all_bins==2, :);
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(all_bin2_pairs(:,4), set_name);
    set_bin2_pairs = all_bin2_pairs(in_set, :);
    if size(set_bin2_pairs,1) < n_bin2_per_set
        error('%s has insufficient bin-2 pairs (need %d, found %d).', ...
            set_name, n_bin2_per_set, size(set_bin2_pairs,1));
    end
    idx = randperm(size(set_bin2_pairs,1), n_bin2_per_set);
    experimental_pairs_l2 = [experimental_pairs_l2; set_bin2_pairs(idx,:)];
end
fprintf('Selected %d Bin 1 pairs and %d Bin 2 pairs.\n', ...
    size(experimental_pairs_l1,1), size(experimental_pairs_l2,1));

%% SELECT INSTRUCTION PAIR
%--------------------------------------------------------------------------
fprintf('\nSetting aside 1 L1 pair for instruction...\n');
% Find all L1 images that were NOT selected for the experiment
used_l1_imgs = [string(experimental_pairs_l1(:,1)); string(experimental_pairs_l1(:,2))];
is_l1_avail = ~ismember(string(all_bin1_pairs(:,1)), used_l1_imgs) & ...
              ~ismember(string(all_bin1_pairs(:,2)), used_l1_imgs);
remaining_l1_pool = all_bin1_pairs(is_l1_avail, :);
assert(height(remaining_l1_pool) >= 1, 'No remaining L1 pairs to use for instruction!');
% Select one
idx_instr_l1 = randi(height(remaining_l1_pool));
instr_pair_l1 = remaining_l1_pool(idx_instr_l1, :);
fprintf('... 1 L1 pair set aside.\n');

%% DEFINE POOL OF ALL REMAINING PAIRS
%--------------------------------------------------------------------------
fprintf('\nDefining pool of all remaining available pairs...\n');
% --- Create a complete list of ALL "used" images so far ---
used_exp_imgs = [string(experimental_pairs_l1(:,1)); string(experimental_pairs_l1(:,2)); ...
                 string(experimental_pairs_l2(:,1)); string(experimental_pairs_l2(:,2))];
instr_l1_imgs = [string(instr_pair_l1(:,1)); string(instr_pair_l1(:,2))];
used_all_images = unique([used_exp_imgs; instr_l1_imgs]);

% A pair is available ONLY IF neither its 'a' nor 'b' image is in the used list
is_pair_available = ~ismember(string(master_pair_list(:,1)), used_all_images) & ...
                    ~ismember(string(master_pair_list(:,2)), used_all_images);

available_pool = master_pair_list(is_pair_available, :);
fprintf('Found %d total remaining pairs available for foils/practice.\n', height(available_pool));

%% SET ASIDE INSTRUCTION/PRACTICE FOILS
%--------------------------------------------------------------------------
n_instr_foil = 1;
n_prac_pairs = 25; % 15 for lure + 10 for repeat = 25 total pairs
n_total_setaside = n_instr_foil + n_prac_pairs; % 26 pairs

fprintf('Setting aside %d foil pairs from available pool for instruction/practice...\n', n_total_setaside);
assert(height(available_pool) >= n_total_setaside, ...
    'Not enough remaining foils for practice! Need %d, have %d', ...
    n_total_setaside, height(available_pool));

% Select and REMOVE these pairs from the available_pool
idx_setaside = randperm(height(available_pool), n_total_setaside);
setaside_pairs = available_pool(idx_setaside, :);
available_pool(idx_setaside, :) = []; % This is the critical step
fprintf('... %d pairs set aside.\n', n_total_setaside);

% Split them into their buckets
instr_pair_foil = setaside_pairs(1:n_instr_foil, :);
prac_pairs = setaside_pairs(n_instr_foil+1 : end, :); % One bucket for all 25

%% SELECT MAIN FOILS
%--------------------------------------------------------------------------
fprintf('\nAssigning all remaining pairs as main foils...\n');
% Whatever is left in the available_pool is now the main foil pool
main_foils = available_pool;
n_foils = size(main_foils, 1); 
fprintf('Assigned %d pairs as main foils.\n', n_foils);

%% COPY EXPERIMENTAL PAIRS (L1 and L2)
%--------------------------------------------------------------------------
fprintf('\nCopying %d experimental pairs (%d L1, %d L2)...\n', ...
    n_experimental_pairs_total, n_bin1_pairs, n_bin2_pairs);
file_idx = 1;
% --- Copy L1 pairs ---
for i = 1:size(experimental_pairs_l1, 1)
    copyfile(experimental_pairs_l1{i,1}, fullfile(out_dir, sprintf('mst_%03d_A_l1.png', file_idx)));
    copyfile(experimental_pairs_l1{i,2}, fullfile(out_dir, sprintf('mst_%03d_B_l1.png', file_idx)));
    file_idx = file_idx + 1;
end
% --- Copy L2 pairs ---
for i = 1:size(experimental_pairs_l2, 1)
    copyfile(experimental_pairs_l2{i,1}, fullfile(out_dir, sprintf('mst_%03d_A_l2.png', file_idx)));
    copyfile(experimental_pairs_l2{i,2}, fullfile(out_dir, sprintf('mst_%03d_B_l2.png', file_idx)));
    file_idx = file_idx + 1;
end
%% COPY MAIN FOILS
%--------------------------------------------------------------------------
fprintf('Copying %d main foil pairs...\n', n_foils);
for i = 1:n_foils
    source_foil_a = main_foils{i,1};
    source_foil_b = main_foils{i,2};
    copyfile(source_foil_a, fullfile(out_dir, sprintf('mst_%03d_A_foil.png', i)));
    copyfile(source_foil_b, fullfile(out_dir, sprintf('mst_%03d_B_foil.png', i)));
end

%% COPY INSTRUCTION/PRACTICE STIMULI
%--------------------------------------------------------------------------
fprintf('\nCopying instruction and practice files...\n');

% --- Instruction pairs (1 L1, 1 Foil) ---
copyfile(instr_pair_l1{1,1}, fullfile(out_dir, 'instr_A.png'));
copyfile(instr_pair_l1{1,2}, fullfile(out_dir, 'instr_B.png'));
copyfile(instr_pair_foil{1,1}, fullfile(out_dir, 'instr_A_foil.png'));
copyfile(instr_pair_foil{1,2}, fullfile(out_dir, 'instr_B_foil.png'));

% --- Practice pairs (25 pairs, indexed 01-25) ---
for i = 1:height(prac_pairs)
    copyfile(prac_pairs{i,1}, fullfile(out_dir, sprintf('prac_%02d_A.png', i)));
    copyfile(prac_pairs{i,2}, fullfile(out_dir, sprintf('prac_%02d_B.png', i)));
end

fprintf('...copied 2 instruction pairs and %d practice pairs.\n', height(prac_pairs));

%% SUMMARY
%--------------------------------------------------------------------------
fprintf('\n------------------------------------------------------------\n');
fprintf('Final output folder: %s\n', out_dir);
fprintf('  -> Experimental pairs: %d total\n', n_experimental_pairs_total);
fprintf('  -> Main foils: %d pairs (all remaining pairs)\n', n_foils);
fprintf('  -> Instruction pairs: 2 (1 L1, 1 Foil)\n');
fprintf('  -> Practice pairs: %d\n', height(prac_pairs));
fprintf('------------------------------------------------------------\n');