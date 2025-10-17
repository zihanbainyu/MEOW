%==========================================================================
%                  Generate stimulus set
%==========================================================================
% author: Zihan Bai

% note:
%   load standard MST stimulus sets (Kirwan & Stark, 2007)
%   randomly, equally select pairs with lure bin = 2 for experimental use (for now)
%   randomly select practice pairs and foils
%   randomly select experimental foils
%   copy and rename into stim_selected folder
clear; 
clc;
rng('shuffle');

%% SETUP
%--------------------------------------------------------------------------
n_experimental_pairs = 120;   % main experimental pairs
n_foils = 640;                % total number of unrelated foils
n_practice_pairs = 20;        % practice experimental pairs
n_practice_foils = 40;        % practice foils (a images only)

target_lure_bin = 2;          % select only pairs from this lure bin

base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/stim_norm/');
out_dir = fullfile(base_dir, 'stimulus/stim_selected');

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


%% REPORT: PAIRS PER LURE BIN PER SET
%--------------------------------------------------------------------------
fprintf('\nPairs per lure bin for each set:\n');
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(master_pair_list(:,4), set_name);
    lure_bins = all_bins(in_set);
    fprintf('\n%s:\n', set_name);
    for b = 1:5
        fprintf('  Bin %d: %d pairs\n', b, sum(lure_bins==b));
    end
    fprintf('  --> Total with lure bin = %d: %d pairs\n', ...
            target_lure_bin, sum(lure_bins==target_lure_bin));
end


%% SELECT EXPERIMENTAL PAIRS (BIN=2)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d total experimental pairs (≈20 per set)...\n', n_experimental_pairs);

pairs_per_set = floor(n_experimental_pairs / numel(mst_set_folders));
experimental_pairs = [];

for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(master_pair_list(:,4), set_name);
    in_bin = all_bins == target_lure_bin;
    set_bin_pairs = master_pair_list(in_set & in_bin, :);

    if size(set_bin_pairs,1) < pairs_per_set
        error('%s has insufficient bin-%d pairs.', set_name, target_lure_bin);
    end

    idx = randperm(size(set_bin_pairs,1), pairs_per_set);
    experimental_pairs = [experimental_pairs; set_bin_pairs(idx,:)];
end

fprintf('Selected %d experimental pairs total.\n', size(experimental_pairs,1));


%% SELECT PRACTICE PAIRS AND FOILS (ALSO BIN=2)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d practice pairs and %d practice foils (from remaining bin-2 pairs)...\n', ...
        n_practice_pairs, n_practice_foils);

% remaining bin-2 pairs (exclude experimental)
bin2_pairs = master_pair_list(all_bins==target_lure_bin, :);
used_targets = string(experimental_pairs(:,1));
remaining_bin2 = bin2_pairs(~ismember(string(bin2_pairs(:,1)), used_targets), :);

% select practice pairs
prac_idx = randperm(size(remaining_bin2,1), n_practice_pairs);
practice_pairs = remaining_bin2(prac_idx,:);
remaining_bin2(prac_idx,:) = [];

% select practice foils (a images only)
foil_idx = randperm(size(remaining_bin2,1), n_practice_foils);
practice_foils = remaining_bin2(foil_idx,:);

fprintf('Remaining bin-2 pool after practice selection: %d pairs.\n', size(remaining_bin2,1));


%% SELECT MAIN FOILS (A-IMAGES ONLY, BALANCED ACROSS SETS)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d main foils (balanced across sets)...\n', n_foils);

foils_per_set = floor(n_foils / numel(mst_set_folders));
main_foils = [];

used_all = [experimental_pairs(:,1); practice_pairs(:,1)];

for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(master_pair_list(:,4), set_name);
    available = master_pair_list(in_set & ~ismember(master_pair_list(:,1), used_all), :);

    n_this = min(foils_per_set, size(available,1));
    idx = randperm(size(available,1), n_this);
    main_foils = [main_foils; available(idx,:)];
end

% fill up remainder if needed
while size(main_foils,1) < n_foils
    extra_idx = randi(size(master_pair_list,1));
    if ~ismember(master_pair_list(extra_idx,1), [used_all; main_foils(:,1)])
        main_foils = [main_foils; master_pair_list(extra_idx,:)];
    end
end

fprintf('Selected %d main foils total.\n', size(main_foils,1));


%% COPY EXPERIMENTAL PAIRS
%--------------------------------------------------------------------------
fprintf('\nCopying %d experimental pairs...\n', n_experimental_pairs);
for i = 1:n_experimental_pairs
    copyfile(experimental_pairs{i,1}, fullfile(out_dir, sprintf('mst_%03d_targ_l2.png', i)));
    copyfile(experimental_pairs{i,2}, fullfile(out_dir, sprintf('mst_%03d_lure_l2.png', i)));
end


%% COPY MAIN FOILS (ONLY "A" IMAGES)
%--------------------------------------------------------------------------
fprintf('Copying %d main foils (a-images only)...\n', n_foils);
for i = 1:n_foils
    source_foil = main_foils{i,1};
    copyfile(source_foil, fullfile(out_dir, sprintf('mst_%03d_foil.png', i)));
end


%% COPY PRACTICE PAIRS & FOILS
%--------------------------------------------------------------------------
practice_dir = fullfile(out_dir, 'practice');
if ~exist(practice_dir, 'dir'), mkdir(practice_dir); end

fprintf('\nCopying %d practice pairs...\n', n_practice_pairs);
for i = 1:n_practice_pairs
    copyfile(practice_pairs{i,1}, fullfile(practice_dir, sprintf('prac_%03d_targ_l2.png', i)));
    copyfile(practice_pairs{i,2}, fullfile(practice_dir, sprintf('prac_%03d_lure_l2.png', i)));
end

fprintf('Copying %d practice foils (a-images only)...\n', n_practice_foils);
for i = 1:n_practice_foils
    source_foil = practice_foils{i,1};
    copyfile(source_foil, fullfile(practice_dir, sprintf('prac_foil_%03d.png', i)));
end


%% SUMMARY
%--------------------------------------------------------------------------
fprintf('\nDone!\n');
fprintf('Final output folder: %s\n', out_dir);
fprintf('  -> Experimental pairs: %d (bin=%d)\n', n_experimental_pairs, target_lure_bin);
fprintf('  -> Main foils: %d (a-images only)\n', n_foils);
fprintf('  -> Practice pairs: %d (bin=%d)\n', n_practice_pairs, target_lure_bin);
fprintf('  -> Practice foils: %d (a-images only)\n', n_practice_foils);
fprintf('------------------------------------------------------------\n');
