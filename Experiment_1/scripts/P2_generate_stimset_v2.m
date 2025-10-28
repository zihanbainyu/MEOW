%==========================================================================
%                  Generate stimulus set
%==========================================================================
% author: Zihan Bai
% note:
%   load standard MST stimulus sets (Kirwan & Stark, 2007)
%   randomly, equally select pairs with lure bin = 1 for experimental use
%   randomly select practice pairs and foils
%   randomly select experimental foils (MUST EXCLUDE BIN 1 AND BIN 2 PAIRS)
%   copy and rename into stim_selected folder
clear; 
clc;
rng('shuffle');
%% SETUP
%--------------------------------------------------------------------------
n_experimental_pairs_total = 360; % 180 from Bin 1, 180 from Bin 2
n_bin1_pairs = 180;
n_bin2_pairs = 180;
n_foils = 640;                % total number of unrelated foils
n_practice_pairs = 10;         % practice experimental pairs
n_practice_foils = 20;        % practice foils (a images only)

base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/stim_norm/');
out_dir = fullfile(base_dir, 'stimulus/stim');
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

% This is the *correct* way to do it: balance by set *within* each bin
n_bin1_per_set = floor(n_bin1_pairs / numel(mst_set_folders)); % 180/6 = 30
n_bin2_per_set = floor(n_bin2_pairs / numel(mst_set_folders)); % 180/6 = 30

experimental_pairs_l1 = [];
experimental_pairs_l2 = [];

% --- Get Bin 1 Pairs (Balanced by Set) ---
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

% --- Get Bin 2 Pairs (Balanced by Set) ---
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


%% SELECT MAIN FOILS (A-IMAGES ONLY, EXCLUDING BIN 1, BIN 2, AND EXPERIMENTAL PAIRS)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d main foils (balanced across sets, EXCLUDING Bin 1, Bin 2, and Experimental pairs)...\n', n_foils);
foils_per_set = floor(n_foils / numel(mst_set_folders));
main_foils = [];

% --- Create a complete list of ALL experimental images (A and B) ---
used_exp_imgs = [string(experimental_pairs_l1(:,1)); string(experimental_pairs_l1(:,2)); ...
                 string(experimental_pairs_l2(:,1)); string(experimental_pairs_l2(:,2))];
used_all_images = unique(used_exp_imgs);

fprintf('Identifying main foils. %d experimental images are off-limits.\n', numel(used_all_images));

% Identify pairs NOT in Lure Bin 1 or Lure Bin 2 (i.e., Bin 3, 4, or 5)
not_bin1_or_2 = all_bins > 2; 
all_candidate_pairs = master_pair_list(not_bin1_or_2, :);

% --- Create the *truly* available pool ---
% A pair is available ONLY IF *neither* its 'a' nor 'b' image is in the used list
is_pair_available = ~ismember(string(all_candidate_pairs(:,1)), used_all_images) & ...
                    ~ismember(string(all_candidate_pairs(:,2)), used_all_images);
                    
all_available_foils = all_candidate_pairs(is_pair_available, :);

fprintf('Found %d *truly available* foil pairs (Bins 3-5, not used elsewhere).\n', ...
        size(all_available_foils,1));

% Check if there are enough pairs remaining for the main foils
if size(all_available_foils, 1) < n_foils
    error('INSUFFICIENT_FOILS: Not enough pairs remain in Bins 3-5 to select %d unused foils.', n_foils);
end

% --- Select foils, balanced across sets ---
for s = 1:numel(mst_set_folders)
    set_name = mst_set_folders{s};
    in_set = strcmp(all_available_foils(:,4), set_name);
    
    available_in_set = all_available_foils(in_set, :);
    
    n_this = min(foils_per_set, size(available_in_set,1));
    if (n_this == 0)
        fprintf('Warning: No available foils found for %s.\n', set_name);
        continue;
    end
    
    idx = randperm(size(available_in_set,1), n_this);
    main_foils = [main_foils; available_in_set(idx,:)];
    
    % Remove these selected pairs from the available pool
    all_available_foils(ismember(all_available_foils(:,1), available_in_set(idx,1)), :) = [];
end

% fill up remainder if needed
fprintf('Selected %d foils. Need %d more to reach %d.\n', ...
        size(main_foils,1), n_foils - size(main_foils,1), n_foils);
        
n_needed = n_foils - size(main_foils,1);
if n_needed > 0
    if size(all_available_foils, 1) < n_needed
        error('INSUFFICIENT_FOILS: Not enough foils left for fill-up.');
    end
    
    % Grab from the remaining shuffled pool
    idx = randperm(size(all_available_foils,1), n_needed);
    main_foils = [main_foils; all_available_foils(idx,:)];
end

fprintf('Selected %d main foils total.\n', size(main_foils,1));


%% SELECT PRACTICE PAIRS AND FOILS (ALSO BIN=1)
%--------------------------------------------------------------------------
fprintf('\nSelecting %d practice pairs and %d practice foils (from remaining bin-1 pairs)...\n', ...
        n_practice_pairs, n_practice_foils);

% Identify remaining bin-1 pairs (exclude experimental)
all_bin1_pairs = master_pair_list(all_bins==1, :);
used_l1_targets = string(experimental_pairs_l1(:,1));
remaining_bin1 = all_bin1_pairs(~ismember(string(all_bin1_pairs(:,1)), used_l1_targets), :);

% Check if we have enough left for practice
if size(remaining_bin1,1) < (n_practice_pairs + n_practice_foils)
    error('Insufficient *remaining* Bin 1 pairs for practice items (need %d, found %d).', ...
          (n_practice_pairs + n_practice_foils), size(remaining_bin1,1));
end

% select practice pairs
prac_idx = randperm(size(remaining_bin1,1), n_practice_pairs);
practice_pairs = remaining_bin1(prac_idx,:);
remaining_bin1(prac_idx,:) = []; % Remove them from the pool

% select practice foils (a images only)
foil_idx = randperm(size(remaining_bin1,1), n_practice_foils);
practice_foils = remaining_bin1(foil_idx,:);
fprintf('Remaining bin-1 pool after practice selection: %d pairs.\n', size(remaining_bin1,1));

%% COPY EXPERIMENTAL PAIRS (L1 and L2)
%--------------------------------------------------------------------------
fprintf('\nCopying %d experimental pairs (%d L1, %d L2)...\n', ...
        n_experimental_pairs_total, n_bin1_pairs, n_bin2_pairs);

file_idx = 1; % Global counter for filenames

% --- Copy L1 pairs ---
for i = 1:size(experimental_pairs_l1, 1)
    copyfile(experimental_pairs_l1{i,1}, fullfile(out_dir, sprintf('mst_%03d_targ_l1.png', file_idx)));
    copyfile(experimental_pairs_l1{i,2}, fullfile(out_dir, sprintf('mst_%03d_lure_l1.png', file_idx)));
    file_idx = file_idx + 1;
end

% --- Copy L2 pairs ---
for i = 1:size(experimental_pairs_l2, 1)
    copyfile(experimental_pairs_l2{i,1}, fullfile(out_dir, sprintf('mst_%03d_targ_l2.png', file_idx)));
    copyfile(experimental_pairs_l2{i,2}, fullfile(out_dir, sprintf('mst_%03d_lure_l2.png', file_idx)));
    file_idx = file_idx + 1;
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
    copyfile(practice_pairs{i,1}, fullfile(practice_dir, sprintf('prac_%03d_targ_l1.png', i)));
    copyfile(practice_pairs{i,2}, fullfile(practice_dir, sprintf('prac_%03d_lure_l1.png', i)));
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
fprintf('  -> Experimental pairs: %d total\n', n_experimental_pairs_total);
fprintf('     -> %d from Bin 1\n', n_bin1_pairs);
fprintf('     -> %d from Bin 2\n', n_bin2_pairs);
fprintf('  -> Main foils: %d (a-images only, EXCLUDING Bin 1, 2, and all used pairs)\n', n_foils);
fprintf('  -> Practice pairs: %d (bin=1)\n', n_practice_pairs);
fprintf('  -> Practice foils: %d (a-images only, from remaining Bin 1)\n', n_practice_foils);
fprintf('------------------------------------------------------------\n');