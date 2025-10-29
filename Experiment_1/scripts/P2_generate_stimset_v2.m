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

base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/stim_norm/');
out_dir = fullfile(base_dir, 'stimulus/stim_1028');
mst_set_folders = {'Set 1','Set 2','Set 3','Set 4','Set 5','Set 6'};
lure_bin_files  = {'Set1_bins.txt','Set2_bins.txt','Set3_bins.txt',...
                   'Set4_bins.txt','Set5_bins.txt','Set6_bins.txt'};
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
%% LOAD ALL STIMULUS PAIRS AND LURE-BIN LABELS115
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


%% SELECT MAIN FOILS (ALL REMAINING PAIRS)
%--------------------------------------------------------------------------
fprintf('\nSelecting main foils from ALL remaining pairs...\n');

% --- Create a complete list of ALL experimental images (A and B) ---
used_exp_imgs = [string(experimental_pairs_l1(:,1)); string(experimental_pairs_l1(:,2)); ...
                 string(experimental_pairs_l2(:,1)); string(experimental_pairs_l2(:,2))];
used_all_images = unique(used_exp_imgs);

% --- Create the *truly* available pool (ANY pair not used for experiment) ---
% A pair is available ONLY IF *neither* its 'a' nor 'b' image is in the used list
is_pair_available = ~ismember(string(master_pair_list(:,1)), used_all_images) & ...
                    ~ismember(string(master_pair_list(:,2)), used_all_images);
                    
main_foils = master_pair_list(is_pair_available, :);
n_foils = size(main_foils, 1); % Get the final count for logging

fprintf('Assigned all %d remaining pairs as main foils.\n', n_foils);


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

%% COPY MAIN FOILS (AS PAIRS: A-image=foil, B-image=foil_lure)
%--------------------------------------------------------------------------
fprintf('Copying %d main foil pairs (a-images as _foil, b-images as _foil_lure)...\n', n_foils);
for i = 1:n_foils
    source_foil_a = main_foils{i,1};
    source_foil_b = main_foils{i,2};
    copyfile(source_foil_a, fullfile(out_dir, sprintf('mst_%03d_foil.png', i)));
    copyfile(source_foil_b, fullfile(out_dir, sprintf('mst_%03d_foil_lure.png', i)));
end

%% SUMMARY
%--------------------------------------------------------------------------
fprintf('\nDone!\n');
fprintf('Final output folder: %s\n', out_dir);
fprintf('  -> Experimental pairs: %d total\n', n_experimental_pairs_total);
fprintf('     -> %d from Bin 1\n', n_bin1_pairs);
fprintf('     -> %d from Bin 2\n', n_bin2_pairs);
fprintf('  -> Main foils: %d pairs (all remaining pairs)\n', n_foils);
fprintf('------------------------------------------------------------\n');