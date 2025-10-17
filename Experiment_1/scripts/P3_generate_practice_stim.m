%==========================================================================
%                  Create Practice and Main Stimulus Sets
%==========================================================================
% author: Zihan Bai
%==========================================================================

clear; 
clc;
rng('shuffle');

%% DEFINE PARAMETERS AND DIRECTORIES
%--------------------------------------------------------------------------
% --- Define how many pairs you want for practice blocks ---
% We need a few for each task (1-back repeats, lures; 2-back repeats, lures)
% 20 pairs is a very safe number.
n_practice_pairs = 20; 

% --- Directories ---
base_dir = '..';
in_dir = fullfile(base_dir, 'stimulus/stim_norm/');
out_dir_base = fullfile(base_dir, 'stimulus/');

% The new output directories this script will create
out_dir_practice = fullfile(out_dir_base, 'stim_final');
out_dir_main = fullfile(out_dir_base, 'stim_main');

mst_set_folders = {'Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5', 'Set 6'};

%% CREATE OUTPUT DIRECTORIES
%--------------------------------------------------------------------------
if exist(out_dir_practice, 'dir') || exist(out_dir_main, 'dir')
    error('Output directories already exist. Please delete them before running this one-time script to ensure a clean partition.');
end
mkdir(out_dir_practice);
mkdir(out_dir_main);

%% LOAD ALL AVAILABLE MST PAIRS
%--------------------------------------------------------------------------
fprintf('loading all available pairs from original source...\n');
master_pair_list = [];
for i = 1:numel(mst_set_folders)
    for j = 1:192 % Assuming 192 pairs per set
        % We will store relative paths
        path_a_rel = fullfile(mst_set_folders{i}, sprintf('%03da.png', j));
        path_b_rel = fullfile(mst_set_folders{i}, sprintf('%03db.png', j));
        
        if exist(fullfile(in_dir, path_a_rel), 'file') && exist(fullfile(in_dir, path_b_rel), 'file')
            master_pair_list = [master_pair_list; {path_a_rel, path_b_rel}];
        end
    end
end
fprintf('found %d total pairs.\n', size(master_pair_list, 1));

%% PARTITION INTO PRACTICE AND MAIN SETS
%--------------------------------------------------------------------------
fprintf('partitioning into practice and main experiment sets...\n');
shuffled_indices = randperm(size(master_pair_list, 1));
master_pair_list = master_pair_list(shuffled_indices, :);

practice_pairs = master_pair_list(1:n_practice_pairs, :);
main_exp_pairs = master_pair_list(n_practice_pairs+1 : end, :);
fprintf('%d pairs selected for practice, %d pairs remaining for main experiment.\n', ...
    size(practice_pairs, 1), size(main_exp_pairs, 1));

%% PROCESS AND RENAME PRACTICE PAIRS
%--------------------------------------------------------------------------
fprintf('copying and renaming practice pairs...\n');
for i = 1:n_practice_pairs
    source_path_a = fullfile(in_dir, practice_pairs{i, 1});
    source_path_b = fullfile(in_dir, practice_pairs{i, 2});
    
    % Give them simple, clear names
    new_name_a = sprintf('practice_%02d_a.png', i);
    new_name_b = sprintf('practice_%02d_b.png', i);
    
    dest_path_a = fullfile(out_dir_practice, new_name_a);
    dest_path_b = fullfile(out_dir_practice, new_name_b);
    
    copyfile(source_path_a, dest_path_a);
    copyfile(source_path_b, dest_path_b);
end
fprintf('practice set created in "%s".\n', out_dir_practice);

%% COPY REMAINING MAIN EXPERIMENT PAIRS
%--------------------------------------------------------------------------
fprintf('copying main experiment pairs to new source directory...\n');
for i = 1:size(main_exp_pairs, 1)
    [set_folder, name, ext] = fileparts(main_exp_pairs{i,1});
    
    % Create the new Set subfolder if it doesn't exist
    new_set_dir = fullfile(out_dir_main, set_folder);
    if ~exist(new_set_dir, 'dir'), mkdir(new_set_dir); end
    
    % Copy file A
    source_path_a = fullfile(in_dir, main_exp_pairs{i, 1});
    dest_path_a   = fullfile(out_dir_main, main_exp_pairs{i, 1});
    copyfile(source_path_a, dest_path_a);
    
    % Copy file B
    source_path_b = fullfile(in_dir, main_exp_pairs{i, 2});
    dest_path_b   = fullfile(out_dir_main, main_exp_pairs{i, 2});
    copyfile(source_path_b, dest_path_b);
end
fprintf('main experiment source created in "%s".\n', out_dir_main);


%% DONE
%--------------------------------------------------------------------------
fprintf('\nSUCCESS! Partitioning is complete.\n');
fprintf('Step 1: Your practice stimuli are ready in "%s".\n', out_dir_practice);
fprintf('Step 2: You should now point your ORIGINAL `generate_stimulus_set.m` script to read from "%s".\n', out_dir_main);