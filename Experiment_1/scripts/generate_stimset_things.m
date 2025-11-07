clear;
clc;
rng('shuffle'); % Seed the random number generator for A/B assignment

%% SETUP
%--------------------------------------------------------------------------
base_dir = '..';
in_dir  = fullfile(base_dir, 'stimulus/main_ready/');
out_dir = fullfile(base_dir, 'stimulus/stim_things_final');

% Create output directory if it doesn't exist
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
    fprintf('Created output directory: %s\n', out_dir);
end

% Get a list of all image files (assuming they are all PNGs)
all_files = dir(fullfile(in_dir, '*.jpg'));
if isempty(all_files)
    error('No PNG files found in the input directory: %s', in_dir);
end

fprintf('Found %d total image files.\n', length(all_files));
if mod(length(all_files), 2) ~= 0
    warning('The total number of files (%d) is odd. Pairing might fail.', length(all_files));
end

%% GROUPING AND PARSING FILES
%--------------------------------------------------------------------------
% A map to store grouped file names based on category and number
paired_groups = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(all_files)
    filename_full = all_files(i).name;
    
    % Remove the extension for clean parsing
    [~, filename_no_ext, ~] = fileparts(filename_full);
    
    % Find all underscore positions
    underscore_indices = strfind(filename_no_ext, '_');
    
    if isempty(underscore_indices)
        fprintf('Skipping file: %s (No underscore found)\n', filename_full);
        continue;
    end
    
    % Find the position of the *LAST* underscore
    last_underscore_index = underscore_indices(end);
    
    % The shared key is everything BEFORE the last underscore
    % e.g., 'apple_07b' -> 'apple_07'
    group_key = filename_no_ext(1:last_underscore_index-1);
    
    if ~isempty(group_key)
        % Store the full file path in the map under its key
        full_path = fullfile(in_dir, filename_full);
        if isKey(paired_groups, group_key)
            current_list = paired_groups(group_key);
            paired_groups(group_key) = [current_list; full_path];
        else
            paired_groups(group_key) = {full_path};
        end
    else
        fprintf('Skipping file: %s (Group key is empty)\n', filename_full);
    end
end

num_groups = paired_groups.Count;
fprintf('Found **%d** potential pairs based on filename parsing.\n', num_groups);

%% RENAME AND COPY FILES
%--------------------------------------------------------------------------
keys = paired_groups.keys;
pair_counter = 0; % Counter for the new sequential filename (001, 002, ...)

for k = 1:num_groups
    key = keys{k};
    group_files = paired_groups(key);
    
    if length(group_files) == 2
        pair_counter = pair_counter + 1;
        
        % Randomly assign A and B
        assignment = randperm(2); 
        file_A = group_files{assignment(1)};
        file_B = group_files{assignment(2)};
        
        % Construct the new filenames
        new_name_A = sprintf('things_%03d_A.jpg', pair_counter);
        new_name_B = sprintf('things_%03d_B.jpg', pair_counter);
        
        new_path_A = fullfile(out_dir, new_name_A);
        new_path_B = fullfile(out_dir, new_name_B);
        
        % Copy and rename the files
        [status_A, ~] = copyfile(file_A, new_path_A);
        [status_B, ~] = copyfile(file_B, new_path_B);
        
        if status_A == 0 || status_B == 0
            warning('Error copying pair %d (%s). Check file permissions.', pair_counter, key);
        end

    elseif length(group_files) ~= 2
        fprintf('**Skipping key %s: Found %d files, expected 2 for a pair.**\n', key, length(group_files));
    end
end

fprintf('A total of **%d** pairs (should be 360) were processed and saved to:\n', pair_counter);
fprintf('%s\n', out_dir);