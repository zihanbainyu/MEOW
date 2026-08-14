%==========================================================================
%                  MATCH LUMINANCE (in place, on the final pool)
%==========================================================================
% Author: Zihan Bai
%
% Run order: P1 (normalize -> stim_norm) -> P2 (select -> stim_pool) -> P3.
% P3 matches the mean luminance & contrast of the images subjects actually see
% (the selected pool) and OVERWRITES them in place, so A_subject_setup loads the
% matched files. No resizing here -- P1 already sets the canvas size (300x300).
%
% NOTE: re-running P2 re-copies fresh, UNMATCHED images from stim_norm into
% stim_pool, so P3 must be run again after every P2.
%==========================================================================
clear;
clc;

%% SETUP
%--------------------------------------------------------------------------
base_dir   = '..';
shine_path = fullfile(base_dir, 'scripts/functions/SHINE_toolbox');
pool_dir   = fullfile(base_dir, 'stimulus/stim_pool');   % matched IN PLACE
%--------------------------------------------------------------------------
addpath(shine_path);
fprintf('Pool folder (matched in place): %s\n\n', pool_dir);


%% 1. GATHER ALL POOL IMAGES
%--------------------------------------------------------------------------
fprintf('Finding all .png images in the pool...\n');
image_files = dir(fullfile(pool_dir, '*.png'));         % stim_pool is flat
if isempty(image_files)
    error('No .png files found in %s. Did you run P2_generate_stimset.m first?', pool_dir);
end
all_source_paths = fullfile({image_files.folder}, {image_files.name})';
fprintf('Found %d pool images to match.\n', numel(all_source_paths));


%% 2. BATCH LUMINANCE PROCESSING & STATISTICS
%--------------------------------------------------------------------------
fprintf('Loading image data for batch processing...\n');
all_images_raw = cell(numel(all_source_paths), 1);
for i = 1:numel(all_source_paths)
    all_images_raw{i} = imread(all_source_paths{i});    % read only; file unchanged
end

% --- Calculate statistics BEFORE matching ---
stats_before = get_luminance_stats(all_images_raw);

% --- Use SHINE to match luminance and contrast ---
fprintf('Matching luminance and contrast across all %d images...\n', numel(all_images_raw));
all_images_lum_matched = lumMatch(all_images_raw);      % in memory
fprintf('Luminance matching complete.\n\n');

% --- Calculate statistics AFTER matching ---
stats_after = get_luminance_stats(all_images_lum_matched);


%% 3. REPORT LUMINANCE STATISTICS
%--------------------------------------------------------------------------
fprintf('--- LUMINANCE STATISTICS (0-255 scale) ---\n');
fprintf('BEFORE Matching:\n');
fprintf('  Mean of means: %.2f\n', mean(stats_before.means));
fprintf('  Std Dev of means: %.2f\n', std(stats_before.means));
fprintf('  Range of means: [%.2f to %.2f]\n', min(stats_before.means), max(stats_before.means));
fprintf('\n');
fprintf('AFTER Matching:\n');
fprintf('  Mean of means: %.2f\n', mean(stats_after.means));
fprintf('  Std Dev of means: %.2f (should be close to 0)\n', std(stats_after.means));
fprintf('  Range of means: [%.2f to %.2f]\n', min(stats_after.means), max(stats_after.means));
fprintf('------------------------------------------\n\n');


%% 4. SAVE MATCHED IMAGES (overwrite the pool in place)
%--------------------------------------------------------------------------
fprintf('Overwriting pool images with luminance-matched versions...\n');
for i = 1:numel(all_source_paths)
    imwrite(all_images_lum_matched{i}, all_source_paths{i});   % same path, PNG

    if mod(i, 100) == 0 || i == numel(all_source_paths)
        fprintf('  Saved %d / %d images.\n', i, numel(all_source_paths));
    end
end

fprintf('\nAll %d pool images matched and saved in place:\n%s\n', ...
    numel(all_source_paths), pool_dir);


%% HELPER FUNCTIONS
%--------------------------------------------------------------------------
function stats = get_luminance_stats(image_cell_array)
    num_images = numel(image_cell_array);
    means = zeros(num_images, 1);
    stds = zeros(num_images, 1);
    for i = 1:num_images
        gray_img = rgb2gray(image_cell_array{i});
        means(i) = mean(gray_img(:));
        stds(i) = std(double(gray_img(:)));
    end
    stats.means = means;
    stats.stds = stds;
end
