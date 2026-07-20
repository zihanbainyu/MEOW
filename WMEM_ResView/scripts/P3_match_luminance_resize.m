%==========================================================================
%                  MATCH LUMINANCE AND RESIZE
%==========================================================================
% Author: Zihan Bai
%==========================================================================
clear;
clc;
rng('shuffle');

%% SETUP
%--------------------------------------------------------------------------
base_dir   = '..'; 
shine_path = fullfile(base_dir, 'scripts/functions/SHINE_toolbox');
in_dir     = fullfile(base_dir, 'stimulus/stim_selected'); 
out_dir    = fullfile(base_dir, 'stimulus/stim_processed');
p.target_size      = [500, 500];
%--------------------------------------------------------------------------
addpath(shine_path);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fprintf('Input folder: %s\n', in_dir);
fprintf('Output folder: %s\n\n', out_dir);


%% 1. GATHER ALL IMAGE FILES
%--------------------------------------------------------------------------
fprintf('Finding all .png images in the input directory...\n');
image_files = dir(fullfile(in_dir, '**', '*.png'));
if isempty(image_files)
    error('No .png files found in %s. Did you run generate_stim_set.m first?', in_dir);
end
all_source_paths = fullfile({image_files.folder}, {image_files.name})';
fprintf('Found %d unique images to process.\n', numel(all_source_paths));


%% 2. BATCH LUMINANCE PROCESSING & STATISTICS
%--------------------------------------------------------------------------
fprintf('Loading image data for batch processing...\n');
all_images_raw = cell(numel(all_source_paths), 1);
for i = 1:numel(all_source_paths)
    % This step only READS data into a variable. The original file is not changed.
    all_images_raw{i} = imread(all_source_paths{i});
end

% --- Calculate statistics BEFORE matching ---
stats_before = get_luminance_stats(all_images_raw);

% --- Use SHINE to match luminance and contrast ---
fprintf('Matching luminance and contrast across all %d images...\n', numel(all_images_raw));
% This operation happens entirely in MATLAB's memory.
all_images_lum_matched = lumMatch(all_images_raw);
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


%% 4. RESIZE AND SAVE PROCESSED IMAGES
%--------------------------------------------------------------------------
fprintf('Resizing and saving processed images to %dx%d...\n', p.target_size(1), p.target_size(2));

for i = 1:numel(all_source_paths)
    source_path = all_source_paths{i};
    img_data = all_images_lum_matched{i}; % Use the processed data from memory
    
    [~, filename, ext] = fileparts(source_path);
    dest_path = fullfile(out_dir, [filename, ext]); % Define path for the NEW file
    
    % This function SAVES A NEW FILE to the destination path.
    resize_and_save_image(img_data, dest_path, p.target_size);
    
    if mod(i, 100) == 0 || i == numel(all_source_paths)
        fprintf('  Processed and saved %d / %d images.\n', i, numel(all_source_paths));
    end
end

fprintf('\nAll images have been processed and saved to:\n%s\n', out_dir);


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

function resize_and_save_image(img_data, dest_path, target_size)
    % This function takes image data from memory and saves it to a new file.
    resized_img = imresize(img_data, [target_size(2), target_size(1)], 'bicubic');
    imwrite(resized_img, dest_path); % Creates the new file on disk
end