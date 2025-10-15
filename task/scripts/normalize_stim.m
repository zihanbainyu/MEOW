%--------------------------------------------------------------------------
%                       IMAGE NORMALIZATION SCRIPT
%--------------------------------------------------------------------------

%% 1. SETUP PARAMETERS
canvas_size = 400;                        % final canvas size
canvas_color = [255, 255, 255];           % white background
root_directory = fullfile(pwd, 'Sets1_6'); % parent folder containing Set 1...Set 6
output_directory = 'norm_stim';           % main output folder
file_types = {'*.jpg'};                   % file types

%--------------------------------------------------------------------------
%% 2. LOOP OVER ALL SETS
fprintf('Starting image normalization process...\n');

% Create main output directory if it doesn't exist
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Loop over Set 1...Set 6
for set_idx = 1:6
    input_directory = fullfile(root_directory, sprintf('Set %d', set_idx));
    set_output_directory = fullfile(output_directory, sprintf('Set %d', set_idx));

    % Create subfolder for this set
    if ~exist(set_output_directory, 'dir')
        mkdir(set_output_directory);
    end

    % Gather images
    image_files = [];
    for i = 1:length(file_types)
        files = dir(fullfile(input_directory, file_types{i}));
        image_files = [image_files; files]; %#ok<AGROW>
    end

    if isempty(image_files)
        fprintf('No images found in %s. Skipping...\n', input_directory);
        continue;
    end

    fprintf('Found %d images in %s.\n', length(image_files), input_directory);

    % Prepare canvas
    canvas = uint8(zeros(canvas_size, canvas_size, 3));
    canvas(:,:,1) = canvas_color(1);
    canvas(:,:,2) = canvas_color(2);
    canvas(:,:,3) = canvas_color(3);

    % Process each image
    for k = 1:length(image_files)
        current_filename = image_files(k).name;
        fprintf('Processing (%d/%d) in Set %d: %s\n', ...
                k, length(image_files), set_idx, current_filename);

        try
            [img, ~, alpha] = imread(fullfile(input_directory, current_filename));

            % Convert grayscale → RGB
            if size(img,3) == 1
                img = repmat(img, [1 1 3]);
            end

            % Scale with aspect ratio
            scale_factor = min(canvas_size / size(img,1), canvas_size / size(img,2));
            new_height = round(size(img,1) * scale_factor);
            new_width  = round(size(img,2) * scale_factor);
            resized_img = imresize(img, [new_height, new_width]);

            % Resize alpha if present
            if ~isempty(alpha)
                resized_alpha = imresize(alpha, [new_height, new_width]);
            end

            % Center position
            y_pos = floor((canvas_size - new_height) / 2) + 1;
            x_pos = floor((canvas_size - new_width) / 2) + 1;

            % Place onto canvas
            final_image = canvas;
            final_image(y_pos:(y_pos+new_height-1), x_pos:(x_pos+new_width-1), :) = resized_img;

            % Transparency handling
            final_alpha = [];
            if ~isempty(alpha)
                final_alpha = uint8(255*ones(canvas_size, canvas_size));
                final_alpha(y_pos:(y_pos+new_height-1), x_pos:(x_pos+new_width-1)) = resized_alpha;
            end

            % Save output
            [~, name, ~] = fileparts(current_filename);
            out_path = fullfile(set_output_directory, [name '.png']); % standardize to PNG
            if ~isempty(final_alpha)
                imwrite(final_image, out_path, 'Alpha', final_alpha);
            else
                imwrite(final_image, out_path);
            end

        catch ME
            fprintf(2, '!!! Error with %s: %s\n', current_filename, ME.message);
            continue;
        end
    end
end

fprintf('\nNormalization complete!\nAll sets saved inside "%s".\n', output_directory);
