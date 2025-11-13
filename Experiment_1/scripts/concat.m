clear;
clc;

% Configuration
subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id);
results_dir = fullfile(base_dir, 'data', subj_folder);
nBlocks = 4;

% Initialize the concatenated struct
in_concat = struct();
in_concat.xPos = cell(1, nBlocks);
in_concat.yPos = cell(1, nBlocks);
in_concat.pupilArea = cell(1, nBlocks);
in_concat.sampleRate = cell(1, nBlocks);
in_concat.startInds = cell(1, nBlocks);
in_concat.trialTypes = cell(1, nBlocks);

% Load each block and populate the concatenated struct
fprintf('Loading and concatenating blocks...\n');
for b = 1:nBlocks
    model_input_file = fullfile(results_dir, sprintf('%03d_2_pcdm_input.mat', subj_id, b));
    
    if ~exist(model_input_file, 'file')
        warning('Block %d input file not found: %s', b, model_input_file);
        continue;
    end
    
    % Load the block
    load(model_input_file, 'in');
    
    % Concatenate each field
    in_concat.xPos{b} = in.xPos{b};
    in_concat.yPos{b} = in.yPos{b};
    in_concat.pupilArea{b} = in.pupilArea{b};
    in_concat.sampleRate{b} = in.sampleRate{b};
    in_concat.startInds{b} = in.startInds{b};
    in_concat.trialTypes{b} = in.trialTypes{b};
    
    fprintf('  Block %d loaded and added.\n', b);
end

% Rename for clarity
in = in_concat;

% Save the concatenated input
output_file = fullfile(results_dir, sprintf('%03d_2_pcdm_input_all_blocks.mat', subj_id));
save(output_file, 'in');

fprintf('\nConcatenation complete!\n');
fprintf('Saved to: %s\n', output_file);
fprintf('Structure: 1x1 struct with 6 fields\n');
fprintf('Each field: 1x%d cell array\n', nBlocks);

% Display summary
fprintf('\n--- Summary ---\n');
for b = 1:nBlocks
    if ~isempty(in.xPos{b})
        fprintf('Block %d: %d samples, %d trials\n', b, length(in.xPos{b}), size(in.startInds{b}, 1));
    else
        fprintf('Block %d: EMPTY\n', b);
    end
end