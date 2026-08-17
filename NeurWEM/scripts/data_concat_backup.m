% One-off: rebuild sub103_concat.mat.
% The post-task MST was run separately, which left sub103_concat.mat holding
% only results_mst. This concatenates the 1-back / 2-back block files into
% results_1_back_all / results_2_back_all exactly as main.m does, and keeps
% the existing results_mst.

subj_id = 101;
nBlocks = 4;
results_dir = fullfile('..', 'data', sprintf('sub%03d', subj_id));
concat_file = fullfile(results_dir, sprintf('sub%03d_concat.mat', subj_id));

% consolidate block files (same logic as consolidate_data in main.m)
results_1_back_all = consolidate_data(results_dir, subj_id, nBlocks, '1_back');
results_2_back_all = consolidate_data(results_dir, subj_id, nBlocks, '2_back');

% preserve the MST that was run separately
loaded = load(concat_file, 'final_data_output');
if isfield(loaded, 'final_data_output') && isfield(loaded.final_data_output, 'results_mst')
    results_mst = loaded.final_data_output.results_mst;
else
    mst_file = fullfile(results_dir, sprintf('sub%03d_mst.mat', subj_id));
    tmp = load(mst_file, 'results_mst');
    results_mst = tmp.results_mst;
end

% assemble final_data_output as in main.m
final_data_output = struct();
final_data_output.subj_id = subj_id;
final_data_output.results_1_back_all = results_1_back_all;
final_data_output.results_2_back_all = results_2_back_all;
final_data_output.results_mst = results_mst;

save(concat_file, 'final_data_output');
fprintf('Rebuilt %s\n', concat_file);
fprintf('  results_1_back_all: %d rows\n', height(results_1_back_all));
fprintf('  results_2_back_all: %d rows\n', height(results_2_back_all));
fprintf('  results_mst: %d rows\n', height(results_mst));

function full_data_table = consolidate_data(results_dir, subj_id, nBlocks, task_name)
    full_data_table = table();
    for b = 1:nBlocks
        block_filename = sprintf('sub%03d_%s_b%d.mat', subj_id, task_name, b);
        block_filepath = fullfile(results_dir, block_filename);
        if exist(block_filepath, 'file')
            loaded_data = load(block_filepath);
            if isfield(loaded_data, 'results_1_back')
                full_data_table = [full_data_table; loaded_data.results_1_back];
            elseif isfield(loaded_data, 'results_2_back')
                full_data_table = [full_data_table; loaded_data.results_2_back];
            else
                warning('Variable not found in file: %s.', block_filename);
            end
        else
            warning('Data file not found: %s.', block_filename);
        end
    end
end