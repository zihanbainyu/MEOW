clear; clc;

base_dir = '..';
data_dir = fullfile(base_dir, 'data');
setup_dir = fullfile(base_dir, 'subj_setup');
out_dir = fullfile(data_dir, 'rec_mod');

subject_list = [602 603 604 605 606 607];

fprintf('Processing %d subjects...\n\n', numel(subject_list));

for s = 1:numel(subject_list)

    subj_id = subject_list(s);
    subj_folder = fullfile(data_dir, sprintf('sub%03d', subj_id));
    
    rec_filename = fullfile(subj_folder, sprintf('sub%03d_concat.mat', subj_id));
    load(rec_filename, 'final_data_output');

    is_identity_A   = final_data_output.results_2_back_all.identity == "A";
    is_goal_AN = final_data_output.results_2_back_all.goal == "A-N";
    is_condition_ok = ismember(final_data_output.results_2_back_all.condition, ["compared", "isolated"]);
    
    % Combine them to find the indices
    mask_2b_target = is_identity_A & is_goal_AN & is_condition_ok;
    
    % EXTRACT the stim_ids. These are the specific IDs we want to find in 1-back.
    target_ids = final_data_output.results_2_back_all.stim_id(mask_2b_target);

    target_ids_B = replace(target_ids, "_A_", "_B_");

    % 2. Stack them in the same column (A list on top, B list on bottom)
    all_target_ids = [target_ids; target_ids_B];

    rows_to_remove = ismember(final_data_output.results_recognition.stim_id, all_target_ids);

    % 2. Overwrite the variable, keeping only the rows that are NOT (~) in that list
    test = final_data_output.results_recognition(~rows_to_remove, :);

    new_indices = find(test.trial_type == "new");

    % 2. Determine how many to remove (half of them)
    num_to_remove = 80;
    
    % 3. Randomly shuffle the indices and pick the first 'num_to_remove'
    shuffled_indices = new_indices(randperm(length(new_indices)));
    indices_to_delete = shuffled_indices(1:num_to_remove);
    
    % 4. Remove those specific rows from the table
    test(indices_to_delete, :) = [];
    final_data_output.results_recognition = test;

    out_filename = fullfile(out_dir, sprintf('sub%03d_concat_mod.mat', subj_id));

    save(out_filename, 'final_data_output');

end



fprintf('Done.\n');
