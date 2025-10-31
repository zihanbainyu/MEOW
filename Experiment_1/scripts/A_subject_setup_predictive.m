%==========================================================================
%         Generate Subject-Specific Task Variables
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
% note: you must run this for each subject before running the task
%
% MODIFIED: 2025-10-31 to implement proper 2-back test generator
%==========================================================================

clear;
clc;
rng('shuffle');

subj_id_str = input('please enter subject ID (e.g., 101): ', 's');
if isempty(subj_id_str)
    error('subject ID cannot be empty.');
end
subj_id = str2double(subj_id_str);
p.subj_id = subj_id;
encoding_schedule_all = table();
test_schedule_all = table();

% response counters
p.counts.encoding = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);
p.counts.test = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);

% directory
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_1028/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end

%% ========================================================================
%  SECTION 1: SUBJECT ID AND PARAMETERS
%  ========================================================================
% # of stimulus pairs per condition (total across all blocks)
p.nComparison = 120;
p.nIsolated_Both = 120;
p.nNovel = 120;
p.nTotalPairs = p.nComparison + p.nIsolated_Both + p.nNovel; % 360 total pairs

% # of blocks in experiment
p.nBlocks = 4; % 90 pairs per block

% keyboard mappings
p.keys.same = 'j';
p.keys.diff = 'k';
p.keys.quit = 'escape';

% timing parameters in seconds
p.timing.image_dur = 1.5;           % stimulus presentation
p.timing.fix_dur = 0.75;            % base fixation
p.timing.fix_jitter = 0.25;         % jitter range: ±0.25s (so 0.5 to 1.0s total)

%% ========================================================================
%  SECTION 2: LOAD STIMULI (L1 and L2 BALANCED)
%  ========================================================================

output_filename = fullfile(p.setup_dir, sprintf('sub%03d_setup.mat', subj_id));

if exist(output_filename, 'file')
    overwrite = input(sprintf('setup for subject %03d exists. overwrite? (y/n): ', subj_id), 's');
    if ~strcmpi(overwrite, 'y')
        fprintf('aborted.\n');
        return;
    end
end

fprintf('loading stimuli (L1 and L2 pools)...\n');

% --- Load L1 Pairs (Bin 1) ---
all_targ_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_targ_l1.png'));
all_lure_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_lure_l1.png'));
assert(numel(all_targ_files_l1) == 180, 'expected 180 target L1 files!');
targ_names_l1 = string({all_targ_files_l1.name}');
lure_names_l1 = string({all_lure_files_l1.name}');
master_pair_list_l1 = table(sort(targ_names_l1), sort(lure_names_l1), ...
    'VariableNames', {'Img_A', 'Img_B'});

% --- Load L2 Pairs (Bin 2) ---
all_targ_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_targ_l2.png'));
all_lure_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_lure_l2.png'));
assert(numel(all_targ_files_l2) == 180, 'expected 180 target L2 files!');
targ_names_l2 = string({all_targ_files_l2.name}');
lure_names_l2 = string({all_lure_files_l2.name}');
master_pair_list_l2 = table(sort(targ_names_l2), sort(lure_names_l2), ...
    'VariableNames', {'Img_A', 'Img_B'});

% Load foil pairs
all_foil_files = dir(fullfile(p.stim_dir, 'mst_*_foil.png'));
all_foil_lure_files = dir(fullfile(p.stim_dir, 'mst_*_foil_lure.png'));

all_foil_pairs = table('Size', [numel(all_foil_files), 2], ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'foil', 'foil_lure'});

all_foil_pairs.foil = string({all_foil_files.name}');
all_foil_pairs.foil_lure = string({all_foil_lure_files.name}');

% Shuffle them
all_foil_pairs = all_foil_pairs(randperm(height(all_foil_pairs)), :);

fprintf('Found %d L1 pairs, %d L2 pairs, and %d foils.\n', ...
    height(master_pair_list_l1), height(master_pair_list_l2), height(all_foil_pairs));

%% ========================================================================
%  SECTION 3: COUNTERBALANCE TARGET-LURE ASSIGNMENT (STRATIFIED)
%  ========================================================================

n_cond_l1 = height(master_pair_list_l1) / 3; % 180 / 3 = 60
n_cond_l2 = height(master_pair_list_l2) / 3; % 180 / 3 = 60

% --- Split L1 Pool (180 pairs -> 60/60/60) ---
fprintf('Splitting L1 pool...\n');
shuffled_list_l1 = master_pair_list_l1(randperm(height(master_pair_list_l1)), :);
final_list_l1 = table('Size', size(shuffled_list_l1), ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Target', 'Lure'});
% Randomly assign which image is target vs lure
for i_pair = 1:height(shuffled_list_l1)
    if rand() > 0.5
        final_list_l1(i_pair, :) = shuffled_list_l1(i_pair, :);
    else
        final_list_l1(i_pair, :) = [shuffled_list_l1(i_pair, 2), shuffled_list_l1(i_pair, 1)];
    end
end
% Split this pool into 3 equal chunks
comp_l1 = final_list_l1(1:n_cond_l1, :);
iso_l1  = final_list_l1(n_cond_l1+1 : 2*n_cond_l1, :);
nov_l1  = final_list_l1(2*n_cond_l1+1 : end, :);

% --- Split L2 Pool (180 pairs -> 60/60/60) ---
fprintf('Splitting L2 pool...\n');
shuffled_list_l2 = master_pair_list_l2(randperm(height(master_pair_list_l2)), :);
final_list_l2 = table('Size', size(shuffled_list_l2), ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Target', 'Lure'});
% Randomly assign which image is target vs lure
for i_pair = 1:height(shuffled_list_l2)
    if rand() > 0.5
        final_list_l2(i_pair, :) = shuffled_list_l2(i_pair, :);
    else
        final_list_l2(i_pair, :) = [shuffled_list_l2(i_pair, 2), shuffled_list_l2(i_pair, 1)];
    end
end
% Split this pool into 3 equal chunks
comp_l2 = final_list_l2(1:n_cond_l2, :);
iso_l2  = final_list_l2(n_cond_l2+1 : 2*n_cond_l2, :);
nov_l2  = final_list_l2(2*n_cond_l2+1 : end, :);

% --- Combine pools to create final 120-pair condition lists ---
comparison_pairs = [comp_l1; comp_l2];
iso_both_pairs = [iso_l1; iso_l2];
novel_pairs = [nov_l1; nov_l2];

% Shuffle the final lists so L1/L2 pairs are mixed
comparison_pairs = comparison_pairs(randperm(height(comparison_pairs)), :);
iso_both_pairs = iso_both_pairs(randperm(height(iso_both_pairs)), :);
novel_pairs = novel_pairs(randperm(height(novel_pairs)), :);

% store condition assignments
p.stim.comparison_l1 = comp_l1;
p.stim.comparison_l2 = comp_l2;
p.stim.iso_both_l1 = iso_l1;
p.stim.iso_both_l2 = iso_l2;
p.stim.novel_l1 = nov_l1;
p.stim.novel_l2 = nov_l2;

% Store the final combined lists
p.stim.comparison = comparison_pairs;
p.stim.iso_both = iso_both_pairs;
p.stim.novel = novel_pairs;

fprintf('Counterbalanced L1/L2 pairs across 3 conditions.\n');
fprintf('  -> Comparison: %d L1, %d L2 (Total %d)\n', height(comp_l1), height(comp_l2), height(comparison_pairs));
fprintf('  -> Isolated:   %d L1, %d L2 (Total %d)\n', height(iso_l1), height(iso_l2), height(iso_both_pairs));
fprintf('  -> Novel:      %d L1, %d L2 (Total %d)\n', height(nov_l1), height(nov_l2), height(novel_pairs));

%% ========================================================================
%  SECTION 4: PARTITION STIMULI INTO BLOCKS (STRATIFIED)
%  ========================================================================

partition_idx = @(N, nblocks) arrayfun(@(k) ...
    ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), ...
    1:nblocks, 'UniformOutput', false);

n_comp_l1 = height(p.stim.comparison_l1); % 60
n_iso_l1 = height(p.stim.iso_both_l1);   % 60
n_nov_l1 = height(p.stim.novel_l1);     % 60

n_comp_l2 = height(p.stim.comparison_l2); % 60
n_iso_l2 = height(p.stim.iso_both_l2);   % 60
n_nov_l2 = height(p.stim.novel_l2);     % 60

% --- Partition L1 pools into 4 blocks (15 pairs each) ---
p.block_indices.comp_l1 = partition_idx(n_comp_l1, p.nBlocks);
p.block_indices.iso_l1  = partition_idx(n_iso_l1, p.nBlocks);
p.block_indices.nov_l1  = partition_idx(n_nov_l1, p.nBlocks);

% --- Partition L2 pools into 4 blocks (15 pairs each) ---
p.block_indices.comp_l2 = partition_idx(n_comp_l2, p.nBlocks);
p.block_indices.iso_l2  = partition_idx(n_iso_l2, p.nBlocks);
p.block_indices.nov_l2  = partition_idx(n_nov_l2, p.nBlocks);

fprintf('Partitioned L1 and L2 pools into %d blocks.\n', p.nBlocks);
fprintf('  -> Each block gets: 15 L1-Comp, 15 L2-Comp, 15 L1-Iso, etc.\n');

%% ========================================================================
%  SECTION 5: BUILD EACH BLOCK
%  ========================================================================

% CRITICAL: Initialize persistent foil tracker BEFORE block loop
all_foils_remain = all_foil_pairs;

for b = 1:p.nBlocks
    fprintf('\n=== BLOCK %d ===\n', b);
   
    % =====================================================================
    % 5A: SELECT BLOCK-SPECIFIC PAIRS (STRATIFIED)
    % =====================================================================

    comp_l1_b = p.stim.comparison_l1(p.block_indices.comp_l1{b}, :);
    comp_l2_b = p.stim.comparison_l2(p.block_indices.comp_l2{b}, :);
    comparison_pairs_b = [comp_l1_b; comp_l2_b];
    comparison_pairs_b = comparison_pairs_b(randperm(height(comparison_pairs_b)), :);

    iso_l1_b = p.stim.iso_both_l1(p.block_indices.iso_l1{b}, :);
    iso_l2_b = p.stim.iso_both_l2(p.block_indices.iso_l2{b}, :);
    iso_both_pairs_b = [iso_l1_b; iso_l2_b];
    iso_both_pairs_b = iso_both_pairs_b(randperm(height(iso_both_pairs_b)), :);

    nov_l1_b = p.stim.novel_l1(p.block_indices.nov_l1{b}, :);
    nov_l2_b = p.stim.novel_l2(p.block_indices.nov_l2{b}, :);
    novel_pairs_b = [nov_l1_b; nov_l2_b];
    novel_pairs_b = novel_pairs_b(randperm(height(novel_pairs_b)), :);

    nComp_b = height(comparison_pairs_b);
    nIso_b = height(iso_both_pairs_b);
    nNov_b = height(novel_pairs_b);

    fprintf('  Loaded %d Comp, %d Iso, %d Nov pairs\n', nComp_b, nIso_b, nNov_b);

    % =====================================================================
    % 5B: BUILD ENCODING PHASE
    % =====================================================================

    comp_miniblocks = {};
    repeat_miniblocks = {};
    iso_single_trials = {};

    % Bucket 1: comparison (C-C')
    for i = 1:nComp_b
        comp_miniblocks{end+1} = { ...
            comparison_pairs_b.Target(i), "comparison", "target", "new", "none"; ...
            comparison_pairs_b.Lure(i), "comparison", "lure", "lure", "k"};
    end

    % Bucket 2: foil repeats (R-R)
    n_foil_repeats = nComp_b;
    
    if height(all_foils_remain) < n_foil_repeats
        error('Not enough foils for encoding in Block %d', b);
    end
    repeat_foil_pairs_b = all_foils_remain(1:n_foil_repeats, :);
    all_foils_remain(1:n_foil_repeats, :) = []; % Remove used rows

    for i = 1:n_foil_repeats
        foil_item = repeat_foil_pairs_b.foil(i);
        repeat_miniblocks{end+1} = { ...
            foil_item, "foil_repeat", "filler", "new", "none"; ...
            foil_item, "foil_repeat", "filler", "repeat", "j"};
    end

    % Bucket 3: isolated items (I and I')
    for i = 1:nIso_b
        iso_single_trials{end+1} = { ...
            iso_both_pairs_b.Target(i), "iso_both", "target", "new", "none"};
        iso_single_trials{end+1} = { ...
            iso_both_pairs_b.Lure(i), "iso_both", "lure", "new", "none"};
    end

    % Combine and shuffle
    all_encoding_blocks = [comp_miniblocks, repeat_miniblocks, iso_single_trials];
    shuffled_indices = randperm(numel(all_encoding_blocks));
    final_ordered_blocks = all_encoding_blocks(shuffled_indices);
    final_encoding_list = vertcat(final_ordered_blocks{:});

    % Convert to table
    encoding_schedule_block = cell2table(final_encoding_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});

    encoding_schedule_block.block = repmat(b, height(encoding_schedule_block), 1);
    encoding_schedule_block.nback_target_id = strings(height(encoding_schedule_block), 1);
    encoding_schedule_block.trial_type_final = encoding_schedule_block.trial_type_designed;

    % Check for 1-back repeats
    for i = 2:height(encoding_schedule_block)
        encoding_schedule_block.nback_target_id(i) = ...
            encoding_schedule_block.stimulus_id(i-1);

        if string(encoding_schedule_block.stimulus_id(i)) == ...
           string(encoding_schedule_block.nback_target_id(i))
            encoding_schedule_block.trial_type_final(i) = "repeat";
            encoding_schedule_block.correct_response(i) = "j";
        elseif encoding_schedule_block.trial_type_designed(i) == "lure"
            encoding_schedule_block.trial_type_final(i) = "lure";
            encoding_schedule_block.correct_response(i) = "k";
        else
            encoding_schedule_block.trial_type_final(i) = "new";
            encoding_schedule_block.correct_response(i) = "none";
        end
    end

    % Cleanup columns
    encoding_schedule_block.item_type = encoding_schedule_block.trial_type_final;
    encoding_schedule_block = removevars(encoding_schedule_block, ...
        {'trial_type_designed', 'nback_target_id', 'trial_type_final'});

    % =====================================================================
    % 5C: ENSURE TARGET ENCODED BEFORE LURE (ISOLATED CONDITION)
    % =====================================================================

    for i = 1:height(iso_both_pairs_b)
        tgt = iso_both_pairs_b.Target(i);
        lur = iso_both_pairs_b.Lure(i);

        idx_tgt_trial = find(encoding_schedule_block.stimulus_id == tgt, 1);
        idx_lur_trial = find(encoding_schedule_block.stimulus_id == lur, 1);

        if ~isempty(idx_tgt_trial) && ~isempty(idx_lur_trial) && idx_lur_trial < idx_tgt_trial
            tmp = encoding_schedule_block(idx_tgt_trial, :);
            encoding_schedule_block(idx_tgt_trial, :) = encoding_schedule_block(idx_lur_trial, :);
            encoding_schedule_block(idx_lur_trial, :) = tmp;
        end
    end

    % Append to full schedule
    encoding_schedule_all = [encoding_schedule_all; encoding_schedule_block];

    % Update counts
    p.counts.encoding.j_presses = p.counts.encoding.j_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'j'));
    p.counts.encoding.k_presses = p.counts.encoding.k_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'k'));
    p.counts.encoding.no_response = p.counts.encoding.no_response + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'none'));

    % =====================================================================
    % 5D: BUILD TEST PHASE (2-BACK WITH OVERLAPPING GOALS)
    % =====================================================================

    fprintf('  Building 2-back test sequence...\n');

    % --- 5D-1: Prepare Goal List ---

    block_pairs = [comparison_pairs_b; iso_both_pairs_b; novel_pairs_b];
    block_conditions = [repmat("comparison", nComp_b, 1); ...
                        repmat("iso_both", nIso_b, 1); ...
                        repmat("novel", nNov_b, 1)];

    goal_list = table('Size', [90, 5], ...
        'VariableTypes', {'string', 'string', 'string', 'string', 'string'}, ...
        'VariableNames', {'O_item', 'O_lure', 'condition', 'goal_type', 'X_item'});

    goal_list.O_item = block_pairs.Target;
    goal_list.O_lure = block_pairs.Lure;
    goal_list.condition = block_conditions;

    % Assign goal types (1/3 each)
    idx_shuf = randperm(90);
    goal_list.goal_type(idx_shuf(1:30)) = "lure";
    goal_list.goal_type(idx_shuf(31:60)) = "same";
    goal_list.goal_type(idx_shuf(61:90)) = "new";

    % Assign X_item for same/lure (deterministic)
    for i = 1:90
        if goal_list.goal_type(i) == "lure"
            goal_list.X_item(i) = goal_list.O_lure(i);
        elseif goal_list.goal_type(i) == "same"
            goal_list.X_item(i) = goal_list.O_item(i);
        else
            goal_list.X_item(i) = ""; % Assigned dynamically
        end
    end

    % Shuffle goal order
    goal_list = goal_list(randperm(height(goal_list)), :);

    % --- 5D-2: Build Sequence ---

    sequence = cell(300, 5);
    row_idx = 1;

    % Init with 2 junk
    if height(all_foils_remain) < 2
        error('Need foils for test init in Block %d', b);
    end
    junk1 = all_foils_remain(1,:); all_foils_remain(1,:) = [];
    junk2 = all_foils_remain(1,:); all_foils_remain(1,:) = [];

    sequence(row_idx,:) = {junk1.foil, "init_junk", "filler", "new", "none"};
    row_idx = row_idx + 1;
    sequence(row_idx,:) = {junk2.foil, "init_junk", "filler", "new", "none"};
    row_idx = row_idx + 1;

    active_goals = [];
    goals_started = false(height(goal_list), 1);
    goal_pointer = 1;

    while goal_pointer <= height(goal_list) || ~isempty(active_goals)

        % Check if any goal needs completion (at N-2)
        goals_to_complete = [];
        for k = 1:size(active_goals, 1)
            goal_idx = active_goals(k, 1);
            goal_pos = active_goals(k, 2);

            if row_idx - goal_pos == 2
                goals_to_complete(end+1) = k;
            end
        end

        if ~isempty(goals_to_complete)
            % Complete first ready goal
            k = goals_to_complete(1);
            goal_idx = active_goals(k, 1);
            goal = goal_list(goal_idx, :);

            X_type = goal.goal_type;

            if strcmp(X_type, "new")
                % Dynamically select next available goal's O_item
                active_goal_indices = active_goals(:, 1);

                next_unstarted = [];
                for check_idx = 1:height(goal_list)
                    if ~goals_started(check_idx) && ...
                       ~ismember(check_idx, active_goal_indices)
                        next_unstarted = check_idx;
                        break;
                    end
                end

                if ~isempty(next_unstarted)
                    X_item = goal_list.O_item(next_unstarted);
                else
                    X_item = goal_list.O_lure(randi(height(goal_list)));
                end
                X_resp = "none";
                X_role = "target";
            else
                % Predetermined for same/lure
                X_item = goal.X_item;
                if strcmp(X_type, "lure")
                    X_resp = "k";
                    X_role = "lure";
                else
                    X_resp = "j";
                    X_role = "target";
                end
            end

            sequence(row_idx,:) = {X_item, goal.condition, X_role, X_type, X_resp};
            row_idx = row_idx + 1;

            active_goals(k, :) = [];

            % If "new" goal, X_item can start its own goal
            if strcmp(X_type, "new")
                X_goal_idx = find(strcmp(goal_list.O_item, X_item), 1);
                if ~isempty(X_goal_idx) && ~goals_started(X_goal_idx)
                    goals_started(X_goal_idx) = true;
                    active_goals(end+1, :) = [X_goal_idx, row_idx - 1];
                end
            end

        else
            % Start new goal
            while goal_pointer <= height(goal_list) && goals_started(goal_pointer)
                goal_pointer = goal_pointer + 1;
            end

            if goal_pointer <= height(goal_list)
                goal = goal_list(goal_pointer, :);

                sequence(row_idx,:) = {goal.O_item, goal.condition, "target", "original", "none"};

                goals_started(goal_pointer) = true;
                active_goals(end+1, :) = [goal_pointer, row_idx];

                row_idx = row_idx + 1;
                goal_pointer = goal_pointer + 1;
            else
                break;
            end
        end
    end

    % End junk
    if height(all_foils_remain) >= 5
        for jj = 1:5
            sequence(row_idx,:) = {all_foils_remain.foil(jj), "end_junk", "filler", "new", "none"};
            row_idx = row_idx + 1;
        end
        all_foils_remain(1:5,:) = [];
    end

    % Convert to table
    test_schedule_block = cell2table(sequence(1:row_idx-1, :), ...
        'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type', 'correct_response'});
    
    test_schedule_block.block = repmat(b, height(test_schedule_block), 1);

    % Append to full schedule
    test_schedule_all = [test_schedule_all; test_schedule_block];

    % Update counts
    p.counts.test.j_presses = p.counts.test.j_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'k'));
    p.counts.test.no_response = p.counts.test.no_response + ...
        sum(strcmp(test_schedule_block.correct_response, 'none'));

    fprintf('  Generated %d encoding trials, %d test trials\n', ...
            height(encoding_schedule_block), height(test_schedule_block));

    % =====================================================================
    % SANITY CHECK FOR THIS BLOCK
    % =====================================================================

    fprintf('\n  === SANITY CHECK (Block %d) ===\n', b);

    % Check response distribution
    j_resp = sum(strcmp(test_schedule_block.correct_response, 'j'));
    k_resp = sum(strcmp(test_schedule_block.correct_response, 'k'));
    none_resp = sum(strcmp(test_schedule_block.correct_response, 'none'));
    fprintf('    Responses: j=%d, k=%d, none=%d\n', j_resp, k_resp, none_resp);

    % Check unique O-items
    targets_original = test_schedule_block(strcmp(test_schedule_block.role, "target") & ...
                                           strcmp(test_schedule_block.trial_type, "original"), :);
    unique_targets = unique(targets_original.stimulus_id);
    fprintf('    Unique O-items started: %d (expected: 90)\n', length(unique_targets));

    % Check for duplicates
    [unique_vals, ~, idx] = unique(targets_original.stimulus_id);
    counts = accumarray(idx, 1);
    duplicates = unique_vals(counts > 1);
    if isempty(duplicates)
        fprintf('    ✓ No duplicate goal starts\n');
    else
        fprintf('    ✗ Found %d duplicates!\n', length(duplicates));
    end

    % Verify 2-back structure
    errors = 0;
    for i = 1:height(goal_list)
        O_item = goal_list.O_item(i);
        goal_type = goal_list.goal_type(i);

        O_idx = find(strcmp(test_schedule_block.stimulus_id, O_item) & ...
                     strcmp(test_schedule_block.role, "target") & ...
                     strcmp(test_schedule_block.trial_type, "original"), 1);

        if ~isempty(O_idx)
            expected_X_idx = O_idx + 2;
            if expected_X_idx <= height(test_schedule_block)
                actual_X_type = test_schedule_block.trial_type(expected_X_idx);
                if ~strcmp(actual_X_type, goal_type)
                    errors = errors + 1;
                end
            end
        end
    end
    
    if errors == 0
        fprintf('    ✓ All goals have correct 2-back structure\n');
    else
        fprintf('    ✗ Found %d 2-back structure errors\n', errors);
    end
    
    fprintf('  === END SANITY CHECK ===\n\n');
end

%% ========================================================================
%  SECTION 6: SAVE OUTPUT
%  ========================================================================

p.encoding_schedule = encoding_schedule_all;
p.test_schedule = test_schedule_all;

fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Total encoding trials: %d\n', height(encoding_schedule_all));
fprintf('  j presses: %d, k presses: %d, no response: %d\n', ...
        p.counts.encoding.j_presses, p.counts.encoding.k_presses, p.counts.encoding.no_response);
fprintf('Total test trials: %d\n', height(test_schedule_all));
fprintf('  j presses: %d, k presses: %d, no response: %d\n', ...
        p.counts.test.j_presses, p.counts.test.k_presses, p.counts.test.no_response);

save(output_filename, 'p');
fprintf('\nSetup saved to: %s\n', output_filename);
fprintf('Done!\n');