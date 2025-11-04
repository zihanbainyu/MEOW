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
sequence_1_back = table();
sequence_2_back = table();

% response counters
p.counts.oneback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);
p.counts.twoback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);

% directory
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end

%% P1: Setup
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

%% P2: Load stimuli
output_filename = fullfile(p.setup_dir, sprintf('sub%03d_setup.mat', subj_id));

if exist(output_filename, 'file')
    overwrite = input(sprintf('setup for subject %03d exists. overwrite? (y/n): ', subj_id), 's');
    if ~strcmpi(overwrite, 'y')
        fprintf('aborted.\n');
        return;
    end
end

fprintf('loading stimuli...\n');

% --- Load L1 Pairs (Bin 1) ---
all_A_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_A_l1.png'));
all_B_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_B_l1.png'));
A_names_l1 = string({all_A_files_l1.name}');
B_names_l1 = string({all_B_files_l1.name}');
master_pair_list_l1 = table(sort(A_names_l1), sort(B_names_l1), ...
    'VariableNames', {'A', 'B'});

% --- Load L2 Pairs (Bin 2) ---
all_A_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_A_l2.png'));
all_B_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_B_l2.png'));
A_names_l2 = string({all_A_files_l2.name}');
B_names_l2 = string({all_B_files_l2.name}');
master_pair_list_l2 = table(sort(A_names_l2), sort(B_names_l2), ...
    'VariableNames', {'A', 'B'});

% Load foil pairs
all_foil_A_files = dir(fullfile(p.stim_dir, 'mst_*_A_foil.png'));
all_foil_B_files = dir(fullfile(p.stim_dir, 'mst_*_B_foil.png'));

all_foil_pairs = table('Size', [numel(all_foil_A_files), 2], ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'A_foil', 'B_foil'});

all_foil_pairs.A_foil = string({all_foil_A_files.name}');
all_foil_pairs.B_foil = string({all_foil_B_files.name}');

% Shuffle them
all_foil_pairs = all_foil_pairs(randperm(height(all_foil_pairs)), :);

fprintf('Found %d L1 pairs, %d L2 pairs, and %d foil pairs.\n', ...
    height(master_pair_list_l1), height(master_pair_list_l2), height(all_foil_pairs));

%% P3: Stimuli assignment
n_cond_l1 = height(master_pair_list_l1) / 3; % 180 / 3 = 60
n_cond_l2 = height(master_pair_list_l2) / 3; % 180 / 3 = 60

% --- Split L1 Pool (180 pairs -> 60/60/60) ---
final_list_l1 = master_pair_list_l1(randperm(height(master_pair_list_l1)), :);

% Split this pool into 3 equal chunks
comp_l1 = final_list_l1(1:n_cond_l1, :);
iso_l1  = final_list_l1(n_cond_l1+1 : 2*n_cond_l1, :);
nov_l1  = final_list_l1(2*n_cond_l1+1 : end, :);

% --- Split L2 Pool (180 pairs -> 60/60/60) ---
final_list_l2 = master_pair_list_l2(randperm(height(master_pair_list_l2)), :);

% Split this pool into 3 equal chunks
comp_l2 = final_list_l2(1:n_cond_l2, :);
iso_l2  = final_list_l2(n_cond_l2+1 : 2*n_cond_l2, :);
nov_l2  = final_list_l2(2*n_cond_l2+1 : end, :);

% --- Combine pools to create final 120-pair condition lists ---
comp_pairs = [comp_l1; comp_l2];
iso_pairs = [iso_l1; iso_l2];
novel_pairs = [nov_l1; nov_l2];

% Shuffle the final lists so L1/L2 pairs are mixed
comp_pairs = comp_pairs(randperm(height(comp_pairs)), :);
iso_pairs = iso_pairs(randperm(height(iso_pairs)), :);
novel_pairs = novel_pairs(randperm(height(novel_pairs)), :);

% store condition assignments
p.stim.comp_l1 = comp_l1;
p.stim.comp_l2 = comp_l2;
p.stim.iso_l1 = iso_l1;
p.stim.iso_l2 = iso_l2;
p.stim.novel_l1 = nov_l1;
p.stim.novel_l2 = nov_l2;

% Store the final combined lists
p.stim.compared = comp_pairs;
p.stim.isolated = iso_pairs;
p.stim.novel = novel_pairs;

fprintf('  -> Compared: %d L1, %d L2 (Total %d)\n', height(comp_l1), height(comp_l2), height(comp_pairs));
fprintf('  -> Isolated:   %d L1, %d L2 (Total %d)\n', height(iso_l1), height(iso_l2), height(iso_pairs));
fprintf('  -> Novel:      %d L1, %d L2 (Total %d)\n', height(nov_l1), height(nov_l2), height(novel_pairs));

%% P4: Split stimuli into blocks
partition_idx = @(N, nblocks) arrayfun(@(k) ...
    ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), ...
    1:nblocks, 'UniformOutput', false);

n_comp_l1 = height(p.stim.comp_l1); % 60
n_iso_l1 = height(p.stim.iso_l1);   % 60
n_nov_l1 = height(p.stim.novel_l1);     % 60

n_comp_l2 = height(p.stim.comp_l2); % 60
n_iso_l2 = height(p.stim.iso_l2);   % 60
n_nov_l2 = height(p.stim.novel_l2);     % 60

% --- Partition L1 pools into 4 blocks (15 pairs each) ---
p.block_indices.comp_l1 = partition_idx(n_comp_l1, p.nBlocks);
p.block_indices.iso_l1  = partition_idx(n_iso_l1, p.nBlocks);
p.block_indices.nov_l1  = partition_idx(n_nov_l1, p.nBlocks);

% --- Partition L2 pools into 4 blocks (15 pairs each) ---
p.block_indices.comp_l2 = partition_idx(n_comp_l2, p.nBlocks);
p.block_indices.iso_l2  = partition_idx(n_iso_l2, p.nBlocks);
p.block_indices.nov_l2  = partition_idx(n_nov_l2, p.nBlocks);

fprintf('  -> Each block gets: 15 L1-Comp, 15 L2-Comp, 15 L1-Iso, etc.\n');

%% P5: Build sequence
all_foils_remain = all_foil_pairs;

for b = 1:p.nBlocks
    fprintf('\nBlock %d \n', b);
   
    %%% Assign block-specific stimuli
    comp_l1_b = p.stim.comp_l1(p.block_indices.comp_l1{b}, :);
    comp_l2_b = p.stim.comp_l2(p.block_indices.comp_l2{b}, :);
    comp_pairs_b = [comp_l1_b; comp_l2_b];
    comp_pairs_b = comp_pairs_b(randperm(height(comp_pairs_b)), :);

    iso_l1_b = p.stim.iso_l1(p.block_indices.iso_l1{b}, :);
    iso_l2_b = p.stim.iso_l2(p.block_indices.iso_l2{b}, :);
    iso_pairs_b = [iso_l1_b; iso_l2_b];
    iso_pairs_b = iso_pairs_b(randperm(height(iso_pairs_b)), :);

    nov_l1_b = p.stim.novel_l1(p.block_indices.nov_l1{b}, :);
    nov_l2_b = p.stim.novel_l2(p.block_indices.nov_l2{b}, :);
    novel_pairs_b = [nov_l1_b; nov_l2_b];
    novel_pairs_b = novel_pairs_b(randperm(height(novel_pairs_b)), :);

    nComp_b = height(comp_pairs_b);
    nIso_b = height(iso_pairs_b);
    nNov_b = height(novel_pairs_b);

    %% P5A: Build sequence 1-back
    comp_miniblocks = {};
    repeat_miniblocks = {};
    iso_trials = {};

    %%%% compared (C-C')
    for i = 1:nComp_b
        comp_miniblocks{end+1} = { ...
            comp_pairs_b.A(i), "compared", "A", "none"; ...
            comp_pairs_b.B(i), "compared", "B", "k"};
    end

    %%%% foil repeats (R-R)
    n_foil_repeats = nComp_b;
    
    if height(all_foils_remain) < n_foil_repeats
        error('Not enough foils for 1-back in Block %d', b);
    end
    repeat_foil_pairs = all_foils_remain(1:n_foil_repeats, :);
    all_foils_remain(1:n_foil_repeats, :) = [];

    for i = 1:n_foil_repeats
        foil_item = repeat_foil_pairs.A_foil(i);
        repeat_miniblocks{end+1} = { ...
            foil_item, "repeat", "A", "none"; ...
            foil_item, "repeat", "A", "j"};
    end

    %%%% isolated items (I-I')
    for i = 1:nIso_b
        iso_trials{end+1} = { ...
            iso_pairs_b.A(i), "isolated", "A", "none"};
        iso_trials{end+1} = { ...
            iso_pairs_b.B(i), "isolated", "B", "none"};
    end

    % Combine and shuffle
    miniblocks_1_back = [comp_miniblocks, repeat_miniblocks, iso_trials];
    shuffled_indices = randperm(numel(miniblocks_1_back));
    final_miniblocks = miniblocks_1_back(shuffled_indices);
    final_1_back_list = vertcat(final_miniblocks{:});

    % Convert to table
    sequence_1_back_block = cell2table(final_1_back_list, ...
        'VariableNames', {'stim_id', 'condition', 'identity', 'corr_resp'});

    sequence_1_back_block.block = repmat(b, height(sequence_1_back_block), 1);

    %%%% Ensure A appears earlier than B for Isolated item
    for i = 1:height(iso_pairs_b)
        tgt = iso_pairs_b.A(i);
        lur = iso_pairs_b.B(i);

        idx_tgt_trial = find(sequence_1_back_block.stim_id == tgt, 1);
        idx_lur_trial = find(sequence_1_back_block.stim_id == lur, 1);

        if ~isempty(idx_tgt_trial) && ~isempty(idx_lur_trial) && idx_lur_trial < idx_tgt_trial
            tmp = sequence_1_back_block(idx_tgt_trial, :);
            sequence_1_back_block(idx_tgt_trial, :) = sequence_1_back_block(idx_lur_trial, :);
            sequence_1_back_block(idx_lur_trial, :) = tmp;
        end
    end

    % Append to full schedule
    sequence_1_back = [sequence_1_back; sequence_1_back_block];

    % Update counts
    p.counts.oneback.resp_same = p.counts.oneback.resp_same + ...
        sum(strcmp(sequence_1_back_block.corr_resp, 'j'));
    p.counts.oneback.resp_similar = p.counts.oneback.resp_similar + ...
        sum(strcmp(sequence_1_back_block.corr_resp, 'k'));
    p.counts.oneback.resp_new = p.counts.oneback.resp_new + ...
        sum(strcmp(sequence_1_back_block.corr_resp, 'none'));

    %% P5B: Build sequence 2-back
    fprintf('  Building 2-back sequence...\n');

    %%% Create goal list
    block_pairs = [comp_pairs_b; iso_pairs_b; novel_pairs_b];
    block_conditions = [repmat("compared", nComp_b, 1); ...
                        repmat("isolated", nIso_b, 1); ...
                        repmat("novel", nNov_b, 1)];
    goal_list = table('Size', [90, 5], ...
        'VariableTypes', {'string', 'string', 'string', 'string', 'string'}, ...
        'VariableNames', {'A', 'B', 'condition', 'goal_type', 'X'});
    goal_list.A = block_pairs.A;
    goal_list.B = block_pairs.B;
    goal_list.condition = block_conditions;
    % assign goal types (1/3 each)
    idx_shuf = randperm(90);
    goal_list.goal_type(idx_shuf(1:30)) = "A-B";
    goal_list.goal_type(idx_shuf(31:60)) = "A-A";
    goal_list.goal_type(idx_shuf(61:90)) = "A-N";
    % assign X for A-A and A-B
    for i = 1:90
        if goal_list.goal_type(i) == "A-B"
            goal_list.X(i) = goal_list.B(i);
        elseif goal_list.goal_type(i) == "A-A"
            goal_list.X(i) = goal_list.A(i);
        else
            goal_list.X(i) = ""; % assigned dynamically
        end
    end
    % shuffle goal order
    goal_list = goal_list(randperm(height(goal_list)), :);
    
    %%% Build sequence
    sequence = cell(300, 5); % Increased size to prevent crash
    row_idx = 1;
    % init with 2 junk
    junk1 = all_foils_remain(1,:); all_foils_remain(1,:) = [];
    junk2 = all_foils_remain(1,:); all_foils_remain(1,:) = [];
    sequence(row_idx,:) = {junk1.A_foil, "init_junk", "J", "JUNK", "none"};
    row_idx = row_idx + 1;
    sequence(row_idx,:) = {junk2.A_foil, "init_junk", "J", "JUNK", "none"};
    row_idx = row_idx + 1;
    active_goals = [];
    goals_started = false(height(goal_list), 1);
    goal_pointer = 1;
    
    while goal_pointer <= height(goal_list) || ~isempty(active_goals)
        % check if any goal needs completion (at N-2)
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
            goal = goal_list(goal_idx, :); % This is the N-2 goal
            X_type = goal.goal_type;      % This is the N-2 goal_type (e.g., "A-N")
            
            if strcmp(X_type, "A-N")
                % A-N GOAL: 'X' is a NEW 'A' item
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
                    X = goal_list.A(next_unstarted); % This is the new 'A' item
                else
                    X = all_foils_remain.A_foil(1); % Grab a foil
                    all_foils_remain(1,:) = [];     % Consume the foil
                end
                
                % --- THIS IS THE FUCKING FIX ---
                % Find the properties of the NEW item 'X'
                X_props = goal_list(strcmp(goal_list.A, X), :);
                if isempty(X_props)
                    % This is a failsafe, e.g., we grabbed a 'B' item
                    log_condition = "junk_foil";
                    log_goal_type = "JUNK"; 
                    X_identity = "J";
                else
                    % Log the NEW item's OWN properties
                    log_condition = X_props.condition(1);
                    log_goal_type = X_props.goal_type(1); % This logs X's OWN goal
                    X_identity = "A";
                end
                X_resp = "none";
                % --- END FIX ---
                
            else
                % A-A or A-B GOAL: 'X' is the 'A' or 'B' from the N-2 pair
                X = goal.X;
                log_condition = goal.condition; % Log N-2's condition
                log_goal_type = X_type;       % Log N-2's goal (it's consistent)
                
                if strcmp(X_type, "A-B")
                    X_resp = "k";
                    X_identity = "B";
                else
                    X_resp = "j";
                    X_identity = "A";
                end
            end
            
            sequence(row_idx,:) = {X, log_condition, X_identity, log_goal_type, X_resp};
            row_idx = row_idx + 1;
            active_goals(k, :) = [];
            
            % if "new" goal, X_item can start its own goal
            if strcmp(X_type, "A-N")
                X_goal_idx = find(strcmp(goal_list.A, X), 1);
                if ~isempty(X_goal_idx) && ~goals_started(X_goal_idx)
                    goals_started(X_goal_idx) = true;
                    active_goals(end+1, :) = [X_goal_idx, row_idx - 1];
                end
            end
        else
            % start new goal
            while goal_pointer <= height(goal_list) && goals_started(goal_pointer)
                goal_pointer = goal_pointer + 1;
            end
            
            if goal_pointer <= height(goal_list)
                goal = goal_list(goal_pointer, :);
                
                % --- THIS IS THE FIX FOR THE UNDEFINED VARIABLE ---
                sequence(row_idx,:) = {goal.A, goal.condition, "A", goal.goal_type, "none"};
                
                goals_started(goal_pointer) = true;
                active_goals(end+1, :) = [goal_pointer, row_idx];
                row_idx = row_idx + 1;
                goal_pointer = goal_pointer + 1;
            else
                break; % All goals started
            end
        end
    end
    
    % end junk
    if height(all_foils_remain) >= 5
        for jj = 1:5
            % --- THIS IS THE FIX FOR THE JUNK LOOP ---
            sequence(row_idx,:) = {all_foils_remain.A_foil(jj), "end_junk", "J", "JUNK", "none"};
            row_idx = row_idx + 1;
        end
        all_foils_remain(1:5,:) = [];
    end
    
    % Convert to table
    sequence_2_back_block = cell2table(sequence(1:row_idx-1, :), ...
        'VariableNames', {'stim_id', 'condition', 'identity', 'goal','corr_resp'});
    
    sequence_2_back_block.block = repmat(b, height(sequence_2_back_block), 1);

    % Append to full schedule
    sequence_2_back = [sequence_2_back; sequence_2_back_block];

    % Update counts
    p.counts.twoback.resp_same = p.counts.twoback.resp_same + ...
        sum(strcmp(sequence_2_back_block.corr_resp, 'j'));
    p.counts.twoback.resp_similar = p.counts.twoback.resp_similar + ...
        sum(strcmp(sequence_2_back_block.corr_resp, 'k'));
    p.counts.twoback.resp_new = p.counts.twoback.resp_new + ...
        sum(strcmp(sequence_2_back_block.corr_resp, 'none'));

    fprintf('  Generated %d 1-back trials, %d 2-back trials\n', ...
            height(sequence_1_back_block), height(sequence_2_back_block));
end


%% P6: Build recognition task
fprintf('\nBuilding final recognition task...\n');

% --- 1. Get 'Old' items (N=240) from 'compared' and 'isolated' only ---
% For each of the 240 pairs, randomly select A or B
n_old_items = p.nComparison + p.nIsolated_Both;

all_study_pairs = [p.stim.compared; p.stim.isolated];
all_study_conds = [repmat("compared", p.nComparison, 1); ...
                   repmat("isolated", p.nIsolated_Both, 1)];

selected_old_items = strings(n_old_items, 1);
selected_identity = strings(n_old_items, 1);

fprintf('  Sampling 240 old items (random A/B) from compared/isolated...\n');
for i = 1:n_old_items
    if rand() > 0.5
        % Select the 'A' item
        selected_old_items(i) = all_study_pairs.A(i);
        selected_identity(i) = "A";
    else
        % Select the 'B' item
        selected_old_items(i) = all_study_pairs.B(i);
        selected_identity(i) = "B";
    end
end

% Build the final 'Old' items table
all_old_items = table(selected_old_items, all_study_conds, selected_identity, ...
    'VariableNames', {'stim_id', 'condition', 'identity'});
all_old_items.trial_type = repmat("old", n_old_items, 1);
all_old_items.correct_response = repmat(p.keys.same, n_old_items, 1); 

% --- 2. Get 240 'New' Foil items ---
n_rec_foils = n_old_items; % 240
assert(height(all_foils_remain) >= n_rec_foils, ...
    'Not enough remaining foils for recognition task! Need %d, have %d', ...
    n_rec_foils, height(all_foils_remain));
    
% Grab 240 'A_foil' images 
new_foils_list = all_foils_remain.A_foil(1:n_rec_foils);
all_foils_remain(1:n_rec_foils, :) = []; % Remove them

all_new_items = table(new_foils_list, repmat("foil", n_rec_foils, 1), ...
    repmat("N/A", n_rec_foils, 1), repmat("new", n_rec_foils, 1), ...
    repmat(p.keys.diff, n_rec_foils, 1), ... 
    'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'correct_response'});
    
% --- 3. Combine, shuffle, and add to 'p' struct ---
sequence_recognition = [all_old_items; all_new_items];
sequence_recognition = sequence_recognition(randperm(height(sequence_recognition)), :);

fprintf('  Generated %d recognition trials (240 old, 240 new).\n', height(sequence_recognition));

% Add to 'p' struct
p.sequence_recognition = sequence_recognition;

%% ========================================================================
%  P7: ADD JITTER & SAVE OUTPUT
%  ========================================================================
fprintf('\nAdding jittered fixations and saving all schedules...\n');

% --- Add jittered fixation durations ---
n_enc_trials = height(sequence_1_back);
n_test_trials = height(sequence_2_back);
n_rec_trials = height(p.sequence_recognition);

sequence_1_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_enc_trials, 1) * 2 - 1) * p.timing.fix_jitter;

sequence_2_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_test_trials, 1) * 2 - 1) * p.timing.fix_jitter;
    
p.sequence_recognition.fix_duration = ...
    p.timing.fix_dur + (rand(n_rec_trials, 1) * 2 - 1) * p.timing.fix_jitter;

% --- Add subj_id to all schedules ---
sequence_1_back.subj_id = repmat(subj_id, n_enc_trials, 1);
sequence_2_back.subj_id = repmat(subj_id, n_test_trials, 1);
p.sequence_recognition.subj_id = repmat(subj_id, n_rec_trials, 1);

% --- Save all schedules to 'p' struct ---
p.sequence_1_back = sequence_1_back;
p.sequence_2_back = sequence_2_back;
p.stim.all_foils_remaining = all_foils_remain; % Save the final unused stack

% --- Final Summary ---
fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Total 1-back trials: %d\n', n_enc_trials);
fprintf('  j presses: %d, k presses: %d, no response: %d\n', ...
        p.counts.oneback.resp_same, p.counts.oneback.resp_similar, p.counts.oneback.resp_new);
fprintf('Total 2-back trials: %d\n', n_test_trials);
fprintf('  j presses: %d, k presses: %d, no response: %d\n', ...
        p.counts.twoback.resp_same, p.counts.twoback.resp_similar, p.counts.twoback.resp_new);
fprintf('Total recognition trials: %d\n', n_rec_trials);
fprintf('  Old: %d, New: %d\n', sum(p.sequence_recognition.trial_type == "old"), ...
        sum(p.sequence_recognition.trial_type == "new"));
fprintf('Foils remaining in stack: %d\n', height(all_foils_remain));

% --- Save the file ---
save(output_filename, 'p');
fprintf('\nSetup saved to: %s\n', output_filename);
