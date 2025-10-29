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
% NOTE: We need 120 pairs *per condition* to balance 40 lure, 40 same, 40 new
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

% Load foil pairs (Section 2)
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
    height(master_pair_list_l1), height(master_pair_list_l2), numel(all_foil_pairs));


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
% Each list now has 60 L1 pairs and 60 L2 pairs, perfectly balanced.
comparison_pairs = [comp_l1; comp_l2];
iso_both_pairs = [iso_l1; iso_l2];
novel_pairs = [nov_l1; nov_l2];

% Shuffle the final lists so L1/L2 pairs are mixed
comparison_pairs = comparison_pairs(randperm(height(comparison_pairs)), :);
iso_both_pairs = iso_both_pairs(randperm(height(iso_both_pairs)), :);
novel_pairs = novel_pairs(randperm(height(novel_pairs)), :);

% store condition assignments
% We also store the separated pools for block partitioning
p.stim.comparison_l1 = comp_l1;
p.stim.comparison_l2 = comp_l2;
p.stim.iso_both_l1 = iso_l1;
p.stim.iso_both_l2 = iso_l2;
p.stim.novel_l1 = nov_l1;
p.stim.novel_l2 = nov_l2;

% Store the final combined lists (for encoding script, etc.)
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
% We partition L1 and L2 pools *separately* to ensure each block
% gets a balanced number of L1 and L2 pairs for each condition.

% helper function: split N items into nblocks approximately equal chunks
partition_idx = @(N, nblocks) arrayfun(@(k) ...
    ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), ...
    1:nblocks, 'UniformOutput', false);

% Total per condition = 120. Per block = 30
% Total L1 per condition = 60. Per block = 15
% Total L2 per condition = 60. Per block = 15

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

% This variable needs to be defined *outside* the loop so it persists
% across blocks.
all_foils_persistent = all_foil_pairs; 

for b = 1:p.nBlocks
    fprintf('--- Block %d ---\n', b);
    
    % =====================================================================
    % 5A: SELECT BLOCK-SPECIFIC PAIRS (STRATIFIED)
    % =====================================================================
    % --- Get Comparison pairs for this block (15 L1 + 15 L2) ---
    comp_l1_b = p.stim.comparison_l1(p.block_indices.comp_l1{b}, :);
    comp_l2_b = p.stim.comparison_l2(p.block_indices.comp_l2{b}, :);
    comparison_pairs_b = [comp_l1_b; comp_l2_b];
    comparison_pairs_b = comparison_pairs_b(randperm(height(comparison_pairs_b)), :); % Shuffle them

    % --- Get Isolated pairs for this block (15 L1 + 15 L2) ---
    iso_l1_b = p.stim.iso_both_l1(p.block_indices.iso_l1{b}, :);
    iso_l2_b = p.stim.iso_both_l2(p.block_indices.iso_l2{b}, :);
    iso_both_pairs_b = [iso_l1_b; iso_l2_b];
    iso_both_pairs_b = iso_both_pairs_b(randperm(height(iso_both_pairs_b)), :); % Shuffle them

    % --- Get Novel pairs for this block (15 L1 + 15 L2) ---
    % NOTE: These are ONLY used in the TEST phase (5D)
    nov_l1_b = p.stim.novel_l1(p.block_indices.nov_l1{b}, :);
    nov_l2_b = p.stim.novel_l2(p.block_indices.nov_l2{b}, :);
    novel_pairs_b = [nov_l1_b; nov_l2_b];
    novel_pairs_b = novel_pairs_b(randperm(height(novel_pairs_b)), :); % Shuffle them
    
    nComp_b = height(comparison_pairs_b); % Should be 30
    nIso_b = height(iso_both_pairs_b); % Should be 30
    nNov_b = height(novel_pairs_b); % Should be 30
    
    fprintf('... Loaded %d Comp, %d Iso, %d Nov pairs (L1/L2 balanced)\n', ...
            nComp_b, nIso_b, nNov_b);

    % =====================================================================
    % 5B: BUILD ENCODING PHASE FOR THIS BLOCK (EFFICIENT & FULLY SHUFFLED)
    % =====================================================================
    
    comp_miniblocks = {};  % 2-trial 'k' blocks
    repeat_miniblocks = {}; % 2-trial 'j' blocks
    iso_single_trials = {}; % 1-trial 'new' blocks

    % --- Bucket 1: comparison condition (C-C') ---
    % 30 blocks, 2 trials each
    for i = 1:nComp_b
        comp_miniblocks{end+1} = { ...
            comparison_pairs_b.Target(i), "comparison", "target", "new", "none"; ...
            comparison_pairs_b.Lure(i), "comparison", "lure", "lure", "k"};
    end

    % --- Bucket 2: foil repeats (R-R) ---
    % 30 blocks, 2 trials each (to balance 'k' presses)
    n_foil_repeats = nComp_b; % 30
    
    % *** FIXED: This now correctly uses the 'all_foils_persistent' table ***
    if height(all_foils_persistent) < n_foil_repeats
        error('Not enough foil pairs for encoding foil_repeat spacers in Block %d', b);
    end
    % Get 30 pairs from the top of the persistent stack
    repeat_foil_pairs_b = all_foils_persistent(1:n_foil_repeats, :);
    all_foils_persistent(1:n_foil_repeats, :) = []; % Remove them from the stack
    
    for i = 1:n_foil_repeats
        % We only use the 'foil' (Target) image for the repeat
        foil_item = repeat_foil_pairs_b.foil(i); 
        repeat_miniblocks{end+1} = { ...
            foil_item, "foil_repeat", "filler", "new", "none"; ...
            foil_item, "foil_repeat", "filler", "repeat", "j"};
    end
    
    % --- Bucket 3: isolated items (I and I') ---
    % 60 blocks, 1 trial each (30 targets + 30 lures)
    % NO dedicated fillers. They are just single 'new' trials.
    for i = 1:nIso_b
        % Add Target
        iso_single_trials{end+1} = { ...
            iso_both_pairs_b.Target(i), "iso_both", "target", "new", "none"};
        % Add Lure
        iso_single_trials{end+1} = { ...
            iso_both_pairs_b.Lure(i), "iso_both", "lure", "new", "none"};
    end
    
    % --- combine ALL blocks and shuffle them together ---
    % We have 30 comp blocks + 30 repeat blocks + 60 iso singles = 120 blocks
    all_encoding_blocks = [comp_miniblocks, repeat_miniblocks, iso_single_trials];
    
    % Shuffle the *order* of the blocks
    shuffled_indices = randperm(numel(all_encoding_blocks));
    final_ordered_blocks = all_encoding_blocks(shuffled_indices);
    
    % Concatenate the blocks into one long trial list
    % Total trials: (30*2) + (30*2) + (60*1) = 180 trials. Perfect.
    final_encoding_list = vertcat(final_ordered_blocks{:});
    
    % convert to table
    encoding_schedule_block = cell2table(final_encoding_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    
    % add block number column
    encoding_schedule_block.block = repmat(b, height(encoding_schedule_block), 1);
    
    % add n-back verification columns
    encoding_schedule_block.nback_target_id = strings(height(encoding_schedule_block), 1);
    encoding_schedule_block.trial_type_final = encoding_schedule_block.trial_type_designed;
    
    % check for 1-back repeats
    % This is CRITICAL: it catches accidental 'j' repeats from the shuffle
    for i = 2:height(encoding_schedule_block)
        encoding_schedule_block.nback_target_id(i) = ...
            encoding_schedule_block.stimulus_id(i-1);
        
        if string(encoding_schedule_block.stimulus_id(i)) == ...
           string(encoding_schedule_block.nback_target_id(i))
            
            % This trial is now a 'repeat', even if not designed as one
            % 'j' response ALWAYS takes priority
            encoding_schedule_block.trial_type_final(i) = "repeat";
            encoding_schedule_block.correct_response(i) = "j";
        
        elseif encoding_schedule_block.trial_type_designed(i) == "lure"
            % This was designed as 'k' and is NOT an accidental repeat.
            encoding_schedule_block.trial_type_final(i) = "lure";
            encoding_schedule_block.correct_response(i) = "k";
        else
            % It's not a repeat, not a lure, so it must be 'new'.
            encoding_schedule_block.trial_type_final(i) = "new";
            encoding_schedule_block.correct_response(i) = "none";
        end
    end

    % --- FUCKING CLEANUP SECTION ---
    % Rename 'trial_type_final' to 'item_type'
    encoding_schedule_block.item_type = encoding_schedule_block.trial_type_final;
    
    % Delete the junk columns
    encoding_schedule_block = removevars(encoding_schedule_block, ...
        {'trial_type_designed', 'nback_target_id', 'trial_type_final'});
    
    % =====================================================================
    % 5C: ENSURE TARGET ENCODED BEFORE LURE (ISOLATED CONDITION ONLY)
    % =====================================================================
    % This logic is still necessary *because* we fully shuffled.
    % We just find the single 'I' row and the single 'I'' row and swap if needed.
    
    all_iso_pairs_block = iso_both_pairs_b;
    
    for i = 1:height(all_iso_pairs_block)
        tgt = all_iso_pairs_block.Target(i);
        lur = all_iso_pairs_block.Lure(i);
        
        % Find the *row index* of the target and lure
        idx_tgt_trial = find(encoding_schedule_block.stimulus_id == tgt, 1);
        idx_lur_trial = find(encoding_schedule_block.stimulus_id == lur, 1);
        
        % swap if lure row appears before target row
        if ~isempty(idx_tgt_trial) && ~isempty(idx_lur_trial) && idx_lur_trial < idx_tgt_trial
            
            % Just swap the two single rows
            tmp = encoding_schedule_block(idx_tgt_trial, :);
            encoding_schedule_block(idx_tgt_trial, :) = encoding_schedule_block(idx_lur_trial, :);
            encoding_schedule_block(idx_lur_trial, :) = tmp;
        end
    end
    
    % append into full schedule
    encoding_schedule_all = [encoding_schedule_all; encoding_schedule_block];
    
    % update response counts based on the *final* verified responses
    p.counts.encoding.j_presses = p.counts.encoding.j_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'j'));
    p.counts.encoding.k_presses = p.counts.encoding.k_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'k'));
    p.counts.encoding.no_response = p.counts.encoding.no_response + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'none'));

        



    % =====================================================================
    % =====================================================================
    %
    % 5D: BUILD TEST PHASE FOR THIS BLOCK (STATE-BASED GENERATOR)
    %
    % =====================================================================
    % =====================================================================
    
    fprintf('... building state-based test sequence for block %d ...\n', b);
    
    
    % --- 5D-1: Build the balanced "To-Do List" (Goal Table) ---
    goal_list = table();

    % --- Create goals for: Comparison ---
    pairs = comparison_pairs_b;
    cond_name = "comparison";
    n = height(pairs); % 30
    n_lure = 10; n_same = 10; n_new = 10; % 30/3
    idx_shuf = randperm(n);
    idx_lure = idx_shuf(1:n_lure);
    idx_same = idx_shuf(n_lure+1 : n_lure+n_same);
    idx_new = idx_shuf(n_lure+n_same+1 : end);
    goals_lure = table(pairs.Target(idx_lure), pairs.Lure(idx_lure), ...
        repmat(string(cond_name), n_lure, 1), repmat("lure", n_lure, 1), repmat("k", n_lure, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_same = table(pairs.Target(idx_same), pairs.Target(idx_same), ...
        repmat(string(cond_name), n_same, 1), repmat("same", n_same, 1), repmat("j", n_same, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_new = table(pairs.Target(idx_new), repmat("F_NEW", n_new, 1), ...
        repmat(string(cond_name), n_new, 1), repmat("new", n_new, 1), repmat("none", n_new, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    comp_goals = vertcat(goals_lure, goals_same, goals_new);

    % --- Create goals for: Isolated ---
    pairs = iso_both_pairs_b;
    cond_name = "isolated_both";
    n = height(pairs); % 30
    n_lure = 10; n_same = 10; n_new = 10; % 30/3
    idx_shuf = randperm(n);
    idx_lure = idx_shuf(1:n_lure);
    idx_same = idx_shuf(n_lure+1 : n_lure+n_same);
    idx_new = idx_shuf(n_lure+n_same+1 : end);
    goals_lure = table(pairs.Target(idx_lure), pairs.Lure(idx_lure), ...
        repmat(string(cond_name), n_lure, 1), repmat("lure", n_lure, 1), repmat("k", n_lure, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_same = table(pairs.Target(idx_same), pairs.Target(idx_same), ...
        repmat(string(cond_name), n_same, 1), repmat("same", n_same, 1), repmat("j", n_same, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_new = table(pairs.Target(idx_new), repmat("F_NEW", n_new, 1), ...
        repmat(string(cond_name), n_new, 1), repmat("new", n_new, 1), repmat("none", n_new, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    iso_goals = vertcat(goals_lure, goals_same, goals_new);

    % --- Create goals for: Novel ---
    pairs = novel_pairs_b;
    cond_name = "novel";
    n = height(pairs); % 30
    n_lure = 10; n_same = 10; n_new = 10; % 30/3
    idx_shuf = randperm(n);
    idx_lure = idx_shuf(1:n_lure);
    idx_same = idx_shuf(n_lure+1 : n_lure+n_same);
    idx_new = idx_shuf(n_lure+n_same+1 : end);
    goals_lure = table(pairs.Target(idx_lure), pairs.Lure(idx_lure), ...
        repmat(string(cond_name), n_lure, 1), repmat("lure", n_lure, 1), repmat("k", n_lure, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_same = table(pairs.Target(idx_same), pairs.Target(idx_same), ...
        repmat(string(cond_name), n_same, 1), repmat("same", n_same, 1), repmat("j", n_same, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    goals_new = table(pairs.Target(idx_new), repmat("F_NEW", n_new, 1), ...
        repmat(string(cond_name), n_new, 1), repmat("new", n_new, 1), repmat("none", n_new, 1), ...
        'VariableNames', {'O_item', 'X_item', 'condition', 'goal_type', 'correct_response'});
    nov_goals = vertcat(goals_lure, goals_same, goals_new);

    % --- Combine all 90 goals into one master list and shuffle ---
    goal_list = [comp_goals; iso_goals; nov_goals];
    goal_list = goal_list(randperm(height(goal_list)), :);
    
    % --- 5D-2: Initialize Generator State ---
    
    final_test_list = {}; 
    
    % --- State-tracking variables ---
    % We need to track the item string, if it's a foil, and its lure
    if height(all_foils_persistent) < 2
        error('Ran out of foil pairs for generator init in Block %d', b);
    end
    junk1_pair = all_foils_persistent(1,:); 
    all_foils_persistent(1,:) = [];
    junk2_pair = all_foils_persistent(1,:); 
    all_foils_persistent(1,:) = [];

    stream_items = {junk1_pair.foil, junk2_pair.foil}; % string
    stream_is_foil = [true, true]; % boolean
    stream_foil_lures = {junk1_pair.foil_lure, junk2_pair.foil_lure}; % string
    stream_responses = {"none", "none"}; % string
    
    % Add junk trials to the final list
    final_test_list(1,:) = {stream_items{1}, "init_junk", "filler", "new", "none"};
    final_test_list(2,:) = {stream_items{2}, "init_junk", "filler", "new", "none"};
    
    % --- Generator control variables ---
    probe_N2_next = false; % Rule 1 flag
    current_goal = [];     % Rule 2/4 flag
    goal_idx = 1;
    
    % --- 5D-3: The Generator Loop ---
    
    while goal_idx <= height(goal_list) || ~isempty(current_goal) || probe_N2_next
        
        % Get current state
        N1_item = stream_items{end};
        N2_item = stream_items{end-1};
        N1_resp = stream_responses{end};
        N2_resp = stream_responses{end-1};
        
        % Get N-2 item's properties (for Rule 1)
        N2_is_foil = stream_is_foil(end-1);
        N2_lure = stream_foil_lures{end-1};
        
        % Rule 3: Check for response streaks to avoid
        avoid_resp = "none";
        if N1_resp == "j" && N2_resp == "j"
            avoid_resp = "j";
        elseif N1_resp == "k" && N2_resp == "k"
            avoid_resp = "k";
        end
        
        % Initialize next trial variables
        next_item = "";
        next_is_foil = false;
        next_foil_lure = "";
        next_condition = "filler_stream";
        next_role = "filler";
        next_trial_type = "new";
        next_response = "none";
        
        % --- APPLY RULES IN ORDER OF PRIORITY ---
        
        % --- Priority 1: Rule 1 (Probe Filler) ---
        if probe_N2_next
            probe_N2_next = false; % Reset flag
            
            % Sanity check: Rule 1 should only fire on a foil
            if ~N2_is_foil
                error('CRITICAL ERROR: Rule 1 (Probe) triggered on a non-foil item!');
            end
            
            prob = rand();
            
            % Choice 1: 'j' response (probe with F)
            if prob < 0.33 && avoid_resp ~= "j"
                next_item = N2_item;
                next_is_foil = true;
                next_foil_lure = N2_lure;
                next_response = "j";
                next_trial_type = "repeat";
                
            % Choice 2: 'k' response (probe with F')
            elseif prob < 0.66 && avoid_resp ~= "k"
                next_item = N2_lure;
                next_is_foil = true;
                next_foil_lure = N2_lure; % Track the original lure
                next_response = "k";
                next_trial_type = "lure";
                
            % Choice 3: 'new' response (probe with F_new)
            else
                if height(all_foils_persistent) < 1, error('Ran out of foil pairs for Rule 1 probe in Block %d', b), end;
                new_F_pair = all_foils_persistent(1,:); all_foils_persistent(1,:) = [];
                
                next_item = new_F_pair.foil;
                next_is_foil = true;
                next_foil_lure = new_F_pair.foil_lure;
                next_response = "none";
                next_trial_type = "new";
            end
            
            next_condition = "filler_probe";
            next_role = "filler";

        % --- Priority 2: Rule 2 (Complete Triplet) ---
        % Check: State is O-F, where F is a foil
        elseif ~isempty(current_goal) && string(N2_item) == string(current_goal.O_item) && (startsWith(N1_item, "mst_") && endsWith(N1_item, "_foil.png"))
            
            % Check Rule 3: Is streak-breaking needed?
            if string(current_goal.correct_response) == avoid_resp
                % DELAY GOAL: Insert a filler trial
                if height(all_foils_persistent) < 1, error('Ran out of foil pairs for Rule 3 delay in Block %d', b), end;
                new_F_pair = all_foils_persistent(1,:); all_foils_persistent(1,:) = [];
                
                next_item = new_F_pair.foil;
                next_is_foil = true;
                next_foil_lure = new_F_pair.foil_lure;
                next_response = "none";
                next_trial_type = "new";
                next_condition = "filler_delay";
                % We *do not* clear current_goal. We try again next trial.
            
            else
                % COMPLETE GOAL: Append the X item
                next_item = current_goal.X_item;
                next_is_foil = false; % This is an O or O' item
                next_foil_lure = "";
                
                % Check if the goal is a 'new' type (O-F-F_NEW)
                if string(next_item) == "F_NEW"
                    if height(all_foils_persistent) < 1, error('Ran out of foil pairs for F_NEW completion in Block %d', b), end;
                    new_F_pair = all_foils_persistent(1,:); all_foils_persistent(1,:) = [];
                    
                    next_item = new_F_pair.foil;
                    next_is_foil = true;
                    next_foil_lure = new_F_pair.foil_lure;
                end
                
                next_condition = current_goal.condition;
                if current_goal.goal_type == "same", next_role = "target"; else, next_role = "lure"; end
                next_trial_type = current_goal.goal_type;
                next_response = current_goal.correct_response;
                
                probe_N2_next = true; % SET RULE 1 FLAG
                current_goal = [];    % Clear sub-goal
            end
            
        % --- Priority 4: Rule 4 (Continue/Start New) ---
        elseif ~isempty(current_goal)
            % Case 4a: Continue Triplet (State is ...O)
            % We just placed 'O', now we MUST place 'F'
            if height(all_foils_persistent) < 1, error('Ran out of foil pairs for Rule 4 filler in Block %d', b), end;
            new_F_pair = all_foils_persistent(1,:); all_foils_persistent(1,:) = [];

            next_item = new_F_pair.foil;
            next_is_foil = true;
            next_foil_lure = new_F_pair.foil_lure;
            next_response = "none";
            next_trial_type = "new";
            next_condition = "triplet_filler";
            % current_goal remains active
            
        else
            % Case 4b: Start New Triplet (Neutral State)
            if goal_idx > height(goal_list)
                % We finished all goals, but probe_N2_next might be true
                if ~probe_N2_next
                    break; % Generator is done
                end
                % If probe_N2_next is true, loop will run 1 more time
            else
                % Get next goal from the list
                current_goal = goal_list(goal_idx, :);
                goal_idx = goal_idx + 1;
                
                next_item = current_goal.O_item;
                next_is_foil = false;
                next_foil_lure = "";
                next_response = "none";
                next_trial_type = "new"; % Always 'new' against N-2 junk
                next_condition = current_goal.condition;
                next_role = "target";
            end
        end
        
        % --- Append the chosen trial to stream and list ---
        if ~isempty(next_item)
            stream_items{end+1} = next_item;
            stream_is_foil(end+1) = next_is_foil;
            stream_foil_lures{end+1} = next_foil_lure;
            stream_responses{end+1} = next_response;
            
            final_test_list(end+1,:) = {next_item, next_condition, next_role, ...
                                        next_trial_type, next_response};
        end
    end
    
    % --- 5D-4: Convert to Table and Final Verification ---
    
    % Convert cell array to table (skipping the first 2 junk trials)
    test_schedule_block = cell2table(final_test_list(3:end,:), ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    
    % add block number column
    test_schedule_block.block = repmat(b, height(test_schedule_block), 1);
    
    % add n-back verification columns (we delete these later)
    test_schedule_block.nback_target_id = strings(height(test_schedule_block), 1);
    test_schedule_block.trial_type_final = test_schedule_block.trial_type_designed;
    
    % --- FINAL N-BACK VERIFICATION PASS ---
    % Re-build the full stream including the 2 junk items for N-back
    full_stream_stims = [string(final_test_list(1,1)); 
                         string(final_test_list(2,1)); 
                         test_schedule_block.stimulus_id];
                     
    for i = 1:height(test_schedule_block)
        current_trial_idx = i + 2; % Offset by 2 junk trials
        n2_target_idx = current_trial_idx - 2;
        
        n2_target_stim = full_stream_stims(n2_target_idx);
        current_stim = full_stream_stims(current_trial_idx);
        
        test_schedule_block.nback_target_id(i) = n2_target_stim;
        
        if current_stim == n2_target_stim
            % It's a 'j' press, regardless of what we designed
            test_schedule_block.trial_type_final(i) = "repeat";
            test_schedule_block.correct_response(i) = "j";
        elseif test_schedule_block.trial_type_designed(i) == "lure"
            % It's a 'k' press (we trust the design)
            test_schedule_block.trial_type_final(i) = "lure";
            test_schedule_block.correct_response(i) = "k";
        else
            % It's a 'new' trial
            test_schedule_block.trial_type_final(i) = "new";
            test_schedule_block.correct_response(i) = "none";
        end
    end
    
    % --- FUCKING CLEANUP SECTION ---
    % Rename 'trial_type_final' to 'item_type'
    test_schedule_block.item_type = test_schedule_block.trial_type_final;
    
    % Delete the junk columns
    test_schedule_block = removevars(test_schedule_block, ...
        {'trial_type_designed', 'nback_target_id', 'trial_type_final'});
    
    
    % append into full schedule
    test_schedule_all = [test_schedule_all; test_schedule_block]; %#ok<AGROW>
    
    % update response counts based on the *final* verified responses
    p.counts.test.j_presses = p.counts.test.j_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'k'));
    p.counts.test.no_response = p.counts.test.no_response + ...
        sum(strcmp(test_schedule_block.correct_response, 'none'));
        
end % --- END OF BLOCK LOOP ---

%% ========================================================================
%  SECTION 6: GENERATE JITTERED FIXATION DURATIONS
%  ========================================================================
% create random fixation durations for each trial
% base duration ± jitter: 0.75 ± 0.25 = [0.5, 1.0] seconds

n_enc_trials = height(encoding_schedule_all);
n_test_trials = height(test_schedule_all);

% generate jittered durations using uniform distribution
encoding_schedule_all.fix_duration = ...
    p.timing.fix_dur + (rand(n_enc_trials, 1) * 2 - 1) * p.timing.fix_jitter;

test_schedule_all.fix_duration = ...
    p.timing.fix_dur + (rand(n_test_trials, 1) * 2 - 1) * p.timing.fix_jitter;

%% ========================================================================
%  SECTION 7: SAVE OUTPUT
%  ========================================================================

p.stim.all_foils_remaining = all_foils_persistent; % Save unused foils
p.encoding_schedule = encoding_schedule_all;
p.test_schedule = test_schedule_all;

save(output_filename, 'p');

fprintf('========================================\n');
fprintf('Setup for subject %03d complete.\n', subj_id);
fprintf('Saved to: %s\n', output_filename);
fprintf('========================================\n');
fprintf('Encoding task: %d trials total (%d j, %d k, %d new)\n', ...
    n_enc_trials, p.counts.encoding.j_presses, p.counts.encoding.k_presses, p.counts.encoding.no_response);
fprintf('Test task: %d trials total (%d j, %d k, %d new)\n', ...
    n_test_trials, p.counts.test.j_presses, p.counts.test.k_presses, p.counts.test.no_response);

% --- END OF SCRIPT ---