%==========================================================================
%         Generate Subject-Specific Task Variables
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
% note: you must run this for each subject before running the task
%
% MODIFIED: 2025-10-28 to implement state-based 2-back test generator
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
assert(numel(all_targ_files_l1) == 140, 'expected 180 target L1 files!');
targ_names_l1 = string({all_targ_files_l1.name}');
lure_names_l1 = string({all_lure_files_l1.name}');
master_pair_list_l1 = table(sort(targ_names_l1), sort(lure_names_l1), ...
    'VariableNames', {'Img_A', 'Img_B'});

% --- Load L2 Pairs (Bin 2) ---
all_targ_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_targ_l2.png'));
all_lure_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_lure_l2.png'));
assert(numel(all_targ_files_l2) == 140, 'expected 180 target L2 files!');
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

for b = 1:p.nBlocks
    fprintf('Block %d\n', b);
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
    n_foil_repeats = nComp_b;

    if isempty(all_foils) || numel(all_foils) < n_foil_repeats
        error('Not enough foils for foil_repeat spacers');
    end
    repeat_foils = all_foils(1:n_foil_repeats);
    all_foils(1:n_foil_repeats) = [];

    for i = 1:n_foil_repeats
        repeat_miniblocks{end+1} = { ...
            repeat_foils(i), "foil_repeat", "filler", "new", "none"; ...
            repeat_foils(i), "foil_repeat", "filler", "repeat", "j"};
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

    % recode variable names
    encoding_schedule_block.item_type = encoding_schedule_block.trial_type_final;

    % delete the junk columns
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
    % 5D: BUILD 2-BACK TEST PHASE FOR THIS BLOCK
    % =====================================================================

    fprintf('... Generating 2-back test sequence\n');

    test_schedule_block = generate_2back_test(comparison_pairs_b, ...
        iso_both_pairs_b, novel_pairs_b, all_foils, b);

    % Append to full test schedule
    test_schedule_all = [test_schedule_all; test_schedule_block];

    % Update test response counts
    p.counts.test.j_presses = p.counts.test.j_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'k'));
    p.counts.test.no_response = p.counts.test.no_response + ...
        sum(strcmp(test_schedule_block.correct_response, 'none'));

end



% =========================================================================
% 2-BACK TEST PHASE GENERATOR (PART B)
% Trial-by-trial state-aware sequence builder
% =========================================================================

function test_schedule_block = generate_2back_test(comparison_pairs_b, iso_both_pairs_b, novel_pairs_b, all_foils, block_num)
    
    % =====================================================================
    % STEP 1: BUILD THE TO-DO LIST (90 GOALS)
    % =====================================================================
    % Each goal represents ONE pair to test with ONE response type
    % Goal structure: O-F-X where O and X come from the pair
    
    todo_list = {};  % List of goals (will shrink as we consume them)
    
    conditions = {'comparison', 'iso_both', 'novel'};
    probe_types = {'lure', 'repeat', 'new'};
    pairs_by_cond = {comparison_pairs_b, iso_both_pairs_b, novel_pairs_b};
    
    for c = 1:3
        cond_name = conditions{c};
        pairs = pairs_by_cond{c};
        
        % Shuffle pairs for this condition
        pairs = pairs(randperm(height(pairs)), :);
        
        pair_idx = 1;
        
        for p = 1:3
            probe_type = probe_types{p};
            
            % 10 goals of this type (10 pairs × 3 probe types = 30 per condition)
            for rep = 1:10
                goal = struct();
                goal.condition = cond_name;
                goal.probe_type = probe_type;
                
                % Assign stimuli for this pair
                goal.target = pairs.Target(pair_idx);
                goal.lure = pairs.Lure(pair_idx);
                
                % O is ALWAYS the target (must appear first!)
                goal.O = goal.target;
                goal.O_role = 'target';
                
                % X depends on probe_type
                switch probe_type
                    case 'lure'
                        % X is the lure
                        goal.X = goal.lure;
                        goal.X_role = 'lure';
                        goal.X_response = 'k';
                        
                    case 'repeat'
                        % X is the SAME as O
                        goal.X = goal.O;
                        goal.X_role = goal.O_role;
                        goal.X_response = 'j';
                        
                    case 'new'
                        % X is a NEW foil
                        goal.X = '';  % Will be assigned from foil pool later
                        goal.X_role = 'foil';
                        goal.X_response = 'j';
                end
                
                goal.O_placed = false;
                goal.F_placed = false;
                goal.X_placed = false;
                
                todo_list{end+1} = goal;
                pair_idx = pair_idx + 1;
            end
        end
    end
    
    % SHUFFLE THE TO-DO LIST (CRITICAL!)
    todo_indices = randperm(length(todo_list));
    todo_list = todo_list(todo_indices);
    
    fprintf('  ... Built shuffled to-do list: %d goals\n', length(todo_list));
    
    % =====================================================================
    % STEP 2: INITIALIZE GENERATOR STATE
    % =====================================================================
    
    state = struct();
    state.sequence = {};  % The growing trial list
    state.N_minus_1 = '';  % Item at N-1
    state.N_minus_2 = '';  % Item at N-2
    state.current_goal = [];  % Current goal being worked on
    state.last_two_responses = {'', ''};  % Track last 2 REQUIRED responses
    state.foil_pool = all_foils;  % Foil pool for drawing new items
    state.just_completed_triplet = false;  % Flag for Rule 1
    
    % =====================================================================
    % STEP 3: TRIAL-BY-TRIAL GENERATOR LOOP
    % =====================================================================
    
    max_trials = 300;  % Safety limit
    trial_count = 0;
    
    while trial_count < max_trials
        trial_count = trial_count + 1;
        
        % Check if we're done (to-do list is empty)
        if isempty(todo_list)
            fprintf('  ... All 90 pairs tested after %d trials\n', trial_count - 1);
            break;
        end
        
        % ================================================================
        % APPLY THE 4-RULE HIERARCHY
        % ================================================================
        
        next_item = struct();
        next_item.stimulus_id = '';
        next_item.condition = '';
        next_item.role = '';
        next_item.correct_response = 'none';
        next_item.rule_applied = '';
        
        % -----------------------------------------------------------------
        % RULE 1: PROBE THE FILLER (Priority 1)
        % -----------------------------------------------------------------
        % Check if a triplet (O-F-X) just completed on previous trial
        if state.just_completed_triplet
            
            % The filler F is now at N-2
            F_item = state.N_minus_2;
            
            % Randomly choose probe type
            probe_choice = randi(3);
            
            if probe_choice == 1
                % Probe with F itself (repeat -> 'j')
                next_item.stimulus_id = F_item;
                next_item.role = 'filler_probe_repeat';
                next_item.correct_response = 'j';
                
            elseif probe_choice == 2
                % Probe with F's lure (if it has one -> 'k')
                % Find F's lure
                F_lure = '';
                all_pairs_temp = [comparison_pairs_b; iso_both_pairs_b; novel_pairs_b];
                
                % Check if F is in Target column
                idx_temp = find(strcmp(all_pairs_temp.Target, F_item), 1);
                if ~isempty(idx_temp)
                    F_lure = all_pairs_temp.Lure(idx_temp);
                else
                    % Check if F is in Lure column
                    idx_temp = find(strcmp(all_pairs_temp.Lure, F_item), 1);
                    if ~isempty(idx_temp)
                        F_lure = all_pairs_temp.Target(idx_temp);
                    end
                end
                
                if ~isempty(F_lure)
                    next_item.stimulus_id = F_lure;
                    next_item.role = 'filler_probe_lure';
                    next_item.correct_response = 'k';
                else
                    % No lure available, use new foil instead
                    if isempty(state.foil_pool)
                        error('Ran out of foils!');
                    end
                    next_item.stimulus_id = state.foil_pool(1);
                    state.foil_pool(1) = [];
                    next_item.role = 'filler_probe_foil';
                    next_item.correct_response = 'j';
                end
                
            else
                % Probe with new foil (-> 'j')
                if isempty(state.foil_pool)
                    error('Ran out of foils!');
                end
                next_item.stimulus_id = state.foil_pool(1);
                state.foil_pool(1) = [];
                next_item.role = 'filler_probe_foil';
                next_item.correct_response = 'j';
            end
            
            next_item.condition = 'probe';
            next_item.rule_applied = 'Rule1_ProbeFiller';
            
            % Clear the flag
            state.just_completed_triplet = false;
            
        % -----------------------------------------------------------------
        % RULE 2: COMPLETE A TRIPLET (Priority 2)
        % -----------------------------------------------------------------
        elseif ~isempty(state.current_goal) && ...
               state.current_goal.O_placed && state.current_goal.F_placed && ...
               ~state.current_goal.X_placed
            
            % State is O-F, now place X to complete the triplet
            goal = state.current_goal;
            
            % If X is 'new', assign a foil now
            if strcmp(goal.probe_type, 'new') && isempty(goal.X)
                if isempty(state.foil_pool)
                    error('Ran out of foils!');
                end
                goal.X = state.foil_pool(1);
                state.foil_pool(1) = [];
            end
            
            next_item.stimulus_id = goal.X;
            next_item.condition = goal.condition;
            next_item.role = goal.X_role;
            next_item.correct_response = goal.X_response;
            next_item.rule_applied = 'Rule2_CompleteTriplet';
            
            % Mark X as placed
            goal.X_placed = true;
            state.current_goal = goal;
            
            % Set flag to trigger Rule 1 on next trial
            state.just_completed_triplet = true;
            
            % Clear current goal (it's done, don't need it anymore)
            state.current_goal = [];
            
        % -----------------------------------------------------------------
        % RULE 3: BREAK RESPONSE STREAKS (Priority 3)
        % -----------------------------------------------------------------
        elseif ~isempty(state.current_goal) && ...
               state.current_goal.O_placed && ~state.current_goal.F_placed
            
            % We're at O-?, need to place F
            goal = state.current_goal;
            
            % Check if placing X (in 2 trials) would create a 3-in-a-row streak
            next_response = goal.X_response;
            
            if ~isempty(state.last_two_responses{1}) && ~isempty(state.last_two_responses{2}) && ...
               strcmp(state.last_two_responses{1}, state.last_two_responses{2}) && ...
               strcmp(state.last_two_responses{2}, next_response)
                
                % ABORT this goal, place a streak-breaker foil trial
                if isempty(state.foil_pool)
                    error('Ran out of foils!');
                end
                next_item.stimulus_id = state.foil_pool(1);
                state.foil_pool(1) = [];
                next_item.condition = 'streak_breaker';
                next_item.role = 'foil';
                next_item.correct_response = 'j';
                next_item.rule_applied = 'Rule3_BreakStreak';
                
                % Reset current goal and put it BACK on the to-do list
                goal.O_placed = false;
                goal.F_placed = false;
                goal.X_placed = false;
                todo_list{end+1} = goal;  % Add back to end of list
                state.current_goal = [];
                
            else
                % No streak, place F normally
                if isempty(state.foil_pool)
                    error('Ran out of foils!');
                end
                next_item.stimulus_id = state.foil_pool(1);
                state.foil_pool(1) = [];
                next_item.condition = goal.condition;
                next_item.role = 'filler';
                next_item.correct_response = 'none';
                next_item.rule_applied = 'Rule4_PlaceFiller';
                
                goal.F_placed = true;
                state.current_goal = goal;
            end
            
        % -----------------------------------------------------------------
        % RULE 4: CONTINUE/START NEW (Priority 4)
        % -----------------------------------------------------------------
        else
            % No current goal, pick a new one from the to-do list
            
            if isempty(todo_list)
                % Should never hit this
                break;
            end
            
            % Pop the first goal from the list
            goal = todo_list{1};
            todo_list(1) = [];  % REMOVE IT FROM THE LIST
            
            % Place O
            next_item.stimulus_id = goal.O;
            next_item.condition = goal.condition;
            next_item.role = goal.O_role;
            next_item.correct_response = 'none';
            next_item.rule_applied = 'Rule4_StartGoal';
            
            % Update goal state
            goal.O_placed = true;
            state.current_goal = goal;
        end
        
        % ================================================================
        % APPEND TRIAL TO SEQUENCE
        % ================================================================
        
        next_item.trial_num = trial_count;
        next_item.block = block_num;
        next_item.N_minus_1 = state.N_minus_1;
        next_item.N_minus_2 = state.N_minus_2;
        
        state.sequence{end+1} = next_item;
        
        % Update state
        state.N_minus_2 = state.N_minus_1;
        state.N_minus_1 = next_item.stimulus_id;
        
        % Update response history (only for REQUIRED responses, not 'none')
        if ~strcmp(next_item.correct_response, 'none')
            state.last_two_responses{1} = state.last_two_responses{2};
            state.last_two_responses{2} = next_item.correct_response;
        end
    end
    
    % =====================================================================
    % STEP 4: CONVERT TO TABLE
    % =====================================================================
    
    if isempty(state.sequence)
        error('Generated empty sequence! Check logic.');
    end
    
    test_schedule_block = struct2table(vertcat(state.sequence{:}));
    
    % Reorder columns to match encoding format
    test_schedule_block = test_schedule_block(:, {'trial_num', 'stimulus_id', ...
        'condition', 'role', 'correct_response', 'block', 'rule_applied', ...
        'N_minus_1', 'N_minus_2'});
    
    fprintf('  ... Generated %d trials for block %d\n', height(test_schedule_block), block_num);
    
    % Count responses
    n_j = sum(strcmp(test_schedule_block.correct_response, 'j'));
    n_k = sum(strcmp(test_schedule_block.correct_response, 'k'));
    n_none = sum(strcmp(test_schedule_block.correct_response, 'none'));
    fprintf('  ... Responses: %d j, %d k, %d none\n', n_j, n_k, n_none);
end% =========================================================================





% =====================================================================
% 5D: BUILD TEST PHASE (DETERMINISTIC SEQUENTIAL CHAIN)
% =====================================================================

fprintf('... building deterministic test sequence for block %d ...\n', b);

% --- 5D-1: Build the Block's Item Pool and Goal List ---

block_pairs = [comparison_pairs_b; iso_both_pairs_b; novel_pairs_b];
block_conditions = [repmat("comparison", nComp_b, 1); ...
                    repmat("isolated_both", nIso_b, 1); ...
                    repmat("novel", nNov_b, 1)];

% Goal list: each row is one 2-back comparison to test
goal_list = table('Size', [90, 6], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'string', 'logical'}, ...
    'VariableNames', {'O_item', 'O_lure', 'condition', 'goal_type', 'X_item', 'used'});
    
goal_list.O_item = block_pairs.Target;
goal_list.O_lure = block_pairs.Lure;
goal_list.condition = block_conditions;
goal_list.used = false(90, 1);

% Assign goal types (1/3 lure, 1/3 same, 1/3 new)
idx_shuf = randperm(90);
goal_list.goal_type(idx_shuf(1:30)) = "lure";
goal_list.goal_type(idx_shuf(31:60)) = "same";
goal_list.goal_type(idx_shuf(61:90)) = "new";

% --- 5D-2: Assign X_item Deterministically ---
% Strategy: For "new" goals, X_item should be the O_item of another unused goal
% For "same" and "lure", X_item is forced

% Separate goals by type
new_goals_idx = find(goal_list.goal_type == "new");
same_lure_goals_idx = find(goal_list.goal_type ~= "new");

% Shuffle new goals to create random order
new_goals_shuffled = new_goals_idx(randperm(length(new_goals_idx)));

% For each "new" goal, assign X_item as the O_item of the next "new" goal
for i = 1:length(new_goals_shuffled)
    goal_idx = new_goals_shuffled(i);
    
    if i < length(new_goals_shuffled)
        % Link to next new goal
        next_goal_idx = new_goals_shuffled(i+1);
        goal_list.X_item(goal_idx) = goal_list.O_item(next_goal_idx);
    else
        % Last new goal: link to first new goal (creates loop)
        first_goal_idx = new_goals_shuffled(1);
        goal_list.X_item(goal_idx) = goal_list.O_item(first_goal_idx);
    end
end

% For same/lure goals, X is forced
for i = 1:length(same_lure_goals_idx)
    goal_idx = same_lure_goals_idx(i);
    if goal_list.goal_type(goal_idx) == "lure"
        goal_list.X_item(goal_idx) = goal_list.O_lure(goal_idx);
    else % "same"
        goal_list.X_item(goal_idx) = goal_list.O_item(goal_idx);
    end
end

fprintf('  -> 90 goals: %d same, %d lure, %d new\n', ...
    sum(goal_list.goal_type == "same"), ...
    sum(goal_list.goal_type == "lure"), ...
    sum(goal_list.goal_type == "new"));

% --- 5D-3: Create Processing Order ---
% Interleave new goals with same/lure goals for variety

processing_order = zeros(90, 1);
new_goals_remaining = new_goals_shuffled;
same_lure_remaining = same_lure_goals_idx(randperm(length(same_lure_goals_idx)));

order_idx = 1;
while ~isempty(new_goals_remaining) || ~isempty(same_lure_remaining)
    % Add 2-3 new goals
    n_new = min(randi([2, 3]), length(new_goals_remaining));
    for k = 1:n_new
        processing_order(order_idx) = new_goals_remaining(1);
        new_goals_remaining(1) = [];
        order_idx = order_idx + 1;
    end
    
    % Add 1 same/lure goal
    if ~isempty(same_lure_remaining)
        processing_order(order_idx) = same_lure_remaining(1);
        same_lure_remaining(1) = [];
        order_idx = order_idx + 1;
    end
end

% --- 5D-4: Build Sequence Sequentially ---

sequence = cell(500, 5);
trial_idx = 1;

% Initialize with 2 junk trials
if height(all_foils_remain) < 2
    error('Ran out of foil pairs for init in Block %d', b);
end
junk1 = all_foils_remain(1,:); all_foils_remain(1,:) = [];
junk2 = all_foils_remain(1,:); all_foils_remain(1,:) = [];

sequence(trial_idx,:) = {junk1.foil, "init_junk", "filler", "new", "none"};
trial_idx = trial_idx + 1;
sequence(trial_idx,:) = {junk2.foil, "init_junk", "filler", "new", "none"};
trial_idx = trial_idx + 1;

buffer = {junk1.foil, junk2.foil}; % Track last 2 items

% Process goals in order
for g = 1:90
    goal_idx = processing_order(g);
    goal = goal_list(goal_idx, :);
    
    if goal.used
        continue; % Skip if already used
    end
    
    % --- Trial 1: Append O_item (start goal) ---
    O_item = goal.O_item;
    O_cond = goal.condition;
    
    sequence(trial_idx,:) = {O_item, O_cond, "target", "new", "none"};
    trial_idx = trial_idx + 1;
    buffer{end+1} = O_item;
    
    % --- Trial 2: Append filler (another O_item if available) ---
    % Look ahead for next unused goal
    next_unused_idx = find(~goal_list.used & processing_order > goal_idx, 1);
    if ~isempty(next_unused_idx)
        filler_goal_idx = processing_order(next_unused_idx);
        filler_item = goal_list.O_item(filler_goal_idx);
    else
        % Use a random O_item that won't interfere
        filler_item = goal_list.O_item(randi(height(goal_list)));
    end
    
    sequence(trial_idx,:) = {filler_item, "filler", "filler", "new", "none"};
    trial_idx = trial_idx + 1;
    buffer{end+1} = filler_item;
    
    % --- Trial 3: Append X_item (complete goal) ---
    X_item = goal.X_item;
    X_type = goal.goal_type;
    
    if X_type == "lure"
        X_resp = "k";
        X_role = "lure";
    elseif X_type == "same"
        X_resp = "j";
        X_role = "target";
    else % "new"
        X_resp = "none";
        X_role = "target";
    end
    
    sequence(trial_idx,:) = {X_item, O_cond, X_role, X_type, X_resp};
    trial_idx = trial_idx + 1;
    buffer{end+1} = X_item;
    
    % Mark goal as used
    goal_list.used(goal_idx) = true;
    
    % --- Trial 4: Filler Probe (make filler task-relevant) ---
    filler_at_N2 = string(buffer{end-2}); % The filler from 2 trials ago
    
    % Find if filler is an O_item with a goal
    filler_goal_idx = find(goal_list.O_item == filler_at_N2, 1);
    
    prob = rand();
    if prob < 0.33 && ~isempty(filler_goal_idx) && ~goal_list.used(filler_goal_idx)
        % Probe with same (filler itself)
        probe_item = filler_at_N2;
        probe_resp = "j";
        probe_type = "repeat";
    elseif prob < 0.66 && ~isempty(filler_goal_idx)
        % Probe with lure
        probe_item = goal_list.O_lure(filler_goal_idx);
        probe_resp = "k";
        probe_type = "lure";
    else
        % Probe with new (look ahead for next unused goal)
        next_unused = find(~goal_list.used, 1);
        if ~isempty(next_unused)
            probe_item = goal_list.O_item(next_unused);
        else
            probe_item = goal_list.O_lure(randi(height(goal_list)));
        end
        probe_resp = "none";
        probe_type = "new";
    end
    
    sequence(trial_idx,:) = {probe_item, "filler_probe", "filler", probe_type, probe_resp};
    trial_idx = trial_idx + 1;
    buffer{end+1} = probe_item;
end

% Add end junk
if height(all_foils_remain) >= 5
    for jj = 1:5
        sequence(trial_idx,:) = {all_foils_remain.foil(jj), "end_junk", "filler", "new", "none"};
        trial_idx = trial_idx + 1;
    end
    all_foils_remain(1:5,:) = [];
end

% Convert to table and trim
final_test_list = cell2table(sequence(1:trial_idx-1, :), ...
    'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type', 'correct_response'});

fprintf('... Generated %d test trials for block %d\n', height(final_test_list), b);
fprintf('... Goals used: %d / 90\n', sum(goal_list.used));


%% ----------

% =====================================================================
% 5D: BUILD TEST PHASE (PROPER 2-BACK WITH GOAL REMOVAL)
% =====================================================================

fprintf('... building 2-back test sequence for block %d ...\n', b);

% --- 5D-1: Prepare Goal List ---

block_pairs = [comparison_pairs_b; iso_both_pairs_b; novel_pairs_b];
block_conditions = [repmat("comparison", nComp_b, 1); ...
                    repmat("isolated_both", nIso_b, 1); ...
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

% Assign X_item based on goal_type
% For "new" goals, we need to carefully select X_items from OTHER goals
% that won't cause conflicts

% First, assign X for same/lure (these are deterministic)
for i = 1:90
    if goal_list.goal_type(i) == "lure"
        goal_list.X_item(i) = goal_list.O_lure(i);
    elseif goal_list.goal_type(i) == "same"
        goal_list.X_item(i) = goal_list.O_item(i);
    else
        goal_list.X_item(i) = ""; % Placeholder for "new" goals
    end
end

% Now assign X for "new" goals from a pool of O_items
% Strategy: Create pairs of "new" goals where each serves as the other's X_item
new_goal_indices = find(goal_list.goal_type == "new");
new_goals_shuffled = new_goal_indices(randperm(length(new_goal_indices)));

for i = 1:length(new_goals_shuffled)
    idx = new_goals_shuffled(i);
    
    % Link to next new goal (circular)
    if i < length(new_goals_shuffled)
        next_idx = new_goals_shuffled(i+1);
    else
        next_idx = new_goals_shuffled(1);
    end
    
    goal_list.X_item(idx) = goal_list.O_item(next_idx);
end

% Shuffle entire list to randomize processing order
goal_list = goal_list(randperm(height(goal_list)), :);

% Mark which goals should NOT start (because they're used as X_items for "new" goals)
goals_used_as_X = false(height(goal_list), 1);
for i = 1:height(goal_list)
    if goal_list.goal_type(i) == "new"
        % Find which goal's O_item matches this X_item
        X_val = goal_list.X_item(i);
        matching_goal = find(strcmp(goal_list.O_item, X_val), 1);
        if ~isempty(matching_goal)
            goals_used_as_X(matching_goal) = true;
        end
    end
end

fprintf('  -> %d goals will start as O_items\n', sum(~goals_used_as_X));
fprintf('  -> %d goals used only as X_items (no separate start)\n', sum(goals_used_as_X));

fprintf('  -> Testing %d goals\n', height(goal_list));

% --- 5D-2: Build Sequence with Goal Tracking ---

sequence = cell(300, 5);
row_idx = 1;

% Init with 2 junk
if height(all_foils_remain) < 2
    error('Need foils for init');
end
junk1 = all_foils_remain(1,:); all_foils_remain(1,:) = [];
junk2 = all_foils_remain(1,:); all_foils_remain(1,:) = [];

sequence(row_idx,:) = {junk1.foil, "init_junk", "filler", "new", "none"};
row_idx = row_idx + 1;
sequence(row_idx,:) = {junk2.foil, "init_junk", "filler", "new", "none"};
row_idx = row_idx + 1;

% Track active goals (goals that have started but not completed)
% Each entry: [goal_index, position_where_O_appeared]
active_goals = [];

goal_pointer = 1; % Which goal to start next

while goal_pointer <= height(goal_list) || ~isempty(active_goals)
    
    % Check if any goal at position N-2 needs completion
    goals_to_complete = [];
    for k = 1:size(active_goals, 1)
        goal_idx = active_goals(k, 1);
        goal_pos = active_goals(k, 2);
        
        % If this goal's O appeared 2 trials ago, complete it now
        if row_idx - goal_pos == 2
            goals_to_complete(end+1) = k;
        end
    end
    
    if ~isempty(goals_to_complete)
        % Complete the first goal that's ready
        k = goals_to_complete(1);
        goal_idx = active_goals(k, 1);
        goal = goal_list(goal_idx, :);
        
        % Append X_item
        X_type = goal.goal_type;
        if X_type == "lure"
            X_resp = "k";
            X_role = "lure";
        elseif X_type == "same"
            X_resp = "j";
            X_role = "target";
        else % "new"
            X_resp = "none";
            X_role = "target";
        end
        
        sequence(row_idx,:) = {goal.X_item, goal.condition, X_role, X_type, X_resp};
        row_idx = row_idx + 1;
        
        % Remove this goal from active list (it's completed!)
        active_goals(k, :) = [];
        
    else
        % Start a new goal (only if it's not marked as used_as_X)
        while goal_pointer <= height(goal_list) && goals_used_as_X(goal_pointer)
            goal_pointer = goal_pointer + 1; % Skip goals used as X_items
        end
        
        if goal_pointer <= height(goal_list)
            goal = goal_list(goal_pointer, :);
            
            % Append O_item
            sequence(row_idx,:) = {goal.O_item, goal.condition, "target", "new", "none"};
            
            % Track this goal as active
            active_goals(end+1, :) = [goal_pointer, row_idx];
            
            row_idx = row_idx + 1;
            goal_pointer = goal_pointer + 1;
        else
            % No more goals to start
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

% Convert to table and trim
final_test_list = cell2table(sequence(1:row_idx-1, :), ...
    'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type', 'correct_response'});

fprintf('... Generated %d test trials\n', height(final_test_list));
fprintf('... All 90 goals tested\n');

% --- 5D-3: SANITY CHECK ---
fprintf('\n=== SANITY CHECK ===\n');

% Check 1: Count unique O_items that appear as "target" with "new" response
targets_new = final_test_list(strcmp(final_test_list.role, "target") & ...
                               strcmp(final_test_list.trial_type, "new"), :);
unique_targets = unique(targets_new.stimulus_id);
fprintf('✓ Unique O-items that started goals: %d (expected: 90)\n', length(unique_targets));

% Check 2: Verify 2-back structure for each goal
fprintf('\n--- Checking 2-back structure for each goal ---\n');
errors = 0;
for i = 1:height(goal_list)
    O_item = goal_list.O_item(i);
    X_item = goal_list.X_item(i);
    goal_type = goal_list.goal_type(i);
    
    % Find where O_item appears as target/new
    O_idx = find(strcmp(final_test_list.stimulus_id, O_item) & ...
                 strcmp(final_test_list.role, "target") & ...
                 strcmp(final_test_list.trial_type, "new"), 1);
    
    if isempty(O_idx)
        fprintf('  ERROR: Goal %d O_item (%s) never appeared!\n', i, O_item);
        errors = errors + 1;
        continue;
    end
    
    % X_item should appear at O_idx + 2
    expected_X_idx = O_idx + 2;
    
    if expected_X_idx > height(final_test_list)
        fprintf('  ERROR: Goal %d X_item would be out of bounds!\n', i);
        errors = errors + 1;
        continue;
    end
    
    actual_X_item = final_test_list.stimulus_id(expected_X_idx);
    actual_X_type = final_test_list.trial_type(expected_X_idx);
    
    if ~strcmp(actual_X_item, X_item)
        fprintf('  ERROR: Goal %d - Expected X=%s at pos %d, got %s\n', ...
                i, X_item, expected_X_idx, actual_X_item);
        errors = errors + 1;
    elseif ~strcmp(actual_X_type, goal_type)
        fprintf('  ERROR: Goal %d - Expected type=%s, got %s\n', ...
                i, goal_type, actual_X_type);
        errors = errors + 1;
    end
end

if errors == 0
    fprintf('✓ All 90 goals have correct 2-back structure!\n');
else
    fprintf('✗ Found %d errors in 2-back structure\n', errors);
end

% Check 3: No O_item should appear as target/new more than once
fprintf('\n--- Checking for duplicate goal starts ---\n');
[unique_vals, ~, idx] = unique(targets_new.stimulus_id);
counts = accumarray(idx, 1);
duplicates = unique_vals(counts > 1);
if isempty(duplicates)
    fprintf('✓ No O-item started multiple goals\n');
else
    fprintf('✗ Found %d O-items that started multiple goals:\n', length(duplicates));
    for d = 1:length(duplicates)
        fprintf('  - %s appeared %d times\n', duplicates(d), counts(find(strcmp(unique_vals, duplicates(d)))));
    end
end

% Check 4: Show first 20 trials for manual inspection
fprintf('\n--- First 20 trials (manual inspection) ---\n');
fprintf('%-5s %-25s %-20s %-10s %-15s %-10s\n', 'Idx', 'Stimulus', 'Condition', 'Role', 'Trial Type', 'Response');
fprintf('%s\n', repmat('-', 1, 95));
for i = 1:min(20, height(final_test_list))
    fprintf('%-5d %-25s %-20s %-10s %-15s %-10s\n', i, ...
            final_test_list.stimulus_id(i), ...
            final_test_list.condition(i), ...
            final_test_list.role(i), ...
            final_test_list.trial_type(i), ...
            final_test_list.correct_response(i));
end

fprintf('\n=== END SANITY CHECK ===\n\n');