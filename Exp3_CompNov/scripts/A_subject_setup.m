%==========================================================================
%         Generate Subject-Specific Task Variables 
%==========================================================================
% Author: Zihan Bai, zihan.bai@nyu.edu, Michelmann Lab at NYU
% note: you must run this for each subject before running the task
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
sequence_1_back = [];
sequence_2_back = [];
goal_list_full = [];

% response counters
p.counts.oneback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);
p.counts.twoback = struct('resp_same', 0, 'resp_similar', 0, 'resp_new', 0);

% directory
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end
p.subj_dir = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
if ~exist(p.subj_dir, 'dir'), mkdir(p.subj_dir); end

%% P1: Setup
p.nComparison = 120;
p.nNovel = 120;
p.nTotalPairs = p.nComparison + p.nNovel; % 240 total pairs

% # of blocks in experiment
% 4 blocks -> 30 compared + 30 novel pairs per block, which divides exactly
% into 10/10/10 goal types per condition (no remainder to spread around).
% Runs come out ~6.2 min of 1-back and ~5.2 min of 2-back.
p.nBlocks = 4;

% keyboard mappings
p.keys.same = 'j';
p.keys.diff = 'k';
p.keys.quit = 'escape';

p.timing.image_dur = 1.5;           % stimulus presentation
p.timing.fix_dur = 1.5;             % base fixation (mean ITI)
p.timing.fix_jitter = 0.75;         % uniform ±0.75s -> ITI 0.75 to 2.25s

% Sequence aperiodicity control (see P5B). For each block we build
% n_candidates complete 2-back sequences and keep the one whose
% compared-trial train is least rhythmic, rather than stopping at the first
% acceptable one (which tends to land just under threshold). period_thresh
% is the 95th percentile of a 500-permutation null for blocks of this
% length, and is used only to warn if even the best candidate is rhythmic.
p.seq.n_candidates = 30;
p.seq.period_thresh = 11;
p.seq.n_perm = 500;      % permutations for the built-in validation
p.seq.alpha = 0.05;      % a block fails if p < alpha

% Fixation-only lead-in at the start of every phase (1-back, 2-back and
% recognition), replacing the junk trials that used to open each block.
% Honoured by C_run_1_back.m, D_run_2_back.m and E_run_recognition.m,
% which mark its onset with an Eyelink 'LEAD_IN_ONSET' message.
p.timing.block_lead_in = 10;        % seconds before the first trial of a run

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
% Pool is still split into thirds so that per-condition counts match the
% 3-condition design; the third chunk (formerly ISOLATED) goes unused.
n_cond_l1 = height(master_pair_list_l1) / 3; % 180 / 3 = 60
n_cond_l2 = height(master_pair_list_l2) / 3; % 180 / 3 = 60

% --- Split L1 Pool (180 pairs -> 60 comp / 60 nov / 60 unused) ---
final_list_l1 = master_pair_list_l1(randperm(height(master_pair_list_l1)), :);

comp_l1 = final_list_l1(1:n_cond_l1, :);
nov_l1  = final_list_l1(n_cond_l1+1 : 2*n_cond_l1, :);

% --- Split L2 Pool (180 pairs -> 60 comp / 60 nov / 60 unused) ---
final_list_l2 = master_pair_list_l2(randperm(height(master_pair_list_l2)), :);

comp_l2 = final_list_l2(1:n_cond_l2, :);
nov_l2  = final_list_l2(n_cond_l2+1 : 2*n_cond_l2, :);

% --- Combine pools to create final 120-pair condition lists ---
comp_pairs = [comp_l1; comp_l2];
novel_pairs = [nov_l1; nov_l2];

% Shuffle the final lists so L1/L2 pairs are mixed
comp_pairs = comp_pairs(randperm(height(comp_pairs)), :);
novel_pairs = novel_pairs(randperm(height(novel_pairs)), :);

% store condition assignments
p.stim.comp_l1 = comp_l1;
p.stim.comp_l2 = comp_l2;
p.stim.novel_l1 = nov_l1;
p.stim.novel_l2 = nov_l2;

% Store the final combined lists
p.stim.compared = comp_pairs;
p.stim.novel = novel_pairs;

fprintf('  -> Compared: %d L1, %d L2 (Total %d)\n', height(comp_l1), height(comp_l2), height(comp_pairs));
fprintf('  -> Novel:    %d L1, %d L2 (Total %d)\n', height(nov_l1), height(nov_l2), height(novel_pairs));

%% P4: Split stimuli into blocks
partition_idx = @(N, nblocks) arrayfun(@(k) ...
    ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), ...
    1:nblocks, 'UniformOutput', false);

n_comp_l1 = height(p.stim.comp_l1); % 60
n_nov_l1 = height(p.stim.novel_l1); % 60

n_comp_l2 = height(p.stim.comp_l2); % 60
n_nov_l2 = height(p.stim.novel_l2); % 60

% --- Partition L1 pools into blocks ---
p.block_indices.comp_l1 = partition_idx(n_comp_l1, p.nBlocks);
p.block_indices.nov_l1  = partition_idx(n_nov_l1, p.nBlocks);

% --- Partition L2 pools into blocks ---
p.block_indices.comp_l2 = partition_idx(n_comp_l2, p.nBlocks);
p.block_indices.nov_l2  = partition_idx(n_nov_l2, p.nBlocks);

fprintf('  -> Each block gets: %d L1-Comp, %d L2-Comp, %d L1-Nov, %d L2-Nov\n', ...
    numel(p.block_indices.comp_l1{1}), numel(p.block_indices.comp_l2{1}), ...
    numel(p.block_indices.nov_l1{1}),  numel(p.block_indices.nov_l2{1}));

% --- Goal-type templates ---
% Per-block pair counts need not be divisible by 3 (e.g. 40 pairs/block at
% nBlocks=3), so goal types are dealt from a session-long cyclic template
% rather than assigned within each block. Cycling A-B / A-A / A-N over the
% full 120-pair pool guarantees exactly 40 of each type per condition
% across the session, while spreading the remainder across blocks as
% 14/13/13, 13/14/13, 13/13/14.
%
% The cycle starts on A-B deliberately and is NOT shifted: the leftover
% goal goes to A-B in block 1, so the goal type that feeds the recognition
% test is over-weighted toward the earliest, highest-quality block.
goal_type_cycle = repmat(["A-B"; "A-A"; "A-N"], p.nComparison/3, 1);
p.goal_cycle.compared = goal_type_cycle;
p.goal_cycle.novel    = goal_type_cycle;

p.block_indices.goal_comp = partition_idx(p.nComparison, p.nBlocks);
p.block_indices.goal_nov  = partition_idx(p.nNovel, p.nBlocks);

%% P5: Build sequence
all_foils_remain = all_foil_pairs;

for b = 1:p.nBlocks
    fprintf('\nBlock %d \n', b);

    %%% Assign block-specific stimuli
    comp_l1_b = p.stim.comp_l1(p.block_indices.comp_l1{b}, :);
    comp_l2_b = p.stim.comp_l2(p.block_indices.comp_l2{b}, :);
    comp_pairs_b = [comp_l1_b; comp_l2_b];
    comp_pairs_b = comp_pairs_b(randperm(height(comp_pairs_b)), :);

    nov_l1_b = p.stim.novel_l1(p.block_indices.nov_l1{b}, :);
    nov_l2_b = p.stim.novel_l2(p.block_indices.nov_l2{b}, :);
    novel_pairs_b = [nov_l1_b; nov_l2_b];
    novel_pairs_b = novel_pairs_b(randperm(height(novel_pairs_b)), :);

    nComp_b = height(comp_pairs_b);
    nNov_b = height(novel_pairs_b);

    %% P5A: Build sequence 1-back
    comp_miniblocks = {};
    repeat_miniblocks = {};

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

    % Combine and shuffle
    miniblocks_1_back = [comp_miniblocks, repeat_miniblocks];
    shuffled_indices = randperm(numel(miniblocks_1_back));
    final_miniblocks = miniblocks_1_back(shuffled_indices);
    final_1_back_list = vertcat(final_miniblocks{:});

    % Convert to table
    sequence_1_back_block = cell2table(final_1_back_list, ...
        'VariableNames', {'stim_id', 'condition', 'identity', 'corr_resp'});

    sequence_1_back_block.block = repmat(b, height(sequence_1_back_block), 1);

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
    block_pairs = [comp_pairs_b; novel_pairs_b];
    block_conditions = [repmat("compared", nComp_b, 1); ...
                        repmat("novel", nNov_b, 1)];
    nGoals_b = nComp_b + nNov_b; % 60
    goal_list = table('Size', [nGoals_b, 5], ...
        'VariableTypes', {'string', 'string', 'string', 'string', 'string'}, ...
        'VariableNames', {'A', 'B', 'condition', 'goal_type', 'X'});
    goal_list.A = block_pairs.A;
    goal_list.B = block_pairs.B;
    goal_list.condition = block_conditions;

    % assign goal types (as close to 1/3 each PER CONDITION as the block
    % size allows; exactly 1/3 each per condition across the session)
    goal_types_comp = p.goal_cycle.compared(p.block_indices.goal_comp{b});
    goal_types_nov  = p.goal_cycle.novel(p.block_indices.goal_nov{b});

    fprintf('  Goal types this block -- compared %d/%d/%d, novel %d/%d/%d (A-B/A-A/A-N)\n', ...
        sum(goal_types_comp == "A-B"), sum(goal_types_comp == "A-A"), sum(goal_types_comp == "A-N"), ...
        sum(goal_types_nov  == "A-B"), sum(goal_types_nov  == "A-A"), sum(goal_types_nov  == "A-N"));

    % Get indices for each condition
    idx_comp = 1:nComp_b;
    idx_nov = nComp_b+1 : nGoals_b;

    % Assign for 'compared' (shuffled so goal type is unrelated to pair)
    goal_list.goal_type(idx_comp) = goal_types_comp(randperm(nComp_b));

    % Assign for 'novel'
    goal_list.goal_type(idx_nov) = goal_types_nov(randperm(nNov_b));

    % assign X for A-A and A-B
    for i = 1:nGoals_b
        if goal_list.goal_type(i) == "A-B"
            goal_list.X(i) = goal_list.B(i);
        elseif goal_list.goal_type(i) == "A-A"
            goal_list.X(i) = goal_list.A(i);
        else
            goal_list.X(i) = ""; % assigned dynamically
        end
    end

    % ---- Search for an aperiodic compared-trial sequence ---------------
    % Condition is assigned independently of serial position, so *on
    % average* compared trials land aperiodically -- but any single shuffle
    % can still come out rhythmic, and a periodic trial-of-interest
    % sequence risks aliasing with drift and physiological noise. So build
    % n_candidates complete sequences and keep the one whose compared-trial
    % train has the flattest spectrum.
    goal_list_pool = goal_list;         % pre-shuffle copy, reused per attempt
    foils_snapshot = all_foils_remain;  % rewind point for discarded attempts
    best_stat = Inf; best = struct();

    for attempt = 1:p.seq.n_candidates
    all_foils_remain = foils_snapshot;

    % shuffle goal order
    goal_list = goal_list_pool(randperm(height(goal_list_pool)), :);

    % Keep A-N goals out of the last few slots. An A-N goal completes by
    % pulling in the *next unstarted* goal as its "new" item, so an A-N
    % sitting at the very end of the list has nothing left to chain to and
    % has to burn a filler foil; it also leaves a single straggler goal in
    % the buffer, which forces a padding trial. Swapping those A-N goals
    % back into the body of the list removes both sources of junk.
    n_tail = 3;
    tail_slots = height(goal_list) - n_tail + 1 : height(goal_list);
    is_AN_tail = find(goal_list.goal_type(tail_slots) == "A-N");
    if ~isempty(is_AN_tail)
        body_non_AN = find(goal_list.goal_type(1:end-n_tail) ~= "A-N");
        picks = body_non_AN(randperm(numel(body_non_AN), numel(is_AN_tail)));
        rows_tail = tail_slots(is_AN_tail);
        tmp = goal_list(rows_tail, :);
        goal_list(rows_tail, :) = goal_list(picks, :);
        goal_list(picks, :) = tmp;
    end

    %%% Build sequence
    % No lead-in junk trials: the block opens with p.timing.block_lead_in
    % seconds of blank screen instead, which serves the same purpose (let
    % the BOLD/pupil signal settle before the first real trial) without
    % spending stimuli. The first two trials are goal-starting 'A' items
    % with corr_resp "none", which is the correct response anyway since
    % nothing has appeared 2 back yet.
    sequence = cell(300, 5); % Oversized to prevent crash
    row_idx = 1;
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

                % Find the properties of the NEW item 'X'
                X_props = goal_list(strcmp(goal_list.A, X), :);
                if isempty(X_props)
                    log_condition = "junk_foil";
                    log_goal_type = "JUNK";
                    X_identity = "J"; % 'J' for Junk
                else
                    % Log the NEW item's OWN properties
                    log_condition = X_props.condition(1);
                    log_goal_type = X_props.goal_type(1); % This logs X's OWN goal
                    X_identity = "A";
                end
                X_resp = "none";

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
            % Can we start a new goal?
            while goal_pointer <= height(goal_list) && goals_started(goal_pointer)
                goal_pointer = goal_pointer + 1; % find next unstarted
            end

            if goal_pointer <= height(goal_list)
                % YES: Start the next goal
                goal = goal_list(goal_pointer, :);
                sequence(row_idx,:) = {goal.A, goal.condition, "A", goal.goal_type, "none"};
                goals_started(goal_pointer) = true;
                active_goals(end+1, :) = [goal_pointer, row_idx];
                row_idx = row_idx + 1;
                goal_pointer = goal_pointer + 1;
            else
                % NO: All goals are started, but buffer is not empty.
                % We MUST add a junk trial to advance time.
                if isempty(active_goals)
                    % Only hit if all goals started AND the buffer is empty.
                    break;
                end

                % Add a junk trial to pad
                if height(all_foils_remain) < 1
                    error('Out of foils for junk padding in block %d', b);
                end
                junk = all_foils_remain(1,:); all_foils_remain(1,:) = [];
                sequence(row_idx,:) = {junk.A_foil, "padding_junk", "J", "JUNK", "none"};
                row_idx = row_idx + 1; % Advance time
            end
        end
    end

    % ---- Score this attempt: how periodic is the compared-trial train? --
    % Binary indicator of compared trials over the block; peak / median
    % power across all non-DC frequencies of its DFT. A random order scores
    % ~6-7; a rhythmic one scores well above.
    cond_col = string(sequence(1:row_idx-1, 2));
    stat = seq_peak_ratio(cond_col == "compared");

    if stat < best_stat
        best_stat = stat;
        best.sequence = sequence;
        best.row_idx = row_idx;
        best.goal_list = goal_list;
        best.foils = all_foils_remain;
        best.indicator = (cond_col == "compared");
    end
    end % attempt

    % Take the flattest candidate
    sequence = best.sequence;
    row_idx = best.row_idx;
    goal_list = best.goal_list;
    all_foils_remain = best.foils;

    % ---- Validate the retained sequence against a permutation null -----
    % Same trials, random order, p.seq.n_perm times. If the real sequence
    % scores no higher than the shuffles, it carries no dominant rhythm.
    null_stats = zeros(p.seq.n_perm, 1);
    for it = 1:p.seq.n_perm
        null_stats(it) = seq_peak_ratio(best.indicator(randperm(numel(best.indicator))));
    end
    p.seq.block_stat(b)   = best_stat;
    p.seq.block_p(b)      = mean(null_stats >= best_stat);
    p.seq.block_null95(b) = prctile(null_stats, 95);
    p.seq.block_nullmean(b) = mean(null_stats);
    p.seq.block_pass(b)   = p.seq.block_p(b) >= p.seq.alpha;

    if p.seq.block_pass(b), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  Aperiodicity: %.2f (null M=%.2f, 95th=%.2f, p=%.3f) -> %s\n', ...
        best_stat, mean(null_stats), prctile(null_stats, 95), p.seq.block_p(b), verdict);
    if ~p.seq.block_pass(b)
        warning(['Block %d compared-trial sequence is periodic ' ...
                 '(stat %.2f, p=%.3f) despite %d candidates.'], ...
                 b, best_stat, p.seq.block_p(b), p.seq.n_candidates);
    end

    goal_list_b = goal_list; % Make a copy for the full list
    goal_list_b.block = repmat(b, height(goal_list_b), 1);
    goal_list_full = [goal_list_full; goal_list_b];

    % No trailing junk trials: the loop above already runs until every goal
    % has both started AND completed, so the block ends on a real trial.

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


%% ========================================================================
%  P6: BUILD FINAL RECOGNITION TASK
%  ========================================================================

fprintf('\nBuilding final recognition task...\n');

% Test the A-B goals from BOTH conditions: these are the pairs where the
% subject actually saw A and then B in the 2-back, so either item can serve
% as the old probe. Exactly nComparison/3 = 40 A-B goals per condition
% across the session -> 40 compared + 40 novel.
tested_goals = goal_list_full(goal_list_full.goal_type == "A-B", :);

n_tested = height(tested_goals);
n_expected = p.nComparison / 3;
n_old_comp = sum(tested_goals.condition == "compared");
n_old_nov  = sum(tested_goals.condition == "novel");
assert(n_old_comp == n_expected && n_old_nov == n_expected, ...
    'Expected %d compared + %d novel old items, got %d + %d', ...
    n_expected, n_expected, n_old_comp, n_old_nov);

% --- Pick A or B from each pair, half of each ---
% Which pairs give A is random; only the counts are fixed, so each
% condition contributes an equal number of A and B probes.
selected_items = strings(n_tested, 1);
selected_identity = strings(n_tested, 1);

conds = unique(tested_goals.condition);
for c = 1:numel(conds)
    idx_c = find(tested_goals.condition == conds(c));
    n_c = numel(idx_c);

    lbl = repmat("B", n_c, 1);
    lbl(1:floor(n_c/2)) = "A";
    lbl = lbl(randperm(n_c));

    for k = 1:n_c
        r = idx_c(k);
        selected_identity(r) = lbl(k);
        if lbl(k) == "A"
            selected_items(r) = tested_goals.A(r);
        else
            selected_items(r) = tested_goals.B(r);
        end
    end
end

all_old_items = table(selected_items, tested_goals.condition, ...
    selected_identity, repmat("old", n_tested, 1), ...
    repmat(p.keys.same, n_tested, 1), ...
    'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'corr_resp'});
all_old_items.trial_type = repmat("old", n_tested, 1);
all_old_items.corr_resp = repmat(p.keys.same, n_tested, 1);

% --- Pool 2b: The 'Foil' Pool (matched 1:1 to the old-item pool -> 80) ---
n_rec_foils = n_tested;
assert(height(all_foils_remain) >= n_rec_foils, ...
    'Not enough remaining foils for recognition task! Need %d, have %d', ...
    n_rec_foils, height(all_foils_remain));

new_foils_list = all_foils_remain.A_foil(1:n_rec_foils);
all_foils_remain(1:n_rec_foils, :) = []; % Remove them

all_new_foils = table(new_foils_list, repmat("foil", n_rec_foils, 1), ...
    repmat("N", n_rec_foils, 1), repmat("new", n_rec_foils, 1), ...
    repmat(p.keys.diff, n_rec_foils, 1), ...
    'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'corr_resp'});

% --- 3. Combine, shuffle, and add to 'p' struct ---
sequence_recognition = [all_old_items; all_new_foils];
sequence_recognition = sequence_recognition(randperm(height(sequence_recognition)), :);

fprintf('  -> Recognition: %d compared + %d novel old, %d foils (%d trials total)\n', ...
    n_old_comp, n_old_nov, n_rec_foils, height(sequence_recognition));
for c = 1:numel(conds)
    id_c = all_old_items.identity(all_old_items.condition == conds(c));
    fprintf('     %-9s old items: %d A / %d B\n', ...
        conds(c), sum(id_c == "A"), sum(id_c == "B"));
end

%% ========================================================================
%  P7: ADD JITTER & SAVE OUTPUT
%  ========================================================================
fprintf('\nAdding jittered fixations and saving all schedules...\n');

% --- Add jittered fixation durations ---
n_1_back_trials = height(sequence_1_back);
n_2_back_trials = height(sequence_2_back);
n_rec_trials = height(sequence_recognition);

sequence_1_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_1_back_trials, 1) * 2 - 1) * p.timing.fix_jitter;

sequence_2_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_2_back_trials, 1) * 2 - 1) * p.timing.fix_jitter;

sequence_recognition.fix_duration = repmat(0.5, n_rec_trials, 1);

% --- Add subj_id to all schedules ---
sequence_1_back.subj_id = repmat(subj_id, n_1_back_trials, 1);
sequence_2_back.subj_id = repmat(subj_id, n_2_back_trials, 1);
sequence_recognition.subj_id = repmat(subj_id, n_rec_trials, 1);


subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.sequence_1_back = sequence_1_back;
subject_data.sequence_2_back = sequence_2_back;
subject_data.sequence_recognition = sequence_recognition;

%% ========================================================================
%  P8: SEQUENCE VALIDATION SUMMARY
%  ========================================================================
% The schedule is only written if every block's compared-trial sequence
% passed the permutation test, so a saved file is always ready to run.
fprintf('\n===== Aperiodicity check (compared trials, %d permutations) =====\n', p.seq.n_perm);
fprintf('  %-7s %-10s %-11s %-11s %-8s %s\n', 'block', 'observed', 'null mean', 'null 95th', 'p', 'verdict');
for b = 1:p.nBlocks
    if p.seq.block_pass(b), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-7d %-10.2f %-11.2f %-11.2f %-8.3f %s\n', b, p.seq.block_stat(b), ...
        p.seq.block_nullmean(b), p.seq.block_null95(b), p.seq.block_p(b), verdict);
end

if ~all(p.seq.block_pass)
    error(['Aperiodicity check FAILED for block(s) %s -- schedule NOT saved. ' ...
           'Re-run to draw a new sequence, or raise p.seq.n_candidates.'], ...
           mat2str(find(~p.seq.block_pass)));
end
fprintf('  All %d blocks passed (observed %.2f-%.2f, all p >= %.3f).\n', ...
    p.nBlocks, min(p.seq.block_stat), max(p.seq.block_stat), min(p.seq.block_p));

subject_data.parameters = p;   % refresh: p gained the validation results

% --- Save the file ---
save(output_filename, 'subject_data');
fprintf('\nSetup saved to: %s\n', output_filename);

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------
function pk = seq_peak_ratio(ind)
% Peak-to-median power ratio of a binary trial indicator. Higher = more
% rhythmic. DC term dropped; only frequencies up to Nyquist are considered.
ind = double(ind(:));
ind = ind - mean(ind);
N = numel(ind);
P = abs(fft(ind)).^2 / N;
halfband = 2:floor(N/2);
pk = max(P(halfband)) / median(P(halfband));
end
