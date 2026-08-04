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
p.nIsolated_Both = 120;
p.nNovel = 120;
p.nRepeat = 120;   % real bin-4 pairs: A-A in 1-back, B is the MST lure
p.nTotalPairs = p.nComparison + p.nIsolated_Both + p.nNovel + p.nRepeat; % 480 total pairs

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
p.timing.block_lead_in = 10;        % blank-screen s before the first trial (replaces lead-in junk)
p.timing.block_tail    = 10;        % blank-screen s after the last trial (replaces end junk)

% 2-back sequence: build n_candidates full sequences per block and keep the one
% that spends zero foils on padding (=> FIXED trial count) and whose compared-
% trial train has the flattest spectrum (aperiodic; avoids aliasing with fMRI
% drift/physiological noise). Validated against a permutation null.
p.seq.n_candidates = 30;
p.seq.n_perm = 500;
p.seq.alpha = 0.05;

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

% --- Load experimental pairs. The two experimental levels keep the l1/l2
%     machinery names, but now source from lure bins 2 and 3:
%       l1 machinery  <-  mst_*_l2.png  (bin 2)
%       l2 machinery  <-  mst_*_l3.png  (bin 3)
all_A_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_A_l2.png'));
all_B_files_l1 = dir(fullfile(p.stim_dir, 'mst_*_B_l2.png'));
master_pair_list_l1 = table(sort(string({all_A_files_l1.name}')), ...
    sort(string({all_B_files_l1.name}')), 'VariableNames', {'A', 'B'});

all_A_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_A_l3.png'));
all_B_files_l2 = dir(fullfile(p.stim_dir, 'mst_*_B_l3.png'));
master_pair_list_l2 = table(sort(string({all_A_files_l2.name}')), ...
    sort(string({all_B_files_l2.name}')), 'VariableNames', {'A', 'B'});

% --- Load repeat pairs (bin 4, mst_*_l4.png): real A-B pairs. The A item is
%     shown twice (A-A, "same") in the 1-back; the pairmate B is never shown and
%     becomes the unseen lure in the post-task MST. ---
all_A_files_rep = dir(fullfile(p.stim_dir, 'mst_*_A_l4.png'));
all_B_files_rep = dir(fullfile(p.stim_dir, 'mst_*_B_l4.png'));
master_pair_list_repeat = table(sort(string({all_A_files_rep.name}')), ...
    sort(string({all_B_files_rep.name}')), 'VariableNames', {'A', 'B'});
master_pair_list_repeat = master_pair_list_repeat(randperm(height(master_pair_list_repeat)), :);

% Load foil pairs (the A item is used as a single new/foil item)
all_foil_A_files = dir(fullfile(p.stim_dir, 'mst_*_A_foil.png'));
all_foil_B_files = dir(fullfile(p.stim_dir, 'mst_*_B_foil.png'));
all_foil_pairs = table(string({all_foil_A_files.name}'), string({all_foil_B_files.name}'), ...
    'VariableNames', {'A_foil', 'B_foil'});
all_foil_pairs = all_foil_pairs(randperm(height(all_foil_pairs)), :);

fprintf('Found %d l2 pairs, %d l3 pairs, %d repeat(l4) pairs, %d foil pairs.\n', ...
    height(master_pair_list_l1), height(master_pair_list_l2), ...
    height(master_pair_list_repeat), height(all_foil_pairs));

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
p.stim.repeat = master_pair_list_repeat;   % 120 real bin-4 pairs (A-A in 1-back, B is MST lure)

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

% --- Partition repeat pool (120 real bin-4 pairs) into 4 blocks (30 each) ---
p.block_indices.repeat = partition_idx(height(p.stim.repeat), p.nBlocks);

fprintf('  -> Each block gets: 15 L1-Comp, 15 L2-Comp, 15 L1-Iso, etc.\n');

% Deterministic goal-type cycle: exactly 1/3 A-B / A-A / A-N per condition,
% partitioned so every block gets a balanced 10/10/10 share (fixed structure).
goal_type_cycle = repmat(["A-B"; "A-A"; "A-N"], p.nComparison/3, 1);
p.goal_cycle.compared = goal_type_cycle;
p.goal_cycle.isolated = goal_type_cycle;
p.goal_cycle.novel    = goal_type_cycle;
p.block_indices.goal_comp = partition_idx(p.nComparison, p.nBlocks);
p.block_indices.goal_iso  = partition_idx(p.nIsolated_Both, p.nBlocks);
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

    %%%% repeat (R-R): real bin-4 pairs, A shown twice (A-A, "same"). No foils
    %%%% appear in the 1-back; the pairmate B is held out for the post-task MST.
    repeat_pairs_b = p.stim.repeat(p.block_indices.repeat{b}, :);
    for i = 1:height(repeat_pairs_b)
        rep_item = repeat_pairs_b.A(i);
        repeat_miniblocks{end+1} = { ...
            rep_item, "repeat", "A", "none"; ...
            rep_item, "repeat", "A", "j"};
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

    % assign goal types deterministically: exactly 1/3 A-B/A-A/A-N per condition
    % this block (the per-block cycle slice yields a balanced 10/10/10).
    goal_types_comp = p.goal_cycle.compared(p.block_indices.goal_comp{b});
    goal_types_iso  = p.goal_cycle.isolated(p.block_indices.goal_iso{b});
    goal_types_nov  = p.goal_cycle.novel(p.block_indices.goal_nov{b});
    fprintf('  Goal types this block -- each condition %d/%d/%d (A-B/A-A/A-N)\n', ...
        sum(goal_types_comp=="A-B"), sum(goal_types_comp=="A-A"), sum(goal_types_comp=="A-N"));

    idx_comp = 1:nComp_b;
    idx_iso  = nComp_b+1 : nComp_b+nIso_b;
    idx_nov  = nComp_b+nIso_b+1 : 90;

    goal_list.goal_type(idx_comp) = goal_types_comp(randperm(numel(idx_comp)));
    goal_list.goal_type(idx_iso)  = goal_types_iso(randperm(numel(idx_iso)));
    goal_list.goal_type(idx_nov)  = goal_types_nov(randperm(numel(idx_nov)));

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
    
    % ---- Candidate search: build n_candidates full sequences and keep the one
    %      that spends zero foils on padding (=> fixed trial count) and whose
    %      compared-trial train is the most aperiodic. No init/end junk -- the
    %      run brackets each block with block_lead_in / block_tail blank screen.
    goal_list_pool = goal_list;         % pre-shuffle copy, reused per attempt
    foils_snapshot = all_foils_remain;  % rewind point for discarded attempts
    best_stat = Inf; best = struct(); best.foils_used = Inf;

    for attempt = 1:p.seq.n_candidates
    all_foils_remain = foils_snapshot;

    % shuffle goal order
    goal_list = goal_list_pool(randperm(height(goal_list_pool)), :);

    % Keep A-N goals out of the last few slots: an A-N there has no goal left to
    % chain to and would force a padding foil. Swap any such goals into the body.
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

    %%% Build sequence (no init junk; the first two trials are goal-starting A
    %%% items with corr_resp "none", correct since nothing is 2-back yet)
    sequence = cell(300, 5);
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
                
                % --- THIS IS THE FUCKING FIX ---
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
                    % This should only be hit if goal_pointer > 90 AND
                    % the buffer is empty. We are done.
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
    
    % ---- Score this attempt: padding foils (want 0), then aperiodicity ----
    cond_col = string(sequence(1:row_idx-1, 2));
    stat = seq_peak_ratio(cond_col == "compared");
    foils_used = height(foils_snapshot) - height(all_foils_remain);
    if foils_used < best.foils_used || (foils_used == best.foils_used && stat < best_stat)
        best_stat = stat;
        best.foils_used = foils_used;
        best.sequence = sequence;
        best.row_idx = row_idx;
        best.goal_list = goal_list;
        best.foils = all_foils_remain;
        best.indicator = (cond_col == "compared");
    end
    end % attempt

    % Retain the best candidate (zero padding if achieved; flattest within that)
    sequence = best.sequence;
    row_idx = best.row_idx;
    goal_list = best.goal_list;
    all_foils_remain = best.foils;
    p.seq.block_foils_used(b) = best.foils_used;
    if best.foils_used > 0
        warning(['Block %d 2-back still spends %d padding foil(s) after %d candidates; ' ...
                 'raise p.seq.n_candidates for a zero-foil (fixed-length) budget.'], ...
                 b, best.foils_used, p.seq.n_candidates);
    end

    % Validate aperiodicity against a permutation null (same trials, shuffled)
    null_stats = zeros(p.seq.n_perm, 1);
    for it = 1:p.seq.n_perm
        null_stats(it) = seq_peak_ratio(best.indicator(randperm(numel(best.indicator))));
    end
    p.seq.block_stat(b) = best_stat;
    p.seq.block_p(b)    = mean(null_stats >= best_stat);
    p.seq.block_pass(b) = p.seq.block_p(b) >= p.seq.alpha;
    if p.seq.block_pass(b), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  Aperiodicity: %.2f (null M=%.2f, p=%.3f) -> %s\n', ...
        best_stat, mean(null_stats), p.seq.block_p(b), verdict);

    % Record the RETAINED goal_list for the recognition test
    goal_list_b = goal_list;
    goal_list_b.block = repmat(b, height(goal_list_b), 1);
    goal_list_full = [goal_list_full; goal_list_b];

    % Convert to table (no trailing junk: the loop ends on a real trial)
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

% --- First, save the goal_list_full we just built ---
fprintf('\nBuilding final recognition task...\n');

tested_goals = goal_list_full(goal_list_full.goal_type ~= "A-N" & goal_list_full.condition ~= "novel", :);

n_tested = height(tested_goals);
selected_items = strings(n_tested, 1);
selected_identity = strings(n_tested, 1);

for i = 1:n_tested
    if rand() > 0.5
        selected_items(i) = tested_goals.A(i);
        selected_identity(i) = "A";
    else
        selected_items(i) = tested_goals.B(i);
        selected_identity(i) = "B";
    end
end

all_old_items = table(selected_items, tested_goals.condition, ...
    selected_identity, repmat("old", n_tested, 1), ...
    repmat(p.keys.same, n_tested, 1), ...
    'VariableNames', {'stim_id', 'condition', 'identity', 'trial_type', 'corr_resp'});
all_old_items.trial_type = repmat("old", n_tested, 1);
all_old_items.corr_resp = repmat(p.keys.same, n_tested, 1); 

% --- Pool 2b: The 'Foil' Pool (N=160) ---
n_rec_foils = 160;
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

%% ========================================================================
%  P6B: BUILD POST-TASK MST (old / similar / new) -- matches NeurWEM_plt
%  ========================================================================
% Standard MST on the bin-4 repeat pairs. Each repeat pair had its A shown twice
% (A-A) in the 1-back; its B was never shown. Half the repeat pairs are probed
% with the seen item A (old), half with the unseen pairmate B (lure); a matched
% set of fresh foils is added as new. One pair contributes EITHER its A (old) OR
% its B (lure), never both. Responses: j = old, k = similar, withhold = new.
fprintf('\nBuilding post-task MST...\n');
p.mst.n_old = 60; p.mst.n_lure = 60; p.mst.n_new = 60;
p.keys.mst_old = 'j'; p.keys.mst_similar = 'k';

rep = p.stim.repeat(randperm(height(p.stim.repeat)), :);
assert(height(rep) >= p.mst.n_old + p.mst.n_lure, ...
    'MST: need %d repeat pairs, have %d', p.mst.n_old + p.mst.n_lure, height(rep));
old_pairs  = rep(1 : p.mst.n_old, :);
lure_pairs = rep(p.mst.n_old + 1 : p.mst.n_old + p.mst.n_lure, :);

assert(height(all_foils_remain) >= p.mst.n_new, ...
    'MST: need %d new foils, have %d', p.mst.n_new, height(all_foils_remain));
new_items = all_foils_remain.A_foil(1 : p.mst.n_new);
all_foils_remain(1 : p.mst.n_new, :) = [];

mst_old  = table(old_pairs.A,  repmat("old",  p.mst.n_old,  1), ...
    repmat(string(p.keys.mst_old),     p.mst.n_old,  1), ...
    'VariableNames', {'stim_id','trial_type','corr_resp'});
mst_lure = table(lure_pairs.B, repmat("lure", p.mst.n_lure, 1), ...
    repmat(string(p.keys.mst_similar), p.mst.n_lure, 1), ...
    'VariableNames', {'stim_id','trial_type','corr_resp'});
mst_new  = table(new_items,    repmat("new",  p.mst.n_new,  1), ...
    repmat("none",                     p.mst.n_new,  1), ...
    'VariableNames', {'stim_id','trial_type','corr_resp'});

sequence_mst = [mst_old; mst_lure; mst_new];
sequence_mst = sequence_mst(randperm(height(sequence_mst)), :);
fprintf('  -> MST: %d old, %d lure, %d new (%d trials).\n', ...
    p.mst.n_old, p.mst.n_lure, p.mst.n_new, height(sequence_mst));

%% ========================================================================
%  P7: ADD JITTER & SAVE OUTPUT
%  ========================================================================
fprintf('\nAdding jittered fixations and saving all schedules...\n');

% --- Add jittered fixation durations ---
n_1_back_trials = height(sequence_1_back);
n_2_back_trials = height(sequence_2_back);
n_rec_trials = height(sequence_recognition);
n_mst_trials = height(sequence_mst);

sequence_1_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_1_back_trials, 1) * 2 - 1) * p.timing.fix_jitter;

sequence_2_back.fix_duration = ...
    p.timing.fix_dur + (rand(n_2_back_trials, 1) * 2 - 1) * p.timing.fix_jitter;

sequence_recognition.fix_duration = repmat(0.5, n_rec_trials, 1);
sequence_mst.fix_duration = repmat(0.5, n_mst_trials, 1);

% --- Add subj_id to all schedules ---
sequence_1_back.subj_id = repmat(subj_id, n_1_back_trials, 1);
sequence_2_back.subj_id = repmat(subj_id, n_2_back_trials, 1);
sequence_recognition.subj_id = repmat(subj_id, n_rec_trials, 1);
sequence_mst.subj_id = repmat(subj_id, n_mst_trials, 1);


subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.sequence_1_back = sequence_1_back;
subject_data.sequence_2_back = sequence_2_back;
subject_data.sequence_recognition = sequence_recognition;
subject_data.sequence_mst = sequence_mst;

% --- Save the file ---
save(output_filename, 'subject_data');
fprintf('\nSetup saved to: %s\n', output_filename);

%% ========================================================================
%  Local functions
%  ========================================================================
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