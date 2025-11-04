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

% directory
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_processed/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end

%% ========================================================================
%  SECTION 1: SUBJECT ID AND PARAMETERS
%  ========================================================================
% # of stimulus pairs per condition (total across all blocks)
p.nComparison = 40;
p.nIsolated_Both = 40;
p.nNovel = 40;

% test phase trial counts
p.nTest_Lure_Events = 120;
p.nTest_Repeat_Events = 120;

% # of blocks in experiment
p.nBlocks = 4;

% keyboard mappings
p.keys.same = 'j';
p.keys.diff = 'k';
p.keys.quit = 'escape';

% timing parameters in seconds
p.timing.image_dur = 1.5;           % stimulus presentation
p.timing.fix_dur = 0.75;            % base fixation
p.timing.fix_jitter = 0.25;         % jitter range: ±0.25s (so 0.5 to 1.0s total)

%% ========================================================================
%  SECTION 2: LOAD STIMULI
%  ========================================================================

output_filename = fullfile(p.setup_dir, sprintf('sub%03d_setup.mat', subj_id));

if exist(output_filename, 'file')
    overwrite = input(sprintf('setup for subject %03d exists. overwrite? (y/n): ', subj_id), 's');
    if ~strcmpi(overwrite, 'y')
        fprintf('aborted.\n');
        return;
    end
end

fprintf('loading stimuli...\n');

all_targ_files = dir(fullfile(p.stim_dir, 'mst_*_targ_l1.png'));
all_lure_files = dir(fullfile(p.stim_dir, 'mst_*_lure_l1.png'));

assert(numel(all_targ_files) == 120, 'expected 120 target files!');

targ_names = string({all_targ_files.name}');
lure_names = string({all_lure_files.name}');

master_pair_list = table(sort(targ_names), sort(lure_names), ...
                         'VariableNames', {'Img_A', 'Img_B'});

all_foil_files = dir(fullfile(p.stim_dir, 'mst_*_foil.png'));
all_foils = string({all_foil_files.name}');

fprintf('Found %d experimental pairs and %d foils.\n', ...
        height(master_pair_list), numel(all_foils));

%% ========================================================================
%  SECTION 3: COUNTERBALANCE TARGET-LURE ASSIGNMENT
%  ========================================================================

shuffled_pairs = master_pair_list(randperm(height(master_pair_list)), :);

final_pairs = table('Size', size(shuffled_pairs), ...
                    'VariableTypes', {'string', 'string'}, ...
                    'VariableNames', {'Target', 'Lure'});

% randomly assign which image is target vs lure for each pair
for i = 1:height(shuffled_pairs)
    if rand() > 0.5
        final_pairs(i, :) = shuffled_pairs(i, :);
    else
        final_pairs(i, :) = [shuffled_pairs(i, 2), shuffled_pairs(i, 1)];
    end
end

% split pairs into three conditions
idx_start = 1;
idx_end = p.nComparison;
comparison_pairs = final_pairs(idx_start:idx_end, :);

idx_start = idx_end + 1;
idx_end = idx_start + p.nIsolated_Both - 1;
iso_both_pairs = final_pairs(idx_start:idx_end, :);

idx_start = idx_end + 1;
idx_end = idx_start + p.nNovel - 1;
novel_pairs = final_pairs(idx_start:idx_end, :);

% shuffle foils
all_foils = all_foils(randperm(numel(all_foils)));

% store condition assignments
p.stim.comparison = comparison_pairs;
p.stim.iso_both = iso_both_pairs;
p.stim.novel = novel_pairs;


%% ========================================================================
%  SECTION 4: PARTITION STIMULI INTO BLOCKS
%  ========================================================================
% divide each condition's stimuli evenly across the 3 blocks

% helper function: split N items into nblocks approximately equal chunks
partition_idx = @(N, nblocks) arrayfun(@(k) ...
    ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), ...
    1:nblocks, 'UniformOutput', false);

% get index ranges for each block
comp_idx_blocks = partition_idx(p.nComparison, p.nBlocks);
iso_idx_blocks = partition_idx(p.nIsolated_Both, p.nBlocks);
nov_idx_blocks = partition_idx(p.nNovel, p.nBlocks);


%% ========================================================================
%  SECTION 5: BUILD EACH BLOCK
%  ========================================================================

fprintf('Building %d blocks...\n', p.nBlocks);

% tables to accumulate all trials across blocks
encoding_schedule_all = table();
test_schedule_all = table();

% response counters
p.counts.encoding = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);
p.counts.test = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);

for b = 1:p.nBlocks
    fprintf('--- Block %d ---\n', b);
    
    % =====================================================================
    % 5A: SELECT BLOCK-SPECIFIC PAIRS
    % =====================================================================
    
    comparison_pairs_b = p.stim.comparison(comp_idx_blocks{b}, :);
    iso_both_pairs_b = p.stim.iso_both(iso_idx_blocks{b}, :);
    novel_pairs_b = p.stim.novel(nov_idx_blocks{b}, :);
    
    nComp_b = height(comparison_pairs_b);
    nIso_b = height(iso_both_pairs_b);
    nNov_b = height(novel_pairs_b);
    
    % =====================================================================
    % 5B: BUILD ENCODING PHASE FOR THIS BLOCK
    % =====================================================================
    
    encoding_blocks = {};  % experimental trial blocks
    spacer_blocks = {};    % filler trial blocks
    
    % --- comparison condition encoding ---
    % each block: target, then lure (requires 'k' response)
    for i = 1:nComp_b
        encoding_blocks{end+1} = { ...
            comparison_pairs_b.Target(i), "comparison", "target", "new", "none"; ...
            comparison_pairs_b.Lure(i), "comparison", "lure", "lure", "k"};
    end
    
    % --- isolated-both condition encoding ---
    % each stimulus appears in separate block with filler
    iso_event_foils = all_foils(1:2*nIso_b);
    all_foils(1:2*nIso_b) = [];
    
    for i = 1:nIso_b
        % target block
        encoding_blocks{end+1} = { ...
            iso_event_foils(1), "iso_both", "filler", "new", "none"; ...
            iso_both_pairs_b.Target(i), "iso_both", "target", "new", "none"};
        iso_event_foils(1) = [];
        
        % lure block
        encoding_blocks{end+1} = { ...
            iso_event_foils(1), "iso_both", "filler", "new", "none"; ...
            iso_both_pairs_b.Lure(i), "iso_both", "lure", "new", "none"};
        iso_event_foils(1) = [];
    end
    
    % --- spacer blocks: foil repeats (provide 'j' responses) ---
    n_foil_repeats = nComp_b;
    repeat_foils = all_foils(1:n_foil_repeats);
    all_foils(1:n_foil_repeats) = [];
    
    for i = 1:n_foil_repeats
        spacer_blocks{end+1} = { ...
            repeat_foils(i), "foil_repeat", "filler", "new", "none"; ...
            repeat_foils(i), "foil_repeat", "filler", "repeat", "j"};
    end
    
    % --- spacer blocks: new foils (fill remaining slots) ---
    n_new_spacers = numel(encoding_blocks) - numel(spacer_blocks);
    
    new_spacer_foils_1 = all_foils(1:n_new_spacers);
    all_foils(1:n_new_spacers) = [];
    
    new_spacer_foils_2 = all_foils(1:n_new_spacers);
    all_foils(1:n_new_spacers) = [];
    
    for i = 1:n_new_spacers
        spacer_blocks{end+1} = { ...
            new_spacer_foils_1(i), "foil_new", "filler", "new", "none"; ...
            new_spacer_foils_2(i), "foil_new", "filler", "new", "none"};
    end
    
    % --- interleave encoding and spacer blocks ---
    shuffled_encoding = encoding_blocks(randperm(numel(encoding_blocks)));
    shuffled_spacers = spacer_blocks(randperm(numel(spacer_blocks)));
    
    final_ordered_blocks = reshape([shuffled_spacers; shuffled_encoding], 1, []);
    final_encoding_list = vertcat(final_ordered_blocks{:});
    
    % convert to table
    encoding_schedule_block = cell2table(final_encoding_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    
    % add block number column
    encoding_schedule_block.block = repmat(b, height(encoding_schedule_block), 1);
    
    % aadd n-back verification columns
    encoding_schedule_block.nback_target_id = strings(height(encoding_schedule_block), 1);
    encoding_schedule_block.trial_type_final = encoding_schedule_block.trial_type_designed;
    
    % check for 1-back repeats
    for i = 2:height(encoding_schedule_block)
        encoding_schedule_block.nback_target_id(i) = ...
            encoding_schedule_block.stimulus_id(i-1);
        
        if string(encoding_schedule_block.stimulus_id(i)) == ...
           string(encoding_schedule_block.nback_target_id(i))
            encoding_schedule_block.trial_type_final(i) = "repeat";
            encoding_schedule_block.correct_response(i) = "j";
        end
    end
    
    % =====================================================================
    % 5C: ENSURE TARGET ENCODED BEFORE LURE (ISOLATED CONDITION ONLY)
    % =====================================================================
    
    all_iso_pairs_block = iso_both_pairs_b;
    
    for i = 1:height(all_iso_pairs_block)
        tgt = all_iso_pairs_block.Target(i);
        lur = all_iso_pairs_block.Lure(i);
        
        idx_tgt = find(encoding_schedule_block.stimulus_id == tgt, 1);
        idx_lur = find(encoding_schedule_block.stimulus_id == lur, 1);
        
        % swap if lure appears before target
        if ~isempty(idx_tgt) && ~isempty(idx_lur) && idx_lur < idx_tgt
            tmp = encoding_schedule_block(idx_tgt, :);
            encoding_schedule_block(idx_tgt, :) = encoding_schedule_block(idx_lur, :);
            encoding_schedule_block(idx_lur, :) = tmp;
        end
    end
    
    % append into full schedule
    encoding_schedule_all = [encoding_schedule_all; encoding_schedule_block];
    
    % update response counts
    p.counts.encoding.j_presses = p.counts.encoding.j_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'j'));
    p.counts.encoding.k_presses = p.counts.encoding.k_presses + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'k'));
    p.counts.encoding.no_response = p.counts.encoding.no_response + ...
        sum(strcmp(encoding_schedule_block.correct_response, 'none'));
    
    % =====================================================================
    % 5D: BUILD TEST PHASE FOR THIS BLOCK
    % =====================================================================
    
    phase2_crit_blocks = {};  % Lure triads (require 'k')
    phase2_safe_blocks = {};  % Repeat triads (require 'j')
    
    % --- critical triads: target → filler → lure ---
    n_crit_expected = nComp_b + nIso_b + nNov_b;
    
    assert(numel(all_foils) >= n_crit_expected, ...
           'not enough foils for critical triads');
    
    intervening_foils = all_foils(1:n_crit_expected);
    all_foils(1:n_crit_expected) = [];
    
    % comparison condition
    for i = 1:nComp_b
        phase2_crit_blocks{end+1} = { ...
            comparison_pairs_b.Target(i), "comparison", "target", "new", "none"; ...
            intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
            comparison_pairs_b.Lure(i), "comparison", "lure", "lure", "k"};
        intervening_foils(1) = [];
    end
    
    % isolated-both condition
    for i = 1:nIso_b
        phase2_crit_blocks{end+1} = { ...
            iso_both_pairs_b.Target(i), "isolated_both", "target", "new", "none"; ...
            intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
            iso_both_pairs_b.Lure(i), "isolated_both", "lure", "lure", "k"};
        intervening_foils(1) = [];
    end
    
    % novel condition
    for i = 1:nNov_b
        phase2_crit_blocks{end+1} = { ...
            novel_pairs_b.Target(i), "novel", "target", "new", "none"; ...
            intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
            novel_pairs_b.Lure(i), "novel", "lure", "lure", "k"};
        intervening_foils(1) = [];
    end
    
    % --- safe triads: foilA → foilB → foilA (2-back repeat) ---
    % calculate proportional number of repeat triads for this block
    total_pairs_block = nComp_b + nIso_b + nNov_b;
    proportion = total_pairs_block / (p.nComparison + p.nIsolated_Both + p.nNovel);
    n_repeat_triads_block = max(round(p.nTest_Repeat_Events * proportion), 1);
    
    n_needed_foils = 2 * n_repeat_triads_block;
    assert(numel(all_foils) >= n_needed_foils, ...
           'not enough foils for repeat triads');
    
    repeat_foil_A = all_foils(1:n_repeat_triads_block);
    repeat_foil_B = all_foils(n_repeat_triads_block+1 : 2*n_repeat_triads_block);
    all_foils(1:2*n_repeat_triads_block) = [];
    
    for i = 1:n_repeat_triads_block
        phase2_safe_blocks{end+1} = { ...
            repeat_foil_A(i), "foil_repeat", "filler", "new", "none"; ...
            repeat_foil_B(i), "foil_repeat", "filler", "new", "none"; ...
            repeat_foil_A(i), "foil_repeat", "filler", "repeat", "j"};
    end
    
    % --- interleave safe and critical triads ---
    n_crit = numel(phase2_crit_blocks);
    n_safe = numel(phase2_safe_blocks);
    
    shuffled_crit = phase2_crit_blocks(randperm(n_crit));
    shuffled_safe = phase2_safe_blocks(randperm(n_safe));
    
    L = max(n_crit, n_safe);
    interleaved_blocks = {};
    
    for i = 1:L
        if i <= n_safe
            interleaved_blocks{end+1} = shuffled_safe{i};
        end
        if i <= n_crit
            interleaved_blocks{end+1} = shuffled_crit{i};
        end
    end
    
    final_test_list = vertcat(interleaved_blocks{:});
    
    % convert to table
    test_schedule_block = cell2table(final_test_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    
    % add block number column
    test_schedule_block.block = repmat(b, height(test_schedule_block), 1);
    
    % add n-back verification columns
    test_schedule_block.nback_target_id = strings(height(test_schedule_block), 1);
    test_schedule_block.trial_type_final = test_schedule_block.trial_type_designed;
    
    % check for 2-back repeats (starting at trial 3)
    for i = 3:height(test_schedule_block)
        test_schedule_block.nback_target_id(i) = ...
            test_schedule_block.stimulus_id(i-2);
        
        if string(test_schedule_block.stimulus_id(i)) == ...
           string(test_schedule_block.nback_target_id(i))
            test_schedule_block.trial_type_final(i) = "repeat";
            test_schedule_block.correct_response(i) = "j";
        end
    end
    
    % append into full schedule
    test_schedule_all = [test_schedule_all; test_schedule_block]; %#ok<AGROW>
    
    % update response counts
    p.counts.test.j_presses = p.counts.test.j_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + ...
        sum(strcmp(test_schedule_block.correct_response, 'k'));
    p.counts.test.no_response = p.counts.test.no_response + ...
        sum(strcmp(test_schedule_block.correct_response, 'none'));
end

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
%  SECTION 7: SANITY CHECKS
%  ========================================================================

fprintf('\n===== FINAL SANITY CHECKS =====\n');

% --- check 1: target precedes lure in encoding (isolated_both only) ---
violations_enc = [];
all_iso_pairs = p.stim.iso_both;

for i = 1:height(all_iso_pairs)
    tgt = all_iso_pairs.Target(i);
    lur = all_iso_pairs.Lure(i);
    
    idx_tgt = find(encoding_schedule_all.stimulus_id == tgt, 1);
    idx_lur = find(encoding_schedule_all.stimulus_id == lur, 1);
    
    if isempty(idx_tgt) || isempty(idx_lur)
        continue;
    end
    
    if idx_lur < idx_tgt
        violations_enc(end+1, :) = [i, idx_tgt, idx_lur];
    end
end

if isempty(violations_enc)
    fprintf('✓ Encoding schedule OK: All isolated_both targets precede lures.\n');
else
    fprintf('✗ Encoding violations found: %d pairs!\n', size(violations_enc, 1));
    disp(array2table(violations_enc, ...
                     'VariableNames', {'PairIndex', 'TargetRow', 'LureRow'}));
end

% --- check 2: fixation duration range ---
fprintf('\nFixation Duration Range:\n');
fprintf('  Encoding: %.3f to %.3f seconds\n', ...
        min(encoding_schedule_all.fix_duration), ...
        max(encoding_schedule_all.fix_duration));
fprintf('  Test:     %.3f to %.3f seconds\n', ...
        min(test_schedule_all.fix_duration), ...
        max(test_schedule_all.fix_duration));

% --- check 3: # of experimental conditions and response types ---
for b = 1:p.nBlocks
    fprintf('\n--- BLOCK %d ---\n', b);
    
    enc_block = encoding_schedule_all(encoding_schedule_all.block == b, :);
    test_block = test_schedule_all(test_schedule_all.block == b, :);
    
    % ENCODING: Count conditions
    fprintf('Encoding Phase:\n');
    fprintf('  Comparison:    %d trials\n', sum(strcmp(enc_block.condition, 'comparison')));
    fprintf('  Isolated-both: %d trials\n', sum(strcmp(enc_block.condition, 'iso_both')));
    fprintf('  Foil-repeat:   %d trials\n', sum(strcmp(enc_block.condition, 'foil_repeat')));
    fprintf('  Foil-new:      %d trials\n', sum(strcmp(enc_block.condition, 'foil_new')));
    
    % ENCODING: Count responses
    enc_j = sum(strcmp(enc_block.correct_response, 'j'));
    enc_k = sum(strcmp(enc_block.correct_response, 'k'));
    enc_none = sum(strcmp(enc_block.correct_response, 'none'));
    fprintf('  Responses: j=%d, k=%d, none=%d\n', enc_j, enc_k, enc_none);
    
    % TEST: Count conditions
    fprintf('Test Phase:\n');
    fprintf('  Comparison:    %d trials\n', sum(strcmp(test_block.condition, 'comparison')));
    fprintf('  Isolated-both: %d trials\n', sum(strcmp(test_block.condition, 'isolated_both')));
    fprintf('  Novel:         %d trials\n', sum(strcmp(test_block.condition, 'novel')));
    fprintf('  Foil-repeat:   %d trials\n', sum(strcmp(test_block.condition, 'foil_repeat')));
    fprintf('  Unrelated:     %d trials\n', sum(strcmp(test_block.condition, 'unrelated_filler')));
    
    % TEST: Count responses
    test_j = sum(strcmp(test_block.correct_response, 'j'));
    test_k = sum(strcmp(test_block.correct_response, 'k'));
    test_none = sum(strcmp(test_block.correct_response, 'none'));
    fprintf('  Responses: j=%d, k=%d, none=%d\n', test_j, test_k, test_none);
end

fprintf('===== SANITY CHECK COMPLETE =====\n\n');


%% ========================================================================
%  SECTION 8: SAVE SUBJECT DATA
%  ========================================================================

subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.encoding_schedule = encoding_schedule_all;
subject_data.test_schedule = test_schedule_all;

save(output_filename, 'subject_data');

fprintf('setup saved to:\n%s\n\n', output_filename);
fprintf('Summary:\n');
fprintf('  Encoding: j=%d, k=%d, none=%d\n', ...
        p.counts.encoding.j_presses, ...
        p.counts.encoding.k_presses, ...
        p.counts.encoding.no_response);
fprintf('  Test:     j=%d, k=%d, none=%d\n', ...
        p.counts.test.j_presses, ...
        p.counts.test.k_presses, ...
        p.counts.test.no_response);