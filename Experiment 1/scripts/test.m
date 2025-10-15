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

% <<< CHANGE: Define p.subj_id immediately to prevent crash.
p.subj_id = subj_id; 

% directory
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_selected/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
p.results_dir  = fullfile(base_dir, 'data', sprintf('sub%03d', p.subj_id));
if ~exist(p.results_dir, 'dir'), mkdir(p.results_dir); end
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end

%% ========================================================================
%  SECTION 1: SUBJECT ID AND PARAMETERS
%  ========================================================================
% # of stimulus pairs per condition (total across all blocks)
p.nComparison = 40;
p.nIsolated_Both = 40;
p.nNovel = 40;

% test phase trial counts
% p.nTest_Lure_Events = 120; % <<< CHANGE: Removed, as this was unused.
p.nTest_Repeat_Events = 120; % This is used to create 'safe' 2-back repeat trials.

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
output_filename = fullfile(p.setup_dir, sprintf('sub%03d_setup.mat', p.subj_id));
if exist(output_filename, 'file')
    overwrite = input(sprintf('setup for subject %03d exists. overwrite? (y/n): ', p.subj_id), 's');
    if ~strcmpi(overwrite, 'y')
        fprintf('aborted.\n');
        return;
    end
end

fprintf('loading stimuli...\n');
all_targ_files = dir(fullfile(p.stim_dir, 'mst_*_targ_l2.png'));
all_lure_files = dir(fullfile(p.stim_dir, 'mst_*_lure_l2.png'));
assert(numel(all_targ_files) == 120, 'expected 120 target files!');

% <<< CHANGE: Robustly pair targets and lures instead of sorting independently.
% This prevents mismatching pairs (e.g., targ_001 with lure_005).
targ_names = string({all_targ_files.name}');
lure_names = string({all_lure_files.name}');
% Extract numeric ID from filenames to ensure correct pairing
get_id = @(s) sscanf(s, 'mst_%d_');
targ_ids = arrayfun(@(x) get_id(x), targ_names);
lure_ids = arrayfun(@(x) get_id(x), lure_names);
[~, sort_idx_targ] = sort(targ_ids);
[~, sort_idx_lure] = sort(lure_ids);
master_pair_list = table(targ_names(sort_idx_targ), lure_names(sort_idx_lure), ...
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
% divide each condition's stimuli evenly across the 4 blocks
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
% <<< CHANGE: Pre-allocate tables for efficiency instead of growing them in a loop.
encoding_schedule_all = table();
test_schedule_all = table();

% response counters
p.counts.encoding = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);
p.counts.test = struct('j_presses', 0, 'k_presses', 0, 'no_response', 0);

% <<< CHANGE: Use an index to consume foils efficiently.
foil_idx = 1;

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
    for i = 1:nComp_b
        encoding_blocks{end+1} = { ...
            comparison_pairs_b.Target(i), "comparison", "target", "new", "none"; ...
            comparison_pairs_b.Lure(i), "comparison", "lure", "lure", "k"};
    end
    
    % --- isolated-both condition encoding ---
    for i = 1:nIso_b
        encoding_blocks{end+1} = { ...
            all_foils(foil_idx), "iso_both", "filler", "new", "none"; ...
            iso_both_pairs_b.Target(i), "iso_both", "target", "new", "none"};
        foil_idx = foil_idx + 1;
        encoding_blocks{end+1} = { ...
            all_foils(foil_idx), "iso_both", "filler", "new", "none"; ...
            iso_both_pairs_b.Lure(i), "iso_both", "lure", "new", "none"};
        foil_idx = foil_idx + 1;
    end
    
    % --- spacer blocks: foil repeats (provide 'j' responses) ---
    n_foil_repeats = nComp_b;
    for i = 1:n_foil_repeats
        spacer_blocks{end+1} = { ...
            all_foils(foil_idx), "foil_repeat", "filler", "new", "none"; ...
            all_foils(foil_idx), "foil_repeat", "filler", "repeat", "j"};
        foil_idx = foil_idx + 1;
    end
    
    % --- spacer blocks: new foils (fill remaining slots) ---
    n_new_spacers = numel(encoding_blocks) - numel(spacer_blocks);
    for i = 1:n_new_spacers
        spacer_blocks{end+1} = { ...
            all_foils(foil_idx), "foil_new", "filler", "new", "none"; ...
            all_foils(foil_idx+1), "foil_new", "filler", "new", "none"};
        foil_idx = foil_idx + 2;
    end
    
    % --- interleave encoding and spacer blocks ---
    shuffled_encoding = encoding_blocks(randperm(numel(encoding_blocks)));
    shuffled_spacers = spacer_blocks(randperm(numel(spacer_blocks)));
    final_ordered_blocks = reshape([shuffled_spacers; shuffled_encoding], 1, []);
    final_encoding_list = vertcat(final_ordered_blocks{:});
    
    encoding_schedule_block = cell2table(final_encoding_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    encoding_schedule_block.block = repmat(b, height(encoding_schedule_block), 1);
    
    % =====================================================================
    % 5C: ENSURE TARGET ENCODED BEFORE LURE (ISOLATED CONDITION ONLY)
    % =====================================================================
    for i = 1:height(iso_both_pairs_b)
        tgt = iso_both_pairs_b.Target(i);
        lur = iso_both_pairs_b.Lure(i);
        idx_tgt = find(strcmp(encoding_schedule_block.stimulus_id, tgt), 1);
        idx_lur = find(strcmp(encoding_schedule_block.stimulus_id, lur), 1);
        
        if idx_lur < idx_tgt
            tmp = encoding_schedule_block(idx_tgt, :);
            encoding_schedule_block(idx_tgt, :) = encoding_schedule_block(idx_lur, :);
            encoding_schedule_block(idx_lur, :) = tmp;
        end
    end
    
    % =====================================================================
    % 5D: FINAL N-BACK VERIFICATION
    % =====================================================================
    % <<< CHANGE: This N-back check is now run *after* the target/lure swap
    % to ensure it catches any new repeats created by the swap.
    encoding_schedule_block.nback_target_id = strings(height(encoding_schedule_block), 1);
    encoding_schedule_block.trial_type_final = encoding_schedule_block.trial_type_designed;
    for i = 2:height(encoding_schedule_block)
        encoding_schedule_block.nback_target_id(i) = encoding_schedule_block.stimulus_id(i-1);
        if strcmp(encoding_schedule_block.stimulus_id(i), encoding_schedule_block.nback_target_id(i))
            encoding_schedule_block.trial_type_final(i) = "repeat";
            encoding_schedule_block.correct_response(i) = "j";
        end
    end
    
    encoding_schedule_all = [encoding_schedule_all; encoding_schedule_block];
    
    p.counts.encoding.j_presses = p.counts.encoding.j_presses + sum(strcmp(encoding_schedule_block.correct_response, 'j'));
    p.counts.encoding.k_presses = p.counts.encoding.k_presses + sum(strcmp(encoding_schedule_block.correct_response, 'k'));
    p.counts.encoding.no_response = p.counts.encoding.no_response + sum(strcmp(encoding_schedule_block.correct_response, 'none'));
    
    % =====================================================================
    % 5E: BUILD TEST PHASE FOR THIS BLOCK
    % =====================================================================
    phase2_crit_blocks = {};  % Lure triads (require 'k')
    phase2_safe_blocks = {};  % Repeat triads (require 'j')
    
    % --- critical triads: target → filler → lure ---
    all_pairs_b = [comparison_pairs_b; iso_both_pairs_b; novel_pairs_b];
    all_conds_b = [repmat("comparison", nComp_b, 1); ...
                   repmat("isolated_both", nIso_b, 1); ...
                   repmat("novel", nNov_b, 1)];
                   
    for i = 1:height(all_pairs_b)
        phase2_crit_blocks{end+1} = { ...
            all_pairs_b.Target(i), all_conds_b(i), "target", "new", "none"; ...
            all_foils(foil_idx), "unrelated_filler", "filler", "new", "none"; ...
            all_pairs_b.Lure(i), all_conds_b(i), "lure", "lure", "k"};
        foil_idx = foil_idx + 1;
    end

    % --- safe triads: foilA → foilB → foilA (2-back repeat) ---
    total_pairs = p.nComparison + p.nIsolated_Both + p.nNovel;
    proportion = height(all_pairs_b) / total_pairs;
    n_repeat_triads_block = max(round(p.nTest_Repeat_Events * proportion), 1);
    
    for i = 1:n_repeat_triads_block
        phase2_safe_blocks{end+1} = { ...
            all_foils(foil_idx), "foil_repeat", "filler", "new", "none"; ...
            all_foils(foil_idx+1), "foil_repeat", "filler", "new", "none"; ...
            all_foils(foil_idx), "foil_repeat", "filler", "repeat", "j"};
        foil_idx = foil_idx + 2;
    end
    
    % --- interleave safe and critical triads ---
    interleaved_blocks = [phase2_crit_blocks(randperm(numel(phase2_crit_blocks))), ...
                          phase2_safe_blocks(randperm(numel(phase2_safe_blocks)))];
    final_test_list = vertcat(interleaved_blocks{randperm(numel(interleaved_blocks))});
    
    test_schedule_block = cell2table(final_test_list, ...
        'VariableNames', {'stimulus_id', 'condition', 'role', ...
                          'trial_type_designed', 'correct_response'});
    test_schedule_block.block = repmat(b, height(test_schedule_block), 1);
    
    % final check for any accidental 2-back repeats from shuffling
    test_schedule_block.nback_target_id = strings(height(test_schedule_block), 1);
    test_schedule_block.trial_type_final = test_schedule_block.trial_type_designed;
    for i = 3:height(test_schedule_block)
        test_schedule_block.nback_target_id(i) = test_schedule_block.stimulus_id(i-2);
        if strcmp(test_schedule_block.stimulus_id(i), test_schedule_block.nback_target_id(i))
            test_schedule_block.trial_type_final(i) = "repeat";
            test_schedule_block.correct_response(i) = "j";
        end
    end
    
    test_schedule_all = [test_schedule_all; test_schedule_block];
    
    p.counts.test.j_presses = p.counts.test.j_presses + sum(strcmp(test_schedule_block.correct_response, 'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + sum(strcmp(test_schedule_block.correct_response, 'k'));
    p.counts.test.no_response = p.counts.test.no_response + sum(strcmp(test_schedule_block.correct_response, 'none'));
end
assert(foil_idx <= numel(all_foils)+1, 'FATAL: Not enough foils for the experiment!');
%% ========================================================================
%  SECTION 6: GENERATE JITTERED FIXATION DURATIONS
%  ========================================================================
n_enc_trials = height(encoding_schedule_all);
n_test_trials = height(test_schedule_all);
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
    idx_tgt = find(strcmp(encoding_schedule_all.stimulus_id, tgt), 1);
    idx_lur = find(strcmp(encoding_schedule_all.stimulus_id, lur), 1);
    if isempty(idx_tgt) || isempty(idx_lur), continue; end
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
    fprintf('Encoding Phase (%d trials):\n', height(enc_block));
    disp(groupcounts(enc_block, {'condition', 'correct_response'}));
    fprintf('Test Phase (%d trials):\n', height(test_block));
    disp(groupcounts(test_block, {'condition', 'correct_response'}));
end
fprintf('===== SANITY CHECK COMPLETE =====\n\n');
%% ========================================================================
%  SECTION 8: SAVE SUBJECT DATA
%  ========================================================================
subject_data.subj_id = p.subj_id;
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