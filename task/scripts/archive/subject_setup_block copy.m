%==========================================================================
% Generate subject-specific task variables (3-block design)
% Ensures target encoded before lure and block-local testing
%==========================================================================
clear; clc; rng('shuffle');

%% PARAMETERS AND DIRECTORIES
subj_id_str = input('please enter subject ID (e.g., 101): ', 's');
if isempty(subj_id_str), error('subject ID cannot be empty.'); end
subj_id = str2double(subj_id_str);

p.nComparison = 40;
p.nIsolated_Both = 40;
p.nNovel = 40;
p.nTest_Lure_Events = 120;
p.nTest_Repeat_Events = 120;
p.nBlocks = 3;

base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');
p.keys.same = 'j'; p.keys.diff = 'k'; p.keys.quit = 'escape';
p.timing.image_dur = 1.5; p.timing.fix_dur = 0.5;

if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end
output_filename = fullfile(p.setup_dir, sprintf('sub%03d.mat', subj_id));
if exist(output_filename, 'file')
    overwrite = input(sprintf('WARNING: setup for subject %03d exists. overwrite? (y/n): ', subj_id), 's');
    if ~strcmpi(overwrite, 'y'), fprintf('aborted.\n'); return; end
end

%% LOAD STIMULI
fprintf('Loading stimuli...\n');
all_targ_files = dir(fullfile(p.stim_dir, 'mst_*_targ.png'));
all_lure_files = dir(fullfile(p.stim_dir, 'mst_*_lure.png'));
assert(numel(all_targ_files) == 120, 'expected 120 target files!');
targ_names = string({all_targ_files.name}'); lure_names = string({all_lure_files.name}');
master_pair_list = table(sort(targ_names), sort(lure_names), 'VariableNames', {'Img_A','Img_B'});
all_foil_files = dir(fullfile(p.stim_dir, 'mst_*_foil.png'));
all_foils = string({all_foil_files.name}');
fprintf('Found %d experimental pairs and %d foils.\n', height(master_pair_list), numel(all_foils));

%% COUNTERBALANCE PAIRS
shuffled_pairs = master_pair_list(randperm(height(master_pair_list)), :);
final_pairs = table('Size', size(shuffled_pairs), 'VariableTypes', {'string','string'}, ...
    'VariableNames', {'Target','Lure'});

for i = 1:height(shuffled_pairs)
    if rand() > 0.5
        final_pairs(i,:) = shuffled_pairs(i,:);
    else
        final_pairs(i,:) = [shuffled_pairs(i,2), shuffled_pairs(i,1)];
    end
end

% Split into conditions
idx_start = 1; idx_end = p.nComparison; comparison_pairs = final_pairs(idx_start:idx_end,:);
idx_start = idx_end+1; idx_end = idx_start+p.nIsolated_Both-1; iso_both_pairs = final_pairs(idx_start:idx_end,:);
idx_start = idx_end+1; idx_end = idx_start+p.nNovel-1; novel_pairs = final_pairs(idx_start:idx_end,:);
all_foils = all_foils(randperm(numel(all_foils)));

p.stim.comparison = comparison_pairs; p.stim.iso_both = iso_both_pairs; p.stim.novel = novel_pairs;

%% HELPER: partition indices into n blocks
partition_idx = @(N,nblocks) arrayfun(@(k) ((floor((k-1)*N/nblocks)+1):floor(k*N/nblocks)), 1:nblocks, 'UniformOutput', false);
comp_idx_blocks = partition_idx(p.nComparison,p.nBlocks);
iso_idx_blocks  = partition_idx(p.nIsolated_Both,p.nBlocks);
nov_idx_blocks  = partition_idx(p.nNovel,p.nBlocks);

%% CONTAINERS
encoding_schedule_all = table(); test_schedule_all = table();
p.counts.encoding = struct('j_presses',0,'k_presses',0,'no_response',0);
p.counts.test = struct('j_presses',0,'k_presses',0,'no_response',0);

%% BUILD EACH BLOCK
fprintf('Building %d blocks...\n', p.nBlocks);
for b = 1:p.nBlocks
    fprintf('--- Block %d ---\n', b);

    % Select block-specific pairs
    comparison_pairs_b = p.stim.comparison(comp_idx_blocks{b},:);
    iso_both_pairs_b   = p.stim.iso_both(iso_idx_blocks{b},:);
    novel_pairs_b      = p.stim.novel(nov_idx_blocks{b},:);
    nComp_b = height(comparison_pairs_b); nIso_b = height(iso_both_pairs_b); nNov_b = height(novel_pairs_b);

    %% -------------------- ENCODING --------------------
    encoding_blocks = {}; spacer_blocks = {};

    % Comparison encoding
    for i=1:nComp_b
        encoding_blocks{end+1} = {comparison_pairs_b.Target(i),"comparison","target","new","none"; ...
                                  comparison_pairs_b.Lure(i),"comparison","lure","lure","k"};
    end

    % Isolated-both encoding (two separate blocks per pair)
    iso_event_foils = all_foils(1:2*nIso_b); all_foils(1:2*nIso_b) = [];
    for i=1:nIso_b
        encoding_blocks{end+1} = {iso_event_foils(1),"iso_both","filler","new","none"; iso_both_pairs_b.Target(i),"iso_both","target","new","none"}; iso_event_foils(1)=[];
        encoding_blocks{end+1} = {iso_event_foils(1),"iso_both","filler","new","none"; iso_both_pairs_b.Lure(i),"iso_both","lure","new","none"}; iso_event_foils(1)=[];
    end

    % Spacer blocks: repeat foils
    n_foil_repeats = nComp_b; repeat_foils = all_foils(1:n_foil_repeats); all_foils(1:n_foil_repeats) = [];
    for i=1:n_foil_repeats
        spacer_blocks{end+1} = {repeat_foils(i),"foil_repeat","filler","new","none"; repeat_foils(i),"foil_repeat","filler","repeat","j"};
    end

    % Spacer blocks: new foils to match encoding blocks
    n_new_spacers = numel(encoding_blocks)-numel(spacer_blocks);
    new_spacer_foils_1 = all_foils(1:n_new_spacers); all_foils(1:n_new_spacers) = [];
    new_spacer_foils_2 = all_foils(1:n_new_spacers); all_foils(1:n_new_spacers) = [];
    for i=1:n_new_spacers
        spacer_blocks{end+1} = {new_spacer_foils_1(i),"foil_new","filler","new","none"; new_spacer_foils_2(i),"foil_new","filler","new","none"};
    end

    % Interleave encoding and spacer
    shuffled_encoding = encoding_blocks(randperm(numel(encoding_blocks)));
    shuffled_spacers  = spacer_blocks(randperm(numel(spacer_blocks)));
    final_ordered_blocks = reshape([shuffled_spacers; shuffled_encoding],1,[]);
    final_encoding_list = vertcat(final_ordered_blocks{:});

    encoding_schedule_block = cell2table(final_encoding_list, 'VariableNames', {'stimulus_id','condition','role','trial_type_designed','correct_response'});
    
    % *** MODIFICATION START ***
    encoding_schedule_block.block = repmat(b, height(encoding_schedule_block), 1);
    % *** MODIFICATION END ***

    encoding_schedule_block.nback_target_id = strings(height(encoding_schedule_block),1);
    encoding_schedule_block.trial_type_final = encoding_schedule_block.trial_type_designed;
    for i=2:height(encoding_schedule_block)
        encoding_schedule_block.nback_target_id(i) = encoding_schedule_block.stimulus_id(i-1);
        if string(encoding_schedule_block.stimulus_id(i))==string(encoding_schedule_block.nback_target_id(i))
            encoding_schedule_block.trial_type_final(i)="repeat"; encoding_schedule_block.correct_response(i)="j";
        end
    end

    % -------------------- ENSURE TARGET BEFORE LURE --------------------
    all_iso_pairs_block = iso_both_pairs_b;
    for i=1:height(all_iso_pairs_block)
        tgt = all_iso_pairs_block.Target(i); lur = all_iso_pairs_block.Lure(i);
        idx_tgt = find(encoding_schedule_block.stimulus_id==tgt,1);
        idx_lur = find(encoding_schedule_block.stimulus_id==lur,1);
        if ~isempty(idx_tgt) && ~isempty(idx_lur) && idx_lur<idx_tgt
            tmp = encoding_schedule_block(idx_tgt,:); encoding_schedule_block(idx_tgt,:) = encoding_schedule_block(idx_lur,:); encoding_schedule_block(idx_lur,:) = tmp;
        end
    end

    encoding_schedule_all = [encoding_schedule_all; encoding_schedule_block]; %#ok<AGROW>
    p.counts.encoding.j_presses = p.counts.encoding.j_presses + sum(strcmp(encoding_schedule_block.correct_response,'j'));
    p.counts.encoding.k_presses = p.counts.encoding.k_presses + sum(strcmp(encoding_schedule_block.correct_response,'k'));
    p.counts.encoding.no_response = p.counts.encoding.no_response + sum(strcmp(encoding_schedule_block.correct_response,'none'));

    %% -------------------- TEST --------------------
    phase2_crit_blocks = {}; phase2_safe_blocks = {};
    n_crit_expected = nComp_b+nIso_b+nNov_b;
    assert(numel(all_foils)>=n_crit_expected,'Not enough foils for critical triads!');
    intervening_foils = all_foils(1:n_crit_expected); all_foils(1:n_crit_expected)=[];

    % Comparison crit triads
    for i=1:nComp_b
        phase2_crit_blocks{end+1} = {comparison_pairs_b.Target(i),"comparison","target","new","none"; ...
                                      intervening_foils(1),"unrelated_filler","filler","new","none"; ...
                                      comparison_pairs_b.Lure(i),"comparison","lure","lure","k"}; intervening_foils(1)=[];
    end

    % Isolated-both crit triads
    for i=1:nIso_b
        phase2_crit_blocks{end+1} = {iso_both_pairs_b.Target(i),"isolated_both","target","new","none"; ...
                                      intervening_foils(1),"unrelated_filler","filler","new","none"; ...
                                      iso_both_pairs_b.Lure(i),"isolated_both","lure","lure","k"}; intervening_foils(1)=[];
    end

    % Novel crit triads
    for i=1:nNov_b
        phase2_crit_blocks{end+1} = {novel_pairs_b.Target(i),"novel","target","new","none"; ...
                                      intervening_foils(1),"unrelated_filler","filler","new","none"; ...
                                      novel_pairs_b.Lure(i),"novel","lure","lure","k"}; intervening_foils(1)=[];
    end

    % Safe repeats per block (proportional)
    total_pairs_block = nComp_b+nIso_b+nNov_b;
    proportion = total_pairs_block/(p.nComparison+p.nIsolated_Both+p.nNovel);
    n_repeat_triads_block = max(round(p.nTest_Repeat_Events*proportion),1);
    n_needed_foils = 2*n_repeat_triads_block; assert(numel(all_foils)>=n_needed_foils,'Not enough foils for repeat triads!');
    repeat_foil_A = all_foils(1:n_repeat_triads_block); repeat_foil_B = all_foils(n_repeat_triads_block+1:2*n_repeat_triads_block);
    all_foils(1:2*n_repeat_triads_block) = [];
    for i=1:n_repeat_triads_block
        phase2_safe_blocks{end+1} = {repeat_foil_A(i),"foil_repeat","filler","new","none"; ...
                                      repeat_foil_B(i),"foil_repeat","filler","new","none"; ...
                                      repeat_foil_A(i),"foil_repeat","filler","repeat","j"};
    end

    % Interleave safe and crit triads
    n_crit = numel(phase2_crit_blocks); n_safe = numel(phase2_safe_blocks);
    shuffled_crit = phase2_crit_blocks(randperm(n_crit)); shuffled_safe = phase2_safe_blocks(randperm(n_safe));
    L = max(n_crit,n_safe); interleaved_blocks = {};
    for i=1:L
        if i<=n_safe, interleaved_blocks{end+1} = shuffled_safe{i}; end
        if i<=n_crit, interleaved_blocks{end+1} = shuffled_crit{i}; end
    end
    final_test_list = vertcat(interleaved_blocks{:});
    test_schedule_block = cell2table(final_test_list, 'VariableNames', {'stimulus_id','condition','role','trial_type_designed','correct_response'});
    
    % *** MODIFICATION START ***
    test_schedule_block.block = repmat(b, height(test_schedule_block), 1);
    % *** MODIFICATION END ***

    test_schedule_block.nback_target_id = strings(height(test_schedule_block),1);
    test_schedule_block.trial_type_final = test_schedule_block.trial_type_designed;
    for i=3:height(test_schedule_block)
        test_schedule_block.nback_target_id(i) = test_schedule_block.stimulus_id(i-2);
        if string(test_schedule_block.stimulus_id(i))==string(test_schedule_block.nback_target_id(i))
            test_schedule_block.trial_type_final(i)="repeat"; test_schedule_block.correct_response(i)="j";
        end
    end

    test_schedule_all = [test_schedule_all; test_schedule_block]; %#ok<AGROW>
    p.counts.test.j_presses = p.counts.test.j_presses + sum(strcmp(test_schedule_block.correct_response,'j'));
    p.counts.test.k_presses = p.counts.test.k_presses + sum(strcmp(test_schedule_block.correct_response,'k'));
    p.counts.test.no_response = p.counts.test.no_response + sum(strcmp(test_schedule_block.correct_response,'none'));
end

%% ================= UPDATED FINAL SANITY CHECKS =================
fprintf('\n===== FINAL SANITY CHECKS =====\n');

% 1. Encoding schedule: target precedes lure for isolated_both
violations_enc = [];
all_iso_pairs = p.stim.iso_both;
for i = 1:height(all_iso_pairs)
    tgt = all_iso_pairs.Target(i);
    lur = all_iso_pairs.Lure(i);
    idx_tgt = find(encoding_schedule_all.stimulus_id == tgt,1);
    idx_lur = find(encoding_schedule_all.stimulus_id == lur,1);
    if isempty(idx_tgt) || isempty(idx_lur), continue; end
    if idx_lur < idx_tgt
        violations_enc(end+1,:) = [i, idx_tgt, idx_lur];
    end
end
if isempty(violations_enc)
    fprintf('Encoding schedule OK: All isolated_both targets precede lures.\n');
else
    fprintf('Encoding violations found: %d pairs!\n', size(violations_enc,1));
    disp(array2table(violations_enc, 'VariableNames', {'PairIndex','TargetRow','LureRow'}));
end

% 3. Response balance summary
fprintf('\nEncoding responses: j=%d, k=%d, none=%d\n', p.counts.encoding.j_presses, p.counts.encoding.k_presses, p.counts.encoding.no_response);
fprintf('Test responses:     j=%d, k=%d, none=%d\n', p.counts.test.j_presses, p.counts.test.k_presses, p.counts.test.no_response);

% 4. Total trial counts
fprintf('Total encoding trials: %d (expected >= %d)\n', height(encoding_schedule_all), p.nComparison + p.nIsolated_Both*2 + p.nNovel*2);
fprintf('Total test trials:     %d (expected >= %d)\n', height(test_schedule_all), p.nTest_Lure_Events + p.nTest_Repeat_Events);
fprintf('===== SANITY CHECK COMPLETE =====\n\n');


%% SAVE SUBJECT SETUP
subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.encoding_schedule = encoding_schedule_all;
subject_data.test_schedule = test_schedule_all;
save(output_filename,'subject_data');
fprintf('Setup saved to %s\n', output_filename);
fprintf('Encoding: j=%d, k=%d, none=%d\n', p.counts.encoding.j_presses,p.counts.encoding.k_presses,p.counts.encoding.no_response);
fprintf('Test: j=%d, k=%d, none=%d\n', p.counts.test.j_presses,p.counts.test.k_presses,p.counts.test.no_response);