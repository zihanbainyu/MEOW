%==========================================================================
%                  Generate subject-specific task variables
%==========================================================================
% author: Zihan Bai

% description:
% this script builds a unique, counterbalanced, and randomized trial schedule
% that is BOTH RESPONSE-BALANCED and GUARANTEED TO BE SPACED using a
% direct interleaving method for both phases.
%==========================================================================

clear; 
clc;
rng('shuffle');

%% DEFINE PARAMETERS AND DIRECTORIES
%--------------------------------------------------------------------------
subj_id_str = input('please enter subject ID (e.g., 101): ', 's');
if isempty(subj_id_str), error('subject ID cannot be empty.'); end
subj_id = str2double(subj_id_str);

% define number of experimenal pairs
p.nComparison = 40;
p.nIsolated_Both = 40;
p.nNovel = 40;

% define test phase events
p.nTest_Lure_Events = 120;
p.nTest_Repeat_Events = 120; 

% directories 
base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
p.setup_dir = fullfile(base_dir, 'subj_setup/'); 

% define keypress and timing variables
p.keys.same = 'j'; p.keys.diff = 'k'; p.keys.quit = 'escape';
p.timing.image_dur = 1.5; 
p.timing.fix_dur = 0.5;

%% PREPARE DIRECTORIES & PREVENT OVERWRITE
%--------------------------------------------------------------------------
if ~exist(p.setup_dir, 'dir'), mkdir(p.setup_dir); end
output_filename = fullfile(p.setup_dir, sprintf('sub%03d.mat', subj_id));
if exist(output_filename, 'file')
    overwrite = input(sprintf('WARNING: setup for subject %03d exists. overwrite? (y/n): ', subj_id), 's');
    if ~strcmpi(overwrite, 'y'), fprintf('aborted.\n'); return; end
end

%% LOAD STIMULI FROM FIXED SET
%--------------------------------------------------------------------------
fprintf('loading stimuli from fixed set...\n');
all_targ_files = dir(fullfile(p.stim_dir, 'mst_*_targ.png'));
all_lure_files = dir(fullfile(p.stim_dir, 'mst_*_lure.png'));
assert(numel(all_targ_files) == 120, 'expected 120 target files!');
targ_names=string({all_targ_files.name}'); lure_names=string({all_lure_files.name}');
master_pair_list = table(sort(targ_names), sort(lure_names), 'VariableNames', {'Img_A', 'Img_B'});
all_foil_files = dir(fullfile(p.stim_dir, 'mst_*_foil.png'));
all_foils = string({all_foil_files.name}');
fprintf('found %d experimental pairs and %d foils.\n', height(master_pair_list), numel(all_foils));

%% COUNTERBALANCE AND ASSIGN EXPERIMENTAL PAIRS
%--------------------------------------------------------------------------
fprintf('assigning stimuli to conditions for subject %03d...\n');
shuffled_pairs = master_pair_list(randperm(height(master_pair_list)), :);
final_pairs = table('Size', size(shuffled_pairs), 'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Target', 'Lure'});

for i = 1:height(shuffled_pairs) 
    if rand()>0.5 
        final_pairs(i,:)=shuffled_pairs(i,:); 
    else 
        final_pairs(i,:)=[shuffled_pairs(i,2),shuffled_pairs(i,1)]; 
    end 
end

all_foils = all_foils(randperm(numel(all_foils)));

idx_start=1; idx_end=p.nComparison;
comparison_pairs=final_pairs(idx_start:idx_end,:);

idx_start=idx_end+1; idx_end=idx_start+p.nIsolated_Both-1;   
iso_both_pairs=final_pairs(idx_start:idx_end,:);

idx_start=idx_end+1; idx_end=idx_start+p.nNovel-1;           
novel_pairs=final_pairs(idx_start:idx_end,:);

p.stim.comparison=comparison_pairs; 
p.stim.iso_both=iso_both_pairs; 
p.stim.novel=novel_pairs;


%% BUILD PHASE 1 (ENCODING)
%--------------------------------------------------------------------------
fprintf('building phase 1 with interleaving spacers...\n');

encoding_blocks = {}; % all mnemonically significant blocks
spacer_blocks = {};   % all mnemonically neutral blocks


% comparison encoding events(40 blocks)
for i=1:p.nComparison, encoding_blocks{end+1} = { ...
    comparison_pairs.Target(i), "comparison", "target", "new", "none"; ...
    comparison_pairs.Lure(i),   "comparison", "lure",   "lure", "k" };end

% isolated encoding events (40 blocks)
iso_event_foils = all_foils(1:2*p.nIsolated_Both); 
all_foils(1:2*p.nIsolated_Both) = [];
for i=1:p.nIsolated_Both
    encoding_blocks{end+1} = {iso_event_foils(1),"iso_both",  "filler","new","none"; iso_both_pairs.Target(i),"iso_both",  "target","new","none"}; iso_event_foils(1)=[];
    encoding_blocks{end+1} = {iso_event_foils(1),"iso_both",  "filler","new","none"; iso_both_pairs.Lure(i),"iso_both", "lure","new","none"};   iso_event_foils(1)=[];
end

% 40 spacer blocks will be repeats (providing 40 'j' presses)
n_foil_repeats = 40;
repeat_foils = all_foils(1:n_foil_repeats); all_foils(1:n_foil_repeats) = [];
for i=1:n_foil_repeats, spacer_blocks{end+1} = { ...
    repeat_foils(i), "foil_repeat", "filler", "new", "none"; ...
    repeat_foils(i), "foil_repeat", "filler", "repeat", "j"}; end

% the rest will be new foils to make up the total of 94
n_new_spacers = numel(encoding_blocks) - numel(spacer_blocks); % Should be 54
new_spacer_foils_1 = all_foils(1:n_new_spacers); all_foils(1:n_new_spacers) = [];
new_spacer_foils_2 = all_foils(1:n_new_spacers); all_foils(1:n_new_spacers) = [];
for i=1:n_new_spacers, spacer_blocks{end+1} = { ...
    new_spacer_foils_1(i), "foil_new", "filler", "new", "none"; ...
    new_spacer_foils_2(i), "foil_new", "filler", "new", "none"}; end

expected_blocks = numel(encoding_blocks);
assert(numel(spacer_blocks)==expected_blocks, ...
    sprintf('Block count mismatch! encoding=%d, spacer=%d', numel(encoding_blocks), numel(spacer_blocks)));

% assemble the sequence with guaranteed interleaving
shuffled_encoding = encoding_blocks(randperm(numel(encoding_blocks)));
shuffled_spacers = spacer_blocks(randperm(numel(spacer_blocks)));
final_ordered_blocks = reshape([shuffled_spacers; shuffled_encoding], 1, []);
final_encoding_list = vertcat(final_ordered_blocks{:});

% finalize the table and verify
encoding_schedule = cell2table(final_encoding_list, 'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type_designed', 'correct_response'});
encoding_schedule.nback_target_id = strings(height(encoding_schedule), 1);
encoding_schedule.trial_type_final = encoding_schedule.trial_type_designed;
for i=2:height(encoding_schedule)
    encoding_schedule.nback_target_id(i) = encoding_schedule.stimulus_id(i-1);
    if string(encoding_schedule.stimulus_id(i))==string(encoding_schedule.nback_target_id(i)), encoding_schedule.trial_type_final(i)="repeat"; encoding_schedule.correct_response(i)="j"; end
end

p.counts.encoding.k_presses = sum(strcmp(encoding_schedule.correct_response, "k"));
p.counts.encoding.j_presses = sum(strcmp(encoding_schedule.correct_response, "j"));
p.counts.encoding.no_response = sum(strcmp(encoding_schedule.correct_response, "none"));

fprintf('Phase 1 generated with %d trials.\n', height(encoding_schedule));
fprintf('  Response balance: %d Same(j), %d Diff(k), %d None\n', ...
    p.counts.encoding.j_presses, p.counts.encoding.k_presses, p.counts.encoding.no_response);

%% ENSURE TARGETS ARE ALWAYS ENCODED BEFORE LURES
%--------------------------------------------------------------------------
fprintf('\nchecking encoding order for potential target-lure reversals...\n');

% Combine all pairs and condition labels
all_pair_structs = [...
    table(p.stim.comparison.Target, p.stim.comparison.Lure, repmat("comparison", p.nComparison, 1), ...
          'VariableNames', {'Target','Lure','Condition'}); ...
    table(p.stim.iso_both.Target, p.stim.iso_both.Lure, repmat("iso_both", p.nIsolated_Both, 1), ...
          'VariableNames', {'Target','Lure','Condition'}); ...
    table(p.stim.novel.Target, p.stim.novel.Lure, repmat("novel", p.nNovel, 1), ...
          'VariableNames', {'Target','Lure','Condition'})];

n_swapped = 0;
for i = 1:height(all_pair_structs)
    tgt = all_pair_structs.Target(i);
    lur = all_pair_structs.Lure(i);

    idx_tgt = find(encoding_schedule.stimulus_id == tgt, 1);
    idx_lur = find(encoding_schedule.stimulus_id == lur, 1);

    % only consider pairs where both exist in encoding schedule
    if isempty(idx_tgt) || isempty(idx_lur)
        continue;
    end

    % if lure appears before target, swap their positions
    if idx_lur < idx_tgt
        tmp = encoding_schedule(idx_tgt, :);
        encoding_schedule(idx_tgt, :) = encoding_schedule(idx_lur, :);
        encoding_schedule(idx_lur, :) = tmp;
        n_swapped = n_swapped + 1;
    end
end

fprintf('Reordered %d target-lure pairs to ensure targets are encoded first.\n', n_swapped);
fprintf('Encoding order check complete.\n\n');

%% VERIFY TARGET–LURE ORDER IN FINAL ENCODING SCHEDULE
fprintf('verifying final encoding order integrity...\n');

violations = []; % store info about any remaining problems
for i = 1:height(all_pair_structs)
    tgt = all_pair_structs.Target(i);
    lur = all_pair_structs.Lure(i);
    idx_tgt = find(encoding_schedule.stimulus_id == tgt, 1);
    idx_lur = find(encoding_schedule.stimulus_id == lur, 1);
    
    if isempty(idx_tgt) || isempty(idx_lur)
        continue; % skip if pair wasn’t used in encoding phase
    end

    if idx_lur < idx_tgt
        violations(end+1,:) = [i, idx_tgt, idx_lur]; %#ok<AGROW>
    end
end

if isempty(violations)
    fprintf('All %d pairs verified: targets appear before lures.\n', height(all_pair_structs));
else
    fprintf('⚠️  Found %d remaining order violations after swap!\n', size(violations,1));
    disp(array2table(violations, 'VariableNames', {'PairIndex','TargetRow','LureRow'}));
end


%% BUILD PHASE 2 TEST SCHEDULE
%--------------------------------------------------------------------------
fprintf('building balanced and spaced phase 2 (2-back test) sequence...\n');

phase2_crit_blocks = {}; % critical lure triads (require 'k' presses)
phase2_safe_blocks = {}; % safe repeat triads (require 'j' presses)

% combine all isolated pairs
total_iso_pairs = p.nIsolated_Both;
iso_all_pairs = iso_both_pairs;
iso_all_conds = repmat("isolated_both", p.nIsolated_Both, 1);

% -------------------------------------------------------------
% Create CRITICAL (LURE) triads → respond 'k'
% -------------------------------------------------------------
n_crit_expected = p.nComparison + total_iso_pairs + p.nNovel;
assert(numel(all_foils) >= n_crit_expected, ...
    'Not enough foils to create intervening items for critical triads!');

intervening_foils = all_foils(1:n_crit_expected); 
all_foils(1:n_crit_expected) = [];

% 1. Comparison condition
for i = 1:p.nComparison
    phase2_crit_blocks{end+1} = { ...
        comparison_pairs.Target(i), "comparison", "target", "new", "none"; ...
        intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
        comparison_pairs.Lure(i), "comparison", "lure", "lure", "k"};
    intervening_foils(1) = [];
end

% 2. Isolated-both condition
for i = 1:total_iso_pairs
    phase2_crit_blocks{end+1} = { ...
        iso_all_pairs.Target(i), iso_all_conds(i), "target", "new", "none"; ...
        intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
        iso_all_pairs.Lure(i), iso_all_conds(i), "lure", "lure", "k"};
    intervening_foils(1) = [];
end

% 3. Novel condition
for i = 1:p.nNovel
    phase2_crit_blocks{end+1} = { ...
        novel_pairs.Target(i), "novel", "target", "new", "none"; ...
        intervening_foils(1), "unrelated_filler", "filler", "new", "none"; ...
        novel_pairs.Lure(i), "novel", "lure", "lure", "k"};
    intervening_foils(1) = [];
end

% -------------------------------------------------------------
% Create SAFE (REPEAT) triads → respond 'j'
% -------------------------------------------------------------
n_repeat_triads = p.nTest_Repeat_Events;
n_needed_foils = 2 * n_repeat_triads;
assert(numel(all_foils) >= n_needed_foils, ...
    'Not enough foils for repeat triads!');

repeat_foil_A = all_foils(1:n_repeat_triads); 
repeat_foil_B = all_foils(n_repeat_triads+1:2*n_repeat_triads);
all_foils(1:2*n_repeat_triads) = [];

for i = 1:n_repeat_triads
    phase2_safe_blocks{end+1} = { ...
        repeat_foil_A(i), "foil_repeat", "filler", "new", "none"; ...
        repeat_foil_B(i), "foil_repeat", "filler", "new", "none"; ...
        repeat_foil_A(i), "foil_repeat", "filler", "repeat", "j"};
end

% -------------------------------------------------------------
% Sanity check
% -------------------------------------------------------------
n_crit = numel(phase2_crit_blocks);
n_safe = numel(phase2_safe_blocks);
assert(n_crit == p.nTest_Lure_Events && n_safe == p.nTest_Repeat_Events, ...
    sprintf('Test block count mismatch! crit=%d safe=%d', n_crit, n_safe));

% -------------------------------------------------------------
% Interleave the two block types
% -------------------------------------------------------------
shuffled_crit = phase2_crit_blocks(randperm(n_crit));
shuffled_safe = phase2_safe_blocks(randperm(n_safe));

final_ordered_blocks = reshape([shuffled_safe; shuffled_crit], 1, []);
final_test_list = vertcat(final_ordered_blocks{:});

% -------------------------------------------------------------
% Finalize and verify
% -------------------------------------------------------------
test_schedule = cell2table(final_test_list, ...
    'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type_designed', 'correct_response'});

test_schedule.nback_target_id = strings(height(test_schedule), 1);
test_schedule.trial_type_final = test_schedule.trial_type_designed;

for i = 3:height(test_schedule)
    test_schedule.nback_target_id(i) = test_schedule.stimulus_id(i-2);
    if string(test_schedule.stimulus_id(i)) == string(test_schedule.nback_target_id(i))
        test_schedule.trial_type_final(i) = "repeat";
        test_schedule.correct_response(i) = "j";
    end
end

% Count trial types
p.counts.test.k_presses = sum(strcmp(test_schedule.correct_response, "k"));
p.counts.test.j_presses = sum(strcmp(test_schedule.correct_response, "j"));
p.counts.test.no_response = sum(strcmp(test_schedule.correct_response, "none"));

fprintf('Phase 2 generated with %d trials.\n', height(test_schedule));
fprintf('  Response balance: %d Same(j), %d Diff(k), %d None\n', ...
    p.counts.test.j_presses, p.counts.test.k_presses, p.counts.test.no_response);

%% SAVE FINAL SUBJECT SETUP
%--------------------------------------------------------------------------
subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.encoding_schedule = encoding_schedule;
subject_data.test_schedule = test_schedule;

save(output_filename, 'subject_data');
fprintf('\nSetup file saved to:\n%s\n', output_filename);
