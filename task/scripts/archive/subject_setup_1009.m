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
p.nIsolated_Target = 13;
p.nIsolated_Lure = 13;
p.nIsolated_Both = 14;
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
p.timing.image_dur = 1; 
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

idx_start=idx_end+1; idx_end=idx_start+p.nIsolated_Target-1;
iso_target_pairs=final_pairs(idx_start:idx_end,:);

idx_start=idx_end+1; idx_end=idx_start+p.nIsolated_Lure-1;  
iso_lure_pairs=final_pairs(idx_start:idx_end,:);

idx_start=idx_end+1; idx_end=idx_start+p.nIsolated_Both-1;   
iso_both_pairs=final_pairs(idx_start:idx_end,:);

idx_start=idx_end+1; idx_end=idx_start+p.nNovel-1;           
novel_pairs=final_pairs(idx_start:idx_end,:);

p.stim.comparison=comparison_pairs; 
p.stim.iso_target=iso_target_pairs;
p.stim.iso_lure=iso_lure_pairs; 
p.stim.iso_both=iso_both_pairs; 
p.stim.novel=novel_pairs;


%% BUILD PHASE 1 (ENCODING)
%--------------------------------------------------------------------------
fprintf('building phase 1 with interleaving spacers...\n');

encoding_blocks = {}; % all mnemonically significant blocks
spacer_blocks = {};   % all mnemonically neutral blocks

% create 94 ENCODING blocks (2 trials each)

% comparison encoding events(40 blocks)
for i=1:p.nComparison, encoding_blocks{end+1} = { ...
    comparison_pairs.Target(i), "comparison", "target", "new", "none"; ...
    comparison_pairs.Lure(i),   "comparison", "lure",   "lure", "k" };end

% isolated encoding events (13 T-iso + 13 L-iso + 14 T-both + 14 L-both = 54 blocks)
iso_event_foils = all_foils(1:54); all_foils(1:54) = [];
for i=1:p.nIsolated_Target, encoding_blocks{end+1} = {iso_event_foils(1),"iso_target","filler","new","none"; iso_target_pairs.Target(i),"iso_target","target","new","none"}; iso_event_foils(1)=[]; end
for i=1:p.nIsolated_Lure,   encoding_blocks{end+1} = {iso_event_foils(1),"iso_lure",  "filler","new","none"; iso_lure_pairs.Lure(i),  "iso_lure",  "lure",  "new","none"}; iso_event_foils(1)=[]; end
for i=1:p.nIsolated_Both
    encoding_blocks{end+1} = {iso_event_foils(1),"iso_both",  "filler","new","none"; iso_both_pairs.Target(i),"iso_both",  "target","new","none"}; iso_event_foils(1)=[];
    encoding_blocks{end+1} = {iso_event_foils(1),"iso_both",  "filler","new","none"; iso_both_pairs.Lure(i),"iso_both", "lure","new","none"};   iso_event_foils(1)=[];
end

% create 94 SPACER blocks (2 trials each)

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

assert(numel(encoding_blocks)==94 && numel(spacer_blocks)==94, 'Block count error!');

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

%% BUILD PHASE 2 TEST SCHEDULE
%--------------------------------------------------------------------------
fprintf('building balanced and spaced phase 2 (test) sequence...\n');
phase2_crit_blocks = {}; % T-F-L Lure Triads
phase2_safe_blocks = {}; % F-F-F Repeat Triads

total_iso_pairs = p.nIsolated_Target + p.nIsolated_Lure + p.nIsolated_Both;
iso_all_pairs = [iso_target_pairs; iso_lure_pairs; iso_both_pairs];
iso_all_conds = [repmat("isolated_target",p.nIsolated_Target,1); repmat("isolated_lure",p.nIsolated_Lure,1); repmat("isolated_both",p.nIsolated_Both,1)];

% create 120 Lure Triads ('k' press) ---
intervening_foils = all_foils(1:120); all_foils(1:120) = [];
for i=1:p.nComparison, phase2_crit_blocks{end+1} = {comparison_pairs.Target(i),"comparison","target","new","none"; intervening_foils(1),"unrelated_filler","filler","new","none"; comparison_pairs.Lure(i),"comparison","lure","lure","k"}; intervening_foils(1)=[]; end
for i=1:total_iso_pairs, phase2_crit_blocks{end+1} = {iso_all_pairs.Target(i),iso_all_conds(i),"target","new","none"; intervening_foils(1),"unrelated_filler","filler","new","none"; iso_all_pairs.Lure(i),iso_all_conds(i),"lure","lure","k"}; intervening_foils(1)=[]; end
for i=1:p.nNovel,      phase2_crit_blocks{end+1} = {novel_pairs.Target(i),"novel","target","new","none"; intervening_foils(1),"unrelated_filler","filler","new","none"; novel_pairs.Lure(i),"novel","lure","lure","k"};       intervening_foils(1)=[]; end

% create 120 repeat Triads ('j' press) to balance ---
repeat_foil_A = all_foils(1:120); all_foils(1:120) = [];
repeat_foil_B = all_foils(1:120); all_foils(1:120) = [];
for i=1:p.nTest_Repeat_Events
    phase2_safe_blocks{end+1} = {repeat_foil_A(i),"foil_repeat","filler","new","none"; repeat_foil_B(i),"foil_repeat","filler","new","none"; repeat_foil_A(i),"foil_repeat","filler","repeat","j"};
end

assert(numel(phase2_crit_blocks)==120 && numel(phase2_safe_blocks)==120, 'Test block count error!');

% asemble with interleaving
shuffled_crit = phase2_crit_blocks(randperm(numel(phase2_crit_blocks)));
shuffled_safe = phase2_safe_blocks(randperm(numel(phase2_safe_blocks)));
final_ordered_blocks = reshape([shuffled_safe; shuffled_crit], 1, []);
final_test_list = vertcat(final_ordered_blocks{:});

% finalize and verify
test_schedule = cell2table(final_test_list, 'VariableNames', {'stimulus_id', 'condition', 'role', 'trial_type_designed', 'correct_response'});
test_schedule.nback_target_id = strings(height(test_schedule), 1);
test_schedule.trial_type_final = test_schedule.trial_type_designed;
for i = 3:height(test_schedule)
    test_schedule.nback_target_id(i) = test_schedule.stimulus_id(i-2);
    if string(test_schedule.stimulus_id(i))==string(test_schedule.nback_target_id(i)), test_schedule.trial_type_final(i)="repeat"; test_schedule.correct_response(i)="j"; end
end

p.counts.test.k_presses = sum(strcmp(test_schedule.correct_response, "k"));
p.counts.test.j_presses = sum(strcmp(test_schedule.correct_response, "j"));
p.counts.test.no_response = sum(strcmp(test_schedule.correct_response, "none"));

fprintf('Phase 2 generated with %d trials.\n', height(test_schedule));
fprintf('  Response balance: %d Same(j), %d Diff(k), %d None\n', ...
    p.counts.test.j_presses, p.counts.test.k_presses, p.counts.test.no_response);

%% SAVE THE FINAL SCHEDULE FILE
%--------------------------------------------------------------------------
subject_data.subj_id = subj_id;
subject_data.parameters = p;
subject_data.encoding_schedule = encoding_schedule;
subject_data.test_schedule = test_schedule;

save(output_filename, 'subject_data');
fprintf('\nsetup file saved to:\n%s\n', output_filename);