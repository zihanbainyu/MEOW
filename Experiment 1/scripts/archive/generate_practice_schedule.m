function practice_schedule = generate_practice_schedule(p)
%==========================================================================
%                  Create Practice Encoding Schedule (1-back task)
%==========================================================================
% Author: Zihan Bai, Michelmann Lab @ NYU
%==========================================================================

base_dir = '..';
p.stim_dir = fullfile(base_dir, 'stimulus/stim_final/');
p.setup_dir = fullfile(base_dir, 'subj_setup/');

all_pair_ids = [1:10];  % exclude 11
rng('shuffle');

nPairs = 10;
pair_ids = all_pair_ids(randperm(numel(all_pair_ids), nPairs));

cond_types = repelem(["repeat", "similar", "new"], ceil(nPairs/3));
cond_types = cond_types(1:nPairs);
cond_types = cond_types(randperm(numel(cond_types)));

% preallocate structure
practice_schedule = struct( ...
    'block', {}, ...
    'trial', {}, ...
    'stimulus_id', {}, ...
    'condition', {}, ...
    'pair_id', {}, ...
    'correct_response', {}, ...
    'fix_duration', {});

trial = 0;

% -------------------
% Build trial list
% -------------------
for i = 1:nPairs
    pid = pair_ids(i);
    base_name = sprintf('practice_%02d', pid);

    % --- first item (A) ---
    trial = trial + 1;
    practice_schedule(trial).block = 0;
    practice_schedule(trial).trial = trial;
    practice_schedule(trial).stimulus_id = sprintf('%s_a.png', base_name);
    practice_schedule(trial).condition = 'base';
    practice_schedule(trial).pair_id = pid;
    practice_schedule(trial).correct_response = 'none';
    practice_schedule(trial).fix_duration = 0.75 + rand()*0.25;

    % --- second item depends on condition ---
    trial = trial + 1;
    practice_schedule(trial).block = 0;
    practice_schedule(trial).trial = trial;
    practice_schedule(trial).pair_id = pid;
    practice_schedule(trial).fix_duration = 0.75 + rand()*0.25;

    switch cond_types(i)
        case 'repeat'
            practice_schedule(trial).stimulus_id = sprintf('%s_a.png', base_name);
            practice_schedule(trial).condition = 'repeat';
            practice_schedule(trial).correct_response = 'j';

        case 'similar'
            practice_schedule(trial).stimulus_id = sprintf('%s_b.png', base_name);
            practice_schedule(trial).condition = 'similar';
            practice_schedule(trial).correct_response = 'k';

        case 'new'
            % pick another A image as "new" (unrelated)
            other_ids = setdiff(pair_ids, pid);
            new_pid = other_ids(randi(numel(other_ids)));
            practice_schedule(trial).stimulus_id = sprintf('practice_%02d_a.png', new_pid);
            practice_schedule(trial).condition = 'new';
            practice_schedule(trial).correct_response = 'none';
    end
end

% -------------------
% Convert to table & save
% -------------------
practice_schedule = struct2table(practice_schedule);

fprintf('practice schedule created: %d trials (%d pairs)\n', ...
    height(practice_schedule), nPairs);

save(fullfile(p.setup_dir, 'practice_schedule.mat'), 'practice_schedule');

end
