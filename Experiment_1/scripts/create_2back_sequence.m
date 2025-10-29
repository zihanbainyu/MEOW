function [final_sequence, sequence_log] = create_2back_sequence()
    % ---
    % This script generates a 92-trial, overlapping 2-back sequence
    % with a perfect 30/30/30 split of Repeat/Lure/Foil test events
    % based on 90 stimulus pairs.
    % ---
    
    % 1. DEFINE STIMULI
    % Create 90 pairs (30 of each category)
    comp_pairs = cell(30, 2);
    iso_pairs = cell(30, 2);
    nov_pairs = cell(30, 2);
    
    for i = 1:30
        comp_pairs{i, 1} = sprintf('C_P%d', i);
        comp_pairs{i, 2} = sprintf('C_L%d', i);
        
        iso_pairs{i, 1} = sprintf('I_P%d', i);
        iso_pairs{i, 2} = sprintf('I_L%d', i);
        
        nov_pairs{i, 1} = sprintf('N_P%d', i);
        nov_pairs{i, 2} = sprintf('N_L%d', i);
    end
    
    % Combine all pairs
    all_pairs = [comp_pairs; iso_pairs; nov_pairs];
    
    % Create a lookup map for lures: {P1: L1, L1: P1}
    lure_map = containers.Map;
    stimulus_pool = cell(180, 1); % Our "deck" of 180 stimuli
    
    for i = 1:size(all_pairs, 1) % 1 to 90
        p = all_pairs{i, 1};
        l = all_pairs{i, 2};
        lure_map(p) = l;
        lure_map(l) = p;
        stimulus_pool{(i*2)-1} = p;
        stimulus_pool{i*2} = l;
    end
    
    % Shuffle the stimulus pool
    stimulus_pool = stimulus_pool(randperm(180));

    % 2. DEFINE EVENT ORDER (90 test trials)
    event_types = [repmat({'REPEAT'}, 1, 30), ...
                   repmat({'LURE'}, 1, 30), ...
                   repmat({'FOIL'}, 1, 30)]; % 1x90 cell array

    % Simple non-predictive check: shuffle until no 3-in-a-row
    while true
        event_types = event_types(randperm(90));
        % Create a simple string representation ('R', 'L', 'F')
        s = cellfun(@(c) c(1), event_types);
        if isempty(strfind(s, 'RRR')) && isempty(strfind(s, 'LLL')) && isempty(strfind(s, 'FFF'))
            break;
        end
    end

    % 3. PRIME & BUILD THE SEQUENCE
    % Pre-allocate the final sequence and a human-readable log
    final_sequence = cell(92, 1);
    sequence_log = cell(92, 1);
    
    % Prime with 2 stimuli from the pool (pop from the end)
    s0 = stimulus_pool{end}; stimulus_pool(end) = [];
    s1 = stimulus_pool{end}; stimulus_pool(end) = [];
    
    final_sequence{1} = s0;
    final_sequence{2} = s1;
    
    sequence_log{1} = sprintf('Trial  1: %-6s (Prime 1)', s0);
    sequence_log{2} = sprintf('Trial  2: %-6s (Prime 2)', s1);

    % 4. ITERATE AND BUILD TRIALS 3-92
    for i = 1:90 % This loop generates trials 3 through 92
        event_to_create = event_types{i};
        t_minus_2_stim = final_sequence{i}; % i=1 -> T3, tests T1
        current_trial_num = i + 2;
        
        next_stim = '';
        event_label = '';

        if strcmp(event_to_create, 'REPEAT')
            next_stim = t_minus_2_stim;
            event_label = 'REPEAT';

        elseif strcmp(event_to_create, 'LURE')
            % Find the paired lure for the t-2 stimulus
            next_stim = lure_map(t_minus_2_stim);
            event_label = 'LURE';

        elseif strcmp(event_to_create, 'FOIL')
            lure_of_t_minus_2 = lure_map(t_minus_2_stim);
            
            % Find a NEW stimulus from the pool that is a valid foil
            found = false;
            for j = 1:length(stimulus_pool)
                candidate = stimulus_pool{j};
                if ~strcmp(candidate, t_minus_2_stim) && ~strcmp(candidate, lure_of_t_minus_2)
                    next_stim = candidate;
                    stimulus_pool(j) = []; % Remove from pool
                    found = true;
                    break;
                end
            end
            
            if ~found
                % Fallback: Pool is empty. Pick a random stim from all 180
                % that is still a valid foil for this specific trial.
                all_stims = keys(lure_map);
                invalid_idx = strcmp(all_stims, t_minus_2_stim) | strcmp(all_stims, lure_of_t_minus_2);
                available_foils = all_stims(~invalid_idx);
                next_stim = available_foils{randi(length(available_foils))};
                % Warning for the user
                if i == 90
                     warning('Foil pool exhausted. Re-using stimuli for last few foil trials.')
                end
            end
            event_label = 'FOIL';
        end
        
        % Add the new stimulus to the sequence
        final_sequence{current_trial_num} = next_stim;
        
        % Log the event
        log_entry = sprintf('Trial %2d: %-6s (Test against T-%-2d / %-6s) -> %-6s (Response: %c)', ...
                            current_trial_num, next_stim, current_trial_num-2, t_minus_2_stim, event_label, event_label(1));
        sequence_log{current_trial_num} = log_entry;
    end
    return final_sequence;
end
