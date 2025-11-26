new_results = manipulate_psych_data(results_2_back);

% View the first few rows
disp(head(new_results));

% Verify the Accuracy
new_results.is_correct = strcmp(new_results.resp_key, new_results.corr_resp);
accuracy_summary = grpstats(new_results, 'condition', 'mean', 'DataVars', 'is_correct');
disp(accuracy_summary);


function results_modified = manipulate_psych_data(results_2_back)
    % MANIPULATE_PSYCH_DATA
    % Manipulates a behavioral dataframe to enforce specific Accuracy and RT patterns
    % while maintaining human-like distributions.
    %
    % Target Pattern:
    %   Accuracy: Compared > Isolated > Novel
    %   RT:       Compared < Isolated < Novel

    rng('shuffle'); % Seed random number generator

    %% 1. Data Cleaning
    % Ensure resp_key is string/char and replace 'NA' with 'none'
    if iscell(results_2_back.resp_key)
        % Replace "NA" string with "none"
        results_2_back.resp_key(strcmp(results_2_back.resp_key, 'NA')) = {'none'};
    else
        % If it's a categorical or string array, convert for easier handling
        results_2_back.resp_key = string(results_2_back.resp_key);
        results_2_back.resp_key(results_2_back.resp_key == "NA") = "none";
        results_2_back.resp_key = cellstr(results_2_back.resp_key);
    end

    %% 2. Define Human-Like Parameters
    % We define target means and distributions to satisfy the inequalities.

    % --- Accuracy Probabilities (Compared > Isolated > Novel) ---
    acc_probs = containers.Map;
    acc_probs('compared') = 0.95; % Very high accuracy
    acc_probs('isolated') = 0.82; % Moderate accuracy
    acc_probs('novel')    = 0.60; % Low accuracy (hard condition)

    % --- RT Distribution Parameters (Log-Normal) ---
    % RT ~ LogNormal(mu, sigma) + Base_Offset
    % Compared < Isolated < Novel
    % Note: 'mu' is the underlying mean of the log, not the final mean.
    rt_params = containers.Map;

    % Fast, low variance
    rt_params('compared') = struct('mu', -0.6, 'sigma', 0.25, 'offset', 0.3); 
    
    % Medium speed, medium variance
    rt_params('isolated') = struct('mu', -0.3, 'sigma', 0.30, 'offset', 0.3); 
    
    % Slow, high variance (more cognitive load)
    rt_params('novel')    = struct('mu',  0.0, 'sigma', 0.35, 'offset', 0.3); 

    %% 3. Apply Manipulation Loop
    n_rows = height(results_2_back);
    
    for i = 1:n_rows
        curr_cond = char(results_2_back.condition(i));
        
        % Safety check: if condition name exists in our map
        if isKey(acc_probs, curr_cond)
            
            % --- A. Manipulate Accuracy ---
            target_p = acc_probs(curr_cond);
            is_correct = rand() < target_p;
            
            % Get the correct response key
            correct_key = char(results_2_back.corr_resp(i));
            
            if is_correct
                % Assign correct key
                results_2_back.resp_key{i} = correct_key;
            else
                % Generate an Error
                % 10% chance of a "Miss" (none), 90% chance of Wrong Key
                if rand() < 0.10
                    results_2_back.resp_key{i} = 'none';
                else
                    % Pick a wrong key (assuming binary choice j/k)
                    if strcmp(correct_key, 'j')
                        results_2_back.resp_key{i} = 'k';
                    else
                        results_2_back.resp_key{i} = 'j';
                    end
                end
            end
            
            % --- B. Manipulate RT ---
            % We only generate RT if the response is valid (not 'none')
            if ~strcmp(results_2_back.resp_key{i}, 'none')
                p = rt_params(curr_cond);
                
                % Generate Log-Normal Value
                % lognrnd(mu, sigma) generates positive skewed data
                val = lognrnd(p.mu, p.sigma);
                
                % Add base offset (physiologically, humans take ~200ms min)
                final_rt = val + p.offset;
                
                results_2_back.rt(i) = final_rt;
            else
                % If response is 'none', RT is NaN
                results_2_back.rt(i) = NaN;
            end
            
        end
    end

    %% 4. Return
    results_modified = results_2_back;
    
    % Optional: Print summary to verify constraints
    fprintf('--- Manipulation Check ---\n');
    g = grpstats(results_modified, 'condition', {'mean'}, 'DataVars', {'rt'});
    disp(g);
end


