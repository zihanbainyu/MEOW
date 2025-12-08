clear; clc; close all;

% config
sub_list = [501, 601, 602, 603, 604, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 617];
tasks = [1, 2];
task_lbls = {'1_back', '2_back'};
base_dir = '..';
out_dir = fullfile(base_dir, 'data', 'eye_movement_data');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% containers
[all_fix, all_sac, all_blk, all_sum] = deal({});

% loop subjects
for s_idx = 1:length(sub_list)
    sid = sub_list(s_idx);
    s_fldr = sprintf('sub%03d', sid);
    f_path = fullfile(base_dir, 'data', s_fldr);
    
    % trial accumulators reset per subject, persist per task across blocks
    acc_tr_t1 = 0; % 1-back accumulator
    acc_tr_t2 = 0; % 2-back accumulator
    
    % loop blocks
    for b = 1:4
        
        % loop tasks
        for t_idx = 1:length(tasks)
            tsk = tasks(t_idx);
            t_lbl = task_lbls{t_idx};
            fn = fullfile(f_path, sprintf('%03d_%01d_%01d.asc', sid, tsk, b));
            
            if ~isfile(fn), continue; end
            fid = fopen(fn);
            
            % dense parsing
            [f_dat, s_dat, b_dat] = deal(cell(2000,1)); 
            [n_f, n_s, n_b, blk_max_raw_tr] = deal(0); % blk_max_raw_tr tracks max trial in this file
            [curr_tr, st_t, st_id] = deal(0, NaN, 'NA');
            
            while ~feof(fid)
                ln = fgetl(fid);
                if startsWith(ln, 'MSG')
                    if contains(ln, 'TRIALID')
                        raw_tr = sscanf(ln, 'MSG %*d TRIALID %d');
                        if raw_tr > 0
                            blk_max_raw_tr = max(blk_max_raw_tr, raw_tr); % update max raw trial for this file
                            
                            % apply task-specific accumulation
                            if tsk == 1 
                                curr_tr = raw_tr + acc_tr_t1;
                            else % tsk == 2
                                curr_tr = raw_tr + acc_tr_t2;
                            end
                            
                            all_sum{end+1,1} = {sid, b, t_lbl, curr_tr, 'NA', NaN, 0, 0, 0};
                        end
                    elseif contains(ln, 'STIM_ONSET') && curr_tr > 0
                        tmp = textscan(ln, '%*s %f %*s %s');
                        st_t = tmp{1}; st_id = tmp{2}{1};
                        all_sum{end,1}{5} = st_id; all_sum{end,1}{6} = st_t;
                    end
                elseif curr_tr > 0
                    if startsWith(ln, 'EFIX')
                        d = sscanf(ln, 'EFIX %c %d %d %d %f %f %d');
                        if numel(d) == 7
                            n_f = n_f+1;
                            f_dat{n_f} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7)};
                            all_sum{end,1}{7} = all_sum{end,1}{7} + 1;
                        end
                    elseif startsWith(ln, 'ESACC')
                        d = sscanf(strrep(ln, '.', 'NaN'), 'ESACC %c %d %d %d %f %f %f %f'); 
                        if numel(d) == 8
                            n_s = n_s+1;
                            s_dat{n_s} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4), d(5), d(6), d(7), d(8)};
                            all_sum{end,1}{8} = all_sum{end,1}{8} + 1;
                        end
                    elseif startsWith(ln, 'EBLINK')
                        d = sscanf(ln, 'EBLINK %c %d %d %d');
                        if numel(d) == 4
                            n_b = n_b+1;
                            b_dat{n_b} = {sid, b, t_lbl, curr_tr, st_id, char(d(1)), d(2), d(3), d(4)};
                            all_sum{end,1}{9} = all_sum{end,1}{9} + 1;
                        end
                    end
                end
            end
            fclose(fid);
            
            % update task-specific accumulator for next block
            if tsk == 1
                acc_tr_t1 = acc_tr_t1 + blk_max_raw_tr;
            else % tsk == 2
                acc_tr_t2 = acc_tr_t2 + blk_max_raw_tr;
            end
            
            % append to master
            if n_f > 0, all_fix = [all_fix; f_dat(1:n_f)]; end
            if n_s > 0, all_sac = [all_sac; s_dat(1:n_s)]; end
            if n_b > 0, all_blk = [all_blk; b_dat(1:n_b)]; end
        end
    end

    % Define variable names for saving
    v_fix = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','x','y','pupil'};
    v_sac = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur','sx','sy','ex','ey'};
    v_blk = {'subj_id','block','task','trial_id','stim_id','eye','onset','offset','dur'};
    v_sum = {'subj_id','block','task','trial_id','stim_id','onset_time','n_fix','n_sac','n_blink'};
    
    % convert and save for the current subject
    if ~isempty(all_fix)
        T = cell2table(vertcat(all_fix{:}), 'Var', v_fix);
        writetable(T, fullfile(out_dir, sprintf('sub%03d_fixations.csv', sid)));
    end
    
    if ~isempty(all_sac)
        T = cell2table(vertcat(all_sac{:}), 'Var', v_sac);
        writetable(T, fullfile(out_dir, sprintf('sub%03d_saccades.csv', sid)));
    end
    
    if ~isempty(all_blk)
        T = cell2table(vertcat(all_blk{:}), 'Var', v_blk);
        writetable(T, fullfile(out_dir, sprintf('sub%03d_blinks.csv', sid)));
    end
    
    if ~isempty(all_sum)
        T = cell2table(vertcat(all_sum{:}), 'Var', v_sum);
        writetable(T, fullfile(out_dir, sprintf('sub%03d_summary.csv', sid)));
    end
    
    [all_fix, all_sac, all_blk, all_sum] = deal({});
    fprintf('sub%03d done. 1-back trials accumulated up to %d, 2-back trials up to %d.\n', sid, acc_tr_t1, acc_tr_t2);
end


% define file types to stack
ftypes = {'fixations.csv', 'saccades.csv', 'blinks.csv', 'summary.csv'};

for i = 1:length(ftypes)
    pattern = ['*' ftypes{i}];
    flist = dir(fullfile(out_dir, pattern));
    
    if isempty(flist), continue; end
    
    alltables = cell(length(flist), 1);
    for k = 1:length(flist)
        alltables{k} = readtable(fullfile(out_dir, flist(k).name), 'PreserveVariableNames', true);
    end
    
    Tmaster = vertcat(alltables{:});
    
    writetable(Tmaster, fullfile(out_dir, ['group_' ftypes{i}]));
end