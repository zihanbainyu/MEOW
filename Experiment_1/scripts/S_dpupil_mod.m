clear; clc; close all; dbstop if error;
fprintf('--- processing subj\n');

% setup
subj_id = 611; base_dir = '..'; 
r_dir = fullfile(base_dir, 'data', sprintf('sub%03d', subj_id));
C = struct('W_CM',60, 'W_PX',1920, 'H_PX',1080, 'D_CM',100, 'SR',1000, 'TRIAL_DUR',1.5, 'BASELINE_DUR',0.5);
load(fullfile(r_dir, sprintf('sub%03d_concat.mat', subj_id)), 'final_data_output');
in = struct('xPos',{{}},'yPos',{{}},'pupilArea',{{}},'sampleRate',{{}},'startInds',{{}},'trialTypes',{{}}, 'baselineInds',{{}});

for r = 1:4
    fprintf('run %d: ', r);
    % 1. load/parse asc
    asc_f = fullfile(r_dir, sprintf('%03d_2_%01d.asc', subj_id, r));
    mat_f = fullfile(r_dir, sprintf('%03d_2_%01d_raw_v2.mat', subj_id, r));
    if isfile(mat_f), load(mat_f, 'raw'); fprintf('loaded raw. '); 
    else, raw = parse_asc(asc_f); save(mat_f, 'raw'); fprintf('parsed asc. '); end
    
    % 2. behavior
    behav = final_data_output.results_2_back_all(final_data_output.results_2_back_all.block==r,:);
    behav.resp_key = cellstr(behav.resp_key); behav.resp_key(strcmp(behav.resp_key,'NA'))={'none'};
    behav.correct = strcmp(cellstr(behav.corr_resp), behav.resp_key);
    n_t = height(behav);
    
    % 3. dva conversion
    [x_d, y_d] = px2dva(raw.xp, raw.yp, C);

    in.xPos{r} = x_d; in.yPos{r} = y_d; in.pupilArea{r} = raw.pp; in.sampleRate{r} = C.SR;
    
    % 4. epoching (baseline + trial)
    inds = nan(n_t, 2); base_inds = nan(n_t, 2); types = zeros(1, n_t);
    ev_msg = raw.ev.msg; ev_ts = raw.ev.ts;
    
    for i = 1:n_t
        tid = find(contains(ev_msg, sprintf('TRIALID %d', i)), 1);
        if isempty(tid), continue; end
        
        % Find STIM_ONSET
        onset_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'STIM_ONSET'), 1);
        if isempty(onset_idx), continue; end
        t_onset = ev_ts(tid + onset_idx - 1);
        
        % Find TRIAL_RESULT for this trial
        result_idx = find(contains(ev_msg(tid:min(tid+50,end)), 'TRIAL_RESULT'), 1);
        if isempty(result_idx), continue; end
        t_trial_end = ev_ts(tid + result_idx - 1);
        
        % Baseline: 500ms before stimulus onset
        t_baseline_start = t_onset - (C.BASELINE_DUR * 1000);
        
        % Find sample indices for baseline
        idx_base_s = find(raw.ts >= t_baseline_start, 1);
        idx_base_e = find(raw.ts < t_onset, 1, 'last');
        
        % Find sample indices for trial (from STIM_ONSET to TRIAL_RESULT)
        idx_trial_s = find(raw.ts >= t_onset, 1);
        idx_trial_e = find(raw.ts <= t_trial_end, 1, 'last');
        
        % Validate baseline
        if ~isempty(idx_base_s) && ~isempty(idx_base_e) && idx_base_e > idx_base_s
            base_inds(i,:) = [idx_base_s, idx_base_e];
        end
        
        % Validate trial
        if ~isempty(idx_trial_s) && ~isempty(idx_trial_e) && idx_trial_e > idx_trial_s
            inds(i,:) = [idx_trial_s, idx_trial_e];
        end
        
        % Trial types (correct trials only)
        if behav.correct(i)
            c = behav.condition(i); g = behav.goal(i); cr = behav.corr_resp(i);
            if     c=="compared" && g=="A-B" && cr=="k", types(i)=1;  % primary compared
            elseif c=="isolated" && g=="A-B" && cr=="k", types(i)=2;  % primary isolated
            elseif c=="novel"    && g=="A-B" && cr=="k", types(i)=3;  % primary novel
            elseif c=="compared" && g=="A-A" && cr=="j", types(i)=4;  % secondary compared
            elseif c=="isolated" && g=="A-A" && cr=="j", types(i)=5;  % secondary isolated
            elseif c=="novel"    && g=="A-A" && cr=="j", types(i)=6;  % secondary novel
            end
        end
    end
    % 
    % % exclude trials with bad epochs
    bad = any(isnan(inds), 2) | any(isnan(base_inds), 2);
    in.startInds{r} = inds(~bad,:); 
    in.baselineInds{r} = base_inds(~bad,:);
    in.trialTypes{r} = types(~bad);

    fprintf('done (%d trials, %d excluded).\n', sum(~bad), sum(bad));
end
save(fullfile(r_dir, sprintf('sub%03d_input.mat', subj_id)), 'in');

% model fit
fprintf('running model fit...\n');
pcdm_dir = fullfile(base_dir, 'scripts/PCDM-main');
addpath(genpath(pcdm_dir));
[d, f] = fitModel(in);
save(fullfile(r_dir, sprintf('sub%03d_model_pred.mat', subj_id)), 'd', 'f');
fprintf('--- all done ---\n');

% functions
function r = parse_asc(f)
    fid = fopen(f); N=5e6; ts=zeros(N,1); xp=nan(N,1); yp=nan(N,1); pp=nan(N,1);
    ev_ts=zeros(5000,1); ev_msg=cell(5000,1); c=0; ec=0;
    while ~feof(fid)
        l = fgetl(fid); if isempty(l), continue; end
        if l(1)>='0' && l(1)<='9' % sample
            c=c+1; d = sscanf(l, '%d %f %f %d');
            if numel(d)==4, ts(c)=d(1); xp(c)=d(2); yp(c)=d(3); pp(c)=d(4);
            else, ts(c)=sscanf(l,'%d',1); end % blink/partial
        elseif startsWith(l, 'MSG') % event
            ec=ec+1; dat = textscan(l, 'MSG %d %s', 'Delimiter', ''); % quick parse
            ev_ts(ec)=dat{1}; ev_msg{ec}=l;
        end
    end
    fclose(fid);
    r.ts=ts(1:c); r.xp=xp(1:c); r.yp=yp(1:c); r.pp=pp(1:c);
    r.ev.ts=ev_ts(1:ec); r.ev.msg=ev_msg(1:ec);
end

function [xd, yd] = px2dva(x, y, C)
    px2cm = C.W_PX / C.W_CM; cx = C.W_PX/2; cy = C.H_PX/2;
    xc = nan(size(x)); yc = nan(size(y));
    ok = (x>0 & x<C.W_PX & y>0 & y<C.H_PX); % valid filter
    xc(ok) = (x(ok)-cx)/px2cm; yc(ok) = (y(ok)-cy)/px2cm;
    xd = atan(xc/C.D_CM) * (180/pi); yd = atan(yc/C.D_CM) * (180/pi);
end