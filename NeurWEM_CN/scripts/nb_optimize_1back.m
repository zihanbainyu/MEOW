function [seqB_opt, midfix_new, report] = nb_optimize_1back(seqB, midfix_old, tm, nTry, TR, hp)
% NB_OPTIMIZE_1BACK  Best-of-N reordering + re-jitter of ONE 1-back block to
% maximise fMRI design efficiency, WITHOUT changing the experiment.
%
%   [seqB_opt, midfix_new, report] = nb_optimize_1back(seqB, midfix_old, tm, nTry, TR, hp)
%
% Constraints preserved (so the design stays valid):
%   * miniblock pairs stay intact and adjacent (compared/repeat = 2-trial
%     units; isolated = singletons) -- units are reordered, never split;
%   * isolated "A" is restored before its "B" (does not affect efficiency);
%   * counts, conditions, stimuli, timing model and run length are unchanged.
%
% Optimises the mean *relative* efficiency (vs the incoming design) across
% three contrasts: Similar-Same (pattern separation), Similar-New, Same-New.
%
%   seqB       : one block's rows; columns stim_id, condition, identity,
%                corr_resp, fix_duration (block/subj_id optional, carried along)
%   midfix_old : incoming mid-run break position (trial index) for this block
%   tm         : timing struct (see nb_timeline_1back)
%   nTry       : number of random valid candidate designs to score
%   TR, hp     : TR and high-pass cutoff (s) for the efficiency model
%
%   report.eff_old / eff_new : 1x3 efficiencies (old vs chosen), same contrasts
%   report.pct_gain          : 1x3 percent improvement
%   report.eff_samples       : nTry x 3 efficiencies of every candidate (for plots)
%   report.con               : contrast names
%
% Michelmann Lab @ NYU.
    C = [0 -1 1 0; 0 0 1 -1; 0 1 0 -1];          % Similar-Same, Similar-New, Same-New
    con = {'Similar-Same','Similar-New','Same-New'};

    cond = string(seqB.condition);
    resp = string(seqB.corr_resp);
    n = height(seqB);

    % --- reconstruct miniblock units as row-index groups (pairs stay intact) ---
    units = {}; i = 1;
    while i <= n
        if cond(i) == "isolated"
            units{end+1} = i; i = i + 1;                       %#ok<AGROW> singleton
        elseif (cond(i)=="compared" || cond(i)=="repeat") && resp(i)=="none"
            units{end+1} = [i i+1]; i = i + 2;                 %#ok<AGROW> 2-trial pair
        else
            error('nb_optimize_1back:badrow', ...
                'row %d starts on %s/%s -- expected a unit start', i, cond(i), resp(i));
        end
    end

    labAll = local_label(cond, resp);                          % label per original row

    % --- reference efficiency of the incoming design ---
    onset0 = nb_timeline_1back(seqB.fix_duration, midfix_old, tm);
    run0   = onset0(end) + tm.image_dur + tm.block_tail;
    eff_ref = nb_design_efficiency(onset0, labAll, run0, TR, hp, C);

    % --- best-of-N random valid search ---
    eff_samples = zeros(nTry, 3);
    best_obj = -inf; best = struct('ord',[],'fixd',[],'mf',[],'eff',[]);
    for t = 1:nTry
        up  = units(randperm(numel(units)));      % shuffle whole units
        ord = [up{:}]';                            % resulting row order
        fixd = tm.fix_dur + (rand(n,1)*2 - 1)*tm.fix_jitter;   % fresh jitter, U(0.5,1.0)
        usz = cellfun(@numel, up); cum = cumsum(usz); bnd = cum(1:end-1);
        cand = find(bnd >= 0.4*n & bnd <= 0.6*n);  % break near midpoint, on a boundary
        if isempty(cand), [~,cand] = min(abs(bnd - n/2)); end
        mf = bnd(cand(randi(numel(cand))));
        lab = labAll(ord);
        onset = nb_timeline_1back(fixd, mf, tm);
        rend  = onset(end) + tm.image_dur + tm.block_tail;
        eff = nb_design_efficiency(onset, lab, rend, TR, hp, C);
        eff_samples(t,:) = eff';
        obj = mean(eff ./ eff_ref);                % avg relative improvement
        if obj > best_obj
            best_obj = obj;
            best.ord = ord; best.fixd = fixd; best.mf = mf; best.eff = eff;
        end
    end

    % --- assemble optimised block; restore isolated A-before-B ---
    seqB_opt = seqB(best.ord, :);
    seqB_opt.fix_duration = best.fixd;
    seqB_opt = local_fix_iso_order(seqB_opt);
    midfix_new = best.mf;

    report = struct();
    report.eff_old     = eff_ref';
    report.eff_new     = best.eff';
    report.pct_gain    = 100*(best.eff' - eff_ref')./eff_ref';
    report.eff_samples = eff_samples;
    report.con         = {con};
end

function lab = local_label(cond, resp)
    cond = string(cond); resp = string(resp);
    lab = zeros(numel(cond),1);
    lab(cond=="isolated")                               = 4;   % New (singletons)
    lab(cond=="repeat"   & resp=="1")                   = 2;   % Same (repeat, 2nd trial)
    lab(cond=="compared" & resp=="2")                   = 3;   % Similar (compared, 2nd trial)
    lab((cond=="repeat"|cond=="compared") & resp=="none") = 1; % Encoding (first presentation)
    assert(all(lab>0), 'nb_optimize_1back: unlabeled row(s)');
end

function s = local_fix_iso_order(s)
% Ensure the isolated target "A" precedes its lure "B". Both are New events,
% so swapping their stimulus content leaves onsets/efficiency unchanged; this
% only restores the experimental encoding order (mirrors A_subject_setup.m).
    cond  = string(s.condition);
    iso   = find(cond=="isolated");
    if isempty(iso), return; end
    stim  = string(s.stim_id);
    ident = string(s.identity);
    key   = regexprep(stim(iso), '_(A|B)_', '_');    % pair key (drop the A/B tag)
    uk    = unique(key);
    for j = 1:numel(uk)
        grp = iso(key == uk(j));
        if numel(grp) ~= 2, continue; end
        rowA = grp(ident(grp)=="A");
        rowB = grp(ident(grp)=="B");
        if isempty(rowA) || isempty(rowB), continue; end
        if rowB < rowA                                % B currently earlier -> swap content
            s.stim_id([rowA rowB])  = s.stim_id([rowB rowA]);
            s.identity([rowA rowB]) = s.identity([rowB rowA]);
        end
    end
end
