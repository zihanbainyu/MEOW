%==========================================================================
%   A0_build_stim_pool.m  --  build the locked stimulus pool
%==========================================================================
% Author: Zihan Bai, Michelmann Lab at NYU
%
% Copies a fixed, fully deterministic stimulus set into ../stimulus/stim_pool:
%
%   - ALL L1 + L2 pairs (A + B): 180 per level = 360 pairs. A_subject_setup.m
%     splits each level into thirds -> compared / novel / repeat (60 each).
%   - Foils, A-image only, ids 1..N_FOIL. Foils are a separate object set,
%     only ever shown as single new/filler items (never lures), so no B_foil
%     is copied. N_FOIL = the exact number drawn per subject: 1-back fillers
%     (120) + recognition-new (80) + MST-new (60) = 260; 2-back padding is 0.
%   - practice pairs (prac_*) and legacy instruction images (instr_*).
%
% The selection is deterministic (all L-pairs, foils 1..N_FOIL) -- no RNG, so
% re-running reproduces the identical pool. It clears stim_pool first, then
% rebuilds. stim_final (the master set) is never touched. The full-screen
% ins_*.png slides live in ../stimulus, not here.
%
% Run from the scripts/ directory:   >> A0_build_stim_pool
%==========================================================================
function A0_build_stim_pool()

    N_FOIL = 260;   % exact foils drawn per subject (120 filler + 80 rec + 60 mst)
    SRC = fullfile('..', 'stimulus', 'stim_final');
    DST = fullfile('..', 'stimulus', 'stim_pool');

    assert(isfolder(SRC), 'source not found: %s', SRC);
    if ~isfolder(DST), mkdir(DST); end

    % --- clean rebuild: clear previously-copied pool files ---
    old = [dir(fullfile(DST, '*.png')); dir(fullfile(DST, '*.csv'))];
    for k = 1:numel(old), delete(fullfile(DST, old(k).name)); end
    if ~isempty(old)
        fprintf('  cleared %d existing files from stim_pool\n', numel(old));
    end

    manifest = {'kind,id'};

    % --- copy ALL L1 / L2 pairs (A + B) ---
    for lvl = ["l1", "l2"]
        A = dir(fullfile(SRC, sprintf('mst_*_A_%s.png', lvl)));
        for k = 1:numel(A)
            nameA = A(k).name;
            stem  = extractBefore(nameA, sprintf('_A_%s.png', lvl));  % e.g. 'mst_017'
            nameB = sprintf('%s_B_%s.png', stem, lvl);
            copyfile(fullfile(SRC, nameA), fullfile(DST, nameA));
            copyfile(fullfile(SRC, nameB), fullfile(DST, nameB));
            id = regexp(stem, '\d+', 'match', 'once');
            manifest{end+1} = sprintf('%s,%s', lvl, id); %#ok<AGROW>
        end
        fprintf('  %s: %d pairs copied\n', lvl, numel(A));
    end

    % --- copy foil A-images, ids 1..N_FOIL (no B) ---
    F = dir(fullfile(SRC, 'mst_*_A_foil.png'));
    fids = zeros(numel(F), 1);
    for k = 1:numel(F)
        fids(k) = str2double(regexp(F(k).name, '\d+', 'match', 'once'));
    end
    keep = find(fids >= 1 & fids <= N_FOIL);
    assert(numel(keep) == N_FOIL, ...
        'foils: found %d ids in 1..%d, expected %d', numel(keep), N_FOIL, N_FOIL);
    for k = keep'
        copyfile(fullfile(SRC, F(k).name), fullfile(DST, F(k).name));
        manifest{end+1} = sprintf('foil,%d', fids(k)); %#ok<AGROW>
    end
    fprintf('  foil: %d copied (ids 1-%d, A only)\n', N_FOIL, N_FOIL);

    % --- copy practice pairs + legacy instruction images ---
    n_support = 0;
    S = dir(fullfile(SRC, '*.png'));
    for k = 1:numel(S)
        nm = S(k).name;
        if startsWith(nm, 'prac_') || startsWith(nm, 'instr_')
            copyfile(fullfile(SRC, nm), fullfile(DST, nm));
            n_support = n_support + 1;
        end
    end
    fprintf('  support assets (prac/instr): %d copied\n', n_support);

    % --- write manifest (kept items only) ---
    fid = fopen(fullfile(DST, 'POOL_MANIFEST.csv'), 'w');
    fprintf(fid, '# stim_pool manifest  all L-pairs kept  foils 1-%d (A only)\n', N_FOIL);
    fprintf(fid, '%s\n', manifest{:});
    fclose(fid);

    fprintf('Pool built at %s  (manifest: POOL_MANIFEST.csv)\n', DST);
end
