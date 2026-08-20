function make_image_similarity(stim_dir, out_csv)
% MAKE_IMAGE_SIMILARITY  Objective A-B image similarity for every experimental
% pair (mnemonic-similarity bins l1/l2): Pearson correlation between the
% grayscale pixels of the A and B exemplars. Writes a lookup CSV with columns
% base_id, bin, pix_corr. Called automatically by repro_controls if missing.
if nargin < 1 || isempty(stim_dir), stim_dir = fullfile('..','stimulus','stim_final'); end
if nargin < 2 || isempty(out_csv),  out_csv  = fullfile('..','results','image_similarity_AB.csv'); end

A = dir(fullfile(stim_dir, 'mst_*_A_l*.png'));
rows = cell(0,3);
for i = 1:numel(A)
    b = fullfile(stim_dir, strrep(A(i).name, '_A_l', '_B_l'));
    if ~isfile(b), continue; end
    tok = regexp(A(i).name, '(mst_\d+)_A_(l\d)', 'tokens', 'once');
    if isempty(tok), continue; end
    va = double(togray(imread(fullfile(stim_dir, A(i).name))));
    vb = double(togray(imread(b)));
    if ~isequal(size(va), size(vb)), vb = double(imresize(togray(imread(b)), size(va))); end
    rows(end+1,:) = {tok{1}, tok{2}, corr(va(:), vb(:))}; %#ok<AGROW>
end
writetable(cell2table(rows, 'VariableNames', {'base_id','bin','pix_corr'}), out_csv);
fprintf('wrote %s (%d pairs)\n', out_csv, size(rows,1));
end

function g = togray(im)
if ndims(im) == 3, g = rgb2gray(im); else, g = im; end
end
