%==========================================================================
%   prac_example_images.m
%==========================================================================
% Return the exact same/similar example images used in the practice, so the
% post-practice recap (run via run_text_instructions with <<IMG:...>> markers)
% shows the very objects the subject just saw.
%
% The indices below mirror the practice schedules:
%   1-back (gen_1_back_practice): same = A(1) shown twice; similar = A(2)->B(2)
%   2-back (gen_2_back_practice): same = A(10) shown twice; similar = A(13)->B(13)
% Keep these in sync if the practice sequences change.
%
% Returns a struct: imgs.same    = {fileA, fileA}   (two identical images)
%                   imgs.similar = {fileA, fileB}   (a lure pair)
%==========================================================================
function imgs = prac_example_images(p, task)

    A = dir(fullfile(p.stim_dir, 'prac_*_A.png'));
    B = dir(fullfile(p.stim_dir, 'prac_*_B.png'));
    A = string(sort({A.name}));   % zero-padded names -> alphabetical == numeric
    B = string(sort({B.name}));

    switch lower(task)
        case '1back'
            i_same = 1;   i_sim = 2;
        case '2back'
            i_same = 10;  i_sim = 13;
        otherwise
            error('prac_example_images: unknown task "%s".', task);
    end

    need = max(i_same, i_sim);
    if numel(A) < need || numel(B) < i_sim
        error(['prac_example_images: not enough prac images in %s ' ...
               '(need A>=%d, B>=%d; found A=%d, B=%d).'], ...
               p.stim_dir, need, i_sim, numel(A), numel(B));
    end

    imgs.same    = {char(A(i_same)), char(A(i_same))};   % same object twice
    imgs.similar = {char(A(i_sim)),  char(B(i_sim))};    % lure pair (A then B)
end
