function [onset, run_end] = nb_timeline_1back(fixd, midfix_after, tm)
% NB_TIMELINE_1BACK  Image onsets (s) on the run clock for a 1-back run.
%   [onset, run_end] = nb_timeline_1back(fixd, midfix_after, tm)
%   fixd         : per-trial jittered fixation duration (s), in run order
%   midfix_after : trial index after which the 6 s mid-run break is inserted
%   tm           : timing struct (block_lead_in, image_dur, block_midfix, block_tail)
%   onset        : image onset time (s) of each trial, relative to first trigger
%   run_end      : total run duration (s) incl. lead-in, break, and tail
%
% Michelmann Lab @ NYU. Used by nb_optimize_1back / design-efficiency scoring.
    n = numel(fixd);
    onset = zeros(n,1);
    t = tm.block_lead_in;                 % settle period before the first trial
    for i = 1:n
        t = t + fixd(i);                  % jittered fixation before image i
        onset(i) = t;
        t = t + tm.image_dur;             % image on screen
        if i == midfix_after
            t = t + tm.block_midfix;      % mid-run fixation break
        end
    end
    run_end = onset(end) + tm.image_dur + tm.block_tail;
end
