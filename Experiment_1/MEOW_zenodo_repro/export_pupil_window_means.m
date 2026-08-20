function export_pupil_window_means()
% EXPORT_PUPIL_WINDOW_MEANS  Save the per-subject 1.0-1.5 s pupil window means
% that repro_pupil needs, WITHOUT modifying any existing script.
%
% Usage:
%   1) run  S_pupil_cluster           % performs the cluster-based window selection;
%                                      % leaves pup_lure_corr/incorr/cond, pup_targ_cond
%                                      % in the workspace
%   2) run  export_pupil_window_means % writes results/pupil_window_means.mat
%
% It reads the variables from the base workspace, so S_pupil_cluster.m is left
% exactly as it is.

DATA = getenv('MEOW_DATA'); if isempty(DATA), DATA = fullfile('..','results'); end
lure_corr   = evalin('base', 'pup_lure_corr');
lure_incorr = evalin('base', 'pup_lure_incorr');
lure_cond   = evalin('base', 'pup_lure_cond');
targ_cond   = evalin('base', 'pup_targ_cond');
save(fullfile(DATA,'pupil_window_means.mat'), 'lure_corr','lure_incorr','lure_cond','targ_cond');
fprintf('saved %s\n', fullfile(DATA,'pupil_window_means.mat'));
end
