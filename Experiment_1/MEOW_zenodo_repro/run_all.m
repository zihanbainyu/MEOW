% run_all.m
% Reproduce every statistic reported in the manuscript, from the group-level data.
% Place this zenodo/ folder inside Experiment_1/ (next to results/, data/, stimulus/),
% or edit the three paths below to point at the deposited data.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'));

setenv('MEOW_DATA', fullfile(here,'..','results'));                                        % *.mat + image_similarity_AB.csv
setenv('MEOW_STIM', fullfile(here,'..','stimulus','stim_final'));                          % A/B images (for image similarity)
setenv('MEOW_MW',   fullfile(here,'..','data','eye_movement_data','group_eye_movement_combined.mat'));

% Optional: Bayes factors (reported for the null-control ANOVAs and entropy).
% addpath(genpath(fullfile(here,'..','toolbox','bayesFactor-master')));

repro_behavioral();   % Figure 2, Supplementary Figure 1
repro_pupil();        % Figure 3, Supplementary Figure 3   (needs pupil_window_means.mat; see README)
repro_gaze();         % Figure 4d-f
repro_cumulative();   % Figure 5
repro_entropy();      % Supplementary Information
repro_controls();     % Self-selection control (image similarity + lure bin)

fprintf('\n==== done. Compare printed values against the manuscript. ====\n');
