% run_all.m
% reproduce every statistic in the manuscript from the consolidated group data.
% layout: scripts/ next to data/ (edit MEOW_DATA below if different).

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'), genpath(fullfile(here,'lib','bayesFactor-master')));
setenv('MEOW_DATA', fullfile(here,'..','data'));

repro_behavior();   % figure 2, supplementary figure 1
repro_pupil();      % figure 3, supplementary figure 3
repro_gaze();       % figure 4, figure 5, supplementary entropy, controls

fprintf('\ndone\n');
