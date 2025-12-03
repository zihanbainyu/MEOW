clear;
clc;
close all;
dbstop if error;

subj_id = 501;
base_dir = '..';
subj_folder = sprintf('sub%03d', subj_id); 
results_dir = fullfile(base_dir, 'data', subj_folder);
tool_dir = fullfile('..', base_dir, 'analysis_toolbox', 'PCDM'); 
addpath(genpath(tool_dir));

% --- 2. load the 'in' struct we just made ---
input_file = fullfile(results_dir, sprintf('sub%03d_2_1_for_tool.mat', subj_id));
if ~exist(input_file, 'file'), error('run script 2 (translation layer) first.'); end
fprintf('loading "in" struct from: %s\n', input_file);
load(input_file, 'in');

% --- 3. RUN THE GODDAMN THING ---
fprintf('calling dataAnalysis(in)...\n');
[d, f] = fitModel(in);
fprintf('...model fitting complete. "d" and "f" structs are created.\n');

output_file = fullfile(results_dir, sprintf('sub%03d_2_1_MODEL_FITS.mat', subj_id));
fprintf('saving final "d" and "f" structs to: %s\n', output_file);
save(output_file, 'd', 'f', 'op'); % 'op' is created inside fitModel



