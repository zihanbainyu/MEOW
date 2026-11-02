% ============================================================
% BrainNet Viewer 2.0-beta batch render — NeurWEM 2-back vs 1-back FC
% NOTE: 2.0 dropped BrainNet_MapCfg. It renders from a saved Scene JSON via
%       BrainNetRender(scenePath, outputPath).
%
% STEPS:
%   1) In the GUI, load each variant (surface + node + edge), set the look you
%      want, then click "Save Scene" and save the three .json files below into
%      THIS folder:
%        - top-40 : nodes Schaefer200_Yeo7.node   + edges NeurWEM_top40.edge
%        - hub    : nodes Schaefer200_Yeo7_hub.node + NO edges (clear edges)
%        - NBS    : nodes Schaefer200_Yeo7.node    + edges NeurWEM_NBS.edge
%   2) Set BNV below to your BrainNet Viewer folder, then run this script.
% ============================================================

BNV = '/Users/bai/Downloads/BrainNet-Viewer-MATLAB-2.0.0-beta.1/src/matlab/BrainNetViewer';
addpath(genpath(BNV));
here = fileparts(mfilename('fullpath'));

scenes = {
  'NeurWEM_top40.scene.json', 'render_top40.png';
  'NeurWEM_hub.scene.json',   'render_hub.png';
  'NeurWEM_NBS.scene.json',   'render_NBS.png';
};

for i = 1:size(scenes,1)
  sp = fullfile(here, scenes{i,1});
  op = fullfile(here, scenes{i,2});
  if exist(sp,'file')~=2
    fprintf('SKIP (save this scene from the GUI first): %s\n', sp); continue;
  end
  BrainNetRender(sp, op);
  fprintf('wrote %s\n', op);
end
disp('done');
