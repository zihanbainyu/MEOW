%% Example spatial correlation visualization
% Specific example trial
base_dir = '..'; 
out_dir = fullfile(base_dir, 'data', 'eye_movement_data');
stim_dir = fullfile(base_dir, 'stimulus', 'stim_final');
load(fullfile(out_dir, 'group_eye_movement_m.mat'));
results_comp = reinstat_res.compared;

% Find the specific trial: subject 501, tr_1b=598, tr_2b=613
example = results_comp(results_comp.subj_id==501 & results_comp.tr_1b==598 & results_comp.tr_2b==613, :);

% Get the relevant info
sid = example.subj_id;
tr_1b = example.tr_1b;
tr_2b = example.tr_2b;

fprintf('Example trial: subject %d, 1-back trial %d, 2-back trial %d, correlation r=%.3f\n', ...
    sid, tr_1b, tr_2b, example.match_score);

% Get fixations for this trial pair
fix_1b = Mw(Mw.subj_id==sid & Mw.trial_id==tr_1b, :);
fix_2b = Mw(Mw.subj_id==sid & Mw.trial_id==tr_2b, :);

% Get stimulus filename
stim_name = fix_1b.stim_id{1}; 
img_path = fullfile(stim_dir, stim_name);

% Create fixation maps
spatial_params = reinstat_res.spatial_params;
map_1b = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);

% Create visualization
figure('color','w','position',[100 100 1600 400]);

subplot(1,5,1);
imagesc(map_1b); axis image; colormap('hot'); colorbar;
title('Encoding (1-back)', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
set(gca, 'YDir', 'normal');

subplot(1,5,2);
imagesc(map_2b); axis image; colormap('hot'); colorbar;
title('Retrieval (2-back)', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
set(gca, 'YDir', 'normal');

subplot(1,5,3);
if isfile(img_path)
    img = imread(img_path);
    imshow(img); hold on;
    % Plot encoding fixations in red
    scatter(fix_1b.x, fix_1b.y, fix_1b.dur/3, 'r', 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w', 'LineWidth', 1);
    % Plot retrieval fixations in blue
    scatter(fix_2b.x, fix_2b.y, fix_2b.dur/3, 'b', 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w', 'LineWidth', 1);
    title('Fixations on Image', 'FontSize', 14);
    legend({'Encoding', 'Retrieval'}, 'Location', 'best');
    hold off;
else
    text(0.5, 0.5, 'Image not found', 'HorizontalAlignment', 'center');
    title('Stimulus Image', 'FontSize', 14);
    fprintf('Warning: Image not found at %s\n', img_path);
end

subplot(1,5,4);
imagesc(map_1b - map_2b); axis image; colormap('jet'); colorbar;
title('Difference', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
set(gca, 'YDir', 'normal');

subplot(1,5,5);
hold on;
plot(map_1b(:), map_2b(:), 'k.', 'MarkerSize', 10);
xlabel('Encoding fixation density', 'FontSize', 12);
ylabel('Retrieval fixation density', 'FontSize', 12);
title(sprintf('Spatial Correlation\nr = %.3f', example.match_score), 'FontSize', 14);
lsline;
axis square; grid on;
hold off;

sgtitle(sprintf('Example Gaze Reinstatement (Subject %d, Stimulus %s)', sid, stim_name), ...
    'FontSize', 16, 'FontWeight', 'bold');

print(gcf, 'Example_Gaze_Reinstatement.pdf', '-dpdf', '-r300');

function map = create_fixation_map(x, y, dur, params)
    in_roi = x >= params.roi_x(1) & x <= params.roi_x(2) & y >= params.roi_y(1) & y <= params.roi_y(2);
    x = x(in_roi); y = y(in_roi); dur = dur(in_roi);
    
    if isempty(x), map = zeros(params.grid_y, params.grid_x); return; end
    
    x_bins = linspace(params.roi_x(1), params.roi_x(2), params.grid_x+1);
    y_bins = linspace(params.roi_y(1), params.roi_y(2), params.grid_y+1);
    map = zeros(params.grid_y, params.grid_x);
    
    for i = 1:length(x)
        x_idx = find(x(i) >= x_bins(1:end-1) & x(i) < x_bins(2:end), 1);
        y_idx = find(y(i) >= y_bins(1:end-1) & y(i) < y_bins(2:end), 1);
        if ~isempty(x_idx) && ~isempty(y_idx)
            map(y_idx, x_idx) = map(y_idx, x_idx) + dur(i);
        end
    end
    
    if params.sigma > 0, map = imgaussfilt(map, params.sigma); end
    if sum(map(:)) > 0, map = map / sum(map(:)); end
end