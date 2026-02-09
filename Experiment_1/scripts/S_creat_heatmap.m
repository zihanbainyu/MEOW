%% example spatial correlation visualization
base_dir = '..'; 
out_dir = fullfile(base_dir, 'data', 'eye_movement_data');
stim_dir = fullfile(base_dir, 'stimulus', 'stim_final');
res_dir = fullfile(base_dir, 'results');
load(fullfile(out_dir, 'group_eye_movement_m.mat'));
load(fullfile(res_dir, 'gaze_reinstat_res.mat'));
results_comp = reinstat_res.compared;

example = results_comp(results_comp.subj_id==602 & results_comp.tr_1b==213 & results_comp.tr_2b==167, :);
sid = example.subj_id; tr_1b = example.tr_1b; tr_2b = example.tr_2b;
fix_1b = Mw(Mw.subj_id==sid & Mw.trial_id==tr_1b & strcmp(Mw.task,'1_back'), :);
fix_2b = Mw(Mw.subj_id==sid & Mw.trial_id==tr_2b & strcmp(Mw.task,'2_back'), :);
stim_name = fix_1b.stim_id{1}; 
img_path = fullfile(stim_dir, stim_name);
spatial_params = reinstat_res.spatial_params;
map_1b = create_fixation_map(fix_1b.x, fix_1b.y, fix_1b.dur, spatial_params);
map_2b = create_fixation_map(fix_2b.x, fix_2b.y, fix_2b.dur, spatial_params);

mismatch_tr_1b = 474;
fix_1b_mis = Mw(Mw.subj_id==sid & Mw.trial_id==mismatch_tr_1b & strcmp(Mw.task,'1_back'), :);
fix_2b_mis = Mw(Mw.subj_id==sid & Mw.trial_id==tr_2b & strcmp(Mw.task,'2_back'), :); 
stim_name_mis = fix_1b_mis.stim_id{1};
img_path_mis = fullfile(stim_dir, stim_name_mis);
map_1b_mis = create_fixation_map(fix_1b_mis.x, fix_1b_mis.y, fix_1b_mis.dur, spatial_params);
map_2b_mis = create_fixation_map(fix_2b_mis.x, fix_2b_mis.y, fix_2b_mis.dur, spatial_params);
valid = ~isnan(map_1b_mis(:)) & ~isnan(map_2b_mis(:));
r_mismatch = corr(map_1b_mis(valid), map_2b_mis(valid));

all_maps = cat(3, map_1b, map_2b, map_1b_mis, map_2b_mis);
clim_heatmaps = [0, max(all_maps(:))];
diff_match = map_1b - map_2b;
diff_mismatch = map_1b_mis - map_2b_mis;
max_abs_diff = max(abs([diff_match(:); diff_mismatch(:)]));
clim_diff = [-max_abs_diff, max_abs_diff];

figure('color','w','position',[100 100 1500 600]);
subplot(2,5,1);
imagesc(map_1b); axis image; colormap('hot'); colorbar;
caxis(clim_heatmaps);
title('1-Back Encoding', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,2);
imagesc(map_2b); axis image; colormap('hot'); colorbar;
caxis(clim_heatmaps);
title('2-Back Retrieval', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,3);
img = imread(img_path);
imshow(img); hold on;
x_1b_img = (fix_1b.x - spatial_params.roi_x(1)) / (spatial_params.roi_x(2) - spatial_params.roi_x(1)) * size(img, 2);
y_1b_img = (fix_1b.y - spatial_params.roi_y(1)) / (spatial_params.roi_y(2) - spatial_params.roi_y(1)) * size(img, 1);
x_2b_img = (fix_2b.x - spatial_params.roi_x(1)) / (spatial_params.roi_x(2) - spatial_params.roi_x(1)) * size(img, 2);
y_2b_img = (fix_2b.y - spatial_params.roi_y(1)) / (spatial_params.roi_y(2) - spatial_params.roi_y(1)) * size(img, 1);
scatter(x_1b_img, y_1b_img, fix_1b.dur*0.3, 'r', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(x_2b_img, y_2b_img, fix_2b.dur*0.3, 'b', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
title('Fixations Overlay', 'FontSize', 14);
legend({'1-Back Encoding', '2-Back Retrieval'}, 'Location', 'best');
hold off;
subplot(2,5,4);
imagesc(diff_match); axis image; colormap('jet'); colorbar;
caxis(clim_diff);
title('Difference', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,5); hold on;
plot(map_1b(:), map_2b(:), 'k.', 'MarkerSize', 10);
xlabel('1-Back', 'FontSize', 12);
ylabel('2-Back', 'FontSize', 12);
title(sprintf('Spatial Correlation\nr = %.3f', example.match_score), 'FontSize', 14);
xlim([0 0.04]); ylim([0 0.04]); yticks(0:0.01:0.04); lsline; axis square; grid on; hold off;
subplot(2,5,6);
imagesc(map_1b_mis); axis image; colormap('hot'); colorbar;
caxis(clim_heatmaps);
title('1-Back Encoding (Mismatch)', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,7);
imagesc(map_2b_mis); axis image; colormap('hot'); colorbar;
caxis(clim_heatmaps);
title('2-Back Retrieval (Same as Match)', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,8);
img_mis = imread(img_path_mis);
imshow(img_mis); hold on;
x_1b_mis_img = (fix_1b_mis.x - spatial_params.roi_x(1)) / (spatial_params.roi_x(2) - spatial_params.roi_x(1)) * size(img_mis, 2);
y_1b_mis_img = (fix_1b_mis.y - spatial_params.roi_y(1)) / (spatial_params.roi_y(2) - spatial_params.roi_y(1)) * size(img_mis, 1);
x_2b_mis_img = (fix_2b_mis.x - spatial_params.roi_x(1)) / (spatial_params.roi_x(2) - spatial_params.roi_x(1)) * size(img_mis, 2);
y_2b_mis_img = (fix_2b_mis.y - spatial_params.roi_y(1)) / (spatial_params.roi_y(2) - spatial_params.roi_y(1)) * size(img_mis, 1);
scatter(x_1b_mis_img, y_1b_mis_img, fix_1b_mis.dur*0.3, 'r', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(x_2b_mis_img, y_2b_mis_img, fix_2b_mis.dur*0.3, 'b', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
title('Fixations Overlay', 'FontSize', 14);
legend({'1-Back (Mismatch)', '2-Back Retrieval'}, 'Location', 'best'); hold off;
subplot(2,5,9);
imagesc(diff_mismatch); axis image; colormap('jet'); colorbar;
caxis(clim_diff);
title('Difference', 'FontSize', 14);
xlabel('X position'); ylabel('Y position');
subplot(2,5,10); hold on;
plot(map_1b_mis(:), map_2b_mis(:), 'k.', 'MarkerSize', 10);
xlabel('1-Back Fixation Map (Mismatch)', 'FontSize', 12);
ylabel('2-Back Fixation Map', 'FontSize', 12);
title(sprintf('Spatial Correlation\nr = %.3f', r_mismatch), 'FontSize', 14);
xlim([0 0.04]); ylim([0 0.04]); yticks(0:0.01:0.04); lsline; axis square; grid on; hold off;
print(gcf, 'example_gaze_reinstat.pdf', '-dpdf', '-r300');


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