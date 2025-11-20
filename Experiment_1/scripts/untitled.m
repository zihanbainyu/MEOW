

clear; 
clc; 
close all;
base_dir = '..'; 
in_dir = fullfile(base_dir, 'stimulus/stim_final');
out_dir = fullfile(base_dir, 'stimulus/roi_masks');% Output for binary masks
check_dir = fullfile(out_dir, 'roi_visual_check');  % Output for QA images

if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if ~exist(check_dir, 'dir'), mkdir(check_dir); end

allImgList = dir(fullfile(in_dir, 'mst_*_l*.png')); 

% parameters
blurKernel = 15;      % smooths difference
diffThreshold = 50;   % intensity diff needed to count as "change" (0-255)
dilateSize = 15;      % expansion of the ROI (in pixels) to add buffer
minBlobSize = 100;    % remove tiny noise specks (pixels)

fprintf('Found %d total images. Starting pairwise comparison...\n', length(allImgList));

% interate over A and find their match B
for i = 1:length(allImgList)
    
    fnameA = allImgList(i).name;
    
    if contains(fnameA, '_A_')
        fnameB = strrep(fnameA, '_A_', '_B_');
       
        pathA = fullfile(in_dir, fnameA);
        pathB = fullfile(in_dir, fnameB);
        
        if ~isfile(pathB)
            warning('match not found for: %s (expected %s)', fnameA, fnameB);
            continue;
        end
        
        imgA = imread(pathA);
        imgB = imread(pathB);
        
        % grey scale 
        if size(imgA, 3) == 3
            grayA = rgb2gray(imgA);
            grayB = rgb2gray(imgB);
        else
            grayA = imgA;
            grayB = imgB;
        end
        
        % image alignment
        try
            [optimizer, metric] = imregconfig('monomodal');
            % Register B to A
            tform = imregtform(grayB, grayA, 'translation', optimizer, metric);
            grayB_reg = imwarp(grayB, tform, 'OutputView', imref2d(size(grayA)));
        catch
            grayB_reg = grayB; % Fallback
        end
        
        % calculate diff
        % diffImg = imabsdiff(grayA, grayB_reg);
        
        % smoothing
        diffBlurred = imgaussfilt(diffImg, blurKernel/3);
        
        % thresholding to create binary masks
        binaryMask = diffBlurred > diffThreshold;
        
        % morphological cleanup
        binaryMask = bwareaopen(binaryMask, minBlobSize); % remove noise
        se = strel('disk', dilateSize);
        finalMask = imdilate(binaryMask, se);
        finalMask = imfill(finalMask, 'holes'); % fill small internal holes
        
        % save
        
        [~, name, ~] = fileparts(fnameA);
        out_name = strrep(name, '_A', '');

        imwrite(finalMask, fullfile(out_dir, [out_name '_droi.png']));
        
        h = figure('Visible', 'off'); 
        imshow(imgA); hold on;
        visboundaries(finalMask, 'Color', 'r', 'LineWidth', 3);
        title(['Diagnostic ROI: ' out_name], 'Interpreter', 'none');
        
        saveas(h, fullfile(check_dir, [out_name '_check.jpg']));
        close(h);
        
    end
end

fprintf('\nDone! Check the "%s" folder to validate your ROIs.\n', checkDir);