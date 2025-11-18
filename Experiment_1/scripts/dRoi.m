clear; 
clc; 
close all;
base_dir = '..'; 
in_dir = fullfile(base_dir, 'stimulus/stim_final');
out_dir = fullfile(base_dir, 'stimulus/roi_masks');
check_dir = fullfile(out_dir, 'roi_visual_check');  
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if ~exist(check_dir, 'dir'), mkdir(check_dir); end

% parameters
blurKernel = 15;      
diffThreshold = 100;   
dilateSize = 15;      
minBlobSize = 100;    

allImgList = dir(fullfile(in_dir, 'mst_*_l*.png')); 
A_files = allImgList(contains({allImgList.name}, '_A_'));
num_iterations = length(A_files);
fprintf('Processing %d A-B pairs, mapping dROI to B image...\n', num_iterations);

for i = 1:num_iterations
    fnameA = A_files(i).name;
    fnameB = strrep(fnameA, '_A_', '_B_');
    pathA = fullfile(in_dir, fnameA);
    pathB = fullfile(in_dir, fnameB);
    if ~isfile(pathB)
        warning('match not found for: %s (expected %s)', fnameA, fnameB);
        continue;
    end
    imgA = imread(pathA);
    imgB = imread(pathB);
    if size(imgA, 3) == 3
        grayA = rgb2gray(imgA);
        grayB = rgb2gray(imgB);
    else
        grayA = imgA;
        grayB = imgB;
    end
    grayB_reg = grayB; 

    registration_mask = true(size(grayA)); 
    try
        [optimizer, metric] = imregconfig('monomodal');
        tform = imregtform(grayB, grayA, 'translation', optimizer, metric);
        R_out = imref2d(size(grayA)); 
        [grayB_reg, ~] = imwarp(grayB, tform, 'OutputView', R_out);
        binary_B = true(size(grayB)); 
        [warped_binary_B, ~] = imwarp(binary_B, tform, 'OutputView', R_out);
        registration_mask = warped_binary_B;
    catch
        warning('Registration failed for %s. Proceeding without alignment.', fnameA);
    end

    diffImg_raw = imabsdiff(grayA, grayB_reg);
    diffImg_raw(~registration_mask) = 0; 
    diffDouble = double(diffImg_raw) ./ 255;
    diffSquared = diffDouble .^ 2; 
    diffImg = uint8(diffSquared .* 255); 

    diffBlurred = imgaussfilt(diffImg, blurKernel/3);

    binaryMask = diffBlurred > diffThreshold; 
    binaryMask = bwareaopen(binaryMask, minBlobSize); 

    se = strel('disk', dilateSize);

    finalMask = imdilate(binaryMask, se);
    finalMask = imfill(finalMask, 'holes');
    
    [~, name, ~] = fileparts(fnameA);
    out_name = strrep(name, '_A', '');
    imwrite(finalMask, fullfile(out_dir, [out_name '_droi.png']));
    h = figure('Visible', 'off'); 
    imshow(imgB); hold on; 
    visboundaries(finalMask, 'Color', 'r', 'LineWidth', 3);
    title(['Diagnostic ROI on B: ' out_name], 'Interpreter', 'none');
    saveas(h, fullfile(check_dir, [out_name '_check.jpg']));
    close(h);
end
fprintf('\nDone! Check the "%s" folder to validate your ROIs, which are mapped to the B images.\n', check_dir);