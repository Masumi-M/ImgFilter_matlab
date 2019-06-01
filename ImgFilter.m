%% Info
% @created 2019/5/11 [Masumi Morishige]
% @updated 2019/6/1 Add: Differential Image / Sharpening / Fourier [Masumi Morishige]

%% Reset
clear;
close all;
clc;
save_dir_name = './OutputImage/';

%% Original
appleIMG = imread('apple.jpeg');
IMG_size = size(appleIMG);
[filepath,name,ext] = fileparts('apple.jpeg');

fprintf('>> Info\n');
fprintf('Size: [%d, %d]\n', IMG_size(1,1), IMG_size(1,2));
fprintf('Number of Pixels: %d px\n', IMG_size(1,1)*IMG_size(1,2));
fprintf('Resolution: %d ppi\n', 72);
fprintf('Number of Gradations: %d (%d bit) \n', 2^8, 8);
fprintf('File Extension: %s \n', ext);

figure('Name', 'Original', 'NumberTitle', 'off');
imshow(appleIMG);
saveas(gcf, strcat(save_dir_name, 'original.jpg'));

%% RGB Image
appleIMG_red = appleIMG(:,:,1);
appleIMG_green = appleIMG(:,:,2);
appleIMG_blue = appleIMG(:,:,3);

figure('Name', 'Red', 'NumberTitle', 'off');
imshow(appleIMG_red);
saveas(gcf, strcat(save_dir_name, 'red.jpg'));

figure('Name', 'Green', 'NumberTitle', 'off');
imshow(appleIMG_green);
saveas(gcf, strcat(save_dir_name, 'green.jpg'));

figure('Name', 'Blue', 'NumberTitle', 'off');
imshow(appleIMG_blue);
saveas(gcf, strcat(save_dir_name, 'blue.jpg'));

%% GrayScale
% 1st type
appleIMG_gray_pin_1 = (appleIMG(:,:,1) + appleIMG(:,:,2) + appleIMG(:,:,3))./3;
appleIMG_gray_1 = repmat(appleIMG_gray_pin_1,[1 1 3]);
figure('Name', 'GrayScale1', 'NumberTitle', 'off');
imshow(appleIMG_gray_1);
saveas(gcf, strcat(save_dir_name, 'grayscale1.jpg'));

% 2nd type
weight = [0.299 0.587 0.114];
appleIMG_gray_pin_2 = (weight(1)*appleIMG(:,:,1) + weight(2)*appleIMG(:,:,2) + weight(3)*appleIMG(:,:,3));
appleIMG_gray_2 = repmat(appleIMG_gray_pin_2,[1 1 3]);
figure('Name', 'GrayScale2', 'NumberTitle', 'off');
imshow(appleIMG_gray_2);
saveas(gcf, strcat(save_dir_name, 'grayscale2.jpg'));

% 3rd type
appleIMG_gray_pin_3 = max(appleIMG,[],[(IMG_size(1,1)*IMG_size(1,2)) 3]);
appleIMG_gray_3 = repmat(appleIMG_gray_pin_3,[1 1 3]);
figure('Name', 'GrayScale3', 'NumberTitle', 'off');
imshow(appleIMG_gray_3);
saveas(gcf, strcat(save_dir_name, 'grayscale3.jpg'));

% 4th type
appleIMG_gray_4 = rgb2gray(appleIMG);
figure('Name', 'GrayScale4', 'NumberTitle', 'off');
imshow(appleIMG_gray_4);
saveas(gcf, strcat(save_dir_name, 'grayscale4.jpg'));

%% Affine Translation
% Zoom
zoom_array = [2 0 0;0 2 0; 0 0 1];
zoom_tform = affine2d(zoom_array);
appleIMG_zoom = imwarp(appleIMG, zoom_tform);
figure('Name', 'Zoom', 'NumberTitle', 'off');
imshow(appleIMG_zoom);
saveas(gcf, strcat(save_dir_name, 'zoom.jpg'));

% Translation
% translate_distance = [3000 3000];
% translate_array = [1 0 0;0 1 0; translate_distance(1) translate_distance(2) 1];
% translate_tform = affine2d(translate_array);
% % appleIMG_translate = imtransform(appleIMG, translate_tform, 'XData', [1, IMG_size(1)], 'YData', [1, IMG_size(2)], 'XYScale', 1);
% appleIMG_translate = imwarp(appleIMG, translate_tform, 'XData', [1, IMG_size(1)], 'YData', [1, IMG_size(2)], 'XYScale', 1);
% figure('Name', 'Translate', 'NumberTitle', 'off');
% imshow(appleIMG_translate);
% saveas(gcf, strcat(save_dir_name, 'translate.jpg'));

% Rotation
theta = pi/6;
rotate_array = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
rotate_tform = affine2d(rotate_array);
appleIMG_rotate = imwarp(appleIMG, rotate_tform);
figure('Name', 'Rotate', 'NumberTitle', 'off');
imshow(appleIMG_rotate);
saveas(gcf, strcat(save_dir_name, 'rotation.jpg'));

% Shear
theta2 = pi/6;
shear_array = [1 0 0;tan(theta2) 1 0; 0 0 1];
shear_tform = affine2d(shear_array);
appleIMG_shear = imwarp(appleIMG, shear_tform);
figure('Name', 'Shear', 'NumberTitle', 'off');
imshow(appleIMG_shear);
saveas(gcf, strcat(save_dir_name, 'shear.jpg'));

%% Noise (salt and pepper)
% With Noise(color)
appleIMG_noise = imnoise(appleIMG, 'salt & pepper');
figure('Name', 'Noise', 'NumberTitle', 'off');
imshow(appleIMG_noise);
saveas(gcf, strcat(save_dir_name, 'noise.jpg'));

% With Noise(grayscale)
appleIMG_noise_gray = rgb2gray(appleIMG_noise);
figure('Name', 'Noise(gray)', 'NumberTitle', 'off');
imshow(appleIMG_noise_gray);
saveas(gcf, strcat(save_dir_name, 'noise_gray.jpg'));

%% Moving Average Filter
movingArray = [1 1 1;1 1 1;1 1 1]/9;
appleIMG_noise_moving = filter2(movingArray, appleIMG_noise_gray);
appleIMG_noise_moving_mean = uint8(appleIMG_noise_moving);
figure('Name', 'Moving Filter', 'NumberTitle', 'off');
imshow(appleIMG_noise_moving_mean);
saveas(gcf, strcat(save_dir_name, 'moving_filter.jpg'));

%% Weighted Average Filter
weightedArray = [1 2 1;2 4 2;1 2 1]/16;
appleIMG_noise_weighted = filter2(weightedArray, appleIMG_noise_gray);
appleIMG_noise_weighted_mean = uint8(appleIMG_noise_weighted);
figure('Name', 'Weighted Filter', 'NumberTitle', 'off');
imshow(appleIMG_noise_weighted_mean);
saveas(gcf, strcat(save_dir_name, 'weighted_filter.jpg'));

%% Median Filter
appleIMG_noise_median = medfilt2(appleIMG_noise_gray);
figure('Name', 'Median Filter', 'NumberTitle', 'off');
imshow(appleIMG_noise_median);
saveas(gcf, strcat(save_dir_name, 'median_filter.jpg'));

%% Differential Image
% Sobel Operator


% Prewitt Operator


% Laplacian Operator


%% Sharpening Filter



%% Fourier Transform


%% End of the Script
