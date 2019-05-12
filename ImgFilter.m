%% Info
% @created 2019/5/11 [Masumi Morishige]

%% Reset
clear;
close all;
clc;

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

%% RGB Image
appleIMG_red = appleIMG(:,:,1);
appleIMG_green = appleIMG(:,:,2);
appleIMG_blue = appleIMG(:,:,3);

figure('Name', 'Red', 'NumberTitle', 'off');
imshow(appleIMG_red);
figure('Name', 'Green', 'NumberTitle', 'off');
imshow(appleIMG_green);
figure('Name', 'Blue', 'NumberTitle', 'off');
imshow(appleIMG_blue);

%% GrayScale
weight = [0.299 0.587 0.114];
appleIMG_gray_pin = (weight(1)*appleIMG(:,:,1) + weight(2)*appleIMG(:,:,2) + weight(3)*appleIMG(:,:,3));
appleIMG_gray = repmat(appleIMG_gray_pin,[1 1 3]);

figure('Name', 'GrayScale', 'NumberTitle', 'off');
imshow(appleIMG_gray);

%% Affine Translation
% Zoom
zoom_array = [2 0 0;0 2 0; 0 0 1];
zoom_tform = affine2d(zoom_array);
appleIMG_zoom = imwarp(appleIMG, zoom_tform);
figure('Name', 'Zoom', 'NumberTitle', 'off');
imshow(appleIMG_zoom);

% Translation
translate_distance = [3000 3000];
translate_array = [1 0 0;0 1 0; translate_distance(1) translate_distance(2) 1];
translate_tform = affine2d(translate_array);
% appleIMG_translate = imtransform(appleIMG, translate_tform, 'XData', [1, IMG_size(1)], 'YData', [1, IMG_size(2)], 'XYScale', 1);
appleIMG_translate = imwarp(appleIMG, translate_tform, 'XData', [1, IMG_size(1)], 'YData', [1, IMG_size(2)], 'XYScale', 1);
figure('Name', 'Translate', 'NumberTitle', 'off');
imshow(appleIMG_translate);

% Rotation
theta = pi/6;
rotate_array = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
rotate_tform = affine2d(rotate_array);
appleIMG_rotate = imwarp(appleIMG, rotate_tform);
figure('Name', 'Rotate', 'NumberTitle', 'off');
imshow(appleIMG_rotate);

% Shear
theta2 = pi/6;
shear_array = [1 0 0;tan(theta2) 1 0; 0 0 1];
shear_tform = affine2d(shear_array);
appleIMG_shear = imwarp(appleIMG, shear_tform);
figure('Name', 'Shear', 'NumberTitle', 'off');
imshow(appleIMG_shear);

%% Noise (salt and pepper)
appleIMG_noise = imnoise(appleIMG, 'salt & pepper');
figure('Name', 'Noise', 'NumberTitle', 'off');
imshow(appleIMG_noise);

%% End of the Script
