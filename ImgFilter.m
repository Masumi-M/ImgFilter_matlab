%% Info
% @created 2019/5/11 [Masumi Morishige]

%% Reset
clear;
close all;
clc;

%% Setting
appleIMG = imread('apple.jpeg');
IMG_size = size(appleIMG);

%% Mono
% weight = [0.3 0.3 0.3];   % Simple Average
weight = [0.299 0.587 0.114];   % WithWeight Average
appleIMG_mono_pin = (weight(1)*appleIMG(:,:,1) + weight(2)*appleIMG(:,:,2) + weight(3)*appleIMG(:,:,3));
appleIMG_mono = repmat(appleIMG_mono_pin,[1 1 3]);

figure(1);
imshow(appleIMG);
figure(2);
imshow(appleIMG_mono);

%% Affine Translation
% Zoom
zoom_array = [2 0 0;0 2 0; 0 0 1];
zoom_tform = affine2d(zoom_array);
appleIMG_zoom = imwarp(appleIMG, zoom_tform);

% Translation
translate_distance = [3000 3000];
translate_array = [1 0 0;0 1 0; translate_distance(1) translate_distance(2) 1];
% translate_array = [1 0 translate_distance(1);0 1 translate_distance(2); 0 0 1];
translate_tform = affine2d(translate_array);
% appleIMG_translate = imtransform(appleIMG, translate_tform, 'XData', [1, IMG_size(1)], 'YData', [1, IMG_size(2)], 'XYScale', 1);
appleIMG_translate = imwarp(appleIMG, translate_tform);

% Rotation
theta = pi/6;
rotate_array = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
rotate_tform = affine2d(rotate_array);
appleIMG_rotate = imwarp(appleIMG, rotate_tform);

% Shear
theta2 = pi/6;
shear_array = [1 0 0;tan(theta2) 1 0; 0 0 1];
shear_tform = affine2d(shear_array);
appleIMG_shear = imwarp(appleIMG, shear_tform);

figure(3);
imshow(appleIMG_zoom);

figure(4);
imshow(appleIMG_translate);

figure(5);
imshow(appleIMG_rotate);

figure(6);
imshow(appleIMG_shear);

%% End of the Script
