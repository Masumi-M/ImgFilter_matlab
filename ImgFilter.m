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
sobelArrayX = [-1 0 1; -2 0 2; -1 0 1];
appleIMG_sobelX = uint8(filter2(sobelArrayX, appleIMG_gray_4));
figure('Name', 'Sobel Filter [x]', 'NumberTitle', 'off');
imshow(appleIMG_sobelX);
saveas(gcf, strcat(save_dir_name, 'sobel_filter_x.jpg'));

sobelArrayY = [-1 -2 -1; 0 0 0; 1 2 1];
appleIMG_sobelY = uint8(filter2(sobelArrayY, appleIMG_gray_4));
figure('Name', 'Sobel Filter [y]', 'NumberTitle', 'off');
imshow(appleIMG_sobelY);
saveas(gcf, strcat(save_dir_name, 'sobel_filter_y.jpg'));

% Prewitt Operator
prewittArrayX = [-1 0 1; -1 0 1; -1 0 1];
appleIMG_prewittX = uint8(filter2(prewittArrayX, appleIMG_gray_4));
figure('Name', 'Prewitt Filter [x]', 'NumberTitle', 'off');
imshow(appleIMG_prewittX);
saveas(gcf, strcat(save_dir_name, 'prewitt_filter_x.jpg'));

prewittArrayY = [-1 -1 -1; 0 0 0; 1 1 1];
appleIMG_prewittY = uint8(filter2(prewittArrayY, appleIMG_gray_4));
figure('Name', 'Prewitt Filter [y]', 'NumberTitle', 'off');
imshow(appleIMG_prewittY);
saveas(gcf, strcat(save_dir_name, 'prewitt_filter_y.jpg'));


% Laplacian Operator
laplacianArray = [1 1 1; 1 -8 1; 1 1 1];
appleIMG_laplacian = uint8(filter2(laplacianArray, appleIMG_gray_4));
figure('Name', 'Laplacian Filter', 'NumberTitle', 'off');
imshow(appleIMG_laplacian);
saveas(gcf, strcat(save_dir_name, 'laplacian_filter.jpg'));

%% Sharpening Filter
movingArray = [1 1 1;1 1 1;1 1 1]/9;
appleIMG_moving = filter2(movingArray, appleIMG_gray_4);
appleIMG_moving_mean = uint8(appleIMG_moving);
figure('Name', 'Sharpening Filter [Before]', 'NumberTitle', 'off');
imshow(appleIMG_moving_mean);
saveas(gcf, strcat(save_dir_name, 'sharpening_filter_before.jpg'));

sharpeningArray = [-1 -1 -1; -1 9 -1; -1 -1 -1];
appleIMG_sharpening = uint8(filter2(sharpeningArray, appleIMG_moving_mean));
figure('Name', 'Sharpening Filter [After]', 'NumberTitle', 'off');
imshow(appleIMG_sharpening);
saveas(gcf, strcat(save_dir_name, 'sharpening_filter_after.jpg'));

%% Fourier Transform
% Fast Fourier Transform[FFT]
zoom_out_array = [0.8 0 0;0 0.8 0; 0 0 1];
zoom_out_tform = affine2d(zoom_out_array);
appleIMG_zoom_out = imwarp(appleIMG_gray_4, zoom_out_tform);
appleIMG_fft_before = appleIMG_zoom_out(1:256,101:356,1);
figure('Name', 'FFT [before]', 'NumberTitle', 'off');
imshow(appleIMG_fft_before);
saveas(gcf, strcat(save_dir_name, 'fft_before.jpg'));

appleIMG_fft_after = abs(fft2(appleIMG_fft_before));
appleIMG_fft_after = uint8(appleIMG_fft_after/max(appleIMG_fft_after(:)).^0.35);
figure('Name', 'FFT [after]', 'NumberTitle', 'off');
imshow(appleIMG_fft_after);
saveas(gcf, strcat(save_dir_name, 'fft_after.jpg'));

% Filtering in frequency-area
appleIMG_fft_shift = fftshift(fft2(appleIMG_fft_before));
appleIMG_fft_shift = repmat(appleIMG_fft_shift, [1 1 3]);
% appleIMG_fft_shift_rep = ;
filter_parameter = [5 10 50];
[image_height, image_width] = size(appleIMG_fft_shift(:,:,1));
for i_para = 1:3
    for i_width = 1:image_width
        for i_height = 1:image_height
            r = sqrt((i_width - (image_height/2 + 1))^2 + (i_height - (image_width/2 - 1))^2);
            if r > filter_parameter(i_para)
                appleIMG_fft_shift(i_height, i_width, i_para) = 0;
            end
        end
    end
end

appleIMG_fft_filtering_R5 = ifft2(ifftshift(appleIMG_fft_shift(:,:,1)));
appleIMG_fft_filtering_R5 = uint8(abs(appleIMG_fft_filtering_R5));
figure('Name', 'FFT filtering [R=5]', 'NumberTitle', 'off');
imshow(appleIMG_fft_filtering_R5);
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R5.jpg'));

appleIMG_fft_filtering_R10 = ifft2(ifftshift(appleIMG_fft_shift(:,:,2)));
appleIMG_fft_filtering_R10 = uint8(abs(appleIMG_fft_filtering_R10));
figure('Name', 'FFT filtering [R=10]', 'NumberTitle', 'off');
imshow(appleIMG_fft_filtering_R10);
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R10.jpg'));

appleIMG_fft_filtering_R50 = ifft2(ifftshift(appleIMG_fft_shift(:,:,3)));
appleIMG_fft_filtering_R50 = uint8(abs(appleIMG_fft_filtering_R50));
figure('Name', 'FFT filtering [R=50]', 'NumberTitle', 'off');
imshow(appleIMG_fft_filtering_R50);
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R50.jpg'));

% Inverse Fourier Transform
% R = 5
appleIMG_fft_filtering_R5_fft = fft2(appleIMG_fft_filtering_R5);
appleIMG_fft_filtering_R5_fft_shift = fftshift(appleIMG_fft_filtering_R5_fft);
appleIMG_fft_filtering_R5_fft_abs_shift = fftshift(abs(appleIMG_fft_filtering_R5_fft));
figure('Name', 'FFT filtering (fftshift) [R=5]', 'NumberTitle', 'off');
imshow(uint8(appleIMG_fft_filtering_R5_fft_abs_shift/max(appleIMG_fft_filtering_R5_fft_abs_shift(:)).^0.35));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R5_fftshift.jpg'));

appleIMG_fft_filtering_R5_ifft = ifft2(fftshift(appleIMG_fft_filtering_R5_fft_shift));
figure('Name', 'FFT filtering (ifft) [R=5]', 'NumberTitle', 'off');
imshow(uint8(abs(appleIMG_fft_filtering_R5_ifft)));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R5_ifft.jpg'));

% R = 10
appleIMG_fft_filtering_R10_fft = fft2(appleIMG_fft_filtering_R10);
appleIMG_fft_filtering_R10_fft_shift = fftshift(appleIMG_fft_filtering_R10_fft);
appleIMG_fft_filtering_R10_fft_abs_shift = fftshift(abs(appleIMG_fft_filtering_R10_fft));
figure('Name', 'FFT filtering (fftshift) [R=10]', 'NumberTitle', 'off');
imshow(uint8(appleIMG_fft_filtering_R10_fft_abs_shift/max(appleIMG_fft_filtering_R10_fft_abs_shift(:)).^0.35));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R10_fftshift.jpg'));

appleIMG_fft_filtering_R10_ifft = ifft2(fftshift(appleIMG_fft_filtering_R10_fft_shift));
figure('Name', 'FFT filtering (ifft) [R=10]', 'NumberTitle', 'off');
imshow(uint8(abs(appleIMG_fft_filtering_R10_ifft)));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R10_ifft.jpg'));

% R = 50
appleIMG_fft_filtering_R50_fft = fft2(appleIMG_fft_filtering_R50);
appleIMG_fft_filtering_R50_fft_shift = fftshift(appleIMG_fft_filtering_R50_fft);
appleIMG_fft_filtering_R50_fft_abs_shift = fftshift(abs(appleIMG_fft_filtering_R50_fft));
figure('Name', 'FFT filtering (fftshift) [R=50]', 'NumberTitle', 'off');
imshow(uint8(appleIMG_fft_filtering_R50_fft_abs_shift/max(appleIMG_fft_filtering_R50_fft_abs_shift(:)).^0.35));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R50_fftshift.jpg'));

appleIMG_fft_filtering_R50_ifft = ifft2(fftshift(appleIMG_fft_filtering_R50_fft_shift));
figure('Name', 'FFT filtering (ifft) [R=50]', 'NumberTitle', 'off');
imshow(uint8(abs(appleIMG_fft_filtering_R50_ifft)));
saveas(gcf, strcat(save_dir_name, 'fft_filtering_R50_ifft.jpg'));

%% End of the Script
fprintf('--- End of the Script ---');