clc; 
close all;

format long g

img1 = imread('../toy-car-images-bw/toy_formatted2.png');
img2 = imread('../toy-car-images-bw/toy_formatted3.png');

% Lucas Kanade method
oflow = opticalFlowLK(img1, img2, 3);

% Normal flow
n_flow = computeNormalFlow(img1, img2);



% Show optical flow overlayed on image
figure;
space = 1;
imshow(img1);
hold on;
% quiver(x,y, n_flow(1:space:end, 1:space:end, 1), n_flow(1:space:end, 1:space:end,2), 5, 'color', 'g');
quiver(oflow(1:space:end, 1:space:end, 1), oflow(1:space:end, 1:space:end,2), 5, 'color', 'g');
title('Overlaid Normal Flow');