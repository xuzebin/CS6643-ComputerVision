% Find the corresponding matched point from 2 rectified stereo image pair 
% using normalized cross correlation (NCC), and compute the disparity of
% the two points.
clc; 
close all;
format long;

% Read image stereo pair
% Here we deliberatly swap the left and right images to find match points
% more easily.
L = imread('../res/cubeR.JPG');
R = imread('../res/cubeL.JPG');


% Convert to grayscale
L = rgb2gray(L);
R = rgb2gray(R);

% L and R must have the same size.
if size(L) ~= size(R)
    error('Stereo image pair must have the same length');
end

[row, col] = size(L);

% Define a ratio of image's size divided by window's size to be used for calculating
% proper window size for correlation.
imageWinRatio = 17;

% Window size
winRow = int32(row / imageWinRatio);
winCol = int32(col / imageWinRatio);

% Make sure the window's size is odd
if mod(winRow, 2) == 0
    winRow = winRow + 1;
end
if mod(winCol, 2) == 0
    winCol = winCol + 1;
end

% Half of window size
halfRow = (winRow - 1) / 2;
halfCol = (winCol - 1) / 2;
halfRow = double(halfRow);
halfCol = double(halfCol);

% Estimate a disparity maximum to reduce unnecessary computations.
% This depends greatly on the scene depth of the image.
maxDisparity = winCol * 3;

% Set the number of points to sample, can be more depending on how many
% points you want to calculate.
keyPointsNum = 1;

% Select and draw a point on the left image
figure;
subplot(121); imshow(L);
axis on;
hold on;

keyPoints = zeros(keyPointsNum, 2);
for i=1:keyPointsNum
[samplex, sampley] = ginput(1);
samplex = int16(samplex);
sampley = int16(sampley);
keyPoints(i,:) = [samplex, sampley];

plot(samplex, sampley, 'g.', 'markerSize', 10);
rectangle('Position', [samplex - halfCol, sampley - halfRow, winCol, winRow],...
	'EdgeColor','r', 'LineWidth', 1)
plot([1, col], [sampley, sampley], 'g');

end

% Draw the match point on the right image
subplot(122); imshow(R);
hold on;
axis on;

B = 20.0; % Baseline (mm)
f = 4.15; % Focal length (mm)
pixelSize = 0.00123; % Pixel size: (mm/pixel)
points3D = zeros(keyPointsNum, 3);

% Find matched point using NCC
for i=1:keyPointsNum
    samplex = keyPoints(i, 1);
    sampley = keyPoints(i, 2);
[bestMatchCol, correlation] = ...
    findMatchPointFromStereoPair([samplex, sampley], L, R, [winRow, winCol], maxDisparity);

    plot(bestMatchCol, sampley, 'g.', 'markerSize', 10);
    rectangle('Position', [bestMatchCol - halfCol, sampley - halfRow, winCol, winRow],...
        'EdgeColor','r', 'LineWidth', 1)
    plot([1, col], [sampley, sampley],'g');
    
    % Estimate depth: z=fB/d
    z = f * B / (double(bestMatchCol - samplex) * pixelSize);
    
    % Draw the text beside the key point
%     text(double(bestMatchCol), double(sampley), sprintf('%0.5f', z * 0.1), 'HorizontalAlignment', 'center');
    
    % Calculate x, y coordinates u=fx/z, v=fy/z => x=uz/f, y=vz/f.
%     a = double(f / pixelSize);
%     b = a;
    a = 3368.84;
    b = 3353.62;
    u0 = 1518.8;
    v0 = 1512;
    z = double(z);
    
    xCoord = double(samplex - u0) * z / a;
    yCoord = double(sampley - v0) * z / b;
    points3D(i, 1:3) = [xCoord, yCoord, z];
end

% points3D = dlmread('points3D.txt');
dlmwrite('points3D.txt', points3D);

% Reconstruction
% figure;
% plot3(points3D(:, 1), points3D(:, 2), points3D(:, 3));

% Draw the window uppper and bottom lines
plot([1, col], [sampley - halfRow, sampley - halfRow], 'r--');
plot([1, col], [sampley + halfRow, sampley + halfRow] ,'r--');

% Draw correlation figure
figure; plot(1:col, correlation(1, 1:col));







