% Calculate F matrix using 8 point algorithm
clc;
close all;
format long g

stereoParam = dlmread('stereo_param10.txt');
% stereoParam = dlmread('stereo_param_scene10.txt');
% stereoParam = dlmread('stereo_param_test.txt');
left = stereoParam(:, 1:2);
right = stereoParam(:, 3:4);

sampleCount = size(stereoParam, 1);


% Construct [uu', uv', u, u'v, vv', v, u', v', 1] where (u,v) is the pixel
% cooridnate of left image, and (u',v') is the pixel coordinate right
% image.

% Sovle 8 equations directly
% A = zeros(8, 8);
% for i=1:8
%     uL = stereoParam(i, 1);
%     vL = stereoParam(i, 2);
%     uR = stereoParam(i, 3);
%     vR = stereoParam(i, 4);
%     A(i, :) = [uL*uR, uL*vR, uL, uR*vL, vL*vR, vL, uR, vR];
% end
% mat
% 
% FVec = A\(-ones(8,1));
% F = vec2mat(FVec, 3)';
% F(3,3) = 1;


% Least Square Estimation (over 8 equations)
A = zeros(sampleCount, 9);
for i=1:sampleCount
    uL = stereoParam(i, 1);
    vL = stereoParam(i, 2);
    uR = stereoParam(i, 3);
    vR = stereoParam(i, 4);
    A(i, :) = [uL*uR, uL*vR, uL, uR*vL, vL*vR, vL, uR, vR, 1];
end

[U, S, V] = svd(A);

% Last column of V is the vector associated with the smallest eigenvalue
FVec = V(:, end);
F = reshape(FVec, 3, 3);


% Rank deprivation
[FU, FD, FV] = svd(F);
FD(3, 3) = 0;
F = FU * FD * FV';
FU
FV
FV = FV(:, 3) .* (1 ./ FV(3, 3))%left epipole
FU = FU(:, 3) .* (1 ./ FU(3, 3))%right epipole



% Draw left image with sampled points
% leftImg = imread('res/sceneL.JPG');
leftImg = imread('imL.png');
% leftImg = rgb2gray(leftImg);

subplot(121);imshow(leftImg);
title('left image');
hold on;
plot(stereoParam(:,1), stereoParam(:,2), 'go');
hold off;

% Draw Epipolar lines on left Image
hold on;
epipolarLinesL = zeros(sampleCount, 3);
for i=1:sampleCount
    pl = [stereoParam(i, 3:4), 1] * F;
    pl = pl / norm(pl(1:2));
    epipolarLinesL(i, 1:3) = pl';
end
points = lineToBorderPoints(epipolarLinesL, size(leftImg));
line(points(:, [1,3])', points(:, [2,4])');



% Draw epipole on left image
% plot(FV(1), FV(2), 'ro');



% Draw right image with sample points
% rightImg = imread('res/sceneR.JPG');
rightImg = imread('imR.png');
% rightImg = rgb2gray(rightImg);
subplot(122);imshow(rightImg);
title('right image');
hold on;
plot(stereoParam(:,3), stereoParam(:,4), 'go');

% Draw epipole on right image
% plot(FU(1), FU(2), 'ro');

% Draw Epipolar lines on right Image
hold on;
epipolarLines = zeros(sampleCount, 3);
[m, n] = size(leftImg);
right_epipolar_x = 1:2*m;
for i=1:sampleCount
    pr = F * [stereoParam(i, 1:2), 1]';
    pr = pr / norm(pr(1:2));
    epipolarLines(i, 1:3) = pr;
    
    right_epipolar_y = (-pr(3)-pr(1)*right_epipolar_x)/pr(2);
    plot(right_epipolar_x, right_epipolar_y);

end
% points = lineToBorderPoints(epipolarLines, size(rightImg));
% line(points(:, [1,3])', points(:, [2,4])');


hold off;

truesize;


% Calculate epipole
[u,d] = eigs(F'*F)
uu = u(:, 1)
epipole = uu / uu(3)


% imwrite(rectLeft, 'rectifiedLeft.png');
% imwrite(rectRight, 'rectifiedRight.png');

% figure;
% showMatchedFeatures(leftImg, rightImg, stereoParam(:, 1:2), stereoParam(:, 3:4));
% 


