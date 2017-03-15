% Draw the corresponding epipolar line on the right by selecting a point on left image.
clc;
close all;
format long

% Hardcode the F matrix
F = [5.67727181723446e-06      -0.00104688474793069         0.168428043932593
     0.00104180185458174     -1.57926498477318e-05        -0.304294908613248
     -0.169371539518221         0.301673986785379         0.871401462189994];
 
 
leftImg = imread('../res/imL.png');
rightImg = imread('../res/imR.png');

% Select point on the left image
subplot(121),imshow(leftImg);
hold on;
[x, y] = ginput(1);
plot(x, y, 'go');
title(sprintf('selected point: (%f, %f)', x, y));

% Draw epipolar line on the right image
subplot(122),imshow(rightImg);
pr = F * [x, y, 1]';
pr = pr / norm(pr(1:2));
[m, n] = size(leftImg);
right_epipolar_x = 1:2*m;

right_epipolar_y = (-pr(3)-pr(1)*right_epipolar_x)/pr(2);
hold on;
plot(right_epipolar_x, right_epipolar_y);


% Select point on the right image
% subplot(122),imshow(rightImg);
% hold on;
% [x, y] = ginput(1);
% plot(x, y, 'go');
% title(sprintf('selected point: (%f, %f)', x, y));
% % Draw epipolar lines on the left image
% subplot(121),imshow(leftImg);
% epipolarLinesL = zeros(1, 3);
% pl = [x, y, 1] * F;
% pl = pl / norm(pl(1:2));
% 
% [m, n] = size(rightImg);
% left_epipolar_x = 1:2*m;
% 
% left_epipolar_y = (-pl(3)-pl(1)*left_epipolar_x)/pl(2);
% hold on;
% plot(left_epipolar_x, left_epipolar_y);


% Calculate epipole
[u,d] = eigs(F'*F)
uu = u(:, 1)
epipole= uu / uu(3)
% plot(epipole(1), epipole(2), 'ro');