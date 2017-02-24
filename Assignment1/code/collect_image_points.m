% Collect and measure 12 points and save the image coordianates into file.
clc; 
clear all;
close all;
img = imread('checkerboard.JPG');
imagesc(img);

[x, y] = ginput(12);
points = [x, y];
hold on;
for i = 1:size(points, 1)
    plot(points(i, 1), points(i, 2), 'r.', 'MarkerSize', 15);
    scatter(points(i, 1), points(i, 2), 80, 'b');
end
hold off;
dlmwrite('image_coordinates_12.txt', points);