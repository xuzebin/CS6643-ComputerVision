%Script for estimating camera calibration matrix using least square method and their intrinsic and
%extrinsic parameters as well as reconstruction and plotting.

clc; 
clear all;
close all;

format long g

%Read world coordinates of 12 points from file.
world_points = dlmread('world_coordinates_12.txt');
%Convert to homogenous coordinates
world_points = cart2hom(world_points);
%Read image coordiantes of 12 poitns from file.
image_points = dlmread('image_coordinates_12.txt');

%The form of Q matrix (every two rows):
%[Px, Py, Pz, 1, 0, 0, 0, 0, -uPx, -uPy, -uPz, -u]
%[0, 0, 0, 0, Px, Py, Pz, 1, -vPx, -vPy, -vPz, -v]

%12 measured points
num = size(image_points, 1);
%Q matrix (24x12)
Q = zeros(num*2, 12);
for i = 1:num
    %Odd rows
    Q(i*2-1, 1:4) = world_points(i, 1:4);%Assign world point coordinates
    Q(i*2-1, 9:12) = -image_points(i, 1).*world_points(i, 1:4);%-u*[Px, Py, Pz, 1]
    %Even rows
    Q(i*2, 5:8) = world_points(i, 1:4);%Assign world point coordinates
    Q(i*2, 9:12) = -image_points(i, 2).*world_points(i, 1:4);%-v*[Px, Py, Pz, 1]
end

%Perform SVD of Q
[U, S, V] = svd(Q);
[min_val, min_index] = min(diag(S(1:12, 1:12)));
m = V(1:12, min_index);
%Convert the 12x1 m vector to 3x4 a M matrix
M = vec2mat(m,4);

%Compute the scale factor ? (rho)
rho = 1 / norm(M(3, 1:3));
M = rho * M;

%Compute intrinsic and extrinsic parameters
a1 = M(1,1:3);
a2 = M(2,1:3);
a3 = M(3,1:3);
b = M(1:3, 4);

cross_a1a3 = cross(a1, a3);
cross_a2a3 = cross(a2, a3);

cosTheta = -dot(cross_a1a3, cross_a2a3) / (norm(cross_a1a3) * norm(cross_a2a3));
theta = acos(cosTheta);
skew = radtodeg(theta);

r3 = a3;
u0 = dot(a1, a3);
v0 = dot(a2, a3);

alpha = norm(cross_a1a3) * sin(theta);
beta = norm(cross_a2a3) * sin(theta);

r1 = cross_a2a3 / norm(cross_a2a3);
r2 = cross(r3, r1);
%K matrix
K = [alpha, -alpha * cot(theta), u0;
         0,   beta / sin(theta), v0;
         0,                   0,  1];

%Translation
t = K\b;%faster than inv(K)*b
%Rotation
R = [r1;r2;r3];

%Get Euler Angles in 3 directions
eulerZYX = rotm2eul(R);
eulerZYX = radtodeg(eulerZYX);

%Take the inverse of the extrinsic matrix (world to camera) for dicussion
%in extrinsic parameters.
extrinsicMatrix = [r1, t(1); r2, t(2); r3, t(3); 0, 0, 0, 1];
extrinsicMatrixInv = inv(extrinsicMatrix);

%Print the matrices
M
K
%Print the instrinsic parameters
theta
skew
u0
v0
alpha
beta
%Print the extrinsic parameters
R
eulerZYX
t

%Print the inverse of the extrinsic Matrix (R t) to be used as the camera
%pose
extrinsicMatrixInv



%Image coordinates reconstruction

%Plot points on the original image
figure(1);
img = imread('checkerboard.JPG');
image(img);
hold on;
plot_points(M, 'world_coordinates_96.txt', 'measured_image_coordinates_96.txt');
hold off;

%Plot points on the a white backgroud
figure(2);
axis([0 4032 0 3024]);
axis ij;
hold on;
reconstructed = plot_points(M, 'world_coordinates_96.txt', 'measured_image_coordinates_96.txt');
hold off;

%Save reconstructed image coordinates to file
dlmwrite('reconstructed_image_coordinates_96.txt', reconstructed);


