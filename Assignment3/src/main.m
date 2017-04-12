%% PHOTOMETRIC STEREO ASSIGNMENT
% Estimate albedo, surface normal and reconstruct surface (shape) using at least 3 images under different illuminations. 
% Assume Lambertian surface reflectance is used.
%
% It shows:
% 1. Albedo map
% 2. 2D normal vectors
% 3. xyz components of the surface normal vector
% 4. New shaded image with a new chosen light source
% 5. Height map
% 6. Reconstructed Surface
% 7. Reconstructed mesh
% Note:
% Under each image directory there must be a csv file named 'light_directions.csv' describing the
% mapping of image filename and light source vector.
%
% Author: Zebin Xu (zebinxu@nyu.edu)
%%
clc; 
close all;
format long;

% filePath = '../res/synth-images/';
% filePath = '../res/sphere-images/';
filePath = '../res/dog-images/';

%% Parse csv file and preprocess and store corresponding images and light source vectors
% Read file names and light source directions from a csv file.
% Each row in the file is: imageFileName, x, y, z
% where x, y, z is the light source direction of the corresponding image.
% For example:
% im1.png, 0, 0, 1
% This means that the image named im1.png has a light direction [0 0 1]
lightSrcFile = strcat(filePath, 'light_directions.csv');
table = readtable(lightSrcFile);
num = size(table, 1);

% Check if the directory exists
if( ~exist(filePath, 'dir'))
    error('Directory not found');
end
% Check if there are at least 3 images 
if( num < 3 )
    error('There must be at least 3 images in the directory.');
end

% Create a numx1 cell array
imgs = cell(num, 1);
% Create a matrix for light source vectors
lightSrcDir = zeros(num, 3);

% The vector order in light_directions.csv file must correspond to that of imgFiles.
for i=1:num
    fileName = char(table2array(table(i, 1)));    
    imgs{i} = imread([filePath fileName]);
    % Convert it to grayscale if it is not
    if (size(imgs{i}, 3) > 1)
        imgs{i} = rgb2gray(imgs{i});
    end
    % Normalize the images to [0, 1] and store them
    imgs{i} = im2double(imgs{i});
    % Store light source direction
    lightSrcDir(i, :) = table2array(table(i, 2:4));
end

% Get the size of the image.
[row, col] = size(imgs{1});

%% Compute albedo and surface normal from multiple images and corresponding light source direction vectors
[albedo, normal] = computeAlbedoSurfNorm(imgs, lightSrcDir);

%% Estimated albedo map
figure;
imshow(albedo);
title('Albedo map');

%% Estimated surface normals
% 2D surface normal vectors
figure;
spacing = 6;
[x, y] = meshgrid(1:spacing:row, 1:spacing:col);
quiver(x,y, normal(1:spacing:end, 1:spacing:end, 1),normal(1:spacing:end, 1:spacing:end, 2));
axis tight;
axis square;
title('2D Normal vectors');

% x, y, z components of the surface normal vector
figure;
subplot(1,3,1);imshow(normal(1:end, 1:end, 1));
subplot(1,3,2);imshow(normal(1:end, 1:end, 2));
subplot(1,3,3);imshow(normal(1:end, 1:end, 3));
suptitle('x, y, z components of the surface normal vector');
truesize;

% Normal map
figure;
imshow(normal);

% New shaded image with a new light source direction (1,2,3)
lightSrc = [1;2;3];
lightSrc = lightSrc / norm(lightSrc);
newShaded = zeros(row, col);
for r=1:row
    for c=1:col
        g = [normal(r, c, 1) * albedo(r, c); 
             normal(r, c, 2) * albedo(r, c); 
             normal(r, c, 3) * albedo(r, c)];
        newShaded(r, c) = dot(g, lightSrc);
    end
end

figure;
imshow(newShaded);
title('New shaded image with light source vector [1 2 3]');

%% Compute p, q (p=dz/dx, q=dz/dy)
% Preallocte p and q 
p = zeros(row, col);
q = zeros(row, col);

for r=1:row
    for c=1:col
         % Compute p (dzdx) and q (dzdy)
         if (normal(r, c, 3) == 0)          
            p(r, c) = 0;
            q(r, c) = 0;
         else
            p(r, c) = -normal(r, c, 1) / normal(r, c, 3);
            q(r, c) = -normal(r, c, 2) / normal(r, c, 3);
         end

        % Check if (dp/dy - dq/dx)^2 is small enough
        if (r > 1 && c > 1)
            dpdy = p(r, c) - p(r, c - 1);
            dqdx = q(r, c) - q(r - 1, c);
            err = (dpdy - dqdx).^2;            
            if (err > 1.0)
                disp('(dp/dy - dq/dx)^2 is not small enough');                
            end
        end
    end
end

% Create a mask used to segment the object from background.
% bw = im2bw(albedo, 0.9);
bw = imbinarize(albedo);
mask = imfill(bw, 'holes');
imshow(mask);
title('mask');

%% Reconstruct height surface via integration
% Naive integration (TV-scan)
% z = computeHeightMap(p, q, mask);

% frankotchellappa's method
z = frankotchellappa(p, q);

%% Show reconstructed surface
% Height map
figure;
% Convert height matrix to grayscale image
height = mat2gray(z);
imshow(height);
title('Height map');

% Render surface
[x, y] = meshgrid(1:row, 1:col);
figure;
surf(x, y, z, 'EdgeColor', 'none');
camlight left;
lighting phong;
title('Surface');

% Render the wireframe mesh
figure;
mesh(x, y, z);
title('Wireframe Mesh');
