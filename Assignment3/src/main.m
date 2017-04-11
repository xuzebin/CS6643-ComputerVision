clc; 
close all;
format long;

% filePath = '../res/synth-images/';
% filePath = '../res/sphere-images/';
filePath = '../res/dog-images/';

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

[albedo, normal] = computeAlbedoSurfNorm(imgs, lightSrcDir);

% Albedo map
figure;
imshow(albedo);
title('Albedo map');

% 2D surface normal vectors
figure;
spacing = 6;
[x, y] = meshgrid(1:spacing:row, 1:spacing:col);
quiver(x,y, normal(1:spacing:end, 1:spacing:end, 1),normal(1:spacing:end, 1:spacing:end, 2));
axis tight;
axis square;
title('Normal vectors');

% 3 components (x,y,z) of the normalized surface normal vector
figure;
subplot(1,3,1);imshow(normal(1:end, 1:end, 1));
subplot(1,3,2);imshow(normal(1:end, 1:end, 2));
subplot(1,3,3);imshow(normal(1:end, 1:end, 3));
suptitle('x, y, z components of the normalized surface normal vector');
truesize;

% Get p, q from normal
% Preallocte p and q (p=dz/dx, q=dz/dy)
p = zeros(row, col);
q = zeros(row, col);

for r=1:row
    for c=1:col
         % Compute p (dzdx) and q (dzdy)
         if (normal(r, c, 3) == 0)          
            p(r, c) = 0;
            q(r, c) = 0;
         else
            p(r, c) = normal(r, c, 1) / normal(r, c, 3);
            q(r, c) = normal(r, c, 2) / normal(r, c, 3);
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

% Create a mask used to distinguish the object from background.
bw = im2bw(albedo, 0);
mask = imfill(bw, 'holes');

figure;
imshow(mask);
title('mask');

% Naive integration (TV-scan)
%z = computeHeightMap(p, q, mask);
% frankotchellappa's method
z = frankotchellappa(p, q);

% Since our camera is looking down from positive z position
% We need to reverse the sign here.
z = -z;

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

% Render the mesh
figure;
mesh(x, y, z);
title('Mesh');
