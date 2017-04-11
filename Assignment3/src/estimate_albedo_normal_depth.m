% Estimate albedo, normal, and height given at least 3 images and their light
% source directions using least-square method.
% The input images are assumed to be Lambertian reflectance surfaces,
% and must be the same size.

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

% Create an empty structure arrays with size num
% data: the images data
% dir: the light source direction vector
imgs(num).data = [];
imgs(num).dir = zeros(3, 1);

% The vector order in light_directions.csv file must correspond to that of imgFiles.
for i=1:num
    fileName = char(table2array(table(i, 1)));    
    imgs(i).data = imread([filePath fileName]);
    % Convert it to grayscale if it is not
    if (size(imgs(i).data, 3) > 1)
        imgs(i).data = rgb2gray(imgs(i).data);
    end
    % Normalize the images to [0, 1] and store them
    imgs(i).data = im2double(imgs(i).data);
    % Store light source direction
    imgs(i).dir = table(i, 2:4);
end

% Get the size of the image.
[row, col] = size(imgs(1).data);

%% Calculate Albedo and normal
% Preallocate an albedo map
albedo = zeros(row, col);
% Prealocate a normal map
normals = zeros(row, col, 3);

% Stack multiple light source directions into one vector V
V = [table2array(imgs(1).dir);
     table2array(imgs(2).dir);
     table2array(imgs(3).dir);
     table2array(imgs(4).dir)];
% Normalize vectors
for k=1:size(V)
    V(k, :) = V(k, :) / norm(V(k, :));
end

% Preallocate g
g = zeros(row, col, 3);
% Preallocte p and q (p=dz/dx, q=dz/dy)
p = zeros(row, col);
q = zeros(row, col);

% Calculate albedo and normals
for r=1:row
    for c=1:col
        % Stack the same pixel from different images into a vector i
        i = [imgs(1).data(r, c);
             imgs(2).data(r, c);
             imgs(3).data(r, c);
             imgs(4).data(r, c)];
        % Create a diagonal matrix from the image vector i
        I = diag(i);

        g(r, c, 1:3) = pinv(I * V) * (I * i);

        gvec = [g(r,c,1); g(r,c,2); g(r,c,3)];
        albedo(r, c) = norm(gvec);
        if (albedo(r, c) == 0)
            normals(r, c, 1:3) = [0;0;0];
        else
            normals(r, c, 1:3) = gvec / albedo(r, c);
        end

        if (normals(r, c, 3) == 0)
            p(r, c) = 0;
            q(r, c) = 0;
        else
            p(r, c) = normals(r, c, 1) / normals(r, c, 3);
            q(r, c) = normals(r, c, 2) / normals(r, c, 3);
        end

        % Check if (dp/dy - dq/dx)^2 is small enough
        if (r > 1 && c > 1)
            dpdy = p(r, c) - p(r, c - 1);
            dqdx = q(r, c) - q(r - 1, c);
            error = (dpdy - dqdx).^2;
        end
    end
end

% Save calculated results
dlmwrite('g.txt', g);
dlmwrite('albedo.txt', albedo);
dlmwrite('normals.txt', normals);
dlmwrite('error.txt', error);


% Albedo map
figure;
imshow(albedo);
title('Albedo map');

% 2D surface normal map
figure;
spacing = 6;
[x, y] = meshgrid(1:spacing:row, 1:spacing:col);
quiver(x,y, normals(1:spacing:end, 1:spacing:end, 1),normals(1:spacing:end, 1:spacing:end, 2));
axis tight;
axis square;
title('Normal map');

% 3 components (x,y,z) of the normalized surface normal vector
figure;
subplot(1,3,1);imshow(normals(1:end, 1:end, 1));
subplot(1,3,2);imshow(normals(1:end, 1:end, 2));
subplot(1,3,3);imshow(normals(1:end, 1:end, 3));
suptitle('x, y, z components of the normalized surface normal vector');
truesize;

% New shaded image with a new light source direction (1,2,3)
lightSrc = [1;2;3];
lightSrc = lightSrc / norm(lightSrc);
newShaded = zeros(row, col);
for r=1:row
    for c=1:col
        newShaded(r, c) = dot([g(r,c,1); g(r,c,2); g(r,c,3)], lightSrc);
    end
end

figure;
imshow(newShaded);
title('New shaded image');
 
%% Reconstruct height surface via naive integration
% Preallocate height map
height = zeros(row, col);
% Integrate along the first column (y)
for r=2:row
    height(r, 1) = height(r - 1, 1) + q(r, 1);
end
% Integrate along each row (x)
for r=1:row
    for c=2:col
        height(r, c) = height(r, c - 1) + p(r, c);
    end    
end
% Save height result
dlmwrite('height.txt', height);

% Depth map
figure;
height = mat2gray(-height);
imshow(height);
title('Height map');

% Render surface
[x, y] = meshgrid(1:row, 1:col);
figure;
surf(-x,-y, height, 'EdgeColor', 'none');
camlight left;
lighting phong;
title('Surface');

% Render the mesh
figure;
mesh(-x, -y, height);
title('Mesh');

