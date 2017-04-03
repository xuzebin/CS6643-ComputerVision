% Estimate albedo, normal, and depth given at least 3 images and their light
% source directions.
% The input images are assumed to be Lambertian reflectance surfaces,
% and must be the same size

clc; 
close all;
format long;

% Recalculation may take long, we can read from precalculated results
RECAL = true;

% Read light source directions from a csv file.
% The file has the format: imgFile, x, y, z
% where x, y, z is the light source direction of the imgFile.
% For example:
% im1.png, 0, 0, 1
lightDirs = readtable('../res/dog-images/light_directions.csv');

% Read all images of different illumination from a directory
imgPath = '../res/dog-images/';
imgType = '*.png';
imgFiles = dir([imgPath imgType]);
num = length(imgFiles);

% Check if images exist
if( ~exist(imgPath, 'dir') || num < 1 )
    error('Directory not found or no matching images found.');
end

% Create an empty structure arrays with size num
% data: the images data
% dir: the light source direction vector
imgs(num).data = [];
imgs(num).dir = zeros(3, 1);

% The vector order in light_directions.csv file must correspond to
% that of imgFiles.
for i=1:num
    imgs(i).data = imread([imgPath imgFiles(i).name]);
    % If not grayscaled, convert it to grayscale
    if (size(imgs(i).data, 3) > 1)
        imgs(i).data = rgb2gray(imgs(i).data);
    end
    % Normalize the images to [0, 1]
    imgs(i).data = im2double(imgs(i).data);    
    imgs(i).dir = lightDirs(i, 2:4);
end

% Iterate every point of the image
[row, col] = size(imgs(1).data);


if (RECAL)
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
            normals(r, c, 1:3) = gvec / albedo(r, c);
        end
    end    
    size(normals)
    % Save calculated results
    size(g)
    dlmwrite('g.txt', g);
    dlmwrite('albedo.txt', albedo);
    dlmwrite('normals.txt', normals);

else % Read from precalculated results instead of recalculate it again
    albedo = dlmread('albedo.txt');
    normals = dlmread('normals.txt');
    % Normals is a 3d matrix, needs to be reshaped back to 3d matrix
    imgSize = size(normals, 1);% Assume width equals height
    normals = reshape(normals, imgSize, imgSize, 3);
    
    g = dlmread('g.txt');
    g = reshape(g, imgSize, imgSize, 3);
end


% Albedo map
figure;
imshow(albedo);
title('Albedo map');

% 2D surface normal map
figure;
spacing = 10;
[x, y] = meshgrid(1:spacing:row, 1:spacing:col);
quiver(x,y, normals(1:spacing:end, 1:spacing:end, 1),normals(1:spacing:end, 1:spacing:end, 2));
axis tight;
axis square;
title('Normal map');

% 3 components (x,y,z) of normalized surface normal vector
figure;
subplot(1,3,1);imshow(normals(1:end, 1:end, 1));
subplot(1,3,2);imshow(normals(1:end, 1:end, 2));
subplot(1,3,3);imshow(normals(1:end, 1:end, 3));
title('Surface normal vector in x, y, z');

% New shaded image with a new light source direction (1,2,3)
lightSrc = [1;2;3];
lightSrc = lightSrc / norm(lightSrc);
newShaded = zeros(row, col);
for r=1:row
    for c=1:col
%         newShaded(r, c) = dot(albedo(r, c) * [normals(r, c, 1);normals(r, c, 2);normals(r, c, 3)], lightSrc);
        newShaded(r, c) = dot([g(r,c,1); g(r,c,2); g(r,c,3)], lightSrc);
    end
end

figure;
imshow(newShaded);
title('New shaded image');



% p=dz/dx, q=dz/dy
p = zeros(row, col);
q = zeros(row, col);
for r=1:row
    for c=1:col
%         g = albedo(r, c) * [normals(r, c, 1);normals(r, c, 2);normals(r, c, 3)];   
        if (g(r, c, 3) == 0)
            p(r, c) = 1;
            q(r, c) = 1;
        else
            p(r, c) = g(r, c, 1) / g(r, c, 3);
            q(r, c) = g(r, c, 2) / g(r, c, 3);
            str = [num2str(g(r,c,1)), ' / ', num2str(g(r,c,3)), ' = ',num2str(p(r,c))];
            disp(str);
        end
    end
end


% Integration: Reconstruct height surface

if (RECAL)
    % Preallocate depth map
    depth = zeros(row, col);

    % Integrate along the first column (y)
    for r=2:row
%         g = albedo(r, 1) * [normals(r, 1, 1);normals(r, 1, 2);normals(r, 1, 3)];
%         if (g(r, c, 3) ~= 0)
%             p = g(r, c, 1) / g(r, c, 3);
%             q = g(r, c, 2) / g(r, c, 3);
%         else
%             q = 0;        
%         end
        depth(r, 1) = depth(r - 1, 1) + q(r, 1);
    end
    % Integrate along each row (x)
    for r=1:row
        for c=2:col
%             g = albedo(r, c) * [normals(r, c, 1);normals(r, c, 2);normals(r, c, 3)];
%             if (g(r, c, 3) ~= 0)
%                 p = g(r, c, 1) / g(r, c, 3);
%                 q = g(r, c, 2) / g(r, c, 3);
%             else
%                 p = 0;
%             end
            depth(r, c) = depth(r, c - 1) + p(r, c);
        end    
    end
    % Save depth result
    dlmwrite('depth.txt', depth);
else
    depth = dlmread('depth.txt');
end

% Depth map
figure;
depth = mat2gray(-depth);
imshow(depth);
title('Depth map');

% Render surface
[x, y] = meshgrid(1:row, 1:col);
figure;
surf(x, y, depth, 'EdgeColor', 'none');
camlight left;
lighting phong;
title('Surface');

% Render the mesh
figure;
mesh(x, y, depth);
title('Mesh');

