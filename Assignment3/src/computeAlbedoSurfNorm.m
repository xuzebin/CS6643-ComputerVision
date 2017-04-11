%% Compute albedo and surface normal given images and corresponding light
% source vectors using least-square method.
%
% Usage:     [albedo, normal] = computeAlbedoSurfNorm(images, V)
%
% Arguments: images - an mx1 cell array of images under different illuminations
%                       all the images must be grayscale and have same size
%              V    - a stacked matrix (mx3) of corresponding light source vectors
%                      where m is the number of images.
%
% Returns:   albedo - a 2D matrix with the albedo for each point on the image.
%           normal - a rowxcolx3 matrix with the 3d normal vector for each point.
%
% Author: Zebin Xu (zebinxu@nyu.edu)
%%
function [albedo, normal] = computeAlbedoSurfNorm(images, V)

if (size(images, 1) < 3)
    error('Cannot infer albedo and surface normal from less than 3 images');
end

% Normalize light source vectors
for k=1:size(V)
    V(k, :) = V(k, :) / norm(V(k, :));
end

% Get the size of the image
[row, col] = size(images{1});

% Preallocate an albedo map
albedo = zeros(row, col);
% Prealocate a matrix for surface normal vectors
normal = zeros(row, col, 3);

% Calculate albedo and normal for each point
for r=1:row
    for c=1:col
        % Stack image pixel values for different images into a vector
        i = [images{1}(r, c);images{2}(r, c);images{3}(r, c);images{4}(r, c)];

        g = (V.' * V) \ (V.' * i);
    
        albedo(r, c) = norm(g);
        if (albedo(r, c) > 1)
            albedo(r, c) = 1;
            normal(r, c, 1:3) = g / albedo(r, c);
        elseif (albedo(r, c) == 0)
            albedo(r, c) = 0;
            normal(r, c, 1:3) = g;
        else
            normal(r, c, 1:3) = g / albedo(r, c);
        end
        
    end
end
