% Compute optical flow over 2x2 pixel neighborhoods
clc; 
close all;

format long g

img1 = imread('../toy-car-images-bw/toy_formatted2.png');
img2 = imread('../toy-car-images-bw/toy_formatted3.png');

if size(img1) ~= size(img2)
    error('Image size not equal!')
end
[rows, cols] = size(img1);

figure;
imshow(img1);
figure;
imshow(img2);

% Gaussian filtering
fimg1 = gaussian_filter(img1, 2);
fimg2 = gaussian_filter(img2, 2);
fimg1 = im2double(fimg1);
fimg2 = im2double(fimg2);
figure;
imshow(fimg1);
title('Smoothed Image');

% Temporal gradient
It = fimg2 - fimg1;
figure;
imshow(mat2gray(It));
title('Temporal Gradient');


% Spatial gradients
Ix = zeros(rows, cols);
Iy = zeros(rows, cols);
for i=1:rows-1
    for j=1:cols-1
        Ix(i,j) = fimg1(i,j+1) - fimg1(i,j);
        Iy(i,j) = fimg1(i+1,j) - fimg1(i,j);
    end
end
figure;
imshow(mat2gray(Ix));
title('Spatial Gradient x');
figure;
imshow(mat2gray(Iy));
title('Spatial Gradient y');


% Compute optcial flow over 2x2 pixel neighborhood
o_flow = zeros(rows, cols, 2, 'double');
for i=2:(rows-1)
    for j=2:(cols-1)        
%         A = [Ix(i,j),     Iy(i,j);            
%              Ix(i+1,j),   Iy(i+1,j);
%              Ix(i,j+1),   Iy(i,j+1);
%              Ix(i+1,j+1), Iy(i+1,j+1)];
%         b = -[It(i,j); It(i+1,j); It(i,j+1); It(i+1,j+1)];
        A = [Ix(i-1,j), Iy(i-1,j);            
             Ix(i+1,j), Iy(i+1,j);
             Ix(i,j-1), Iy(i,j-1);
             Ix(i,j+1), Iy(i,j+1);
             Ix(i,j),   Iy(i,j)];
        b = -[It(i-1,j); It(i+1,j); It(i,j-1); It(i,j+1); It(i,j)];        

        o_flow(i,j,:) = pinv(A' * A) * (A' * b);
        
        if isnan(o_flow(i,j,1)) == 1
            disp('zerox');
        end
        if isnan(o_flow(i,j,2)) == 1
            disp('zerox');
        end  
    end
end

figure;
imshow(o_flow(:,:,1));
figure;
imshow(o_flow(:,:,2));

% Show optical flow
figure;
space = 1;
quiver(o_flow(1:space:end,1:space:end,1), o_flow(1:space:end,1:space:end,2), 15, 'color', 'g');
title('Optical Flow');
set(gca,'Ydir','reverse')
axis tight;

% Show overlayed optical flow
figure;
space = 1;
imshow(img1);
hold on;
quiver(o_flow(1:space:end, 1:space:end, 1), o_flow(1:space:end, 1:space:end,2), 15, 'color', 'g');
title('Overlayed Optical Flow');





