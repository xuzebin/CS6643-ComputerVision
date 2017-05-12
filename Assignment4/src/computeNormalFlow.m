% Compute Normal Flow
function n_flow = computeNormalFlow(img1, img2)

if size(img1) ~= size(img2)
    error('Image size not equal!')
end
[rows, cols] = size(img1);


fimg1 = gaussian_filter(img1, 1);
fimg2 = gaussian_filter(img2, 1);
fimg1 = im2double(fimg1);
fimg2 = im2double(fimg2);


% Temporal gradient
It = fimg2 - fimg1;

% Spatial gradients
Ix = zeros(rows, cols);
Iy = zeros(rows, cols);
for i=1:rows-1
    for j=1:cols-1
        Ix(i,j) = fimg1(i,j+1) - fimg1(i,j);
        Iy(i,j) = fimg1(i+1,j) - fimg1(i,j);
    end
end

% Compute normal flow
n_flow = zeros(rows, cols, 2, 'double');
for i=1:rows
    for j=1:cols
        Ix2Iy2 = Ix(i,j)^2 + Iy(i,j)^2;
        n_flow(i,j,1) = -It(i,j) * Ix(i,j) / Ix2Iy2;
        n_flow(i,j,2) = -It(i,j) * Iy(i,j) / Ix2Iy2;
    end
end







