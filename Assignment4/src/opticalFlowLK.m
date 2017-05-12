%% Opitcal flow using Lucas and Kanade method
%%
function OF = opticalFlowLK(I1, I2, kernel_size)

kernel = ones(kernel_size);
ksize = size(kernel, 1);% Kernel size must be square
half_ksize = ksize / 2;
[rows, cols] = size(I1);

fimg1 = gaussian_filter(I1, 2);
fimg2 = gaussian_filter(I2, 2);
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

OF = zeros(rows, cols, 2, 'double');
A = zeros(ksize, 2, 'double');
bb = zeros(ksize, 1, 'double');
for i=(ceil(half_ksize)+1):(rows-ksize+ceil(half_ksize)-1)
    for j=(ceil(half_ksize)+1):(size(I1,2)-ksize+ceil(half_ksize)-1)
        t=1;
        for a=-ceil(half_ksize):ceil(half_ksize)
            for b=-ceil(half_ksize):ceil(half_ksize)
                A(t,1) = Ix(i+a,j+b);
                A(t,2) = Iy(i+a,j+b);
                bb(t) = -It(i+a,j+b);
                t=t+1;
            end
        end
        OF(i,j,:) = inv(A' * A) * A' * bb;
    end
end