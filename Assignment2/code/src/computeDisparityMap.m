% Compute disparity map from 2 rectified stereo image pair using normalized cross correlation (NCC).

% Read image stereo pair
L = imread('../res/cubeL.JPG');
R = imread('../res/cubeR.JPG');

% Convert to grayscale
L = rgb2gray(L);
R = rgb2gray(R);

% % Downsize to speed up computation
L = imresize(L, 0.05);
R = imresize(R, 0.05);

[row, col] = size(L);

% Define a ratio of image's size divided by window's size to be used for calculating
% proper window size for correlation.
imageWinRatio = 17;

% Window size
winRow = int32(row / imageWinRatio);
winCol = int32(col / imageWinRatio);

% Estimate a disparity maximum to reduce unnecessary computations.
% This depends greatly on the scene depth of the image.
maxDisparity = winCol * 3;


% % Half of window size
halfRow = (winRow - 1) / 2;
halfCol = (winCol - 1) / 2;
halfRow = double(halfRow);
halfCol = double(halfCol);

% Prepocessing: image paddding
paddedL = padarray(L, [halfRow, halfCol], 'both');
paddedR = padarray(R, [halfRow, halfCol], 'both');

% Disparity Map
disparityMap = zeros(size(L));

% Normalized cross correlation to compute disparity of two images.
for i=1:row
    for j=1:col
        rowStart = i;
        rowEnd = i + winRow - 1;
        colStart = j;
        colEnd = j + winCol - 1;
        
        % Construct vector w (left block vector)
        w = paddedL(rowStart:rowEnd, colStart:colEnd);
        % Convert the matrix to vector (column-wise).
        w = w(:);
        % Convert to zero-mean vector
        w = double(double(w) - mean(w));
        
        % Compare w with w2 in the same line in right image.
        maxCorrelation = eps(0);
        bestMatchCol = 1;
        for k = j:-1:max(j - maxDisparity, 1)
            colStart2 = k;
            colEnd2 = k + winCol - 1;
            
            % Construct w2 (right block vector)
            w2 = paddedR(rowStart:rowEnd, colStart2:colEnd2);
            w2 = w2(:);
            w2 = double(double(w2) - mean(w2));
            
            if size(w) ~= size(w2)
                error('w and w2 must be the same length')
            end
            
            % Compute normalized cross correlation
            ncc = sum(w.*w2) / sqrt((sum(w.^2).*sum(w2.^2)));

            % Record the best match
            if ncc > maxCorrelation
                maxCorrelation = ncc;
                bestMatchCol = k;
            end
        end
        
        disparityMap(i, j) = j - bestMatchCol;
    end
end


figure;
% Normalize
disparityMap = double(disparityMap) ./ double(max(disparityMap(:)));
imshow(disparityMap)


