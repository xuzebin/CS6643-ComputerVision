function [bestMatch, correlation] = findMatchPointFromStereoPair(leftPoint, imgL, imgR, winSize, maxDisparity)
% Find the matched right point using normalized cross correlation for one point in the left.
%
% The pipeline is optimized to reduce the time complexity of NCC.
% 1. For each pixel(u1,v1) in the left image, find corresponding pixels (windows) in the right image
% starting from (u2=u1, v2=v1) and decreases v2 down to (v1 - DISPARITY_MAX)
% manually set by user (Experimentally, we can find the maximum disparity
% by mannually selecting a few correspongding closest pixels and compute their
% disparities). In this way, we only need to compare DISPARITY_MAX number
% of windows instead of the whole line.
%
% Input: 
% leftPoint = image coordinate of the left point.
% imgL, imgR = left and right rectified images. Must have the same size and be
% grayscale.
% winSize = window size for block matching.
% maxDisparity = an estimated maximum disparity in the image to speed up
% computation.
%
% Output:
% bestMatch = the matched column
% correlation = 1D row vector of the correlation result in one line.

    if nargin == 4
        maxDisparity = intmax('int16'); % Default maximum disparity
    end

    % imgL and imgL must have the same size and be grayscale
    if size(imgL) ~= size(imgL)
        error('Stereo image pair must have the same length');
    end
    if size(imgL, 3) ~= 1
        error('image must be grayscale');
    end
    
    
    imgCol = size(imgL, 2);
    
    winRow = double(winSize(1));
    winCol = double(winSize(2));
    
    % Prepocessing: zero paddding
    paddedL = padarray(imgL, [(winRow - 1) / 2, (winCol - 1) / 2], 'both');
    paddedR = padarray(imgR, [(winRow - 1) / 2, (winCol - 1) / 2], 'both');
   
    % Correlation result of a line in the image
    correlation = zeros(1, imgCol);

    j = leftPoint(1);
    i = leftPoint(2);
    
    rowStart = i;
    rowEnd = i + winRow - 1;
    colStart = j;
    colEnd = j + winCol - 1;

    % Construct vector w (left window vector)
    w = paddedL(rowStart:rowEnd, colStart:colEnd);
    % Convert the matrix to vector (column-wise transformation).
    w = w(:);
    % Convert to zero-mean vector
    w = double(double(w) - mean(w));

    % Compare w with w2 on the same line in right image.
    maxCorrelation = eps(0);
    bestMatch = j;
    
    for k = j:-1:max(j - maxDisparity, 1)
        colStart2 = k;
        colEnd2 = k + winCol - 1;
    
        % Construct w2 (vector of the right window)
        w2 = paddedR(rowStart:rowEnd, colStart2:colEnd2);
        w2 = w2(:);
        w2 = double(double(w2) - mean(w2));

        if size(w) ~= size(w2)
            error('w and w2 must be the same length')
        end

        % Normalized cross correlation
        ncc = sum(w.*w2) / sqrt((sum(w.^2).*sum(w2.^2)));
        % A more concised equation: ncc = dot(w, w2) / (norm(w) * norm(w2));      

        correlation(1, k) = ncc;

        % Find the best match
        if ncc > maxCorrelation
            maxCorrelation = ncc;
            bestMatch = k;
        end
    end
end
