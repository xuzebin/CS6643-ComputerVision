% function compute_disparity(leftImg, rightImg)
% Compute disparity from 2 stereo image pair using correlation
%

L = imread('res/left_cube.JPG');
R = imread('res/right_cube.JPG');
L = rgb2gray(L);
R = rgb2gray(R);


% L = imread('imL.png');
% R = imread('imR.png');
% L = rgb2gray(L);
% R = rgb2gray(R);
L = imresize(L, 0.5);
R = imresize(R, 0.5);
% % 
% disparityRange = [-6 10];
% disparityMap = disparity(L,R,'BlockSize', 15);%, 'DisparityRange',disparityRange);
% figure
% imshow(disparityMap,disparityRange);

% L = rgb2gray(L);
% R = rgb2gray(R);




% arr = ones(7, 7)
% pLeft = padarray(arr, [(bRow-1)/2, (bCol-1)/2], 'both')
% center = floor((block+1)/2)
% for i=1:7
%     leftY = center(2) + i - 1;
%     for j=1:7
%         leftX = center(1) + j - 1;
%         leftBlock = pLeft((leftX-(bRow-1)/2):(leftX+(bRow-1)/2), (leftY-(bCol-1)/2):(leftY+(bCol-1)/2))
%     end
% end



% Padding
% paddedLeft = padarray(L, [(bRow-1)/2, (bCol-1)/2], 'both');
% paddedRight = padarray(R, [(bRow-1)/2, (bCol-1)/2], 'both');
% 
% disparityMap = zeros(row, col);
% center = floor((block+1)/2);
% for i=1:row % scan each horizontal epipolar line
%     currentRow = center(1) + i - 1; 
%     for j=1:col
%         leftCol = center(2) + j - 1;
%         % Extract a block
%         leftBlock = paddedLeft((currentRow-(bRow-1)/2):(currentRow+(bRow-1)/2), (leftCol-(bCol-1)/2):(leftCol+(bCol-1)/2));
%         % Reshape to a vector
%         w = reshape(leftBlock, [bRow * bCol, 1]);
%         % Convert to a zero-mean vector
%         w = double(w - mean(w));
%         
%         maxCorrelation = realmin('double');
%         bestMatch = 1;
%         % Compare this block with blocks on the same line in right image.
%         for line=1:col
%             currentCol = center(2) + line - 1;
%             rightBlock = paddedRight((currentRow-(bRow-1)/2):(currentRow+(bRow-1)/2), (currentCol-(bCol-1)/2):(currentCol+(bCol-1)/2));
%             w2 = reshape(rightBlock, [bRow * bCol, 1]);
%             w2 = double(w2 - mean(w2));
%             % Calculate normalized cross-correlation
%             ncc = dot(w, w2) / (norm(w) * norm(w2));
%             
%             if ncc > maxCorrelation
%                 maxCorrelation = ncc;
%                 bestMatch = line;
%             end
%         end
%         disparityMap(i, j) = (bestMatch - j);% R(i, bestMatch); 
%     end
% end

% Assume L and R has the same size.
[row, col] = size(L);

% Block size
bRow = 55;
bCol = 55;

disparityMap = zeros(row, col);
halfRow = (bRow - 1) / 2;
halfCol = (bCol - 1) / 2;

% for i=1:row % scan each horizontal epipolar line
%     for j=1:col
%         rowStart = max(i - halfRow, 1);
%         rowEnd = min(i + halfRow, row);
%         colStart = max(j - halfCol, 1);
%         colEnd = min(j + halfCol, col);
%         
%         % Construct vector w (left block vector)
%         w = L(rowStart:rowEnd, colStart:colEnd);
%         % Convert the matrix to vector (column-wise).
%         w = w(:);
%         % Convert to zero-mean vector
%         w = double(w - mean(w));
%         
%         % Compare w with w2 in the same line in right image.
%         maxCorrelation = eps(0);
%         bestMatchCol = 1;
%         for k=j:col
%             colStart2 = max(k - halfCol, 1);
%             colEnd2 = min(k + halfCol, col);
%             
%             % Construct w2 (right block vector)
%             w2 = R(rowStart:rowEnd, colStart2:colEnd2);
%             w2 = w2(:);
%             w2 = double(w2 - mean(w2));
% 
%             if size(w) == size(w2)% Ignore blocks with unequal size (at the corners).
%                 % Compute normalized cross correlation
%                 ncc = dot(w, w2) / (norm(w) * norm(w2));
%                 
%                 % Record the best match
%                 if ncc > maxCorrelation
%                     maxCorrelation = ncc;
%                     bestMatchCol = k;
%                 end
%             end
%         end
%         disparityMap(i, j) = bestMatchCol - j;
%     end
% end
% cross_correlation = zeros(1, col);
% for i=389:389 % scan each horizontal epipolar line
%     for j=551:551
%         rowStart = max(i - halfRow, 1);
%         rowEnd = min(i + halfRow, row);
%         colStart = max(j - halfCol, 1);
%         colEnd = min(j + halfCol, col);
%         
%         % Construct vector w (left block vector)
%         w = L(rowStart:rowEnd, colStart:colEnd);
%         % Convert the matrix to vector (column-wise).
%         w = w(:);
%         % Convert to zero-mean vector
%         w = double(double(w) - mean(w));
%         
%         
%         % Compare w with w2 in the same line in right image.
%         maxCorrelation = eps(0);
%         bestMatchCol = 1;
%         for k=j:col
%             colStart2 = max(k - halfCol, 1);
%             colEnd2 = min(k + halfCol, col);
%             
%             % Construct w2 (right block vector)
%             w2 = R(rowStart:rowEnd, colStart2:colEnd2);
%             w2 = w2(:);
%             
%             w2 = double(double(w2) - mean(w2));
%             
%             if size(w) == size(w2)% Ignore blocks with unequal size (at the corners).
%                 % Compute normalized cross correlation
%                 ncc = dot(w, w2) / (norm(w) * norm(w2));                
%                 
%                 cross_correlation(1, k) = ncc;
%                 % Record the best match
%                 if ncc > maxCorrelation
%                     maxCorrelation = ncc;
%                     bestMatchCol = k;
%                 end
%             end
%         end
%     end
% end
% 
% x=1:col;
% 
% figure; plot(x, cross_correlation(1, x));


% figure; imshow(cross_correlation);




% Padding
paddedL = padarray(L, [halfRow, halfCol], 'both');
paddedR = padarray(R, [halfRow, halfCol], 'both');

cross_correlation = zeros(1, col);

figure;
subplot(121); imshow(L);
axis on;
hold on;
[samplex, sampley] = ginput(1);
samplex = int16(samplex);
sampley = int16(sampley);

plot(samplex, sampley, 'g.', 'markerSize', 10);
rectangle('Position', [samplex - halfCol, sampley - halfRow, bCol, bRow],...
	'EdgeColor','r', 'LineWidth', 1)

plot([1, col], [sampley, sampley], 'g');



% Compare w with w2 in the same line in right image.
for i=sampley:sampley % scan each horizontal epipolar line
    for j=samplex:samplex
        rowStart = i;
        rowEnd = i + bRow - 1;
        colStart = j;
        colEnd = j + bCol - 1;
        
        % Construct vector w (left block vector)
        w = paddedL(rowStart:rowEnd, colStart:colEnd);
        % Convert the matrix to vector (column-wise).
        w = w(:);
        % Convert to zero-mean vector
        w = double(double(w) - mean(w));
        
        % Compare w with w2 in the same line in right image.
        maxCorrelation = eps(0);
        bestMatchCol = 1;
        for k=1:col
            colStart2 = k;
            colEnd2 = k + bCol - 1;
            
            % Construct w2 (right block vector)
            w2 = paddedR(rowStart:rowEnd, colStart2:colEnd2);
            w2 = w2(:);
            
            w2 = double(double(w2) - mean(w2));
            if (size(w) ~= size(w2))
                throw(MException('not equal'))
            end
            
            % Compute normalized cross correlation
            ncc = dot(w, w2) / (norm(w) * norm(w2));                

            cross_correlation(1, k) = ncc;
            
            % Record the best match
            if ncc > maxCorrelation
                maxCorrelation = ncc;
                bestMatchCol = k;
            end            
        end
        
        subplot(122); imshow(R);
        hold on;
        axis on;
        plot(bestMatchCol, sampley, 'g.', 'markerSize', 10);
        
        rectangle('Position', [bestMatchCol - halfCol, sampley - halfRow, bCol, bRow],...
            'EdgeColor','r', 'LineWidth', 1)
        plot([1, col], [sampley, sampley] ,'g');
        plot([1, col], [sampley - halfRow, sampley - halfRow], 'r--');
        plot([1, col], [sampley + halfRow, sampley + halfRow] ,'r--');
        x=1:col;
        figure; plot(x, cross_correlation(1, x));
    end
end



disp(bestMatchCol)
disp(test_point)
% TODO estimate depth: z=fB/d
