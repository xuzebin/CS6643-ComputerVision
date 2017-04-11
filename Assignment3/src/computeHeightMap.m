%% Compute height map given p and q using naive integration.
%
% Usage:     z = computeHeightMap(p, q, mask)
%
% Arguments: p    - dz/dx, a 2D matrix giving the partial derivative in x direction.
%            q    - dz/dy, a 2D matrix giving the partial derivative in y direction.
%            mask - a binary image describing the 2d shape of the object, 1
%            indicates the object and 0 indicates the background. With this mask
%            the noise of height map can be reduced significantly.
%
% Returns:   z - a 2D matrix giving the surface heights.
%
% Author: Zebin Xu (zebinxu@nyu.edu)
%%
function z = computeHeightMap(p, q, mask)
    if ~all(size(p) == size(q))
      error('p and q matrices must match');
    end
        
    [row, col] = size(p);
    
    if (nargin < 3)
        mask = ones(row, col);        
    end
    
    % Preallocate height map
    z = zeros(row, col);
    
    % Integrate along the first column (y)
    for r=2:row
        if mask(r, 1) > 0            
            z(r, 1) = z(r - 1, 1) + q(r, 1);
        else
            z(r, 1) = 0;
        end
    end
    % Integrate along each row (x)    
    for r=1:row
        for c=2:col
            if mask(r, c) > 0            
                z(r, c) = z(r, c - 1) + p(r, c); 
            else
                z(r, c) = 0;
            end
        end                     
    end
    