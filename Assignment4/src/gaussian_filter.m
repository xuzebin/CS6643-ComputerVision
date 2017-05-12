%
% Filters an image using a Gaussian filter. The code is meant only 
% as a teaching tool -- it is not optimized for speed (the Gaussian
% is separable, etc. or with matrix/vector operations).
% 
% inputs:  
% -------
% im_in - the input image (grayscale uint8)
% sigma - the variance of the 2D isotropic Gaussian
%
% outputs:  
% --------
% im_out - the filtered image (grayscale uint8)
%
function im_out = gaussian_filter(im_in, sigma)

    %% This function only works for grayscale images
    if (size(im_in,3) > 1)
        error('This function only works for grayscale images.');
    end

    %% Create the gaussian filter

    % Heuristic for the size of the Gaussian filter based on sigma
    n = ceil(sigma*3)*2+1;
    % Allocate space for the output
    gauss_filter = zeros(n);
    
    % Use the same sigma for x and y (isotropic)
    sigma_x = sigma;
    sigma_y = sigma;
    
    % This Gaussian is not rotated
    theta = 0;
   
    % Compute coefficients
    a = (cos(theta)^2 / (2*sigma_x^2)) + (sin(theta)^2 / (2*sigma_y^2));
    b = (-sin(2*theta) / (4*sigma_x^2)) + (sin(2*theta) / (4*sigma_y^2));
    c = (cos(theta)^2 / (2*sigma_y^2)) + (sin(theta)^2 / (2*sigma_x^2));
          
    % These are the array indices
    i=1;
    j=1;
    
    % Loop over size
    for u=linspace(-n/2,n/2,n)
        
        % Loop over size
        for v=linspace(-n/2,n/2,n)
       
            % Compute the Gaussian at this location
            gauss_filter(i,j) = exp(-(a*u^2 + 2*b*u*v + c*v^2));
            
            % Increment the j index
            j=j+1;
            
        end
        
        % Reset the j index
        j=1;
        % Increment the i index
        i=i+1;
        
    end

    %% Apply the filter to the image
    
    % Convert input image to double
    im_in = double(im_in);
    % Get the size of the input image
    [rows, cols] = size(im_in);
    % Create a blank canvas for the output image
    im_out = double(zeros(rows, cols));
    
    % This is how much we need to offset from the center
    half_filt = floor(n/2);
    % This is how much bigger our padded image needs to be
    extra = half_filt*2;
    
    % Create the padded image to deal with boundaries
    padded_im = double(zeros(rows+extra, cols+extra));
    % Place the input image in the center of the padded image
    padded_im(half_filt+1:rows+half_filt, half_filt+1:cols+half_filt) = im_in;
    
    % Loop over all rows of the input image
    for i=1:rows
        
        % Loop over all columns of the input image
        for j=1:cols
            
            % Compute current x and y pixel locations in the padded image
            x = i+half_filt;
            y = j+half_filt;
            
            % Get the neighborhood matrix
            neighborhood = padded_im(x-half_filt:x+half_filt, y-half_filt:y+half_filt);
            % Apply the weights from the Gaussian filter (element wise)
            weighted_neighborhood = gauss_filter.*neighborhood;
            % The value is the sum of the Gaussian weighted neighborhood
            filtered_value = sum(weighted_neighborhood(:));
        
            % This is the new value for the output image
            im_out(i,j) = filtered_value;
       
        end
        
    end
    
    % Convert to range 0-1
    im_out = im_out./max(im_out(:));
    % Convert to uint8
    im_out = uint8(im_out.*255);
    
end