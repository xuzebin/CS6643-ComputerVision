%Plot reconstructed (red dots) and measured (blue circles) points and return the reconstructed image
%coordinates.
%Inputs:
% M: Camera calibration matrix.
% world_coordinates_file: The file storing the world coordinates used to reconstruct the image coordinates using M.
% measured_image_coordiantes_file: The file storing the measured image
% coordinates.
function reconstructed_coords = plot_points(M, world_coordinates_file, measured_image_coordiantes_file)
%Read measured world coordinates.
wc = dlmread(world_coordinates_file);
%Convert to homogenous coordinates
wc = cart2hom(wc);
point_number = size(wc, 1);
reconstructed_coords = zeros(point_number, 2);

%Plot the reconstructed image coordnates with red dots.
for i = 1:point_number
    z = M(3, 1:4) * wc(i, 1:4)';
    u = M(1, 1:4) * wc(i, 1:4)' / z;
    v = M(2, 1:4) * wc(i, 1:4)' / z;
    plot(u, v, 'r.', 'MarkerSize', 12);
    reconstructed_coords(i, :) = [u, v];
end

%Circle the 96 manually measured image coordinates
measured_coords = dlmread(measured_image_coordiantes_file);
plot(measured_coords(:,1), measured_coords(:,2), 'bo');

end