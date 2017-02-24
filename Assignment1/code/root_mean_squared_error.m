measured = dlmread('measured_image_coordinates_96.txt');
measured_sorted = sortrows(measured,1);

estimated = dlmread('reconstructed_image_coordinates_96.txt');
estimated_sorted = sortrows(estimated,1);
%Mean squared error
mse = immse(measured_sorted, estimated_sorted);
%Root mean squared error
rmse = sqrt(mse)