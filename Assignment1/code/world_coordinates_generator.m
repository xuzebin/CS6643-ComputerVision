% Calculate world cooridnates based on the square on the
% checkerboard.

% We use a dual checkerboards perpendicular to each other and glue them on
% the corner of the wall.
% The world coordinate is a right-handed orthonormal coordiantes.
% y-axis is upward, x-axis rightward, and z-axis leftward.
% Square's size is 2.8 cm.
% There is an offset (1.1 cm) in x and z axis.
% We calculate world coordinates of 6x8 points on each checkerboard
% (48x2=96 points in total).

world_coords = zeros(96, 3);%Left checkerboard (1:48), right checkerboard(49:96)
square_size = 2.8;
offset = 1.1;

% Left checkerboard where x = 0, z(1:6), y(1:8).
for y = 1:8
    for z = 1:6
        world_coords((y - 1) * 6 + z, 1:3) = [0.0, y * square_size, z * square_size + offset];
    end
end
% Right checkerboard where z = 0, x(1:7), y(1:9).
for y = 1:8
    for x = 1:6
        world_coords(48 + (y - 1) * 6 + x, 1:3) = [x * square_size + offset, y * square_size, 0.0];
    end
end

dlmwrite('world_coordinates_96.txt', world_coords)

