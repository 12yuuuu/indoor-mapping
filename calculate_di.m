function di = calculate_di(transmitter_coords_fixed, Si)
    % Number of transmitter coordinates and Si points
    num_transmitters = size(transmitter_coords_fixed, 1);
    num_Si_points = size(Si, 1);
    
    % Initialize di matrix
    di = zeros(num_transmitters, num_Si_points);
    
    % Calculate distances between each transmitter coordinate and each Si point
    for i = 1:num_transmitters
        for j = 1:num_Si_points
            di(i, j) = sqrt(sum((transmitter_coords_fixed(i, :) - Si(j, :)).^2));
        end
    end
end