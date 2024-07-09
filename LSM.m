% Fixed transmitter coordinates (commented out)
% transmitter_coords_fixed = [
%     0, 0;
%     30, 50;
%     50, 70;
%     70, 30;
%     90, 50;
%     70, 10
% ];

%{
    Currently, running the program will only display two lines 
    and the positions of 10 receivers on the graph.

    The blue line represents the connection of S (receiver).
    The green line represents the connection of the simulation result.

    If you want to view other icons on the figure, you can find the
    corresponding plot command and uncomment it. :)
%}

% Generate random fixed transmitter coordinates
% This creates a 6x2 matrix of random coordinates within a 100x100 area.
clear;
transmitter_coords_fixed = rand(10, 2) * 100;
transmitter_coords_fixed(1,1) = 0;
transmitter_coords_fixed(1,2) = 0;
% Create a new figure for plotting
figure;
hold on;
% Set the x and y axis limits for the plot
xlim([-60, 140]);
ylim([-40, 140]);

% Add labels to the transmitters (commented out)
% The text function would label each point with P_1, P_2, etc.
% text(transmitter_coords_fixed(:,1), transmitter_coords_fixed(:,2), ...
%     {'P_1'; 'P_2'; 'P_3'; 'P_4'; 'P_5'; 'P_6'}, ...
%     'VerticalAlignment','bottom','HorizontalAlignment','right');

% Plot the fixed transmitter coordinates as red circles
% 'ro' specifies red circles, 'MarkerSize', 4 sets the size, 
% and 'MarkerFaceColor', 'r' fills the circles with red color.
plot(transmitter_coords_fixed(:,1), transmitter_coords_fixed(:,2), ...
     'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
hold off;

% Label the x and y axes, and add a title to the plot
xlabel('X-axis');
ylabel('Y-axis');
title('Indoor Positioning');
grid on;
% Ensure the aspect ratio is equal so that distances are represented accurately
axis equal;

% Calculate Si, which might be the coordinates of some reference points
% based on the transmitter coordinates (the function is assumed to be predefined)
Si = calculate_Si(transmitter_coords_fixed);
disp("Si");
disp(Si);

% Calculate di, which could be the distances between points (function assumed predefined)
di = calculate_di(transmitter_coords_fixed, Si);
disp("di");
disp(di);

% Given x and y values for some known points
% These might represent coordinates of certain transmitters or points of interest
x_values = transmitter_coords_fixed(:,1);
y_values = transmitter_coords_fixed(:,2);

% Calculate the number of known points
N = length(x_values);

% Extract and display the distances of the first transmitter to all Si points
% This creates a 12x1 vector to store the distances
d1 = zeros(12, 1);
for i = 1:12
    d1(i) = di(1, i);
end

disp("d1");
disp(d1);

% Plot circles based on distances di(7, 1) to di(7, 6) (commented out)
% This would create dashed circles around each transmitter position
% hold on;
% theta = linspace(0, 2*pi, 100);
% for i = 1:6
%     x_circle = di(i, 7) * cos(theta) + transmitter_coords_fixed(i, 1);
%     y_circle = di(i, 7) * sin(theta) + transmitter_coords_fixed(i, 2);
%     plot(x_circle, y_circle, 'k:', 'LineWidth', 1); % Black dashed circle with increased thickness
% end
% hold off;

% Initialize a 12xN matrix to store ci values
ci = zeros(12, N);
% Calculate ci values for each pair of points
for i = 1:12
    for j = 1:N
        ci(i, j) = (x_values(j)^2 + y_values(j)^2 - di(j, i).^2 + d1(i).^2);
    end
end
disp("ci");
disp(ci);

% Initialize variables to store sums for each set of ci
sum_xicj = zeros(12, 1);
sum_yicj = zeros(12, 1);

% Calculate matrix elements for later use
sum_xj2 = sum(x_values(2:end).^2);
sum_yj2 = sum(y_values(2:end).^2);
sum_xiyj = sum(x_values(2:end) .* y_values(2:end));

% Calculate vector elements for each set of ci
for i = 1:12
%     disp("strat");
%     disp(x_values(2:end)');
%     disp(ci(i,2:end));
%     disp(x_values(2:end)'.*ci(i,2:end));
%     disp(sum(x_values(2:end).*ci(i,2:end)));
%     disp("end")
    sum_xicj(i) = sum(x_values(2:end)'.*ci(i,2:end));
    sum_yicj(i) = sum(y_values(2:end)'.*ci(i,2:end));
end

% Construct the matrices and vectors for solving linear equations
A = zeros(2, 2, 12);
b = zeros(12, 2);
for i = 1:12
    A(:,:,i) = [sum_xj2, sum_xiyj; sum_xiyj, sum_yj2];
    b(i, :) = [sum_xicj(i); sum_yicj(i)];
end

% Calculate the results for each set of ci
% Using the formula result = 0.5 * inv(A) * b
result = zeros(2, 12);
for i = 1:12
    result(:,i) = 0.5 * inv(A(:,:,i)) * b(i, :)';
end

% Display the results for each set of ci
disp('Results for each set of ci:');
for i = 1:12
    fprintf('Set %d: x = %f, y = %f\n', i, result(1, i), result(2, i));
end

% Plot the resulting points
hold on;
for i = 1:12
    % plot(result(1, i), result(2, i), 'gx', 'MarkerSize', 10, 'LineWidth', 2); % Green crosses for result points
    % text(result(1, i), result(2, i)-2, sprintf('S_r%d', i), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center'); % Label the points
    % Connect all result points with green lines
    plot(result(1, :), result(2, :), 'g-x'); % Green lines connecting result points
    % Plot lines connecting Si and result points
    % plot([Si(i, 1), result(1, i)], [Si(i, 2), result(2, i)], 'r-'); % Red lines
end
hold off;

% Calculate distances (radius of the circles)
% This computes the Euclidean distance between each Si and the corresponding result point
distances = zeros(12, 1);
for i = 1:12
    distances(i) = sqrt(sum((Si(i, :) - result(:, i)').^2));
end
disp(distances);

hold on;