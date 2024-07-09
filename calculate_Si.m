function Si = calculate_Si(transmitter_coords_fixed)
    % Fit a linear regression line
    coefficients = polyfit(transmitter_coords_fixed(:,1), transmitter_coords_fixed(:,2), 4);
    x_limits = xlim;
    x_fit = -10:0.1:x_limits(2)-70;
    y_fit = polyval(coefficients, x_fit);

    % Plot the regression line
    hold on;
    plot(x_fit, y_fit, 'b-', 'LineWidth', 1);

    % Get 12 points on the regression line
    x_fit_sampled = linspace(-10, x_limits(2)-70, 12);
    y_fit_sampled = polyval(coefficients, x_fit_sampled);
    plot(x_fit_sampled,y_fit_sampled,"bo")
    % Calculate Si values
    Si = [x_fit_sampled', y_fit_sampled'];
    
    % Plot Si points
    % plot(Si(:, 1), Si(:, 2), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
%     for i = 1:size(Si, 1)
%         % text(Si(i, 1), Si(i, 2)+1, ['S_{', num2str(i), '}'], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     end
    hold off;
end