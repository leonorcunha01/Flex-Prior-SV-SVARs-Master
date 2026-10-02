function plot_h_comp(array1, array2)
    % Determine the number of variables (columns) in the arrays
    num_vars = size(array1, 2);

    % Create a new figure
    figure;

    % Loop over each variable to create a subplot
    for i = 1:num_vars
        subplot(num_vars, 1, i); % 1 row, num_vars columns, ith subplot
        hold on; % Hold on to plot multiple series on the same subplot

        % Plot the data for the ith variable
        plot(array1(:, i), 'b-', 'DisplayName', 'real');
        plot(array2(:, i), 'r-', 'DisplayName', 'estim');

        % Set the title and labels
        title(['h', num2str(i)]);
        xlabel('time');
        ylabel('value');
        % Show legend and grid
        legend show;
        grid on;

        hold off; % Release the hold for the next subplot
    end
end