function plot_array_comp(array_list, base_name, xvals, y_name)
    % plot_array_comp(array_list, base_name, xvals, y_name)
    % Plots each column of multiple arrays (same size) against xvals in separate subplots.
    %
    % INPUTS:
    %   array_list : cell array containing the arrays to compare (e.g. {A, B, C})
    %   base_name  : string used as base for subplot titles
    %   xvals      : (optional) vector for x-axis values
    %   y_name     : (optional) string for y-axis label


    % Handle optional inputs
    if nargin < 3 || isempty(xvals)
        xvals = 1:size(array_list{1}, 1);
    end
    if nargin < 4 || isempty(y_name)
        y_name = 'value';
    end

    % Determine number of arrays and variables
    num_arrays = numel(array_list);
    num_vars = size(array_list{1}, 2);

    % Colors for plotting (extend automatically if needed)
    colors = lines(num_arrays);

    % Create new figure
    figure;

    for i = 1:num_vars
        subplot(num_vars, 1, i);
        hold on;

        % Plot each array's ith column
        for j = 1:num_arrays
            plot(xvals, array_list{j}(:, i), ...
                'Color', colors(j, :), ...
                'LineWidth', 1, ... 
                'DisplayName', ['Array ', num2str(j)]);
        end

        title([base_name, ' ', num2str(i)]);
        xlabel('time');
        ylabel(y_name);
        legend show;
        grid on;
        hold off;
    end
end

