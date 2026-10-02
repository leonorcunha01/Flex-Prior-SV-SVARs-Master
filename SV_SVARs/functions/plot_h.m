function plot_h(array, xvals, ylab)
    % plot_h(array, xvals, ylab)
    % Plots each column of 'array' against xvals in separate subplots

    % Handle optional inputs
    if nargin < 2 || isempty(xvals)
        xvals = 1:size(array, 1);
    end
    if nargin < 3
        ylab = 'value';
    end

    % Determine number of variables (columns)
    num_vars = size(array, 2);

    % Create a new figure
    figure;

    % Loop over each variable to create a subplot
    for i = 1:num_vars
        subplot(num_vars, 1, i);
        plot(xvals, array(:, i), 'b-', 'LineWidth', 1.2);

        % Title and labels
        title(['h', num2str(i)]);
        ylabel(ylab);
        grid on;

        % % Format x-axis if values are datetime
        % if isdatetime(xvals)
        %     datetick('x', 'keeplimits');
        % end
    end

    % Label x-axis only on the last subplot (optional)
    xlabel('x');
end
