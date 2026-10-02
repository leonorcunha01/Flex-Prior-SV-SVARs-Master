function plot_array_comp_with_bands(array_list, base_name, xvals, y_name, colors)
    % plot_array_comp_with_bands(..., colors)
    %
    % If colors is supplied, it must be num_arrays x 3.
    % Colors are used for estimate lines and darker/lighter shading.

    % Handle optional inputs
    if nargin < 3 || isempty(xvals)
        xvals = 1:size(array_list{1}, 1);
    end
    if nargin < 4 || isempty(y_name)
        y_name = 'value';
    end

    num_arrays = numel(array_list);
    num_vars   = size(array_list{1}, 2);

    % Handle colors
    if nargin < 5 || isempty(colors)
        colors = lines(num_arrays);
    else
        if size(colors,1) ~= num_arrays || size(colors,2) ~= 3
            error('colors must be an n_arrays x 3 RGB matrix');
        end
    end

    figure;

    for i = 1:num_vars
        subplot(num_vars, 1, i);
        hold on;

        for j = 1:num_arrays
            A = array_list{j};

            if ndims(A) == 2
                % No credible bands
                plot(xvals, A(:, i), ...
                    'Color', colors(j,:), ...
                    'LineWidth', 1.4, ...
                    'DisplayName', ['Array ', num2str(j)]);

            elseif ndims(A) == 3 && size(A, 3) == 3
                % Credible band case
                lower = A(:, i, 1);
                est   = A(:, i, 2);
                upper = A(:, i, 3);

                % Lighter shade version of the line color
                shade_color = 0.7 * colors(j,:) + 0.3 * [1 1 1];

                % Shaded credible band
                fill([xvals, fliplr(xvals)], ...
                     [lower', fliplr(upper')], ...
                     shade_color, ...
                     'EdgeColor', 'none', ...
                     'FaceAlpha', 0.2, ...
                     'HandleVisibility', 'off');

                % Estimate line
                plot(xvals, est, ...
                    'Color', colors(j,:), ...
                    'LineWidth', 1.6, ...
                    'DisplayName', ['Array ', num2str(j)]);
            else
                error('Array must be 2D or 3D with size 3 in the third dimension.');
            end
        end

        title([base_name, ' ', num2str(i)]);
        xlabel('time');
        ylabel(y_name);
        legend show;
        grid on;
        hold off;
    end
end
