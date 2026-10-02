function plot_hist_comparison3(X1, X2, X3, nRows, nCols, plotTitles, varNames)
% PLOT_HIST_COMPARISON3 plots smooth density lines for each variable (column)
% in X1, X2, and X3, with a single legend on the right.
%
% Usage:
%   plot_hist_comparison3(X1, X2, X3, nRows, nCols, plotTitles, varNames)
%
% Inputs:
%   X1, X2, X3 : Arrays of size (R x N), where R = draws and N = variables
%   nRows, nCols : Number of subplot rows and columns
%   plotTitles : Cell array of strings, one for each variable/column (subplot titles)
%   varNames : Cell array of 3 strings, names of the datasets for the legend

    [~, N] = size(X1);

    % Colors 
    c1 = [0.4 0.6 0.8]; % soft blue 
    c2 = [0.9 0.4 0.4]; % soft red 
    c3 = [0.5 0.8 0.5]; % soft green    
    
    figure;
    tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Predefine handles for legend
    h1 = []; h2 = []; h3 = [];

    for i = 1:min(N, nRows * nCols)
        nexttile;
        hold on;

        % Kernel density estimates
        [f1, xi1] = ksdensity(X1(:, i));
        [f2, xi2] = ksdensity(X2(:, i));
        [f3, xi3] = ksdensity(X3(:, i));

        % Plot lines
        h1 = plot(xi1, f1, 'Color', c1, 'LineWidth', 1.5);
        h2 = plot(xi2, f2, 'Color', c2, 'LineWidth', 1.5);
        h3 = plot(xi3, f3, 'Color', c3, 'LineWidth', 1.5);

        % Set subplot title
        title(plotTitles{i});
        

        hold off;
    end

    % Add one legend for the entire figure
    lgd = legend([h1, h2, h3], varNames, ...
        'Location', 'best');
    % lgd.Layout.Tile = 'east';
end

