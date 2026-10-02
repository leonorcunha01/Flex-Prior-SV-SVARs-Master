function plot_hist_comparison(X1, X2, nRows, nCols)
% PLOT_HIST_COMPARISON plots histograms for each variable (column) in X1 and X2.
%
% Usage:
%   plot_hist_comparison(X1, X2, nRows, nCols)
%
% Inputs:
%   X1, X2 : Arrays of size (R x N) where R = draws, N = variables
%   nRows, nCols : Number of subplot rows and columns
%
% The function creates a grid of subplots (nRows x nCols) and
% plots the histogram of each variable in X1 and X2 on the same axes.
    
    [~, N] = size(X1);

    % Colors
    c1 = [0.4 0.6 0.8];  % soft blue
    c2 = [0.9 0.4 0.4];  % soft red
    
    figure;
    for i = 1:min(N, nRows * nCols)
        subplot(nRows, nCols, i);
        hold on;
        % Kernel density estimates
        [f1, xi1] = ksdensity(X1(:, i));
        [f2, xi2] = ksdensity(X2(:, i));

        % Plot the lines
        plot(xi1, f1, 'Color', c1, 'LineWidth', 1.5);
        plot(xi2, f2, 'Color', c2, 'LineWidth', 1.5);
        legend({sprintf('Var %d (X1)', i), sprintf('Var %d (X2)', i)}, ...
               'Location', 'best');
        hold off;
    end
end
