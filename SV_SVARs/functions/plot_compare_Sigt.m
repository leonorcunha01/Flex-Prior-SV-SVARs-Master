function plot_compare_Sigt(Sig_mat1, Sig_mat2)
% Ensure the matrices have the same dimensions
assert(isequal(size(Sig_mat1), size(Sig_mat2)), 'Matrices must have the same dimensions');

[T, n, ~] = size(Sig_mat1);

% Determine the layout of subplots
subplotRows = n;
subplotCols = n;

% Create a figure to hold all the subplots
figure;

% Loop over each entry in the n x n matrix
for i = 1:n
    for j = 1:n
        % Calculate the subplot index
        subplotIndex = (i - 1) * subplotCols + j;

        % Create a subplot for the (i,j) entry
        subplot(subplotRows, subplotCols, subplotIndex);

        % Extract the time series for the (i,j) entry from both matrices
        sigij1 = squeeze(Sig_mat1(:, i, j));
        sigij2 = squeeze(Sig_mat2(:, i, j));

        % Plot the time series for both matrices
        plot(1:T, sigij1, 'b', 'DisplayName', 'Matrix 1');
        hold on;
        plot(1:T, sigij2, 'r', 'DisplayName', 'Matrix 2');
        hold off;

        % Add title and labels
        if i == j
            title(sprintf('Variance of variable %d', i));
        else
            title(sprintf('Covariance between variables %d and %d', i, j));
        end
        xlabel('Time');
        ylabel('Value');

        % Add a legend to distinguish between different matrices
        legend('show');

        % Optionally, you can add grid and adjust other plot properties
        grid on;
    end
end

% Adjust layout to prevent overlap
sgtitle('Comparison of Variance and Covariance Over Time');
end
