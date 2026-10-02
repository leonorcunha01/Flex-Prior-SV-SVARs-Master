% This function calculates the variance and covariance for figures
function plot_Sigt(Sig_mat)

[T, n, ~] = size(Sig_mat);


% Create a figure for the combined plot
figure;

% Hold on to plot multiple lines on the same graph
hold on;

% Loop over each entry in the n x n matrix
for i = 1:n
    for j = 1:n
        % Extract the time series for the (i,j) entry
        sigij = squeeze(Sig_mat(:, i, j));

        % Plot the time series with a unique identifier for each line
        plot(1:T, sigij, 'DisplayName', sprintf('sig_(%d,%d)', i, j));
    end
end

% Add title and labels
title('Variance and Covariance Over Time, Combined Plot');
xlabel('Time');
ylabel('Value');

% Add a legend to distinguish between different entries
legend('show');

% Optionally, you can add grid and adjust other plot properties
grid on;

% Release the hold
hold off;


end