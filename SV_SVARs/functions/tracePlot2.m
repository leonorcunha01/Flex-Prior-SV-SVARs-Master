function tracePlot2(parameterEstimates, parameterTrue, rows, cols)
    % Input:
    % parameterEstimates: An array of any dimension where the first dimension
    % represents the iterations of the sampler.

    % Create a new figure
    figure;

    

    % Get the size of the input array
    sizeArray = size(parameterEstimates);
    its = sizeArray(1);

    % Calculate the total number of parameters to plot
    numParameters = prod(sizeArray(2:end));

    
    % Reshape the array into a 2D array where each column represents a parameter
    parameterEstimates2D = reshape(parameterEstimates, sizeArray(1), numParameters);

    % % Ensure the matrices have the same dimensions
    % assert(isequal(size(parameterEstimates2D), size(parameterTrue)), 'Vectors must have the same dimensions');

    % Determine the layout of subplots
    subplotRows = rows;
    subplotCols = cols;

% Loop over each entry in the n x n matrix
for i = 1:subplotRows
    for j = 1:subplotCols
        % Calculate the subplot index
        Index = (i - 1) * subplotCols + j;

        % Create a subplot for the (i,j) entry
        subplot(subplotRows, subplotCols, Index);

        % Extract the time series for the Index for both
        estim = parameterEstimates2D(:, Index);
        true = parameterTrue(Index);

        % Plot the time series for both matrices
        plot(1:its, estim, 'r', 'DisplayName', 'Draw');
        hold on;
        plot(1:its, repelem(true,its), 'b', 'DisplayName', 'True');
        hold off;

        % Add title and labels
        
        title(sprintf('Entry %d and %d', i, j));
        
        xlabel('Iteration');
        ylabel('Value');

        

        % Optionally, you can add grid and adjust other plot properties
        grid on;
    end
end

% Adjust layout to prevent overlap
sgtitle('Trace Plot');
% Add a legend to distinguish between different matrices
        legend('show');
end
