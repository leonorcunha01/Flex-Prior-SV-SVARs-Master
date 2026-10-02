function plot_HD(decomp, time)
%% Plot median or mean Historical Decomposition for variable k
% INPUTS
% decomp : 
% reconstr         : 

[~, n, ~] = size(decomp);


var_names = {'Y1', 'Y2', 'Y3', 'Y4'}; 
shock_names = {'Shock1', 'Shock2', 'Shock3', 'Shock4'};  


for j = 1:n
    figure;
    
    % Extract contributions of each shock to variable j
    contrib = squeeze(decomp(:, j, :));  % T x n matrix
    
    % Plot stacked bars
    bar(time, contrib, 'stacked');
    hold on;
    
    % Add the reconstructed total as a black line
    plot(time, sum(decomp(:, j, :),3), 'k-', 'LineWidth', 1.5);
    
    % Formatting
    title(['Historical Decomposition: ', var_names{j}]);
    xlabel('Time');
    ylabel(var_names{j});
    legend([shock_names, {'Total'}], 'Location', 'bestoutside');
    grid on;
    
    % Aesthetics
    colormap(parula(n)); 
    set(gca, 'FontSize', 12);
end

