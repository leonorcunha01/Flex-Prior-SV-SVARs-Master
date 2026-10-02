function plot_fevd(fevds)

[n, ~, horizon] = size(fevds);
t = 1:horizon;                  % time vector

figure;
% Plot FEVD over time for each variable
for var = 1:n
    % Extract contributions of all shocks to this variable over time
    FEVD_var = zeros(horizon, n);
    for shock = 1:n
        FEVD_var(:,shock) = fevds(var,shock,:);
    end
    
    subplot(n,1,var);
    bar(t-1, FEVD_var, 'stacked');  % stacked bars
    xlabel('Horizon');
    ylabel('FEVD');
    title(['Variable ', num2str(var)]);
    if var == 1
        legend(arrayfun(@(x) ['Shock ', num2str(x)], 1:n, 'UniformOutput', false), 'Location', 'northoutside');
    end
    grid on;
end

end


