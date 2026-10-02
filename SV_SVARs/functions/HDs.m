function [y_hat_tkj, y_hat] = HDs(alpha, B0, p, estT, n, eps)
%% Historical Decomposition from VAR
% INPUTS
% -------------------------------------------------
% alpha, B0 : estimated VAR parameters
% p         : lags
% estT      : number of estimated points
% n         : number of variables
% eps       : T x n matrix of structural shocks

% OUTPUTS
% -------------------------------------------------
% y_hat_tkj : T x n x n matrix, contribution of shock j to variable k
% y_hat     : T x n matrix, sum across shocks (reconstructed series)


%% Compute IRFs
irfR = zeros(n, n,estT-1, size(alpha,1));
for i = 1:size(alpha,1)
    irfR(:,:,:,i) = IRF2(alpha(i,:), B0(i,:), p, estT-2, n);
   
end

    % Quantiles and mean
    Theta = squeeze(median(irfR, 4));


%% Historical Decomposition
y_hat_tkj = zeros(estT, n, n);   % contribution of each shock
y_hat = zeros(estT, n);          % total reconstructed values

for j = 1:n                      % variable
    for k = 1:n                  % shock
        for t = 1:estT
            effect = zeros(t-1, 1);
            for i = 1:t-1
                effect(i) = Theta(j, k, i) * eps(t - i+1, k);
            end
            y_hat_tkj(t, j, k) = sum(effect);
        end
    end
    y_hat(:, j) = sum(y_hat_tkj(:, j, :), 3);
end


%% Plot Historical Decomposition for variable k
% var_names = {'Y1', 'Y2', 'Y3', 'Y4'}; 
% shock_names = {'Shock1', 'Shock2', 'Shock3', 'Shock4'};  
% 
% 
% for j = 1:n
%     figure;
% 
%     % Extract contributions of each shock to variable j
%     contrib = squeeze(y_hat_tkj(:, j, :));  % T x n matrix
% 
%     % Plot stacked bars
%     bar(time, contrib, 'stacked');
%     hold on;
% 
%     % Add the reconstructed total as a black line
%     plot(time, y_hat(:, j), 'k-', 'LineWidth', 1.5);
% 
%     % Formatting
%     title(['Historical Decomposition: ', var_names{j}]);
%     xlabel('Time');
%     ylabel(var_names{j});
%     legend([shock_names, {'Total'}], 'Location', 'bestoutside');
%     grid on;
% 
%     % Aesthetics
%     colormap(parula(n)); 
%     set(gca, 'FontSize', 12);
% end
% 
% 
% 


end % function