function [decompositions, reconstruct] = HD_drawbydraw(alpha_draws, B0_draws, p, estT, n, eps_draws)
%% Historical Decomposition from VAR
% INPUTS
% -------------------------------------------------
% alpha_draws, B0_draws : VAR parameters draws
% p                     : lags
% estT                  : number of estimated points
% n                     : number of variables
% eps_draws             : draws x T x n array of structural shocks draws

% OUTPUTS
% -------------------------------------------------
% y_hat_tkj : T x n x n matrix, contribution of shock j to variable k
% y_hat     : T x n matrix, sum across shocks (reconstructed series)

[draws, ~] = size(B0_draws);


% Historical Decomposition

decompositions = zeros(draws, estT, n, n);
reconstruct = zeros(draws, estT, n);

for d = 1:size(alpha_draws,1)
% disp(d)
    alpha = squeeze(alpha_draws(d,:));
    B0 = squeeze(B0_draws(d,:));
    eps = squeeze(eps_draws(d,:,:));

    Theta = squeeze(IRF2(alpha, B0, p, estT-1, n));  % IRF

    y_hat_tkj = zeros(estT, n, n);   % contribution of each shock
    y_hat = zeros(estT, n);          % total reconstructed values
    
    for j = 1:n                      % variable
        for k = 1:n                  % shock
            for t = 1:estT
                effect = zeros(t, 1);
                for i = 0:(t-1)
                    effect(i+1) = Theta(j, k, i+1) * eps(t - i, k);
                end
                y_hat_tkj(t, j, k) = sum(effect);
            end
        end
        y_hat(:, j) = sum(y_hat_tkj(:, j, :), 3);
    end

    % Store
    decompositions(d,:,:,:) = y_hat_tkj;
    reconstruct(d,:,:) = y_hat;
end









end % function