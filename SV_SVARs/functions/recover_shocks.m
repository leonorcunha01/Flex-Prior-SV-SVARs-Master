function [u_hat, struct_shock] = recover_shocks(data, B0_draws, alpha_draws, p)
% Recover reduced-form residuals and structural shocks
%
%   Inputs:
%       data   - (T x n) matrix of observed data
%       B0     - (n^2) impact matrix
%       alpha  - ((n + n^2*p)) stacked autoregressive coefficients
%       p      - number of lags
%
%   Outputs:
%       u_hat        - reduced-form residuals (T - 4p) x n
%       struct_shock - structural shocks (T - 4p) x n

    % Dimensions
    [T, n] = size(data);
    [draws, ~] = size(B0_draws);


    y = data(4*p + 1 - p:end, :);
    
    struct_shock = NaN(draws,T-4*p,n);
    u_hat = NaN(draws,T-4*p,n);    

    for d = 1:draws
        % Reshape
        alpha = squeeze(alpha_draws(d,:));
        B0 = reshape(B0_draws(d,:), n, n)';
    
        % --- Reconstruct a (constant term) ---
        a_mean = zeros(n, 1);
        ai = 1;
        for i = 1:n
            a_mean(i) = alpha(ai);
            ai = ai + p*n + 1;
        end
    
        % --- Reconstruct A matrices ---
        A_mean = zeros(n, n, p);
        Ai = 2; % start after first a element
        for i = 1:n
            for j = 1:p
                A_mean(i,:,j) = alpha(Ai : Ai + n - 1);
                Ai = Ai + n;
            end
            Ai = Ai + 1; % skip next a
        end
    
        % --- Compute fitted values and residuals ---
        y_fitted = zeros(T - 4*p, n);
          
        for t = 1:(T - 4*p)
            y_hat_t = a_mean;
            for j = 1:p
                y_hat_t = y_hat_t + A_mean(:, :, j) * y(p + t - j, :)';
            end
            y_fitted(t, :) = y_hat_t';
        end
    
        ushock = y(p + 1:end, :) - y_fitted;
    
        % --- Compute structural shocks ---
        sshock = NaN(size(ushock));
        for t = 1:size(ushock, 1)
            sshock(t, :) = (B0 * ushock(t, :)')';
        end

        struct_shock(d,:,:) = sshock;

        u_hat(d,:,:) = ushock;



    end


end % function
