function [u_til_hat, struct_shock] = recover_shocksBraun(data, Btil_draws, alpha_draws, p, k)
% Recover reduced-form residuals and structural shocks from the output of
% Braun and Brüggeman 2022 replication code SVAR_oilmarket function
%
%   Inputs:
%       data           - T x ntil matrix of observed data
%       Btil_draws     - draws x ntil x ntil impact matrix
%       alpha_draws    - draws x (ntil + ntil^2*p) stacked autoregressive coefficients
%       p              - number of lags
%
%   Outputs:
%       u_til_hat      - reduced-form residuals: draws x (T - 4p) x n
%       struct_shock   - structural shocks: draws x (T - 4p) x n

    % Dimensions
    [T, ntil] = size(data);
    n = ntil-k;
    [~, draws] = size(alpha_draws);
    
    struct_shock = NaN(draws, T - p, ntil - k);
    u_til_hat = NaN(draws, T - p, ntil - k);

    for d=1:draws
        
        alpha = alpha_draws(:,d);
        Btil = Btil_draws(:,:,d);

        % Reshape
        Alpha = reshape(alpha, (1+ntil*p), ntil);
        % Btil = reshape(Btil, ntil, ntil);
        B = Btil(1:n,1:n);
    
        % --- Reconstruct nu (constant term) ---
        nu_til = Alpha(1,:);
    
        % --- Reconstruct A matrices ---
        A_til = NaN(ntil, ntil, p);
        rowidx = (1:ntil) + 1;
        for j = 1:p
    
            A_til(:,:,j) = Alpha(rowidx,:);
    
            rowidx = rowidx + ntil;
            
        end
     
    
        % --- Compute fitted values and residuals ---
        y_fitted = NaN(T - p, ntil);
        y = data(1:end, :);
    
        for t = 1:(T - p)
            y_hat_t = nu_til';
    
            for j = 1:p
                y_hat_t = y_hat_t + A_til(:, :, j) * y(p + t - j, :)';
            end
            y_fitted(t, :) = y_hat_t';
        end
    
        u_shock = y(p + 1:end, :) - y_fitted;
    
        % --- Compute structural shocks ---
        sshock = NaN(T - p, ntil - k);
        
        for t = 1:size(sshock, 1)
            eps_t = (B \ u_shock(t, 1:end-k)');
            sshock(t, :) = eps_t';
        end

        struct_shock(d,:,:) = sshock;

        u_til_hat(d,:,:) = u_shock(:,1:n);

    end

end % function
