function fevds = FEVDs_drawbydraw(alpha_draws,B0_draws,p,H,n,h_draws,phi_draws, omega_draws)
%% Forecast Error Variance Decomposition from VAR
% INPUTS
% alpha, B0 : VAR parameters draws
% p         : lags
% H         : horizon
% n         : number of variables
% h         : n vector of last h entries draws
% phi       : n vector of phi draws
% omega

% OUTPUTS
% fevds : n x n x H+1 matrix, contribution of shock j to variable k at t+h

[draws, ~] = size(B0_draws);

fevds = zeros(draws,n,n,H+1);
for d=1:draws

    valpha = squeeze(alpha_draws(d,:));
    vB = squeeze(B0_draws(d,:));
    h = squeeze(h_draws(d,:,:));
    phi = squeeze(phi_draws(d,:));
    omega = squeeze(omega_draws(d,:));

    % Get IRF
    Theta = squeeze(IRF2(valpha,vB,p,H,n));     % n x n x (H+1)
    
    % Estimate h path until t+H (simplified)
    h_path = NaN(n,H+1);
    h_path(:,1) = h;
    for k = 1:n
        for hor = 2:H+1
            h_path(k,hor) = phi(k) * h_path(k,hor-1) + sqrt(omega(k)).*randn(1);
        end
    end
    
    % Get Mean Square Prediction errors
    MSPE_jkh = zeros(n,n,H+1);   % cumulative contributions of shock j to var. k
    MSPE_jh  = zeros(n,H+1);     % overal variance of var. K
    
    for hor = 0:H
        for j = 1:n         % variables
            for k = 1:n         % shocks
                MSPE_elem = zeros(hor+1,1);
                for i = 0:hor
                    theta_jki = Theta(j,k,i+1);
                    eh_khi = exp(h_path(k,hor+1-i));
                    MSPE_elem(i+1) = theta_jki^2 * eh_khi;
                end
                MSPE_jkh(j,k,hor+1) = sum(MSPE_elem);
            end
            MSPE_jh(j,hor+1) = sum(MSPE_jkh(j,:,hor+1));
        end
    end
    
    % Forecast error decompositions
    FEVD = zeros(n,n,H+1);
    for hor = 0:H
        for j = 1:n
            FEVD(j,:,hor+1) = MSPE_jkh(j,:,hor+1) ./ MSPE_jh(j,hor+1);
        end
    end

    fevds(d,:,:,:) = FEVD;


end

end % function
