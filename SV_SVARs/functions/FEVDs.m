function fevds = FEVDs(valpha,vB,p,H,n,h,phi)
%% Forecast Error Variance Decomposition from VAR
% INPUTS
% alpha, B0 : estimated VAR parameters
% p         : lags
% H         : horizon
% n         : number of variables
% h         : n vector of last estimated h entries
% phi       : n vector of estimated phi parameters

% OUTPUTS
% fevds : n x n x H+1 matrix, contribution of shock j to variable k at t+h


% Get IRF
Theta = squeeze(IRF2(valpha,vB,p,H,n));     % n x n x (H+1)

% Estimate h path until t+H (simplified)
h_path = NaN(n,H+1);
h_path(:,1) = h;
for k = 1:n
    for hor = 2:H+1
        h_path(:,hor) = phi(k) * h_path(:,hor-1);
    end
end

% Get Mean Square Prediction errors
MSPE_jkh = zeros(n,n,H+1);   % cumulative contributions of shock j to var. k
MSPE_jh  = zeros(n,H+1);     % overal variance of var. K

for hor = 0:H
    for j = 1:n
        for k = 1:n
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

fevds = FEVD;


end
