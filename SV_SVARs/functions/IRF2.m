function [irfs] = IRF2(alpha_draws,B0_draws,p,h,n)
%% Impulse Response Functions from VAR
% INPUTS
% alpha_draws, B0_draws : VAR parameters draws
% p                     : lags
% n                     : number of variables
% h                     : horizon

% OUTPUTS
% irfs : n x n x H+1 matrix, effect of shock j to variable k at t+h



% Compute IRFs
irfs = NaN(size(alpha_draws,1), n, n,h+1);

for d = 1:size(alpha_draws,1)
    
    valpha = alpha_draws(d,:);
    vB = B0_draws(d,:);

    % Reshape
    Phi = reshape(valpha,p*n+1,n);
    PhiAR = Phi(2:end,:)';
    B =reshape( vB,n,n);
    
    
    % IRFs
    [AA, J] = Companion( PhiAR , 0);
    
    Apoweri = AA;
    Phis = zeros(n,n,h+1);
    Phis(:,:,1) = eye(n);
    for ii = 1:h
        Phis(:,:,ii+1) = J*Apoweri*J';
        Apoweri = Apoweri * AA;
    end
    
    Theta = zeros(n,n,h+1);
    for ii = 0:h
        Theta_ii = Phis(:,:,ii+1)/B;
        Theta(:,:,ii+1) = Theta_ii(1:n,1:n);
    end

    irfs(d,:,:,:) = Theta;
   
end

end

