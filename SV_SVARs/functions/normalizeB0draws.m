function B0_norm = normalizeB0draws(B0)
%Normalize B0 draws by its diagonal means


    [ndraws, total] = size(B0);
    n = sqrt(total);
    
    % Mean across draws
    B0_mean = mean(B0, 1);

    % Indices of diagonal elements
    diag_idx = 1:(n+1):n^2;

    % Build diagonal normalization vector (each diag element repeated n times)
    diag_vals = B0_mean(diag_idx);
    B0_diag = repelem(diag_vals, n);

    % Repeat for all draws
    B0_diag = repmat(B0_diag, ndraws, 1);

    % Normalize
    B0_norm = B0 ./ B0_diag;
end
