%% MC study 6
% random walk h

addpath("functions\utility\")

% Reps number
R = 500;

% Parameters for DGP
T = 400;                        % Sample size
n=2;                            % Nr. variables
p = 2;                          % Nr. lags
k = n*p+1;                      % Nr. coeffs per equation
np = n*p;                       % Nr. AR coefficients per equation
phi = 1;                        % for a RW

% Omega
omega_vec = repelem(0.05,T,n);

% (Inverse) impact matrix
B0 = [3 1.5;                            % Actual B0
      -0.5 2];

B0_norm = B0;                           % Normalized version
for i=1:n
    B0_norm(i,:) = B0(i,:)/B0(i,i);
end    

% VAR coefficients
A = zeros(n, n, p);                                           % Coeff matrices
A(:,:,1) =[0.2, 0.05;                                         % Lag 1
        -0.1, 0.3 ];
A(:,:,2) =[.1 0.2;                                            % Lag 2
          -.15 .3 ];
a = [2 8]';                                                   % intercepts
alpha = [a(1), A(1,:,1), A(1,:,2), a(2), A(2,:,1), A(2,:,2)]; % join them in shape of estimates




%% Parameters for model
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;

% Size dropped observations
szdrop = 4*p;

% Define the quantiles to compute, including the median (0.50)
quantiles = [0.05, 0.16, 0.50, 0.84, 0.95];


% Save parameters
save('simulations\results\mc06_params.mat', 'B0', 'B0_norm', 'alpha', 'a', 'A', ...
    'omega_vec', 'phi', 'n', 'R', 'szdrop', 'T', 'p');
 


%% Preallocation

B0_hs_mean = zeros(R, n*n); B0_hs_quant = zeros(R, n*n, 5);
alpha_hs_mean = zeros(R, k*n); alpha_hs_quant = zeros(R, k*n, 5);
Sig_hs_mean = zeros(R, T-szdrop, n, n); Sig_hs_quant = zeros(R, 5, T-szdrop, n, n); 
h_hs_mean = zeros(R,T-szdrop,n); h_hs_quant = zeros(R,T-szdrop,n,5);

B0_hs_norm_mean = zeros(R, n*n); B0_hs_norm_quant = zeros(R, n*n, 5);

B0_ig_mean = zeros(R, n*n); B0_ig_quant = zeros(R, n*n,5);
alpha_ig_mean = zeros(R, k*n); alpha_ig_quant = zeros(R, k*n,5);
Sig_ig_mean = zeros(R, T-szdrop, n, n); Sig_ig_quant = zeros(R, 5, T-szdrop, n, n); 
h_ig_mean = zeros(R,T-szdrop,n); h_ig_quant = zeros(R,T-szdrop,n,5);

B0_ig_norm_mean = zeros(R, n*n); B0_ig_norm_quant = zeros(R, n*n, 5);

h_real = zeros(R,T-szdrop,n);

effsamp = nsim/thin;


%% Start parallel pool (if not already running)
if isempty(gcp('nocreate'))
    parpool('local', 8);
end

%% Parallel Simulation
parfor r = 1:R
    disp(r)
    rng(123 + r, 'combRecursive'); % reproducible random streams

    % --- Generate h and D ---
    valid = false;  % flag to control redraw
    while ~valid
        h = zeros(T, n);
        D = zeros(n, n, T);    
        for t = 1:T
            uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
            if t == 1
                h(t,:) = phi .* 0 + uh;
            else
                h(t,:) = phi .* h(t-1,:) + uh;
            end
            D(:,:,t) = diag(exp(h(t,:)));
        end
          
        if (all(h(:) <= 6) && all(h(:) >= -6))                                   % check condition
            valid = true;
        end
    end

    % --- Generate data ---
    y = dgp_sim(n,T,p, B0, a, A, D);
    h_real(r,:,:) = h(szdrop+1:end,:);

    % --- Estimate models ---
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

    

    % --- Save estimates (mean and quantiles): HS ---
    B0_hs_mean(r,:) = squeeze(mean(B0_hs));
    alpha_hs_mean(r,:) = squeeze(mean(alpha_hs));
    h_hs_mean(r,:,:) = squeeze(mean(h_hs));
    Sig_hs_mean(r,:,:,:) = Sig_mean_hs;

    B0_hs_quant(r,:,:) = quantile(B0_hs, quantiles)';
    alpha_hs_quant(r,:,:) = quantile(alpha_hs, quantiles)';

    tmp_h_hs_quant = zeros(T-szdrop, n, 5);
    for j = 1:n
        tmp_h_hs_quant(:,j,:) = quantile(squeeze(h_hs(:,:,j)), quantiles)';
    end
    h_hs_quant(r,:,:,:) = tmp_h_hs_quant;
    
    Sig_hs_quant(r,:,:,:,:) = Sig_quant_hs;

    % --- Save estimates (mean and quantiles): IG ---
    B0_ig_mean(r,:) = squeeze(mean(B0_ig));
    alpha_ig_mean(r,:) = squeeze(mean(alpha_ig));
    h_ig_mean(r,:,:) = squeeze(mean(h_ig));
    Sig_ig_mean(r,:,:,:) = Sig_mean_ig;

    B0_ig_quant(r,:,:) = squeeze(quantile(B0_ig, quantiles))';
    alpha_ig_quant(r,:,:) = squeeze(quantile(alpha_ig, quantiles))';

    tmp_h_ig_quant = zeros(T-szdrop, n, 5);
    for j = 1:n
        tmp_h_ig_quant(:,j,:) = quantile(squeeze(h_ig(:,:,j)), quantiles)';
    end
    h_ig_quant(r,:,:,:) = tmp_h_ig_quant;

    Sig_ig_quant(r,:,:,:,:) = Sig_quant_ig;

    % --- Normalize B0 ---
    B0_hs_norm = normalizeB0draws(B0_hs);
    B0_ig_norm = normalizeB0draws(B0_ig);

    B0_hs_norm_mean(r,:) = squeeze(mean(B0_hs_norm));
    B0_hs_norm_quant(r,:,:) = quantile(B0_hs_norm, quantiles)';

    B0_ig_norm_mean(r,:) = squeeze(mean(B0_ig_norm));
    B0_ig_norm_quant(r,:,:) = quantile(B0_ig_norm, quantiles)';
end


%% Save all results
save('simulations\results\mc06_vars.mat', 'B0_hs_mean', 'B0_hs_quant', 'alpha_hs_mean', 'alpha_hs_quant', ...
     'Sig_hs_mean', 'Sig_hs_quant', 'h_hs_mean', 'h_hs_quant',  ...
     'B0_ig_mean', 'B0_ig_quant', 'alpha_ig_mean', 'alpha_ig_quant', ...
     'Sig_ig_mean', 'Sig_ig_quant', 'h_ig_mean', 'h_ig_quant', ...
     'h_real', ...
     'B0_hs_norm_mean', 'B0_hs_norm_quant', 'B0_ig_norm_mean', 'B0_ig_norm_quant');
