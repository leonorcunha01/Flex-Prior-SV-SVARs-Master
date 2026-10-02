% Main auxiliary simulation for illustration of scaling issue

addpath("functions\utility\")

% Reps number
R = 300;

% Parameters for DGP
T = 400;                        % Sample size
n = 3;                          % Nr. variables
p = 2;                          % Nr. lags
k = n*p+1;                      % Nr. coeffs per equation
np = n*p;                       % Nr. AR coefficients per equation
phi = 0.95;                     % Chosen phi

% Omega
omega_vec = repelem(0.05,T,n);

% --- Inverse impact matrix B0 ---
B0 = [3.0   1.2   -0.8;
      0.9   2.0    1.1;
     -0.6   0.7    2.5];

B0_norm = B0;
for i = 1:n
    B0_norm(i, :) = B0(i, :) / B0(i, i);
end

% --- VAR Coefficients ---
A = zeros(n, n, p);
A(:, :, 1) = [ 0.2,  0.05,  0.00;          % Lag 1 
               -0.1,  0.3,   0.10;
                0.05, -0.05, 0.25 ];

A(:, :, 2) = [ 0.1,  0.2,   0.0;           % Lag 2
               -0.15, 0.3,   0.05;
                0.00, 0.10,  0.20 ];

a = [2; 8; 4];                             % intercepts

% Construct alpha vector (flattened)
alpha = [a(1), reshape(A(1,:,1:2), 1, []), ...
          a(2), reshape(A(2,:,1:2), 1, []), ...
          a(3), reshape(A(3,:,1:2), 1, [])];



% --- Parameters for estimation ---
nsim = 100000;
burnin = 10000;
thin = 100;
var_id = 1:3;
p = 2;
szdrop = 4*p;

quantiles = [0.05, 0.16, 0.50, 0.84, 0.95];

save('simulations\results_aux\mcA03_params.mat', 'D', 'B0', 'B0_norm', 'alpha', 'a', 'A',...
    'omega_vec', 'phi', 'n', 'R', 'szdrop', 'T', 'p');

rng(123, 'combRecursive');

%% --- Preallocate storage ---
effsamp = nsim / thin;
nB = n * n;
nA = numel(alpha);

B0_1_mean = zeros(R, nB);  B0_1_quant = zeros(R, nB, 5);
B0_05_mean = zeros(R, nB); B0_05_quant = zeros(R, nB, 5);
B0_01_mean = zeros(R, nB); B0_01_quant = zeros(R, nB, 5);

B0_1_norm_mean = zeros(R, nB);  B0_1_norm_quant = zeros(R, nB, 5);
B0_05_norm_mean = zeros(R, nB); B0_05_norm_quant = zeros(R, nB, 5);
B0_01_norm_mean = zeros(R, nB); B0_01_norm_quant = zeros(R, nB, 5);

alpha_1_mean = zeros(R, nA);  alpha_1_quant = zeros(R, nA, 5);
alpha_05_mean = zeros(R, nA); alpha_05_quant = zeros(R, nA, 5);
alpha_01_mean = zeros(R, nA); alpha_01_quant = zeros(R, nA, 5);

h_1_mean = zeros(R, T-4*p, n);
h_05_mean = zeros(R, T-4*p, n);
h_01_mean = zeros(R, T-4*p, n);

if isempty(gcp('nocreate'))
    parpool('local', 8);
end

%% --- Monte Carlo loop ---
parfor r = 1:R
    disp(r)

    h = zeros(T, n);
    D = zeros(n, n, T);

    for t = 1:T
        uh = mvnrnd(zeros(1, n), diag(omega_vec(t, :)));
        if t == 1
            h(t, :) = phi .* 0 + uh;
        elseif t == 4*T/5
            h(t, :) = phi .* h(t-1, :) + 2*[1 1 1] + uh;
        else
            h(t, :) = phi .* h(t-1, :) + uh;
        end
        D(:, :, t) = diag(exp(h(t, :)));
    end

    y = dgp_sim(n, T, p, B0, a, A, D);

    HC = 1;
    params = [1, 0.5, 0.1];

    for pp = 1:length(params)
        param_val = params(pp);
        [h_out, Sigh_out, B0_out, alpha_out, Sig_mean_out, Sig_quant_out, ...
         phi_out, tau_out, lambda_out] = func_SVAR7(y, HC, nsim, burnin, thin, p, var_id, param_val);

        B0_mean = squeeze(mean(B0_out));
        B0_quant = quantile(B0_out, quantiles)';

        % Normalize B0
        B0_diag = [repelem(B0_mean(1),effsamp)', repelem(B0_mean(1),effsamp)', repelem(B0_mean(1),effsamp)', ...
                   repelem(B0_mean(5),effsamp)', repelem(B0_mean(5),effsamp)', repelem(B0_mean(5),effsamp)', ...
                   repelem(B0_mean(9),effsamp)', repelem(B0_mean(9),effsamp)', repelem(B0_mean(9),effsamp)'];
        B0_norm = B0_out ./ B0_diag;

        B0_norm_mean = squeeze(mean(B0_norm));
        B0_norm_quant = quantile(B0_norm, quantiles)';

        alpha_mean = squeeze(mean(alpha_out));
        alpha_quant = quantile(alpha_out, quantiles)';

        h_mean = squeeze(mean(h_out));

        switch param_val
            case 1
                B0_1_mean(r,:) = B0_mean(:);
                B0_1_quant(r,:,:) = B0_quant;
                B0_1_norm_mean(r,:) = B0_norm_mean(:);
                B0_1_norm_quant(r,:,:) = B0_norm_quant;
                alpha_1_mean(r,:) = alpha_mean(:);
                alpha_1_quant(r,:,:) = alpha_quant;
                h_1_mean(r,:) = h_mean(:);

            case 0.5
                B0_05_mean(r,:) = B0_mean(:);
                B0_05_quant(r,:,:) = B0_quant;
                B0_05_norm_mean(r,:) = B0_norm_mean(:);
                B0_05_norm_quant(r,:,:) = B0_norm_quant;
                alpha_05_mean(r,:) = alpha_mean(:);
                alpha_05_quant(r,:,:) = alpha_quant;
                h_05_mean(r,:) = h_mean(:);

            case 0.1
                B0_01_mean(r,:) = B0_mean(:);
                B0_01_quant(r,:,:) = B0_quant;
                B0_01_norm_mean(r,:) = B0_norm_mean(:);
                B0_01_norm_quant(r,:,:) = B0_norm_quant;
                alpha_01_mean(r,:) = alpha_mean(:);
                alpha_01_quant(r,:,:) = alpha_quant;
                h_01_mean(r,:) = h_mean(:);
        end
    end
end

%% --- Save results ---
save('simulations\results_aux\mcAUX3_vars.mat', ...
    'B0_1_mean', 'B0_1_quant', 'B0_1_norm_mean', 'B0_1_norm_quant', ...
    'B0_05_mean', 'B0_05_quant', 'B0_05_norm_mean', 'B0_05_norm_quant', ...
    'B0_01_mean', 'B0_01_quant', 'B0_01_norm_mean', 'B0_01_norm_quant', ...
    'alpha_1_mean', 'alpha_1_quant', ...
    'alpha_05_mean', 'alpha_05_quant', ...
    'alpha_01_mean', 'alpha_01_quant', ...
    'h_1_mean', 'h_05_mean', 'h_01_mean');

