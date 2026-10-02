%% Main file
% ======================================================================= %
% This is the main script for the Master Thesis 'Flexible Prior
% Distributions for Stochastic Volatility in Bayesian Structural Vector 
% Autoregressions'.
% 
% The file calls the scripts which run the main (6) and auxiliary
% simulations or loads the results, in the case the file exists already and
% the simulations do not need running again. The change is done by
% commenting in/out the script names.
% 
% In the empirical application part, the file calls the script which builds
% the sample, and the ones which estimate the HS and IG models for the 
% new data sample, estimate the Braun & Bruggemann model with the chosen
% restrictions, compute the correlations for labeling the shocks and create
% the VAR tools for the analysis.
% These steps, except the last, are repeated for the 'old' (BB22) sample.
% 
% It further creates and saves the relevant figures and tables.
% 
% ======================================================================= %


%% ====================================================================== %
%  Settings
% ======================================================================= %

clear;
addpath("functions\");

%% ====================================================================== %
%  Simulations
% ======================================================================= %
addpath("simulations\");
addpath("simulations\results\");
addpath("simulations\results_aux\");
addpath("figures\")
addpath("tables\")

% For plotting
variances = {'v_b = 1', 'v_b = 0.5', 'v_b = 0.1'};
entriesB0_n2 = {'(1,1)', '(1,2)', ...
                '(2,1)', '(2,2)'};
entriesB0_n3 = {'(1,1)', '(1,2)', '(1,3)', ...
                '(2,1)', '(2,2)', '(2,3)', ...
                '(3,1)', '(3,2)', '(3,3)'};
entriesAlpha_n2 = {'Int. 1', '(1,1) - Lag 1', '(1,2) - Lag 1', '(1,1) - Lag 2', '(1,2) - Lag 2', ...
                   'Int. 2', '(2,1) - Lag 1', '(2,2) - Lag 1', '(2,1) - Lag 2', '(2,2) - Lag 2'};
entriesAlpha_n3 = {'Int. 1', '(1,1) - Lag 1', '(1,2) - Lag 1', '(1,3) - Lag 1', '(1,1) - Lag 2', '(1,2) - Lag 2', '(1,3) - Lag 2', ...
                   'Int. 2', '(2,1) - Lag 1', '(2,2) - Lag 1', '(2,3) - Lag 1', '(2,1) - Lag 2', '(2,2) - Lag 2', '(2,3) - Lag 2', ...
                   'Int. 3', '(3,1) - Lag 1', '(3,2) - Lag 1', '(3,3) - Lag 1', '(3,1) - Lag 2', '(3,2) - Lag 2', '(3,3) - Lag 2'};


%% Auxiliary simulations, for scaling issue
% MC_aux2_study and MC_aux3_study are the final used to illustrate and test
% the scaling issue
% ======================================================================= %
%                       Auxiliary simulation 2
% This is the alternative one

% MC_aux2_study;                                 % run simulation if needed

load('simulations\results_aux\mcA02_params_var.mat')           % Load parameters
load('simulations\results_aux\mcAUX2_vars.mat')   % Load results

% --- Example h and y ---
rng(102,"twister")

D = zeros(n,n,T);
h = zeros(T,n);
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

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

nsim = 100000; burnin = 10000; thin = 100;
var_id = [1 2 3]; p = 2;
HC = 1;
params = [1, 0.5, 0.1];
h_est = NaN(3, T- 4*p, n); B0_est = NaN(3, n^2); 

for pp = 1:length(params)
    param_val = params(pp);
    [h_out, ~, B0_out, ~] = func_SVAR7(y, HC, nsim, burnin, thin, p, var_id, param_val);
    h_est(pp,:,:) = median(h_out);
    B0_est(pp,:) = median(B0_out);
end

plot_array_comp({h(4*p+1:end,:), squeeze(h_est(1,:,:)), squeeze(h_est(2,:,:)), squeeze(h_est(3,:,:))}, 'h_', 4*p+1:T, 'h')
    savefig('figures\MCaux2\h_comp_examp_aux2.fig');

B0_med = array2table(B0_est', "VariableNames", {'v_b = 1', 'v_b = 0.5', 'v_b = 0.1'}, "RowNames", entriesB0_n3);
    save('tables\B0_examp_aux2.mat', 'B0_med');


% Distributions of estimates across simulations
plot_hist_comparison3(B0_1_quant(:,:,3), B0_05_quant(:,:,3), B0_01_quant(:,:,3), 3, 3, ...                  % B0 estimated
    entriesB0_n3, variances)                
        savefig('figures/MCaux2/MCaux2_hist_B0.fig'); % save
plot_hist_comparison3(B0_1_norm_quant(:,:,3), B0_05_norm_quant(:,:,3), B0_01_norm_quant(:,:,3), 3, 3, ...   % B0 normalized
    entriesB0_n3, variances) 
        savefig('figures/MCaux2/MCaux2_hist_B0norm.fig'); %save
plot_hist_comparison3(alpha_1_quant(:,:,3), alpha_05_quant(:,:,3), alpha_01_quant(:,:,3), 3, 7, ...         % alpha estimated
    entriesAlpha_n3, variances)  
        savefig('figures/MCaux2/MCaux2_hist_alpha.fig'); %save


% ======================================================================= %
%                       Auxiliary simulation 3
% This is the main one

% MC_aux3_study;                                 % run simulation if needed

load('simulations\results_aux\mcA03_params.mat')           % Load parameters
load('simulations\results_aux\mcAUX3_vars.mat')   % Load results

% --- Example h and y generation ---
rng(103,"twister")

D = zeros(n,n,T);
h = zeros(T,n);
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

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

nsim = 100000; burnin = 10000; thin = 100;
var_id = 1:3; p = 2;
HC = 1;
params = [1, 0.5, 0.1];
h_est = NaN(3, T- 4*p, n); B0_est = NaN(3, n^2);

for pp = 1:length(params)
    param_val = params(pp);
    [h_out, ~, B0_out, ~] = func_SVAR7(y, HC, nsim, burnin, thin, p, var_id, param_val);
    h_est(pp,:,:) = median(h_out);
    B0_est(pp,:) = median(B0_out);
end

plot_array_comp({h(4*p+1:end,:), squeeze(h_est(1,:,:)), squeeze(h_est(2,:,:)), squeeze(h_est(3,:,:))}, 'h_', 4*p+1:T, 'h')
    savefig('figures\MCaux3\h_comp_examp_aux3.fig');

B0_med = array2table(B0_est', "VariableNames", {'v_b = 1', 'v_b = 0.5', 'v_b = 0.1'}, "RowNames", entriesB0_n3);
    save('tables\B0_examp_aux3.mat', 'B0_med');


% Distributions of estimates across simulations
plot_hist_comparison3(B0_1_quant(:,:,3), B0_05_quant(:,:,3), B0_01_quant(:,:,3), 3, 3, ...                  % B0 estimated
    entriesB0_n3, variances)                
        savefig('figures/MCaux3/MCaux3_hist_B0.fig'); % save
plot_hist_comparison3(B0_1_norm_quant(:,:,3), B0_05_norm_quant(:,:,3), B0_01_norm_quant(:,:,3), 3, 3, ...   % B0 normalized
    entriesB0_n3, variances) 
        savefig('figures/MCaux3/MCaux3_hist_B0norm.fig'); %save
plot_hist_comparison3(alpha_1_quant(:,:,3), alpha_05_quant(:,:,3), alpha_01_quant(:,:,3), 3, 7, ...         % alpha estimated
    entriesAlpha_n3, variances)  
        savefig('figures/MCaux3/MCaux3_hist_alpha.fig'); %save



%% Main simulations, for comparing HS and IG versions
% MC_01_study, MC_02_study, MC_03_study, MC_04_study, MC_05_study, MC_06_study
clear;
addpath("simulations\");
addpath("figures\")
addpath("tables\")



% ======================================================================= %
%                       Simulation 1
clear;

% MC_01_study;                                  % run simulation if needed
load('simulations\results\mc01_params.mat')             % Load parameters
load('simulations\results\mc01_vars.mat')               % Load results

DGP = 1;                                        % Name DGP for storing plots and tables

% Example run
% h process, with a jump of 1 unit at 4/5 of the sample
rng(31, 'twister')

D = zeros(n,n,T);
h = zeros(T,n);
for t=1:T
        uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
        if t == 1
            h(t,:) = phi .* 0 + uh;
        else
            h(t,:) = phi .* h(t-1,:) + uh;
        end

        D(:,:,t) = diag((exp(h(t,:))));
end

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_hs' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_ig' '_DGP' num2str(DGP) '.fig']));


% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));

% ======================================================================= %
%                       Simulation 2
clear;

% MC_02_study;                                  % run simulation if needed
load('simulations\results\mc02_params.mat')             % Load parameters
load('simulations\results\mc02_vars.mat')               % Load results

DGP = 2;                                        % Name DGP for storing plots and tables


% Example h and y processes
% h process, with a jump of 2 units at 4/5 of the sample
rng(32, 'twister')

D = zeros(n,n,T);
h = zeros(T,n);
for t=1:T
        uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
        if t == 1
            h(t,:) = phi .* 0 + uh;
        elseif t == 4*T/5
            h(t,:) = phi .* h(t-1,:) + [2 2] + uh;
        else
            h(t,:) = phi .* h(t-1,:) + uh;
        end
        D(:,:,t) = diag((exp(h(t,:))));
end

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_hs' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_ig' '_DGP' num2str(DGP) '.fig']));



% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));


% ======================================================================= %
%                       Simulation 3
clear;

% MC_03_study;                                  % run simulation if needed
load('simulations\results\mc03_params.mat')             % Load parameters
load('simulations\results\mc03_vars.mat')               % Load results

DGP = 3;                                        % Name DGP for storing plots and tables


% Example h and y processes
% h process, with extreme value of h at 4/5 of the sample
rng(33, 'twister')

D = zeros(n,n,T);
h = zeros(T,n);
for t=1:T
        uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
        if t == 1
            h(t,:) = phi .* 0 + uh;
        elseif t == 4*T/5
            h(t,:) = [4 4];
        else
            h(t,:) = phi .* h(t-1,:) + uh;
        end
        D(:,:,t) = diag((exp(h(t,:))));
end

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_hs' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_ig' '_DGP' num2str(DGP) '.fig']));


% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

 plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));


% ======================================================================= %
%                       Simulation 4
clear;

% MC_04_study;                                  % run simulation if needed
load('simulations\results\mc04_params.mat')             % Load parameters
load('simulations\results\mc04_vars.mat')               % Load results

DGP = 4;                                        % Name DGP for storing plots and tables


% Example h and y processes
% h process, with 2 extreme points around 4/5 of sample, h1 first
rng(34, 'twister')

D = zeros(n,n,T);
h = zeros(T,n);
for t=1:T
        uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
        if t == 1
            h(t,:) = phi .* 0 + uh;
        elseif t == 4*T/5
            h(t,:) = [4 0];
        elseif t == 4*T/5+4
            h(t,:) = [0 4];
        else
            h(t,:) = phi .* h(t-1,:) + uh;
        end
        D(:,:,t) = diag((exp(h(t,:))));
end

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_hs' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_ig' '_DGP' num2str(DGP) '.fig']));


% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

 plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));

% ======================================================================= %
%                       Simulation 5
clear;

 % MC_05_study;                                  % run simulation if needed
load('simulations\results\mc05_params.mat')             % Load parameters
load('simulations\results\mc05_vars.mat')               % Load results

DGP = 5;                                        % Name DGP for storing plots and tables

% Example h and y processes
% h process, with a jump of 1 unit at 4/5 of the sample
rng(35, 'twister')

D = zeros(n,n,T);
h = zeros(T,n);
for t=1:T
        uh = mvnrnd(zeros(1,n), diag(omega_vec(t,:)));
        if t == 1
            h(t,:) = phi .* 0 + uh;
        elseif t == 4*T/5
            h(t,:) = [0 4];
        elseif t == 4*T/5+4
            h(t,:) = [4 0];
        else
            h(t,:) = phi .* h(t-1,:) + uh;
        end
        D(:,:,t) = diag((exp(h(t,:))));
end

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);

% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_hs' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0_ig' '_DGP' num2str(DGP) '.fig']));

% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));

% ======================================================================= %
%                       Simulation 6
clear;

 % MC_06_study;                                  % run simulation if needed
load('simulations\results\mc06_params.mat')             % Load parameters
load('simulations\results\mc06_vars.mat')               % Load results

DGP = 6;                                        % Name DGP for storing plots and tables

% Example h and y processes
% RW h process, which does not go above 6 or below -6
rng(36, 'twister')

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

% y process
[y] = dgp_sim(n,T,p, B0, a, A, D);


% Estimate models 
nsim = 100000; 
burnin = 10000;
thin = 100;
var_id = [1 2];
p = 2;
    HC = 1;
    [h_hs, Sigh_hs, B0_hs, alpha_hs, Sig_mean_hs, Sig_quant_hs, phi_hs, tau_hs, lambda_hs] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);
    
    HC = 0;
    [h_ig, Sigh_ig, B0_ig, alpha_ig, Sig_mean_ig, Sig_quant_ig, phi_ig] = func_SVAR7(y, HC ,nsim, burnin, thin, p, var_id, 0.1);

tracePlot2(B0_hs, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0' '_DGP' num2str(DGP) '.fig']));
tracePlot2(B0_ig, B0', n, n)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_traceB0' '_DGP' num2str(DGP) '.fig']));


% Analysis of results (Tables and plots of error)
 examine_sims;

% Choose an example r
r = 250;

 % Plot h estimates
% Line colors
colors = [0.0 0.0 0.0;   % black (baseline)
          0.0 0.0 1.0;   % blue
          1.0 0.0 0.0]; % red 

 plot_array_comp_with_bands({squeeze(h_real(r,:,:)), squeeze(h_hs_quant(r,:,:,[1 3 5])), squeeze(h_ig_quant(r,:,:,[1 3 5]))}, 'h', 4*p+1:T, 'h', colors)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\example_h' '_DGP' num2str(DGP) '.fig']));



%% ====================================================================== %
%  Empirical Application
% ======================================================================= %
clear
addpath("Empirical\");
addpath("Empirical\data_raw\")
addpath("Empirical\figures\")
addpath("Empirical\tables\")
addpath("Empirical\estimation\")
addpath("Empirical\Braun_oil_app\")
addpath("Empirical\Braun_oil_app\results")
addpath("Empirical\Braun_oil_app\functions\")
addpath("functions\")
addpath("functions\utility\")

%% Prepare data

data_prep_oil;

%% Load data if prepared already
clear;
load('Empirical\oilmarket_data_extended.mat')       % Main sample (extended)
load('Empirical\BB22_data.mat')                     % Braun & Brüggemann sample
load('Empirical\oil_data.mat')                      % New data (restricted main sample)


%%    % New sample
samplename = 'new';
Dates = oil_data.Date;
startdate = datetime(1988,2,1); enddate = datetime(2025,5,1);

% Estimate HS and IG model
idx_data = and(Dates >=  startdate ,Dates <= enddate );
data = 100*log(table2array(oil_data(idx_data, {'prod', 'rea', 'rpo', 'oecdcrude_sa'})));
Dates = Dates(idx_data);
figure; plot(data)

% Parameters
p = 13;                     % Lags
n = 4;                      % variables
T=size(data,1);             % observations

    estimate_oil;          % estimate model if needed

% Run Braun & Brüggemann models
p = 13; 
y = 100 * log(table2array(oil_data(:, {'prod', 'rea', 'rpo', 'oecdcrude_sa'}))); 
m = [data_BB22.kilian08extended; NaN(size(oil_data,1) - size(data_BB22,1),1)];              % Not used for estimation, but needed for the structure
date = Dates;

    % Runme_Models_Main;              % Run main model estimation file 

% If model is estimated, load
load(sprintf('Empirical/estimation/oil_ig_%s.mat', samplename));
load(sprintf('Empirical/estimation/oil_hs_%s.mat', samplename));

% Labels
 oil_labels;

% Run analysis
 oil_analysis;

%%  BB22 sample - older sample

samplename = 'old';
Dates = data_BB22.date;
startdate = datetime(1978,10-2,1); enddate = datetime(2018,11,1);

% Estimate HS and IG model
idx_data = and(Dates >=  startdate ,Dates <= enddate );
data = 100*log(table2array(data_extended(idx_data, {'Prod', 'WIP', 'RAC', 'inventsa'})));
Dates = Dates(idx_data);
figure; plot(data)

% Parameters
p = 13;                     % Lags
n = 4;                      % variables
T=size(data,1);             % observations

%    estimate_oil;          % estimate model if needed


% Run Braun & Brüggemann models
p = 13; 
y = 100 * log(table2array(data_BB22(:, {'Prod', 'WIP', 'RAC', 'inventsa'}))); 
m = data_BB22.kilian08extended;              % Not used for estimation, but needed for the structure
date = Dates;

%     Runme_Models_Main;              % Run main model estimation file 
 


% If model is estimated, load
load(sprintf('Empirical/estimation/oil_ig_%s.mat', samplename));
load(sprintf('Empirical/estimation/oil_hs_%s.mat', samplename));

% Labels
 oil_labels;








