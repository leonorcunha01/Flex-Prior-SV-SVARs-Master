%% ========================================================================
%%  Estimation Oil Market Data
%  Estimates the model in the HS versionband IG version 
%  with v_b = 0.1.
%  ========================================================================

%% --- Estimation ---

% Parameters
nsim = 1000000;              % Post-burnin iterations
burnin = 50000;              % Burnin iterations
thin = 1000;                 % Thining
var_id = [1 2 3 4];         % All vars


rng(66, 'twister')

 % HS model, v_b = 0.1
HC = 1;
[modelHS.h, modelHS.Sigh, modelHS.B0, modelHS.alpha, modelHS.Sig_mean, modelHS.Sig_quant, modelHS.phi, modelHS.tau, modelHS.lambda] = ...
    func_SVAR7(data, HC ,nsim, burnin, thin, p, var_id, 0.1);

% IG model, v_b = 0.1
HC = 0;
[modelIG.h, modelIG.Sigh, modelIG.B0, modelIG.alpha, modelIG.Sig_mean, modelIG.Sig_quant, modelIG.phi] = ...
    func_SVAR7(data, HC ,nsim, burnin, thin, p, var_id, 0.1);


% Normalized B0
modelHS.B0_norm = normalizeB0draws(modelHS.B0);
modelIG.B0_norm = normalizeB0draws(modelIG.B0);

% Estimates
modelHS.alpha_med = squeeze(median(modelHS.alpha));
modelHS.B0norm_med = squeeze(median(modelHS.B0_norm));
modelIG.alpha_med = squeeze(median(modelIG.alpha));
modelIG.B0norm_med = squeeze(median(modelIG.B0_norm));


% Recover structural shocks
% (use non-normalized B0)
[modelHS.u_hat, modelHS.struct_shock] = recover_shocks(data, ...
    modelHS.B0, modelHS.alpha, p);

[modelIG.u_hat, modelIG.struct_shock] = recover_shocks(data, ...
    modelIG.B0, modelIG.alpha, p);

modelHS.struct_shock_med = squeeze(median(modelHS.struct_shock));
modelIG.struct_shock_med = squeeze(median(modelIG.struct_shock));
modelHS.struct_shock_mean = squeeze(mean(modelHS.struct_shock));
modelIG.struct_shock_mean = squeeze(mean(modelIG.struct_shock));

% Plot
Dates_forplot = Dates(4*p+1:end);
    % Structural Shocks
plot_array(modelHS.struct_shock_mean, Dates_forplot,'struct. shocks', 'Shock')
    saveas(gcf, sprintf('Empirical/figures/struc_shocks_%s_modelHS.fig', samplename));

plot_array(modelIG.struct_shock_mean, Dates_forplot,'struct. shocks', 'Shock')
    saveas(gcf,sprintf('Empirical/figures/struc_shocks_%s_modelIG.fig', samplename));

    % Comparison structural shocks
plot_array_comp({modelHS.struct_shock_mean, modelIG.struct_shock_mean}, 'Struct. Shock', Dates_forplot, 'Shock')
    saveas(gcf,sprintf('Empirical/figures/struct_shocks%scomp.fig', samplename));

 
% Plots: h process comparison
colorsh = [0.0 0.0 1.0;   % blue
            1.0 0.0 0.0]; % red 

quants = [.05 .5 .95];
ModelHS.h_quant = NaN(length(Dates_forplot), n, 3);
ModelIG.h_quant = NaN(length(Dates_forplot), n, 3);
for q = 1:3    
    ModelHS.h_quant(:,:,q) = quantile(modelHS.h,quants(q), 1);
    ModelIG.h_quant(:,:,q) = quantile(modelIG.h,quants(q), 1);
end 
plot_array_comp_with_bands({ModelHS.h_quant, ModelIG.h_quant}, 'h', Dates_forplot', 'h', colorsh)
    saveas(gcf,sprintf('Empirical/figures/h_%scomp.fig', samplename));

% Plots: densities B0 comparison
plot_hist_comparison(modelHS.B0_norm, modelIG.B0_norm, n, n)
    saveas(gcf,sprintf('Empirical/figures/B0_norm_hist_%scomp.fig', samplename));



% Save
save(sprintf('Empirical/estimation/oil_ig_%s.mat', samplename), 'modelIG');
save(sprintf('Empirical/estimation/oil_hs_%s.mat', samplename), 'modelHS');


% Quick convergence check: Trace plots
% B0
tracePlot2(modelHS.B0, median(modelHS.B0), n, n)
    saveas(gcf,sprintf('Empirical/figures/traceB0_%smodelHS.fig', samplename));

tracePlot2(modelIG.B0, mean(modelIG.B0), n, n)
    saveas(gcf,sprintf('Empirical/figures/traceB0_%smodelIG.fig', samplename));

% alpha
k=n*p+1;
tracePlot2(modelHS.alpha(:,2:k), mean(modelHS.alpha(:,2:k)), p, n) %first row
tracePlot2(modelHS.alpha(:,k+2:2*k), mean(modelHS.alpha(:,k+2:2*k)), p, n) %second row
tracePlot2(modelHS.alpha(:,2*k+2:3*k), mean(modelHS.alpha(:,2*k+2:3*k)), p, n) %third row
tracePlot2(modelHS.alpha(:,3*k+2:4*k), mean(modelHS.alpha(:,3*k+2:4*k)), p, n) %fourth row
tracePlot2(modelHS.alpha(:,[1 k+1 2*k+1 3*k+1]), mean(modelHS.alpha(:,[1 k+1 2*k+1 3*k+1])), 1, n) %intercepts

tracePlot2(modelIG.alpha(:,2:k), mean(modelIG.alpha(:,2:k)), p, n) %first row
tracePlot2(modelIG.alpha(:,k+2:2*k), mean(modelIG.alpha(:,k+2:2*k)), p, n) %second row
tracePlot2(modelIG.alpha(:,2*k+2:3*k), mean(modelIG.alpha(:,2*k+2:3*k)), p, n) %third row
tracePlot2(modelIG.alpha(:,3*k+2:4*k), mean(modelIG.alpha(:,3*k+2:4*k)), p, n) %fourth row
tracePlot2(modelIG.alpha(:,[1 k+1 2*k+1 3*k+1]), mean(modelIG.alpha(:,[1 k+1 2*k+1 3*k+1])), 1, n) %intercepts




