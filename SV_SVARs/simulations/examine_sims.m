% File to compute MAE and MSE from raw results from a simulation,
% create tables and plots of these measures and save them.


%% Real quantities

Sig_real = zeros(R, T-4*p, n, n);
for r = 1:R
    Sig_real(r,:,:,:) = construct_Sigt(squeeze(h_real(r,:,:)), B0);
end


%% Error measures
% MSE and MAE, for either version (HS or IG), using median or mean as point
% estimates

%%% Errors
% HS mean and median
B0_hs_mean_err     = B0_hs_norm_mean   - vec(B0_norm')';
alpha_hs_mean_err  = alpha_hs_mean     - alpha;
Sig_hs_mean_err    = Sig_hs_mean       - Sig_real;

B0_hs_med_err      = B0_hs_norm_quant(:,:,3)           - vec(B0_norm')';
alpha_hs_med_err   = alpha_hs_quant(:,:,3)             - alpha;
Sig_hs_med_err     = squeeze(Sig_hs_quant(:,3,:,:,:))  - Sig_real;


% IG mean and median
B0_ig_mean_err     = B0_ig_norm_mean   - vec(B0_norm')';
alpha_ig_mean_err  = alpha_ig_mean     - alpha;
Sig_ig_mean_err    = Sig_ig_mean       - Sig_real;

B0_ig_med_err      = B0_ig_norm_quant(:,:,3)           - vec(B0_norm')';
alpha_ig_med_err   = alpha_ig_quant(:,:,3)             - alpha;
Sig_ig_med_err     = squeeze(Sig_ig_quant(:,3,:,:,:))  - Sig_real;


%%% MAE
% HS
B0_hs_med_MAE_ind   = mean(abs(B0_hs_med_err));
B0_hs_med_MAE       = mean(B0_hs_med_MAE_ind(2:3));
alpha_hs_med_MAE_ind = mean(abs(alpha_hs_med_err));
alpha_hs_med_MAE     = mean(alpha_hs_med_MAE_ind);
Sig_hs_med_MAE_ind   = squeeze(mean(abs(Sig_hs_med_err)));
Sig_hs_med_MAE       = squeeze(mean(Sig_hs_med_MAE_ind));

B0_hs_mean_MAE_ind   = mean(abs(B0_hs_mean_err));
B0_hs_mean_MAE       = mean(B0_hs_mean_MAE_ind(2:3));
alpha_hs_mean_MAE_ind = mean(abs(alpha_hs_mean_err));
alpha_hs_mean_MAE     = mean(alpha_hs_mean_MAE_ind);
Sig_hs_mean_MAE_ind   = squeeze(mean(abs(Sig_hs_mean_err)));
Sig_hs_mean_MAE       = squeeze(mean(Sig_hs_mean_MAE_ind));

% IG
B0_ig_med_MAE_ind   = mean(abs(B0_ig_med_err));
B0_ig_med_MAE       = mean(B0_ig_med_MAE_ind(2:3));
alpha_ig_med_MAE_ind = mean(abs(alpha_ig_med_err));
alpha_ig_med_MAE     = mean(alpha_ig_med_MAE_ind);
Sig_ig_med_MAE_ind   = squeeze(mean(abs(Sig_ig_med_err)));
Sig_ig_med_MAE       = squeeze(mean(Sig_ig_med_MAE_ind));

B0_ig_mean_MAE_ind   = mean(abs(B0_ig_mean_err));
B0_ig_mean_MAE       = mean(B0_ig_mean_MAE_ind(2:3));
alpha_ig_mean_MAE_ind = mean(abs(alpha_ig_mean_err));
alpha_ig_mean_MAE     = mean(alpha_ig_mean_MAE_ind);
Sig_ig_mean_MAE_ind   = squeeze(mean(abs(Sig_ig_mean_err)));
Sig_ig_mean_MAE       = squeeze(mean(Sig_ig_mean_MAE_ind));


%%% MSE
% HS
B0_hs_med_MSE_ind   = mean(B0_hs_med_err.^2);
B0_hs_med_MSE       = mean(B0_hs_med_MSE_ind(2:3));
alpha_hs_med_MSE_ind = mean(alpha_hs_med_err.^2);
alpha_hs_med_MSE     = mean(alpha_hs_med_MSE_ind);
Sig_hs_med_MSE_ind   = squeeze(mean(Sig_hs_med_err.^2));
Sig_hs_med_MSE       = squeeze(mean(Sig_hs_med_MSE_ind));

B0_hs_mean_MSE_ind   = mean(B0_hs_mean_err.^2);
B0_hs_mean_MSE       = mean(B0_hs_mean_MSE_ind(2:3));
alpha_hs_mean_MSE_ind = mean(alpha_hs_mean_err.^2);
alpha_hs_mean_MSE     = mean(alpha_hs_mean_MSE_ind);
Sig_hs_mean_MSE_ind   = squeeze(mean(Sig_hs_mean_err.^2));
Sig_hs_mean_MSE       = squeeze(mean(Sig_hs_mean_MSE_ind));


% IG
B0_ig_med_MSE_ind   = mean(B0_ig_med_err.^2);
B0_ig_med_MSE       = mean(B0_ig_med_MSE_ind(2:3));
alpha_ig_med_MSE_ind = mean(alpha_ig_med_err.^2);
alpha_ig_med_MSE     = mean(alpha_ig_med_MSE_ind);
Sig_ig_med_MSE_ind   = squeeze(mean(Sig_ig_med_err.^2));
Sig_ig_med_MSE       = squeeze(mean(Sig_ig_med_MSE_ind));


B0_ig_mean_MSE_ind   = mean(B0_ig_mean_err.^2);
B0_ig_mean_MSE       = mean(B0_ig_mean_MSE_ind(2:3));
alpha_ig_mean_MSE_ind = mean(alpha_ig_mean_err.^2);
alpha_ig_mean_MSE     = mean(alpha_ig_mean_MSE_ind);
Sig_ig_mean_MSE_ind   = squeeze(mean(Sig_ig_mean_err.^2));
Sig_ig_mean_MSE       = squeeze(mean(Sig_ig_mean_MSE_ind));



%% Tables and Plots for error measures
% === Average error tables ===
suffixes = {'hs_med_MAE','ig_med_MAE','hs_mean_MAE','ig_mean_MAE', ...
            'hs_med_MSE','ig_med_MSE','hs_mean_MSE','ig_mean_MSE'};
eval(['AllTables_DGP' num2str(DGP) ' = struct();']);
for i = 1:length(suffixes)
    suffix = suffixes{i};
    alpha_val = round(eval(['alpha_' suffix]),4);
    B0_val = round(eval(['B0_' suffix]),4);
    Sig_val = round(eval(['Sig_' suffix]),4);
    Sig_entries = Sig_val(:)';
    Tab = table(alpha_val,B0_val,Sig_entries(1),Sig_entries(2),Sig_entries(3),Sig_entries(4), ...
        'VariableNames',{'alpha','B0','Sig11','Sig12','Sig21','Sig22'});
    eval(['AllTables_DGP' num2str(DGP) '.Table_' suffix '_DGP' num2str(DGP) ' = Tab;']);
end
save(['tables\AllTables_DGP' num2str(DGP) '.mat'], ['AllTables_DGP' num2str(DGP)]);


% === Individual error tables, alpha and B0 ===
suffixes = {'hs_med_MAE_ind','ig_med_MAE_ind','hs_mean_MAE_ind','ig_mean_MAE_ind', ...
            'hs_med_MSE_ind','ig_med_MSE_ind','hs_mean_MSE_ind','ig_mean_MSE_ind'};
eval(['AllTables_ind_DGP' num2str(DGP) ' = struct();']);
for i = 1:length(suffixes)
    suffix = suffixes{i};
    alpha_val = round(eval(['alpha_' suffix]),4); alpha_val = alpha_val(:)';
    B0_val = round(eval(['B0_' suffix]),4); B0_val = B0_val(:)';
    alphaNames = strcat('alpha', string(1:length(alpha_val)));
    B0Names = strcat('B0', string(1:length(B0_val)));
    Tab = array2table([alpha_val,B0_val], 'VariableNames', [alphaNames,B0Names]);
    eval(['AllTables_ind_DGP' num2str(DGP) '.Table_' suffix '_DGP' num2str(DGP) ' = Tab;']);
end
save(['tables\AllTables_ind_DGP' num2str(DGP) '.mat'], ['AllTables_ind_DGP' num2str(DGP)]);


% === Sigma individual errors comparison plots ===
suffixes = {'med_MAE_ind','mean_MAE_ind','med_MSE_ind','mean_MSE_ind'};
for i = 1:length(suffixes)
    suffix = suffixes{i};
    Sig_hs_val = eval(['Sig_hs_' suffix]);
    Sig_ig_val = eval(['Sig_ig_' suffix]);
    combined_hs = [Sig_hs_val(:,1,1), Sig_hs_val(:,2,2), Sig_hs_val(:,1,2)];
    combined_ig = [Sig_ig_val(:,1,1), Sig_ig_val(:,2,2), Sig_ig_val(:,1,2)];
    plot_array_comp({combined_hs, combined_ig}, 'sig_', 4*p+1:T);
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\Sig_' suffix '_DGP' num2str(DGP) '_comp.fig']));
end




% Plot MAE Sigt, around jumps
jump = 0.8*T-4*p;                                                           % Time point of jump
idx = jump-10:jump+10;                                                      % Interval to plot
% sigma_(1,1)
vals = [Sig_hs_med_MAE_ind(idx,1,1) Sig_ig_med_MAE_ind(idx,1,1)];
figure; bar(idx+4*p,vals); xticks(idx+4*p); xlabel('Time'); ylabel('MAE'); legend(types)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\Sig11_med_MAE_closeup_DGP' num2str(DGP) '.fig']));

% sigma_(2,2)
vals = [Sig_hs_med_MAE_ind(idx,2,2) Sig_ig_med_MAE_ind(idx,2,2)];
figure; bar(idx+4*p,vals); xticks(idx+4*p); xlabel('Time'); ylabel('MAE'); legend(types)
    saveas(gcf, fullfile('figures\', ['DGP' num2str(DGP) '\Sig22_med_MAE_closeup_DGP' num2str(DGP) '.fig']));


