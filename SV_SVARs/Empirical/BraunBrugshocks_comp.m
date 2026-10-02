%% ================================
%  Compare BraunBrug shocks (supply, demand, invdemand, econactiv)
%  against estimated HS and IG shocks
%  ================================

%% --- Correlation of BH shocks ---
BB_shocks = BraunBrugData.BraunBrugShocks{:,braunNames};
corrBB = corr(BB_shocks);
corrTabBB = array2table(corrBB, ...
    'VariableNames', braunNames, ...
    'RowNames', braunNames);
disp(corrTabBB)

%% Correlations
corr_BraunBrug_HS = corr(table2array(BraunBrug(:, braunNames)), BraunBrug.ShocksHS);
corr_BraunBrug_IG = corr(table2array(BraunBrug(:, braunNames)), BraunBrug.ShocksIG);

disp('Correlation: BraunBrug shocks vs HS estimated shocks');
disp(array2table(corr_BraunBrug_HS, ...
    'VariableNames', {'HS1','HS2','HS3','HS4'}, ...
    'RowNames', braunNames));

disp('Correlation: BraunBrug shocks vs IG estimated shocks');
disp(array2table(corr_BraunBrug_IG, ...
    'VariableNames', {'IG1','IG2','IG3','IG4'}, ...
    'RowNames', braunNames));


% %% --- Auxiliary Regressions (HS and IG) ---
% % Split HS and IG shocks into columns for regression
% shocksMatHS = BraunBrug.ShocksHS;
% shocksMatIG = BraunBrug.ShocksIG;
% 
% for i = 1:4
%     BraunBrug.(['ShockHS' num2str(i)]) = shocksMatHS(:, i);
%     BraunBrug.(['ShockIG' num2str(i)]) = shocksMatIG(:, i);
% end
% 
% % Containers for regression results
% auxReg_BraunBrug_HS = struct();
% auxReg_BraunBrug_IG = struct();
% 
% % Regress each BraunBrug shock on HS and IG shocks
% for v = 1:length(braunNames)
%     shockVar = braunNames{v};
% 
%     % --- HS regression ---
%     modelHS_lm = fitlm(BraunBrug, ...
%         sprintf('%s ~ -1 + ShockHS1 + ShockHS2 + ShockHS3 + ShockHS4', shockVar));
%     auxReg_BraunBrug_HS.(shockVar) = modelHS_lm;
%     fprintf('\nAuxiliary regression (HS): %s\n', shockVar);
%     disp(modelHS_lm);
% 
%     % --- IG regression ---
%     modelIG_lm = fitlm(BraunBrug, ...
%         sprintf('%s ~ -1 + ShockIG1 + ShockIG2 + ShockIG3 + ShockIG4', shockVar));
%     auxReg_BraunBrug_IG.(shockVar) = modelIG_lm;
%     fprintf('\nAuxiliary regression (IG): %s\n', shockVar);
%     disp(modelIG_lm);
% end


%% --- Visual Comparison: BraunBrug vs HS ---
for v = 1:length(braunNames)
    varName = braunNames{v};

    % Multi-line comparison plot
    plot_h_comp(repmat(BraunBrug.(varName), 1, 4)/20, BraunBrug.ShocksHS)

    % % Scatter plots
    % figure;
    % for i = 1:4
    %     subplot(2,2,i)
    %     scatter(BraunBrug.(['ShockHS' num2str(i)]), BraunBrug.(varName))
    %     xlabel(['HS Shock ' num2str(i)])
    %     ylabel(varName)
    %     title([varName ' vs HS Shock ' num2str(i)])
    % end
    % sgtitle(['Scatter plots of ' varName ' (BraunBrug) vs HS shocks'])
end
% 
% 
% %% --- Visual Comparison: BraunBrug vs IG ---
% for v = 1:length(braunNames)
%     varName = braunNames{v};
% 
%     % Multi-line comparison plot
%     plot_h_comp(repmat(BraunBrug.(varName), 1, 4), BraunBrug.ShocksIG)
% 
%     % Scatter plots
%     figure;
%     for i = 1:4
%         subplot(2,2,i)
%         scatter(BraunBrug.(['ShockIG' num2str(i)]), BraunBrug.(varName))
%         xlabel(['IG Shock ' num2str(i)])
%         ylabel(varName)
%         title([varName ' vs IG Shock ' num2str(i)])
%     end
%     sgtitle(['Scatter plots of ' varName ' (BraunBrug) vs IG shocks'])
% end
