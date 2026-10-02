%% ================================
%  Compare Baum shocks (supply, demand, invdemand, econactiv)
%  against estimated HS and IG shocks
%  ================================

varsToCompare = {'supply', 'demand', 'invdemand', 'econactiv'};
shockTypes = {'HS', 'IG'};

%% --- Correlation of Baum shocks ---
Baum_shocks = Baumdata.Baumshocks{:,varsToCompare};
corrBaum = corr(Baum_shocks);
corrTabBaum = array2table(corrBaum, ...
    'VariableNames', varsToCompare, ...
    'RowNames', varsToCompare);
disp(corrTabBaum)

%% --- Correlation tables (Baum shocks × HS/IG shocks) ---
fprintf('\n=== Correlation Tables ===\n')

% HS correlations
corrHS = NaN(length(varsToCompare), 4);
for v = 1:length(varsToCompare)
    for j = 1:4
        corrHS(v, j) = corr(mergedData.(varsToCompare{v}), mergedData.ShocksHS(:, j));
    end
end
corrTableHS = array2table(corrHS, ...
    'VariableNames', {'ShockHS1','ShockHS2','ShockHS3','ShockHS4'}, ...
    'RowNames', varsToCompare);
disp('Correlation between Baum shocks and HS shocks:')
disp(corrTableHS)

% IG correlations
corrIG = NaN(length(varsToCompare), 4);
for v = 1:length(varsToCompare)
    for j = 1:4
        corrIG(v, j) = corr(mergedData.(varsToCompare{v}), mergedData.ShocksIG(:, j));
    end
end
corrTableIG = array2table(corrIG, ...
    'VariableNames', {'ShockIG1','ShockIG2','ShockIG3','ShockIG4'}, ...
    'RowNames', varsToCompare);
disp('Correlation between Baum shocks and IG shocks:')
disp(corrTableIG)

%% --- Loop through all Baum variables and shocks for regressions and plots ---
% for v = 1:length(varsToCompare)
%     varName = varsToCompare{v};
% 
%     for s = 1:length(shockTypes)
%         shockType = shockTypes{s};
%         shockField = ['Shocks' shockType];
% 
%         fprintf('\n=== %s vs %s shocks ===\n', upper(varName), shockType);
% 
%         % --- Auxiliary regression ---
%         shocksMat = mergedData.(shockField);
%         for i = 1:size(shocksMat, 2)
%             mergedData.(['Shock' shockType num2str(i)]) = shocksMat(:, i);
%         end
% 
%         regFormula = sprintf('%s ~ -1 + Shock%s1 + Shock%s2 + Shock%s3 + Shock%s4', ...
%                              varName, shockType, shockType, shockType, shockType);
%         auxReg = fitlm(mergedData, regFormula);
%         disp(auxReg)
% 
%         % --- Visual comparison (line plot) ---
%         plot_h_comp(repmat(mergedData.(varName), 1, 4), mergedData.(shockField))
% 
%         % --- Scatter plots ---
%         figure;
%         for i = 1:4
%             subplot(2, 2, i)
%             scatter(mergedData.(['Shock' shockType num2str(i)]), mergedData.(varName))
%             xlabel(['Shock ' num2str(i)])
%             ylabel(strrep(varName, '_', '\_'))
%             title(sprintf('%s vs Shock %d', strrep(varName, '_', '\_'), i))
%         end
%         sgtitle(sprintf('Scatter plots of %s (Baum19) vs each %s shock', ...
%                         strrep(varName, '_', '\_'), shockType))
%     end
% end

