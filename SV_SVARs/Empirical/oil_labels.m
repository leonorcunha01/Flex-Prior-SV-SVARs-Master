%% ========================================================================
%%  Oil Market Shocks Labeling
%  Examines correlations of structural shocks to Baumeister shocks and
%  BB22 shocks, for HR20 and BH19 restrictions.
%  ========================================================================


%% --- Correlations Baumeister shocks ---

% Load demand & supply shocks
demandBaum = readtable('data_raw/BH2_demand_shocks.xlsx', 'Sheet', 1, 'Range', 'A3:D606');
demandBaum.Properties.VariableNames = {'Date', 'econactiv', 'demand', 'invdemand'};

supplyBaum = readtable('data_raw/BH2_supply_shocks.xlsx', 'Sheet', 1, 'Range', 'A3:B603');
supplyBaum.Properties.VariableNames = {'Date', 'supply'};

% Restrict sample
startDate = datetime('01-Jun-1992', 'InputFormat', 'dd-MMM-yyyy');
endDate   = datetime('01-Feb-2025', 'InputFormat', 'dd-MMM-yyyy');
for t = {'demandBaum', 'supplyBaum'}
    tbl = eval(t{1});
    tbl = tbl(tbl.Date >= startDate & tbl.Date <= endDate, :);
    assignin('base', t{1}, tbl);
end

% Combine Baum shocks
Baumshocks = innerjoin(demandBaum(:, {'Date','demand','invdemand','econactiv'}), supplyBaum, 'Keys', 'Date');


% Estimated shocks (HS and IG)
estimHS = table(Dates(4*p+1:end), modelHS.struct_shock_mean, 'VariableNames', {'Date','ShocksHS'});
estimIG = table(Dates(4*p+1:end), modelIG.struct_shock_mean, 'VariableNames', {'Date','ShocksIG'});
estimShocks = innerjoin(estimHS, estimIG, 'Keys', 'Date');

% Merge all
mergedData = innerjoin(Baumshocks, estimShocks, 'Keys', 'Date');

% --- Compare shocks ---
varsToCompare = {'supply','demand','invdemand','econactiv'};

corrBaumTab = array2table(corr(Baumshocks{:,varsToCompare}),'VariableNames', varsToCompare, 'RowNames', varsToCompare);

disp('Correlation among Baum shocks:')
disp(corrBaumTab)

save(['Empirical\tables\corrBaum_' samplename '.mat'], 'corrBaumTab');


% Preallocate correlation arrays
corrHS = NaN(numel(varsToCompare), 4);
corrIG = NaN(numel(varsToCompare), 4);

for j = 1:4
    for v = 1:numel(varsToCompare)
        corrHS(v,j) = corr(mergedData.(varsToCompare{v}), mergedData.ShocksHS(:,j));
        corrIG(v,j) = corr(mergedData.(varsToCompare{v}), mergedData.ShocksIG(:,j));
    end
end

% Display and save
corrBaumHSTab = array2table(corrHS, 'VariableNames', "ShockHS"+(1:4), 'RowNames', varsToCompare);
disp('Correlation between Baum shocks and HS shocks:')
disp(corrBaumHSTab)
save(['Empirical\tables\corrBaum_HS_' samplename '.mat'], 'corrBaumHSTab');

corrBaumIGTab = array2table(corrIG, 'VariableNames', "ShockIG"+(1:4), 'RowNames', varsToCompare);
disp('Correlation between Baum shocks and IG shocks:')
disp(corrBaumIGTab)
save(['Empirical\tables\corrBaum_IG_' samplename '.mat'], 'corrBaumIGTab');



%% --- Braun–Brüggemann HR20 shocks ---


load(sprintf('Empirical/Braun_oil_app/results/%s_SR_HR20_main.mat', samplename))

BraunBrugHR20.dates = output.input.date;
BraunBrugHR20.alpha = output.alphas;
BraunBrugHR20.vecB = output.vecBs;
BraunBrugHR20.y_til = [output.input.y output.input.m];
[ntilsq, nrep] = size(BraunBrugHR20.vecB);
ntil = sqrt(ntilsq);

BraunBrugHR20.Btil = NaN(ntil,ntil,nrep);
for i = 1:nrep
    Bdraw = reshape(BraunBrugHR20.vecB(:,i), ntil, ntil);
    Bdraw = Bdraw .* (-sign(Bdraw(1,:))); % Braun-Brugg normalization
    BraunBrugHR20.Btil(:,:,i) = Bdraw;
end


[BraunBrugHR20.u_hat, BraunBrugHR20.struct_shock] = recover_shocksBraun( ...
    BraunBrugHR20.y_til, BraunBrugHR20.Btil, BraunBrugHR20.alpha, p, 1);

BraunBrugHR20.struct_shock_mean = squeeze(mean(BraunBrugHR20.struct_shock));

plot_array(BraunBrugHR20.struct_shock_mean, BraunBrugHR20.dates(p+1:end), 'struct. shocks', 'Shock')
    save(sprintf('Empirical/figures/struc_shocks_%sBraunBrugHR20.fig', samplename));

% --- Prepare BraunBrug shocks ---
braunNames = {'supply','flowdem','invdemand','residual'};
BBHR20 = table(BraunBrugHR20.dates(p+1:end), BraunBrugHR20.struct_shock_mean, 'VariableNames', {'Date','BraunBrugShocks'});
for i = 1:numel(braunNames)
    BBHR20.(braunNames{i}) = BBHR20.BraunBrugShocks(:,i);
end

% Add estimated HS/IG shocks
BBHR20 = innerjoin(BBHR20, estimHS, 'Keys', 'Date');
BBHR20 = innerjoin(BBHR20, estimIG, 'Keys', 'Date');


%% --- Correlations BraunBrugg HR20 ---
disp('Correlation among Braun–Brüggemann shocks:')
disp(array2table(corr(BBHR20{:,braunNames}), ...
     'VariableNames', braunNames, 'RowNames', braunNames))

corr_BB_HR20_HS = corr(BBHR20{:,braunNames}, BBHR20.ShocksHS);
corr_BB_HR20_IG = corr(BBHR20{:,braunNames}, BBHR20.ShocksIG);

corr_BB_HR20_HSTab = array2table(corr_BB_HR20_HS, 'VariableNames', "HS"+(1:4), 'RowNames', braunNames);
disp('Correlation: BraunBrug shocks vs HS estimated shocks')
disp(corr_BB_HR20_HSTab)
save(['Empirical\tables\corrBB_HR20_HS_' samplename '.mat'], 'corr_BB_HR20_HSTab');

corr_BB_HR20_IGTab = array2table(corr_BB_HR20_IG, 'VariableNames', "IG"+(1:4), 'RowNames', braunNames);
disp('Correlation: BraunBrug shocks vs IG estimated shocks')
disp(corr_BB_HR20_IGTab)
save(['Empirical\tables\corrBB_HR20_IG_' samplename '.mat'], 'corr_BB_HR20_IGTab');



%% --- Braun–Brüggemann BH19 shocks ---

load(sprintf('Empirical/Braun_oil_app/results/%s_SR_BH19_main.mat', samplename))

BraunBrugBH19.dates = output.input.date;
BraunBrugBH19.alpha = output.alphas;
BraunBrugBH19.vecB  = output.vecBs;
BraunBrugBH19.y_til = [output.input.y output.input.m];
[~, nrep] = size(BraunBrugBH19.vecB);

BraunBrugBH19.Btil = NaN(ntil, ntil, nrep);
for i = 1:nrep
    Bdraw = reshape(BraunBrugBH19.vecB(:, i), ntil, ntil);
    Bdraw = Bdraw .* (-sign(Bdraw(1, :))); % Braun–Brüggemann normalization
    BraunBrugBH19.Btil(:, :, i) = Bdraw;
end

[BraunBrugBH19.u_hat, BraunBrugBH19.struct_shock] = recover_shocksBraun( ...
    BraunBrugBH19.y_til, BraunBrugBH19.Btil, BraunBrugBH19.alpha, p, 1);

BraunBrugBH19.struct_shock_mean = squeeze(mean(BraunBrugBH19.struct_shock));

plot_array(BraunBrugBH19.struct_shock_mean, BraunBrugBH19.dates(p+1:end), 'struct. shocks', 'Shock')
    save(sprintf('Empirical/figures/struc_shocks_%sBraunBrugBH19.fig', samplename));


% --- Prepare Braun–Brüggemann shocks (BH19) ---

braunNames = {'supply', 'flowdem', 'invdemand', 'residual'};
BBBH19 = table(BraunBrugBH19.dates(p+1:end), BraunBrugBH19.struct_shock_mean, 'VariableNames', {'Date', 'BraunBrugShocks'});

for i = 1:numel(braunNames)
    BBBH19.(braunNames{i}) = BBBH19.BraunBrugShocks(:, i);
end

% Add estimated HS/IG shocks
BBBH19 = innerjoin(BBBH19, estimHS, 'Keys', 'Date');
BBBH19 = innerjoin(BBBH19, estimIG, 'Keys', 'Date');

%% --- Correlations Braun–Brüggemann (BH19) ---
disp('Correlation among Braun–Brüggemann shocks (BH19):')
disp(array2table(corr(BBBH19{:, braunNames}), ...
     'VariableNames', braunNames, 'RowNames', braunNames))

corr_BB_BH19_HS = corr(BBBH19{:, braunNames}, BBBH19.ShocksHS);
corr_BB_BH19_IG = corr(BBBH19{:, braunNames}, BBBH19.ShocksIG);

corr_BB_BH19_HSTab = array2table(corr_BB_BH19_HS, 'VariableNames', "HS"+(1:4), 'RowNames', braunNames);
disp('Correlation: BraunBrug BH19 shocks vs HS estimated shocks')
disp(corr_BB_BH19_HSTab)
save(['Empirical\tables\corrBB_BH19_HS_' samplename '.mat'], 'corr_BB_BH19_HSTab');

corr_BB_BH19_IGTab = array2table(corr_BB_BH19_IG, 'VariableNames', "IG"+(1:4), 'RowNames', braunNames);
disp('Correlation: BraunBrug BH19 shocks vs IG estimated shocks')
disp(corr_BB_BH19_IGTab)
save(['Empirical\tables\corrBB_BH19_IG_' samplename '.mat'], 'corr_BB_BH19_IGTab');


