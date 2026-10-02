%% ========================================================================
%  Prepare Oil Market Data
%  Combines multiple sources, constructs oil_data (1988–2025),
%  compares with existing data (1978–2018), and extends the old dataset.
%  ========================================================================

%% ------------------------------------------------------------------------
%  Load existing dataset (Braun & Bruggeman 2022)
%  ------------------------------------------------------------------------
load('Empirical\Braun_oil_app\data\oilmarket_data.mat')  % Contains variable 'data' with fields date, WIP, invent, inventsa, WTI, RAC, Prod, killian08extended.
data_BB22 = data;
clear("data")
save('Empirical\BB22_data.mat', 'data_BB22')

%% ------------------------------------------------------------------------
%  WORLD PRODUCTION OF CRUDE OIL
%  ------------------------------------------------------------------------
worldProdcrude = readtable("Empirical\data_raw\WorldCrudeProd.csv", 'Range', 'G6:H634'); 
worldProdcrude.Properties.VariableNames = {'Date', 'prod'};
worldProdcrude.Date = datetime(worldProdcrude.Date, 'InputFormat', 'yyyy-MM');
worldProdcrude.prod = worldProdcrude.prod / 1000; % convert units

% Plot
figure;
hold on;
plot(worldProdcrude.Date, worldProdcrude.prod, 'LineWidth', 1.2);
plot(data_BB22.date, data_BB22.Prod);
hold off;
xlabel('Date'); ylabel('Production');
title('World Crude Production');
legend({'New sample', 'BB22'})
grid on;
    saveas(gcf, 'Empirical\figures\WorldProd.fig');

%% ------------------------------------------------------------------------
%  INDUSTRIAL PRODUCTION INDEX (BH19)
%  ------------------------------------------------------------------------
indexBH19 = readtable('Empirical/data_raw/BH19IPindex.xlsx', 'Sheet', 1, 'Range', 'A2:B811');
indexBH19.Properties.VariableNames = {'Date', 'rea'};

% Plot comparison
figure;
hold on;
plot(indexBH19.Date, indexBH19.rea, 'LineWidth', 1.2);
plot(data_BB22.date, data_BB22.WIP);
hold off;
xlabel('Date'); ylabel('Index');
title('Baumeister Real Activity Index');
legend({'New sample', 'BB22'})
grid on;
    saveas(gcf, 'Empirical\figures\index.fig');

%% ------------------------------------------------------------------------
%  REAL PRICE OF OIL
%  ------------------------------------------------------------------------
spotOilPrice = readtable('Empirical/data_raw/SpotOilPrice.xls', 'Sheet', 2, 'Range', 'A4:B479');
spotOilPrice.Properties.VariableNames = {'Date', 'price'};
spotOilPrice.Date = dateshift(spotOilPrice.Date, 'start', 'month');

cpi = readtable('Empirical/data_raw/CPIAUCSL.csv', 'Range', 'A2:B632'); 
cpi.Properties.VariableNames = {'Date', 'cpi'};

realPriceOil = innerjoin(spotOilPrice, cpi, 'Keys', 'Date');
realPriceOil.rpo = realPriceOil.price ./ (realPriceOil.cpi / 100);

% Plot comparison
figure;
hold on;
plot(realPriceOil.Date, realPriceOil.rpo);
plot(data_BB22.date, data_BB22.RAC);
hold off;
xlabel('Date'); ylabel('Price');
title('Real Price of Oil');
legend({'New sample', 'BB22'})
grid on;
    saveas(gcf, 'Empirical\figures\realPrice.fig');

%% ------------------------------------------------------------------------
%  INVENTORIES
%  ------------------------------------------------------------------------
USpet = readtable('Empirical\\data_raw\\USpetStocks.csv', 'Range', 'A8:B841'); 
USpet.Properties.VariableNames = {'Date', 'uspet'};
USpet.Date = datetime(USpet.Date, 'InputFormat', 'MMM yyyy');

UScrude = readtable('Empirical/data_raw/UScrudeStocks.csv', 'Range', 'A8:B841'); 
UScrude.Properties.VariableNames = {'Date', 'uscrude'};
UScrude.Date = datetime(UScrude.Date, 'InputFormat', 'MMM yyyy');

OECDpet = readtable('Empirical/data_raw/OECDpetStocks.csv', 'Range', 'BM6:BN550'); 
OECDpet.Properties.VariableNames = {'Date', 'oecdpet'};
OECDpet.Date = datetime(OECDpet.Date, 'InputFormat', 'yyyy-MM');

% Merge inventories
invent = innerjoin(USpet, UScrude, 'Keys', 'Date');
invent = innerjoin(invent, OECDpet, 'Keys', 'Date');

% Compute OECD crude
ratio = invent.oecdpet ./ invent.uspet;
oecdcrude = ratio .* invent.uscrude * 1000;
oecdcrude_table = table(invent.Date, oecdcrude, 'VariableNames', {'Date', 'oecdcrude'});

writetable(oecdcrude_table, 'Empirical/data_raw/oecdcrude_unadjusted.csv'); % Save for adjusting outside

% The seasonal adjustment is made in a small R script, and saved to Empirical/data_raw/oecdcrude_adjusted.csv

% Seasonally adjusted data
oecd_sa = readtable('Empirical/data_raw/oecdcrude_adjusted.csv');
oecd_sa.Date = datetime(oecd_sa.Date, 'InputFormat', 'yyyy-MM-dd');
oecd_sa = oecd_sa(oecd_sa.Date > datetime(1988,1,1), :);
invent = outerjoin(invent, oecdcrude_table(:, {'Date', 'oecdcrude'}), 'Keys', 'Date', 'MergeKeys', true);
invent = outerjoin(invent, oecd_sa(:, {'Date', 'oecdcrude_sa'}), 'Keys', 'Date', 'MergeKeys', true);

% Plot comparison
figure;
hold on;
plot(invent.Date, invent.oecdcrude_sa, 'LineWidth', 1.2);
plot(data_BB22.date, data_BB22.inventsa);
hold off;
xlabel('Date'); ylabel('Stock percent change');
title('OECD Crude Oil Stocks (Seasonally Adjusted)');
legend({'New sample', 'BB22'})
grid on;
    saveas(gcf, 'Empirical\figures\inventories_sa.fig');


%% ------------------------------------------------------------------------
%  JOIN ALL SERIES
%  ------------------------------------------------------------------------
joinedTable = innerjoin(worldProdcrude, indexBH19, 'Keys', 'Date');
joinedTable = innerjoin(joinedTable, realPriceOil, 'Keys', 'Date');
joinedTable = innerjoin(joinedTable, invent, 'Keys', 'Date');
oil_data = joinedTable(:, {'Date', 'prod', 'rea', 'rpo', 'oecdcrude', 'oecdcrude_sa'});
oil_data = oil_data(oil_data.Date > datetime(1988,1,1), :);

save("Empirical\oil_data.mat", "oil_data"); % Jan-1988 to May-2025

%% ------------------------------------------------------------------------
%  CORRELATION: LAST 20 YEARS OF OLD SAMPLE
%  ------------------------------------------------------------------------
last20_start = datetime(1998,10,1);
last20_end   = datetime(2018,10,1);

mask = data_BB22.date >= last20_start & data_BB22.date <= last20_end;

% Align by date
[commonDates, ia_rel, idxNew] = intersect(data_BB22.date(mask), oil_data.Date);

idxOld_fullMask = find(mask);          % indices in full data corresponding to the mask
idxOld = idxOld_fullMask(ia_rel);     % indices in the full 'data' struct

corrTable = table(corr(data_BB22.Prod(idxOld), oil_data.prod(idxNew), 'Rows', 'complete'), ...
    corr(data_BB22.WIP(idxOld), oil_data.rea(idxNew), 'Rows', 'complete'), ...
    corr(data_BB22.RAC(idxOld), oil_data.rpo(idxNew), 'Rows', 'complete'), ...
    corr(data_BB22.inventsa(idxOld), oil_data.oecdcrude_sa(idxNew), 'Rows', 'complete'), ...
    corr(data_BB22.invent(idxOld), oil_data.oecdcrude(idxNew), 'Rows', 'complete'), ...
    'VariableNames', {'corr_Prod', 'corr_WIP', 'corr_RAC', 'corr_InventSA', 'corr_Invent'});

disp('Correlations for overlapping 20-year period (1998–2018):');
disp(corrTable);
save('Empirical\tables\CorrOldNewSample.mat', 'corrTable')

%% Plot all variables, both samples
clear;
load("Empirical\oil_data.mat")
load('Empirical\BB22_data.mat')

% Define the cutoff date: from Braun & Brüggemann replication code
cutoff_date = datetime(1978, 10-2, 1);

% Filter the data
oil_data_filtered = oil_data(oil_data.Date >= cutoff_date, :);
data_BB22_filtered = data_BB22(data_BB22.date >= cutoff_date, :);

% Create the figure
figure;
subplot(2,2,1);
hold on;
plot(oil_data_filtered.Date, oil_data_filtered.prod, 'LineWidth', 1.5, 'DisplayName', 'New sample');
plot(data_BB22_filtered.date, data_BB22_filtered.Prod, 'LineWidth', 1, 'DisplayName', 'BB22 sample');
title('World Oil Production');
legend show;
hold off;

subplot(2,2,2);
hold on;
plot(oil_data_filtered.Date, oil_data_filtered.rea, 'LineWidth', 1.5, 'DisplayName', 'New sample');
plot(data_BB22_filtered.date, data_BB22_filtered.WIP, 'LineWidth', 1, 'DisplayName', 'BB22 sample');
title('Industrial Production');
legend show;
hold off;

subplot(2,2,3);
hold on;
plot(oil_data_filtered.Date, oil_data_filtered.rpo, 'LineWidth', 1.5, 'DisplayName', 'New sample');
plot(data_BB22_filtered.date, data_BB22_filtered.RAC, 'LineWidth', 1, 'DisplayName', 'BB22 sample');
title('Real Price of Oil');
legend show;
hold off;

subplot(2,2,4);
hold on;
plot(oil_data_filtered.Date, oil_data_filtered.oecdcrude_sa, 'LineWidth', 1.5, 'DisplayName', 'New sample');
plot(data_BB22_filtered.date, data_BB22_filtered.inventsa, 'LineWidth', 1, 'DisplayName', 'BB22 sample');
title('OECD oil Inventories, s.a.');
legend show;
hold off;

% Save
saveas(gcf, 'Empirical/figures/samples.fig');
