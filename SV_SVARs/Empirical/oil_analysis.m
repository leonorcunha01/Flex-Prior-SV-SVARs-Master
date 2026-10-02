%% ========================================================================
%%  Analysis Oil Market Data
%  Computes IRFs, FEVDs and HDs.
%  ========================================================================

%% --- IRFs ---
hor = 18;
[modelHS.IRFs] = IRF2(modelHS.alpha, modelHS.B0_norm, p, hor, n);
plot_irf(modelHS.IRFs)
   saveas(gcf,sprintf('Empirical/figures/IRFs_%s_modelHS.fig', samplename));

[modelIG.IRFs] = IRF2(modelIG.alpha, modelIG.B0_norm, p, hor, n);
plot_irf(modelIG.IRFs)
    saveas(gcf,sprintf('Empirical/figures/IRFs_%s_modelIG.fig', samplename));

plot_irf_compare(modelHS.IRFs, modelIG.IRFs, ...
    {'Oil Prod.', 'Ind. Prod.', 'Price of Oil', 'Inventories'}, ...
    {'Supply', 'Econ. Activ.', 'Cons. Demand', 'Inv. Demand'})
    saveas(gcf,sprintf('Empirical/figures/IRFs_%s_comp.fig', samplename));

%% --- FEVDs ---
modelHS.fevd = FEVDs_drawbydraw(modelHS.alpha, modelHS.B0, p, 60, n, squeeze(modelHS.h(:,end,:)), modelHS.phi, squeeze(modelHS.Sigh(:,end,:)));
plot_fevd(squeeze(median(modelHS.fevd)))
    saveas(gcf,sprintf('Empirical/figures/FEVDs_%s_modelHS.fig', samplename));

modelIG.fevd = FEVDs_drawbydraw(modelIG.alpha, modelIG.B0, p, 60, n, squeeze(modelIG.h(:,end,:)), modelIG.phi, squeeze(modelIG.Sigh(:,end,:)));
plot_fevd(squeeze(median(modelIG.fevd)))
    saveas(gcf,sprintf('Empirical/figures/FEVDs_%s_modelIG.fig', samplename));

% Tables



%% --- HDs ---
[modelHS.hd, modelHS.y_hd] = HD_drawbydraw(modelHS.alpha, modelHS.B0, p, T-4*p, n, modelHS.struct_shock);
[modelIG.hd, modelIG.y_hd] = HD_drawbydraw(modelIG.alpha, modelIG.B0, p, T-4*p, n, modelIG.struct_shock);


% Plot and save
Dates_forplot = Dates(4*p+1:end);
close all
plot_HD(squeeze(median(modelHS.hd)), Dates_forplot)
figHDs = findall(0, 'Type', 'figure');
figHDs = flipud(figHDs);
for i = 1:length(figHDs)
    fig = figHDs(i);
    filename = fullfile('Empirical/figures/', sprintf('HD_%d_set_%s_modelHS.fig', i, samplename));
    savefig(fig, filename);
end

close all
plot_HD(squeeze(median(modelIG.hd)), Dates_forplot)
figHDs = findall(0, 'Type', 'figure');
figHDs = flipud(figHDs);
for i = 1:length(figHDs)
    fig = figHDs(i);
    filename = fullfile('Empirical/figures/', sprintf('HD_%d_set_%s_modelIG.fig', i, samplename));
    savefig(fig, filename);
end




