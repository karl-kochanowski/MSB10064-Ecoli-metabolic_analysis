%% 2021-01-02 Karl Kochanowski

% QC of FBA
% 1) compare FBA and 13C fluxes for selected reactions in CCM
% 2) any unexpected secretion products?
% 3) plot all non-zero FBA fluxes
% 4) Pearson correlation coefficient with growth rate
% 5) FVA results

function FBA_QC(res)
%% QC 1: compare FBA and 13C fluxes for selected reactions in CCM
plot_comparison_13C_FBA(res.FBA.flux(res.ixFluxFBA,:),res.data.flux_names,res.data.flux_mean,res.data.flux_std);

%% QC 2: identify secretion products
plot_FBA_secretion_products(res);

%% QC 3: plot all non-zero FBA fluxes
plot_all_FBA_fluxes(res);

%% QC 4: plot Pearson correlation coefficient with growth rate
plot_flux_growth_correlation(res);

%% QC 5: plot FVA results
plot_FVA_results(res)
end