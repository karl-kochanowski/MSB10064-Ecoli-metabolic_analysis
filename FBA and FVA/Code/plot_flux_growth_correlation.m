%% 2021-01-03 Karl Kochanowski

% calculate Pearson correlation between flux and growth rate
% - only consider fluxes that are non-zero in at least 50% of conditions

function plot_flux_growth_correlation(res)

% get FBA-extracted growth rate
mue = res.FBA.flux(7,:);

% get all fluxes except growth rate
ix_fluxes = setdiff(1:length(res.model.rxns),7);
fluxes = res.FBA.flux(ix_fluxes,:);

%% calculate correlation coefficient
rho_all = 1000.*ones(size(fluxes,1),1);
for i = 1:length(rho_all)
    flux_tmp = fluxes(i,:);
    ix_non_zero = find(flux_tmp ~=0);
    
    % only consider fluxes that are non-zero in at least 50% of conditions
    if(length(ix_non_zero) >= 8)
        % if all fluxes are negative (reaction defined in negative
        % direction), flip them
        if(max(flux_tmp(ix_non_zero))<= 0)
            flux_tmp = flux_tmp.*-1;
        end
        % calculate correlation coefficient
        rho_all(i,1) = corr(mue(ix_non_zero)',flux_tmp(ix_non_zero)');
    end
end
ix_relevant = find(rho_all ~= 1000);

%% plot distribution of correlation coefficients
figure('Name','correlation flux with growth rate');
h = histogram(rho_all(ix_relevant),[-1:0.05:1]);
set(gca,'FontSize',12);
xlabel('Pearson correlation coefficient [-]');
ylabel('number in bin');
axis([-1.05 1.05 0 450]);
h.FaceColor = [0 0 0];
h.FaceAlpha = 1;
end