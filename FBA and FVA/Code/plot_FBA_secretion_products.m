%% 2021-01-03 Karl Kochanowski

% identify and plot secretion products

function plot_FBA_secretion_products(res)
% indices of exchange fluxes
ix_secr = [9:332,2584];
flux = res.FBA.flux(ix_secr,:);

% which of these fluxes are non-zero at least once?
flux_max = max(abs(flux),[],2);
flux_max_nonzero = find(flux_max > 0.001);

figure('Name','exchange fluxes');
for i = 1:length(flux_max_nonzero)
    semilogy(i-0.1,abs(flux(flux_max_nonzero(i),1:8)),'ok','MarkerFaceColor','b');
    hold on;
    semilogy(i+0.1,abs(flux(flux_max_nonzero(i),9:16)),'ok','MarkerFaceColor','r');
    flux_name_tmp = res.model.rxns(ix_secr(flux_max_nonzero(i)));
    flux_names(i,1) = strrep(flux_name_tmp,'_','-');
end
set(gca,'FontSize',12,'XTick',[1:length(flux_max_nonzero)],'XTickLabel',flux_names);
xtickangle(-45);
axis([-inf inf 0.001 100]);
xlabel('exchange reactions');
ylabel('absolute flux [mmol/h/gCDW]');

end