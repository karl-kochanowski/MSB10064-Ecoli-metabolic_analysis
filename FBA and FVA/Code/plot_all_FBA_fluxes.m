%% 2021-01-03 Karl Kochanowski

% plot all non-zero fluxes

function plot_all_FBA_fluxes(res)

%% basic visualization of FBA data
%% 1. number of fluxes that are never zero, and number of fluxes that are not
%zero at least in one condition
ix_non_zero = find(res.FBA.flux ~=0);
non_zero_logical = zeros(size(res.FBA.flux));
non_zero_logical(ix_non_zero) = 1;

non_zero_logical_sum = sum(non_zero_logical,2);

nr1 = length(find(non_zero_logical_sum==16));
nr2 = length(find(non_zero_logical_sum));
figure('Name','number of fluxes that are always or at least once non-zero');
bar([nr1,nr2]);
set(gca,'FontSize',16,'XTick',[1 2],'XTickLabel',{'always non-zero','>1 non-zero'});
xtickangle(-45);

%% 2. plot all fluxes
% to simplify things, only consider those fluxes that are non-zero in all
% conditions
ix_consider_flux = find(non_zero_logical_sum==16);

% one visualization issue: fluxes that change direction are difficult to
% plot in log scale - do not show these here
flux_sign = zeros(length(ix_consider_flux),16);
flux_original = res.FBA.flux(ix_consider_flux,:);
flux_sign(find(flux_original>0)) = 1;
flux_sign(find(flux_original<0)) = -1;
flux_sign_sum = sum(flux_sign,2);
ix_consider_flux_adj = find(abs(flux_sign_sum) == 16);

% normalize catabolic and anabolic fluxes separately
flux_cat = flux_original(ix_consider_flux_adj,1:8)./repmat(flux_original(ix_consider_flux_adj,1),1,8);
flux_ana = flux_original(ix_consider_flux_adj,9:end)./repmat(flux_original(ix_consider_flux_adj,9),1,8);

figure('Name','FBA data normalized');
% catabolic limitation
subplot(1,2,1),imagesc(log2(flux_cat),[-2 2]);
colormap(redbluecmap(9));
set(gca,'XTick',[],'YTick',[]);
title('catabolic lim','FontWeight','normal','FontSize',12);
hold on;
colorbar('northoutside');
% anabolic limitation
subplot(1,2,2),imagesc(log2(flux_ana),[-2 2]);
colormap(redbluecmap(9));
hold on;
set(gca,'XTick',[],'YTick',[]);
title('anabolic lim','FontWeight','normal','FontSize',12);
colorbar('northoutside');

end