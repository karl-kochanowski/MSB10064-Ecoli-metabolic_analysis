%% 2021-01-02 Karl Kochanowski

% Quality control: compare FBA flux estimates with 13C-flux data

function plot_comparison_13C_FBA(fluxesFBA,names,fluxesMeasured_mean,fluxesMeasured_std)

%some fluxes are defined in negative direction (e.g. glucose uptake rate)
revSign = ones(length(names),1);
revSign(2,1) =-1;
revSign(25,1) =-1;
revSign(26,1) =-1;
fluxesFBA = fluxesFBA.*repmat(revSign,1,16);

figure('Name','Comparison 13C vs FBA');
mueFBA = fluxesFBA(1,:);
mueData = fluxesMeasured_mean(1,:);

for i=2:size(fluxesFBA,1)
    subplot(3,9,i-1),errorbar(mueData(1:8),fluxesMeasured_mean(i,1:8),fluxesMeasured_std(i,1:8),'-db','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1);
    hold on;
    subplot(3,9,i-1),errorbar(mueData(9:16),fluxesMeasured_mean(i,9:16),fluxesMeasured_std(i,9:16),'-dr','MarkerFaceColor','k','MarkerSize',6,'LineWidth',1);
    
    subplot(3,9,i-1),plot(mueFBA(1:8),fluxesFBA(i,1:8),'ob','MarkerSize',6,'LineWidth',1);
    subplot(3,9,i-1),plot(mueFBA(9:16),fluxesFBA(i,9:16),'or','MarkerSize',6,'LineWidth',1);
    
    
    set(gca,'FontSize',10);
    title(names(i,1),'FontSize',10,'FontWeight','Normal');
    axis([0 1 0 inf]);
end

end