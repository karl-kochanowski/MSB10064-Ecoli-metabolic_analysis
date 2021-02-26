%% Karl Kochanowski

%plot 9 examples of protein-flux regulation analyses

function plot_examples_regulation_analysis(data)

% indices of reactions of interest
ix = [109,272,84,36,13,25,23,188,166];

figure('Name','regulation analysis protein-flux');

for i = 1:length(ix)
    x = data.regulationAnalysis.fluxVSprotein.flux_data(ix(i),:);
    y = data.regulationAnalysis.fluxVSprotein.protein_data(ix(i),:);
    
    minval = min([x,y]);
    maxval = max([x,y]);
   
    rho = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient(ix(i),:);
    rho_SE = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient_SE(ix(i),:);
    
    gene = data.protein.gene(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(ix(i),2));
    
    subplot(3,3,i),...
        line([minval maxval],[minval maxval],'Color','k','LineWidth',1);
        hold on;
        %% add Standard deviation of Fit catabolic limitation
        ub = [minval maxval].*(rho(1)+rho_SE(1).*sqrt(8));
        lb = [minval maxval].*(rho(1)-rho_SE(1).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[0.8 0.8 1],'EdgeColor','none')
        
        %% add standard deviation of Fit anabolic lim
        ub = [minval maxval].*(rho(2)+rho_SE(2).*sqrt(8));
        lb = [minval maxval].*(rho(2)-rho_SE(2).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[1 0.8 0.8],'EdgeColor','none')
        
        
        line([minval maxval],[minval maxval].*rho(1),'Color','b','LineWidth',1);
        line([minval maxval],[minval maxval].*rho(2),'Color','r','LineWidth',1);
        
        plot(x(1:8),y(1:8),'ob','MarkerFaceColor','b');
        plot(x(9:end),y(9:end),'or','MarkerFaceColor','r');
        axis([-inf inf -inf inf]);
        box on;
    set(gca,'FontSize',12);
    t = gene{1};
    t2 = strcat(upper(t(1)),t(2:end));
    title(t2,'FontWeight','Normal','FontSize',16);
    
    %
    r1 = strcat('\rho_P=',num2str(rho(1),2));
    text(minval+abs(minval).*0.1,minval*rho(1),r1,'FontSize',12,'Color','b');
    r2 = strcat('\rho_P=',num2str(rho(2),2));
    text(minval+abs(minval).*0.1,minval*rho(2),r2,'FontSize',12,'Color','r');
    
end

end