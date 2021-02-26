%% 2020-06-18 Karl Kochanowski

% plot examples regulation analysis protein plus enzyme saturation

function plot_examples_regulation_analysis_saturation(data,ix)


figure('Name','regulation analysis protein-flux-metabolites');
runIx = 1;
for i = 1:length(ix)
    
    x = data.regulationAnalysis.enzymeSaturation.flux_data(ix(i),:);
    y = data.regulationAnalysis.enzymeSaturation.protein_data(ix(i),:);
    z = data.regulationAnalysis.enzymeSaturation.metabolite_data(ix(i),:);
    minval = min([x]);
    maxval = max([x]);
    
    rho = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(ix(i),:);
    rho_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(ix(i),:);
    
    gene = data.protein.gene(data.regulationAnalysis.enzymeSaturation.indices(ix(i),2));
    
    % first plot: protein
    subplot(length(ix),3,runIx),...
        hold on;
        %% add Standard deviation of Fit glc limitation
        ub = [minval maxval].*(rho(1)+rho_SE(1).*sqrt(8));
        lb = [minval maxval].*(rho(1)-rho_SE(1).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[0.8 0.8 1],'EdgeColor','none')
        
        %% add standard deviation of Fit glut lim
        ub = [minval maxval].*(rho(2)+rho_SE(2).*sqrt(8));
        lb = [minval maxval].*(rho(2)-rho_SE(2).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[1 0.8 0.8],'EdgeColor','none')
        
        line([minval maxval],[minval maxval],'Color','k','LineWidth',1);
        
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

    %plot 2: metabolite vs flux
    subplot(length(ix),3,runIx+1),...
        
        hold on;
        %% add Standard deviation of Fit glc limitation
        ub = [minval maxval].*(rho(3)+rho_SE(3).*sqrt(8));
        lb = [minval maxval].*(rho(3)-rho_SE(3).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[0.8 0.8 1],'EdgeColor','none')
        
        %% add standard deviation of Fit glut lim
        ub = [minval maxval].*(rho(4)+rho_SE(4).*sqrt(8));
        lb = [minval maxval].*(rho(4)-rho_SE(4).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[1 0.8 0.8],'EdgeColor','none')
        
        line([minval maxval],[minval maxval],'Color','k','LineWidth',1);
        
        line([minval maxval],[minval maxval].*rho(3),'Color','b','LineWidth',1);
        line([minval maxval],[minval maxval].*rho(4),'Color','r','LineWidth',1);
        
        plot(x(1:8),z(1:8),'ob','MarkerFaceColor','b');
        plot(x(9:end),z(9:end),'or','MarkerFaceColor','r');
        axis([-inf inf -inf inf]);
        box on;
    set(gca,'FontSize',12);
    t = gene{1};
    t2 = strcat(upper(t(1)),t(2:end));
    title(t2,'FontWeight','Normal','FontSize',16);
    
    %
    r1 = strcat('\rho_S=',num2str(rho(3),2));
    text(minval+abs(minval).*0.1,minval*rho(3),r1,'FontSize',12,'Color','b');
    r2 = strcat('\rho_S=',num2str(rho(4),2));
    text(minval+abs(minval).*0.1,minval*rho(4),r2,'FontSize',12,'Color','r');    

    %plot 2: protein+metabolite vs flux
    subplot(length(ix),3,runIx+2),...
        
        hold on;
        %% add Standard deviation of Fit glc limitation
        ub = [minval maxval].*(rho(5)+rho_SE(5).*sqrt(8));
        lb = [minval maxval].*(rho(5)-rho_SE(5).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[0.8 0.8 1],'EdgeColor','none')
        
        %% add standard deviation of Fit glut lim
        ub = [minval maxval].*(rho(6)+rho_SE(6).*sqrt(8));
        lb = [minval maxval].*(rho(6)-rho_SE(6).*sqrt(8));
        patch([minval maxval maxval minval],[ub(1) ub(2) lb(2) lb(1)],[1 0.8 0.8],'EdgeColor','none')
        
        line([minval maxval],[minval maxval],'Color','k','LineWidth',1);
        
        line([minval maxval],[minval maxval].*rho(5),'Color','b','LineWidth',1);
        line([minval maxval],[minval maxval].*rho(6),'Color','r','LineWidth',1);
        
        plot(x(1:8),y(1:8)+z(1:8),'ob','MarkerFaceColor','b');
        plot(x(9:end),y(9:end)+z(9:end),'or','MarkerFaceColor','r');
        axis([-inf inf -inf inf]);
        box on;
    set(gca,'FontSize',12);
    t = gene{1};
    t2 = strcat(upper(t(1)),t(2:end));
    title(t2,'FontWeight','Normal','FontSize',16);
    
    %
    r1 = strcat('\rho_P_+_S=',num2str(rho(5),2));
    text(minval+abs(minval).*0.1,minval*rho(5),r1,'FontSize',12,'Color','b');
    r2 = strcat('\rho_P_+_S=',num2str(rho(6),2));
    text(minval+abs(minval).*0.1,minval*rho(6),r2,'FontSize',12,'Color','r');     
    
    runIx = runIx + 3;
end


end
