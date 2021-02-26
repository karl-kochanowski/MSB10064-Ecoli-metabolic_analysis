%% Karl Kochanowski

% plot regulation analysis 2 (protein plus substrate effect) - distribution of regulation coefficients for
% unique metabolic reactions

function plot_regulation_analysis_2_distribution(data)
% plotting order and parameters
ix = [1 3 5 2 4 6];
title_text = {'only proteins','only substrates','both'};
title_text = [title_text,title_text];

figure('Name','distribution of regulation coefficients - protein plus substrate effect');

%% plot regulation coefficients

x = [-1.9:0.2:1.9];
regulation_coefs = data.regulationAnalysis.enzymeSaturation.uniqueReactions_regCoeff;
% for plotting purposes: regulation coefs that lie outside the boundaries:
% set to boundaries
regulation_coefs(find(regulation_coefs<-2)) = -2;
regulation_coefs(find(regulation_coefs>2)) = 2;

for i = 1:length(ix)
    if(i<=3)
        c = 'b';
    else
        c = 'r';
    end
    % get 
    ix_non_nan = ~isnan(regulation_coefs(:,ix(i)));
    % calculate distribution of regulation coefficients
    [n,~] = histcounts(regulation_coefs(ix_non_nan,ix(i)),[-2:0.2:2]);
    
    % calculate fraction of reactions that lie between 0.5 and 1.5
    above_lb = find(regulation_coefs(ix_non_nan,ix(i)) > 0.5);
    below_ub = find(regulation_coefs(ix_non_nan,ix(i)) < 1.5);
    nr_overlap = length(intersect(above_lb,below_ub));
    fraction = 100.*nr_overlap./sum(n);
    
    subplot(2,3,i),...
    
    patch([0.5 1.5 1.5 0.5],[0 0 0.4 0.4],[0.8 0.8 0.8],'EdgeColor','none');
    hold on;
    bar(x,n./sum(n),c);
    axis([-2.1 2.1 0 0.4]);
    set(gca,'FontSize',16);
    set(gca,'XTick',[-2:1:2]);
    title(strcat(title_text{i},'(',num2str(fraction,2),'%)'),'FontWeight','Normal');
    box on;
    
end
end