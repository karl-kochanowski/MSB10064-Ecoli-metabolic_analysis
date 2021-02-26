%% Karl Kochanowski

% plot regulation analysis I - distribution of regulation coefficients for
% unique metabolic reactions

function plot_regulation_analysis_1_distribution(data)

regulation_coefs = data.regulationAnalysis.fluxVSprotein.uniqueReactions_regCoeff;

%% plot distribution of unique reactions
ix_relevant_cat = find(~isnan(regulation_coefs(:,1)));
ix_relevant_ana = find(~isnan(regulation_coefs(:,2)));

x = [-1.9:0.2:2];

% for plotting purposes: regulation coefs that lie outside the boundaries:
% set to boundaries
regulation_coefs(find(regulation_coefs<-2)) = -2;
regulation_coefs(find(regulation_coefs>2)) = 2;


[n_cat,x_cat] = histcounts(regulation_coefs(ix_relevant_cat,1),[-2:0.2:2]);
[n_ana,x_ana] = histcounts(regulation_coefs(ix_relevant_ana,2),[-2:0.2:2]);

%fraction of reactions with regulatory coefficients between 0.5 and 1.5%
% catabolic limitation
it1 = find(regulation_coefs(ix_relevant_cat,1)>0.5);
it2 = find(regulation_coefs(ix_relevant_cat,1)<1.5);
i3 = intersect(it1,it2);
cat_Fraction = length(i3)./sum(n_cat(:));

% anabolic limitation
it1 = find(regulation_coefs(ix_relevant_ana,2)>0.5);
it2 = find(regulation_coefs(ix_relevant_ana,2)<1.5);
i3 = intersect(it1,it2);
ana_Fraction = length(i3)./sum(n_ana(:));

%% plot distributions
figure('name','regulatory analysis - protein-flux - averaged protein reg coefs');
% catabolic limitation
subplot(2,3,1),
patch([0.5 1.5 1.5 0.5],[0 0 0.25 0.25],[0.8 0.8 0.8],'EdgeColor','none');
hold on;
bar(x,n_cat./sum(n_cat),'b');
box on;

axis([-2.1 2.1 0 0.25]);
set(gca,'FontSize',16);
set(gca,'XTick',[-2:1:2]);
hold on;

text(0.8,0.2,strcat(num2str(cat_Fraction.*100,2),'%'),'FontSize',16,'Color','w');
title('catabolic lim','FontWeight','normal','Color','b','FontSize',16);

% anabolic limitation
subplot(2,3,4),
patch([0.5 1.5 1.5 0.5],[0 0 0.25 0.25],[0.8 0.8 0.8],'EdgeColor','none');
hold on;
bar(x,n_ana./sum(n_ana),'r');
box on;

axis([-2.1 2.1 0 0.25]);
set(gca,'FontSize',16);
set(gca,'XTick',[-2:1:2]);

text(0.8,0.2,strcat(num2str(ana_Fraction.*100,2),'%'),'FontSize',16,'Color','w');
title('anabolic lim','FontWeight','normal','Color','r','FontSize',16);


end