%% 2021-01-15 Karl Kochanowski

% extract regulation analysis data 1 (protein only) to generate tables

function data = extract_regulation_analysis_1(data)

%% Table 1: all individual protein-reaction pairs
n = {'Number','Reaction_Name','Protein_Name','Catabolic_Limitation','Catabolic_Limitation_SE','Anabolic_Limitation','Anabolic_Limitation_SE'};

nr_reactions = size(data.regulationAnalysis.fluxVSprotein.regulatory_coefficient,1);
reactions = [1:nr_reactions]';
reaction_name = data.flux.model.rxns(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(:,1));
protein_name = data.protein.gene(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(:,2));
regcoef_cat = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient(:,1);
regcoef_cat_SE = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient_SE(:,1);
regcoef_ana = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient(:,2);
regcoef_ana_SE = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient_SE(:,2);

data.regulationAnalysis.fluxVSprotein.export_all_proteins = table(reactions,reaction_name,protein_name,regcoef_cat,regcoef_cat_SE,regcoef_ana,regcoef_ana_SE,'VariableNames',n);

%% Table 2: all unique reactions
n = {'Number','Reaction_Name','Catabolic_Limitation','Anabolic_Limitation'};
nr_reactions = size(data.regulationAnalysis.fluxVSprotein.uniqueReactions_ix,1);
reactions = [1:nr_reactions]';
reaction_name = data.flux.model.rxns(data.regulationAnalysis.fluxVSprotein.uniqueReactions_ix(:,1));
regcoef_cat = data.regulationAnalysis.fluxVSprotein.uniqueReactions_regCoeff(:,1);
regcoef_ana = data.regulationAnalysis.fluxVSprotein.uniqueReactions_regCoeff(:,2);

data.regulationAnalysis.fluxVSprotein.export_unique_reactions = table(reactions,reaction_name,regcoef_cat,regcoef_ana,'VariableNames',n);

end