%% 2021-01-15 Karl Kochanowski

% extract regulation analysis data 2 (protein and substrates) to generate tables

function data = extract_regulation_analysis_2(data)

%% Table 1: alpha estimate for each metabolite-protein-flux combination

n = {'Reaction_Name','Protein','Metabolite_Name','KEGGID','Mz','alpha'};
reaction_name = data.flux.model.rxns(data.regulationAnalysis.enzymeSaturation.connectionFluxProteinMetAlpha(:,1));
protein_name = data.protein.gene(data.regulationAnalysis.enzymeSaturation.connectionFluxProteinMetAlpha(:,2));
metIx = data.regulationAnalysis.enzymeSaturation.connectionFluxProteinMetAlpha(:,3);
metabolite_name = data.metabolome.raw.AnnotatedIonListNames(metIx);
keggid = data.metabolome.raw.AnnotatedIonListKEGG_ID(metIx);
mz = data.metabolome.raw.annotatedIonMz(metIx);
alpha = data.regulationAnalysis.enzymeSaturation.connectionFluxProteinMetAlpha(:,4);

data.regulationAnalysis.enzymeSaturation.export_alpha_estimates = table(reaction_name,protein_name,metabolite_name,keggid,mz,alpha,'VariableNames',n);

%% Table 2: regulation coefficients for all reaction-protein pairs
n = {'Reaction_Number','Reaction_Name','Protein','catabolic_rho_p','catabolic_rho_p_SE','anabolic_rho_p','anabolic_rho_p_SE','catabolic_rho_s','catabolic_rho_s_SE','anabolic_rho_s','anabolic_rho_s_SE','catabolic_rho_ps','catabolic_rho_ps_SE','anabolic_rho_ps','anabolic_rho_ps_SE'};

nr_reactions = size(data.regulationAnalysis.enzymeSaturation.regulatory_coefficient,1);
reactions = [1:nr_reactions]';
reaction_name = data.flux.model.rxns(data.regulationAnalysis.enzymeSaturation.indices(:,1));
protein_name = data.protein.gene(data.regulationAnalysis.enzymeSaturation.indices(:,2));
cat_r_p = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,1);
cat_r_p_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,1);
ana_r_p = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,2);
ana_r_p_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,2);
cat_r_s = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,3);
cat_r_s_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,3);
ana_r_s = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,4);
ana_r_s_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,4);
cat_r_ps = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,5);
cat_r_ps_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,5);
ana_r_ps = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,6);
ana_r_ps_SE = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(:,6);

data.regulationAnalysis.enzymeSaturation.export_regulation_coefs = table(reactions,reaction_name,protein_name,cat_r_p,cat_r_p_SE,ana_r_p,ana_r_p_SE,cat_r_s,cat_r_s_SE,ana_r_s,ana_r_s_SE,cat_r_ps,cat_r_ps_SE,ana_r_ps,ana_r_ps_SE,'VariableNames',n);

%% Table 3: mean regulation coefficients for all unique reactions
nr_reactions = size(data.regulationAnalysis.enzymeSaturation.uniqueReactions_ix,1);
reactions = [1:nr_reactions]';
reaction_name = data.flux.model.rxns(data.regulationAnalysis.enzymeSaturation.uniqueReactions_ix(:,1));
n = {'Reaction_Number','Reaction_Name','catabolic_rho_p','anabolic_rho_p','catabolic_rho_s','anabolic_rho_s','catabolic_rho_ps','anabolic_rho_ps'};
reg_coefs_unique = data.regulationAnalysis.enzymeSaturation.uniqueReactions_regCoeff;
data.regulationAnalysis.enzymeSaturation.export_regulation_coefs_unique_reactions = table(reactions,reaction_name,reg_coefs_unique(:,1),reg_coefs_unique(:,2),reg_coefs_unique(:,3),reg_coefs_unique(:,4),reg_coefs_unique(:,5),reg_coefs_unique(:,6),'VariableNames',n)
end