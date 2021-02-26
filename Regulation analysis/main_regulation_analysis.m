%% Karl Kochanowski

%this main function calls all scripts to perform regulation analysis
% - step 1: load all data (fluxes, proteomics, metabolomics
% - step 2: regulation analysis I (flux-protein)
% - step 3: regulation analysis II (flux-protein-metabolites)

function data = main_regulation_analysis()
addpath('.\Data');
addpath('.\Code');

%% load flux, proteomics, metabolomics data
% load proteomics
data = load_proteomics();
% load flux data
tmp = load_fluxes();
data.flux = tmp.flux;
% load metabolomics data
tmp = load_metabolomics();
data.metabolome = tmp.metabolome;

%% perform regulation analysis
%% part 1: relate protein and flux changes
data = regulation_analysis_1_proteinOny(data);
% plot distribution of regulation coefficients
plot_regulation_analysis_1_distribution(data);
% plot nine examples 
plot_examples_regulation_analysis(data);

% generate export tables for supplementary files
data = extract_regulation_analysis_1(data);

%% part 2: relate protein, flux, and metabolite changes
data = regulation_analysis_2_protein_and_substrates(data);
% plot distribution of regulation coefficients
plot_regulation_analysis_2_distribution(data);
% plot examples 1: main figure
ix = [8,23,89]; % indices of example proteins in main figure 3C
plot_examples_regulation_analysis_saturation(data,ix);
% plot examples 2: arginine pathway examples
ix = [9,8,18,33]; % indices of example proteins in arginine pathway
plot_examples_regulation_analysis_saturation(data,ix);

% generate export tables for supplementary files
data = extract_regulation_analysis_2(data);
end