%% 2021-01-01 Karl Kochanowski

% main analysis script to estimate genome-scale metabolic fluxes with
% constrained FBA

% part of Kochanowski et al MSB (2021)

% step 1: load model and modify as described in the manuscript
% step 2: load all physiology data needed to constrain the model
% step 3: perform FBA
% step 4: perform FVA
% step 5: plot results and QC

function res = main_analysis_gsFBA

%% add data and code paths
addpath('.\Data');
addpath('.\Code');

%% step 1: load and modify FBA model iJO1366
res = FBA_load_model;
%% step 2: load physiology data (used to constrain model)
res = FBA_load_physiology_data(res);
%% step 3: perform FBA
solver = 'gurobi'; %will be used for all problems
% change Cobra toolbox solver
changeCobraSolver(solver,'all');

res = FBA_perform_FBA_analysis(res);

%% step 4: perform FVA (CAUTION: THIS STEP TAKES A LOT OF TIME!)
flag_all_fluxes = 0; 
% 0: run FVA only on reactions with non-zero flux
% 1: run FVA on all reactions regardless of whether they carry flux or not
res = FBA_perform_FVA_analysis(res,flag_all_fluxes);

%% step 5: plot FBA results
FBA_QC(res);

end