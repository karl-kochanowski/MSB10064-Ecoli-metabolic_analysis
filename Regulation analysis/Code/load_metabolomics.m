%% Karl Kochanowski

% load untargeted metabolomics data

function data = load_metabolomics()

%% load metabolite data
load 'metabolome data.mat';
% save full data series including the metadata (e.g. metabolite KEGG ID
% etc)
data.metabolome = metabolome_FIA;
%sort conditions to match the flux/proteomics data
data.metabolome.sorted.strainIx = [1 2 5 4 3 8 7 6 9 10 16 15 14 13 12 11];
% extract growth rate of metabolomics data
data.metabolome.sorted.mueMet = data.metabolome.sorted.averageMetadata(data.metabolome.sorted.strainIx,2);


end