%% 2021-01-01 Karl Kochanowski

% load E.coli genome-scale FBA model (iJO1366) and modify as described in
% Kochanowski et al MSB

function res = FBA_load_model

% load initial data
load iJO1366.mat;
modelInitial = iJO1366;

% get mapping metabolite KEGG_ID to reaction from original paper
[~,KEGGIDmapping,~] = xlsread('match metabolite IDs to KEGG IDs.xlsx',1,'A2:B1806');
modelInitial.metKEGGID = KEGGIDmapping;

%% modify original model as described in Kochanowski et al
% close glyoxylate shunt (based on 13C flux analysis data)
modelInitial = changeRxnBounds(modelInitial,{'ICL'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'ICL'},0,'u');

% close malic enzymes (based on 13C flux analysis data)
modelInitial = changeRxnBounds(modelInitial,{'ME2'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'ME2'},0,'u');
modelInitial = changeRxnBounds(modelInitial,{'ME1'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'ME1'},0,'u');

% close two bypass reactions in glycolysis
modelInitial = changeRxnBounds(modelInitial,{'F6PA'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'F6PA'},0,'u');
modelInitial = changeRxnBounds(modelInitial,{'FBA3'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'FBA3'},0,'u');

% close exchange metabolite reactions (metabolites not found to be
% secreted)
modelInitial = changeRxnBounds(modelInitial,{'EX_pyr(e)'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'EX_pyr(e)'},0,'u');
modelInitial = changeRxnBounds(modelInitial,{'EX_lac-D(e)'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'EX_lac-D(e)'},0,'u');
modelInitial = changeRxnBounds(modelInitial,{'EX_etoh(e)'},0,'l');
modelInitial = changeRxnBounds(modelInitial,{'EX_etoh(e)'},0,'u');
%
% enable reversible fumarate exchange between periplasm and cytosol through dctA
modelInitial = changeRxnBounds(modelInitial,{'FUMt2_2pp'},-1000,'l');

% add phenylpyruvate secretion
modelInitial = addExchangeRxn(modelInitial,{'phpyr[c]'},0,0);

%output: modified model
res.model = modelInitial;



end