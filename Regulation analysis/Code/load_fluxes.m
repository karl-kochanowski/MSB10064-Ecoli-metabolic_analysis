%% Karl Kochanowski

% 1: load FBA flux data
% 2: interpolate flux data to match growth rates of proteomics/metabolomics

function data = load_fluxes()

%% load fluxome data
load 'FBA flux data.mat';

% growth rate range for interpolation
x_range = [0.2:0.001:1.05];

data.flux.limType = [repmat({'catabolic'},1,length(x_range)),repmat({'anabolic'},1,length(x_range))];
data.flux.fluxValues = res.FBA.flux;
data.flux.model = res.model;
data.flux.mue_interpolated = [x_range,x_range];

%% interpolate fluxes
% only consider fluxes that are non-zero in all titration steps (16
% conditions in total)
mueFlux = data.flux.fluxValues(7,:); %growth rate as inferred by FBA
for i = 1:size(data.flux.fluxValues,1)
    fluxData = data.flux.fluxValues(i,:);
    if(max(abs(fluxData))>0)&&(length(find(fluxData~=0))==16)
        flux_interpol_glc = interp1(mueFlux(1:8),fluxData(1:8),x_range,'linear','extrap');
        flux_interpol_glu = interp1(mueFlux(9:end),fluxData(9:end),x_range,'linear','extrap');     
        flag_flux = 1;
    else
        flux_interpol_glc = 0.*x_range;
        flux_interpol_glu = 0.*x_range;
        flag_flux = 0;
    end
    data.flux.fluxValues_interpolated(i,:)=[flux_interpol_glc,flux_interpol_glu];    
    data.flux.non_zero_flux(i,1) = flag_flux;
end

%% load reaction details from iJO1366 model
[~,FBA_reactions,~]=xlsread('iJO1366_reaction_details.xlsx','B2:B2585');
data.flux.model.reactionDetails = FBA_reactions;


end