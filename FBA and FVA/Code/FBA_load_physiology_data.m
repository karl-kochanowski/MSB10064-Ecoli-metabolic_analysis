%% 2021-01-01 Karl Kochanowski

% load all physiology data needed to constrain the FBA model

function res = FBA_load_physiology_data(res)

%% load physiology, flux ratios and absolute fluxes

%% get major exchange rates
[phys,~,~]=xlsread('flux and physiology data.xlsx','exchange rates I','C5:R8');
[~,physName,~]=xlsread('flux and physiology data.xlsx','exchange rates I','A5:A8');

%% get minor secretion rates
[phys_minor,~,~]=xlsread('flux and physiology data.xlsx','exchange rates II','H7:W14');
[~,physName_minor,~]=xlsread('flux and physiology data.xlsx','exchange rates II','B7:B14');

phys_minor(find(phys_minor<1))=0; %set all super-low secretion rates to zero
phys_minor=phys_minor./1000; %convert into the right unit [mmol/h/gCDW]

% only consider those minor physiology fluxes that at least once reach 10
% micromol/h/gCDW
max_value = max(phys_minor,[],2);
ix_phys_minor_above_threshold = find(max_value > 0.01);

%% get flux ratios
[fluxRatio,~,~]=xlsread('flux and physiology data.xlsx','flux ratios','D5:S13');
[~,fluxRatioName,~]=xlsread('flux and physiology data.xlsx','flux ratios','B5:B13');

%% get absolute fluxes
[absFlux,~,~]=xlsread('flux and physiology data.xlsx','absolute fluxes','E5:T31');
[absFlux_std,~,~]=xlsread('flux and physiology data.xlsx','absolute fluxes','E35:T61');

%manually get sole pyk flux (in E.coli, glucose uptake through the PTS is
%coupled to PEP->PYR conversion as well, and the 13C flux data only report to
%PEP->PYR flux)
rel = absFlux_std(13,:)./absFlux(13,:);
absFlux(13,:) = absFlux(13,:)-absFlux(2,:);
absFlux_std(13,:) = absFlux(13,:).*rel;

[~,absFluxName,~]=xlsread('flux and physiology data.xlsx','absolute fluxes','A5:A31');

%% save data
res.data.physiology = phys;
res.data.physiology_names = physName;
res.data.minor_physiology = phys_minor(ix_phys_minor_above_threshold,:);
res.data.minor_physiology_names = physName_minor(ix_phys_minor_above_threshold);

res.data.fluxRatio = fluxRatio;
res.data.fluxRatio_names = fluxRatioName;
res.data.flux_mean = absFlux;
res.data.flux_std = absFlux_std;
res.data.flux_names = absFluxName;

%% load flux IDs to be able to match FBA and 13C fluxes
[~,fluxIDs,~]=xlsread('match IDs of 13C and FBA flux.xlsx',1,'A2:B28');

for i=1:size(fluxIDs,1)
   tmpIx=find(strcmpi(res.model.rxns,fluxIDs(i,2))); 
   res.ixFluxFBA(i,1)=tmpIx(1,1);   
end




end