%% 2021-01-01 Karl Kochanowski

% perform FBA analysis as described in Kochanowski et al MSB
% for each condition
%   - add constraints based on physiology and 13C
%   - FBA step 1: maximize ATP production rate
%   - FBA step 2: minimize sum of fluxes given the fixed ATP production
%   rate

function res = FBA_perform_FBA_analysis(res)



% loop through each condition
for i=1:16
    i
    % get model
    tmpModel = res.model;
    
    % get physiology constraints (in total 12 constraints
    physiology.names={'Ec_biomass_iJO1366_WT_53p95M','EX_glc(e)','EX_ac(e)','EX_akg(e)','EX_o2(e)'};
    physiology.names=[physiology.names,res.data.minor_physiology_names'];

    %with constrained oxygen uptake rate
    physiology.values = res.data.flux_mean([1,2,23,27,26],i)';
    physiology.values(1,2) = -1.*physiology.values(1,2);
    physiology.values(1,5) = -1.*physiology.values(1,5);
    physiology.values = [physiology.values,res.data.minor_physiology(:,i)'];
   
    if(i>10) %shut down glutamate dehydrogenase in NQ393 strain (where this enzyme was deleted)
        tmpModel = changeRxnBounds(tmpModel,{'GLUDy'},0,'l');
        tmpModel = changeRxnBounds(tmpModel,{'GLUDy'},0,'u');
    
    end
    % add 4 flux ratios as constraints
    tmpModel = addRatioReaction(tmpModel, {'G6PDH2r' 'PFK'}, [res.data.fluxRatio(1,i) 1-res.data.fluxRatio(1,i)]);
    tmpModel = addRatioReaction(tmpModel, {'MDH' 'PPC'}, [res.data.fluxRatio(9,i) 1-res.data.fluxRatio(9,i)]);
    tmpModel = addRatioReaction(tmpModel, {'PGI' 'EDD'}, [res.data.fluxRatio(2,i) 1-res.data.fluxRatio(2,i)]);
    tmpModel = addRatioReaction(tmpModel, {'ENO' 'PPCK'}, [res.data.fluxRatio(8,i) 1-res.data.fluxRatio(8,i)]);

    % run FBA
    FBA_result = runFBA(tmpModel,physiology);
    
    % save output
    if (~isempty(FBA_result.v))
       res.FBA.flux(:,i) = FBA_result.v;
       res.FBA.model{i,1} = FBA_result.model;
    else
        res.FBA.flux(:,i) = zeros(length(tmpModel.rxns),1);
    end    

end


end

%% two-step FBA: maximize ATP production, then minimize sum of fluxes using the first flux solution as upper bound
function res=runFBA(model,constraints)

%set constraints for glc uptake, acetate secretion, aKG secretion, growth
%rate, inferred oxygen uptake rate
for i=1:5
    model = changeRxnBounds(model,constraints.names(i),constraints.values(i)-0.05*abs(constraints.values(i)),'l');
    model = changeRxnBounds(model,constraints.names(i),constraints.values(i)+0.05*abs(constraints.values(i)),'u');   
end
%set constraints for minor secretion products
for i=5:length(constraints.values)
    model = changeRxnBounds(model,constraints.names(i),constraints.values(i)-0.1*abs(constraints.values(i)),'l');
    model = changeRxnBounds(model,constraints.names(i),constraints.values(i)+0.1*abs(constraints.values(i)),'u');   
end

% change objective function flag
model.c(:) = 0; % set all reactions to 0
model.c(748,1) = 1; %optimize for ATP production rate

%step 1: maximize ATP production rate
resTmp = optimizeCbModel(model,'max',1e-6); 

%fix ATP production rate from first FBA
model.ub(748,1)=abs(resTmp.v(748,1));
model.lb(748,1)=abs(resTmp.v(748,1));

%step 2: minimize sum of fluxes (given the physiology and fixed ATP production
%rate)
[resTmp3,irrev] = minimizeModelFlux(model,'min',1e-6);
for i=1+length(model.rxns):length(irrev.rxns)-1
    name = irrev.rxns{i}(1:end-2);
    idx = strcmp(model.rxns,name);
    resTmp3.v(idx)=resTmp3.v(idx)-resTmp3.v(i);
end
resTmp3 = resTmp3.v(1:length(model.rxns));
% these fluxes get flipped, flip them back
resTmp3(164,:) = -resTmp3(164,:);
resTmp3(252,:) = -resTmp3(252,:);

%save output
res.v = resTmp3;
res.model = model;




end