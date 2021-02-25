%% 2020-12-21 Karl Kochanowski

% perform FVA of FBA-estimated metabolic fluxes 

function res = FBA_perform_FVA_analysis(res,flag_all_fluxes)

% get number of experiments
ix = length(res.FBA.model);

%% run FVA
for i = 1:ix
    i
    model = res.FBA.model{i};
    nr_reactions = length(model.rxns);
    
    %options: either run FVA on all reactions, regardless of whether they
    %carry any flux in the given condition (== 1), or only on those
    %reactions that carry flux in a given condition ( ~= 1)
    if(flag_all_fluxes == 1)
        ix_reactions = 1:nr_reactions;
    else
        ix_reactions = find(res.FBA.flux(:,i) ~=0);
    end

    %to avoid numerical issues, relax lower bound of objective function and optimality
    %constraints a tiny bit
    model.lb(748,1)=model.lb(748,1).*0.9999;
    opt = 99.99;
    %opt = 100;
    
    [minFlux maxFlux] = fluxVariability(model,'optPercentage',opt,'osenseStr','max','rxnNameList',model.rxns(ix_reactions),'allowLoops','fastSNP','printLevel',1);
    
    
    res.FBA.minFVA(:,i) = zeros(nr_reactions,1);
    res.FBA.minFVA(:,i) = NaN;
    res.FBA.minFVA(ix_reactions,i) = minFlux;
    
    res.FBA.maxFVA(:,i) = zeros(nr_reactions,1);
    res.FBA.maxFVA(:,i) = NaN;
    res.FBA.maxFVA(ix_reactions,i) = maxFlux;
    
    

    
    
end

end