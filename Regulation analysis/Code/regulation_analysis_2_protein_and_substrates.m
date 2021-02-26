%% Karl Kochanowski

% regulation analysis 2: relate flux, protein, and metabolite abundance changes

% perform regulation analysis as described in Kochanowski et al, MSB
% (2021)

function data = regulation_analysis_2_protein_and_substrates(data)

%% load metabolome growth rates and condition order
mueMetabolome = data.metabolome.sorted.mueMet;
strainIx = data.metabolome.sorted.strainIx;
mueFlux_interpolated = data.flux.mue_interpolated;

%% get fluxes at exactly the same growth rate as the metabolomics data
limType_flux = data.flux.limType;
ix_glc_flux = find(strcmpi(limType_flux,'catabolic'));
ix_glu_flux = find(strcmpi(limType_flux,'anabolic'));

tmp_glc = abs(repmat(mueFlux_interpolated(ix_glc_flux)',1,8)-repmat(mueMetabolome(1:8)',length(ix_glc_flux),1));
[~,ind_glc_overlap] = min(tmp_glc,[],1);
tmp_glu = abs(repmat(mueFlux_interpolated(ix_glu_flux)',1,8)-repmat(mueMetabolome(9:end)',length(ix_glu_flux),1));
[~,ind_glu_overlap] = min(tmp_glu,[],1);
ind_glu_overlap = ind_glu_overlap+length(ix_glu_flux);

%% initialize tables to gather the data
indx_flux_protein_met = [];
regulatory_coefficient = [];
regulatory_coefficient_SE = [];
indx_flux_protein_met_alpha = [];
exitflag = [];

flux_summary = [];
protein_summary = [];
metabolite_summary = [];

%% exclude co-factors such as proton, h2o, ammonia,...
modelIx_excludedMet = [];
KEGGid_excluded={'C00011';'C00080';'C00001';'C00288';'C00009';'C01342';'C00059'};
for i = 1:length(KEGGid_excluded)
    ix_tmp = find(strcmpi(data.flux.model.metKEGGID(:,2),KEGGid_excluded(i)));
    modelIx_excludedMet = [modelIx_excludedMet;ix_tmp];
end

%% loop through all fluxes, perform regulation analysis if fluxes, substrates, and at least one protein were quantified
for i = 1:size(data.flux.fluxValues(:,1))
    
    fluxData = data.flux.fluxValues(i,:);
    fluxData_interpolated = data.flux.fluxValues_interpolated(i,[ind_glc_overlap,ind_glu_overlap]);
    
    %flag: was the flux measured
    flag_flux_measured = data.flux.non_zero_flux(i,1);

    %first selection: only consider reactions were we also measured the
    %flux, and which are not the growth rate flux (i == 7)
    if (flag_flux_measured==1)&&(i~=7)
        %second selection: only consider reactions for which we also
        %measured at least one protein
        genesReaction = data.flux.model.genes(find(data.flux.model.rxnGeneMat(i,:)));
        
        [genes_intersect,ia,~] = intersect(data.protein.bnr,genesReaction);

        %check whether any of the proteins were actually quantified across
        %at least one of the limitations
        flag_protein_quantified = data.protein.non_zero_protein_present(ia);
        ia_original = ia;
        ia = ia(find(flag_protein_quantified));
        %if all criteria are fulfilled: perform regulatory analysis
        %(separately for each limitation)
        if(~isempty(ia))
            proteinData = data.protein.relativeConc_interpolated(ia,[ind_glc_overlap,ind_glu_overlap]);
            
            %third selection: only consider reactions for which I have
            %substrate data as well
            
            %since some fluxes are defined in the negative direction: flip
            %these, and make sure to flip the definition of substrate/product as well
            if(max(fluxData)<0)
                fluxData_interpolated=fluxData_interpolated.*-1;
                ix_associatedSubstrates=find(data.flux.model.S(:,i)>0);
            else
                ix_associatedSubstrates=find(data.flux.model.S(:,i)<0);
            end
            
            %exclude all common metabolites
            ix_associatedSubstrates = setdiff(ix_associatedSubstrates,modelIx_excludedMet);
            
            %find overlap with measured metabolites
            metabolitesReaction = data.flux.model.metKEGGID(ix_associatedSubstrates,2);
            [~,ia_met,ib_met] = intersect(data.metabolome.raw.AnnotatedIonListKEGG_ID,metabolitesReaction);
            
            %list of substrates that have been quantified
            metabolites_intersect_ix = ia_met;
            
            %check: do we capture all substrates?
            nonMeasuredMetabolites_ix = setdiff(ix_associatedSubstrates,ix_associatedSubstrates(ib_met));
            if(isempty(nonMeasuredMetabolites_ix))
                complete_reaction=1;
            else
                complete_reaction=0;
            end
            %if all selection criteria are fulfilled: perform regulation
            %analysis
            %run analysis separately for each flux-protein pair and each
            %metabolite
            
            if(~isempty(metabolites_intersect_ix))&&(complete_reaction==1)
                
                %loop through proteins
                for j=1:size(proteinData,1)
                    % get flux and protein data, normalize
                    normalizedFlux_glc = fluxData_interpolated(1:8)./fluxData_interpolated(1);
                    normalizedFlux_glc_log = log(abs(normalizedFlux_glc));
                    
                    normalizedFlux_glu = fluxData_interpolated(9:end)./fluxData_interpolated(9);
                    normalizedFlux_glu_log = log(abs(normalizedFlux_glu));
                    
                    normalizedProtein_glc = proteinData(j,1:8)./proteinData(j,1);
                    normalizedProtein_glc_log = log(abs(normalizedProtein_glc));
                    
                    normalizedProtein_glu=proteinData(j,9:end)./proteinData(j,9);
                    normalizedProtein_glu_log=log(abs(normalizedProtein_glu));
                    
                    %% estimate alpha's for each substrate
                    normFluxOverNormProtein = [normalizedFlux_glc_log normalizedFlux_glu_log]-[normalizedProtein_glc_log normalizedProtein_glu_log];
                    
                    %fit alpha
                    ixMet = data.metabolome.raw.AnnotatedIonListIndex(metabolites_intersect_ix);
                    metData_mean = data.metabolome.sorted.normalizedAverageIonIntensity(strainIx,ixMet);
                    metData_mean_log = log(metData_mean);
                    
                    %fitting function: constrained least square linear fit
                    %constraints: 0 < alpha < 4
                    l = length(ixMet);
                    
                    ix_nonNan = find(~isnan(normFluxOverNormProtein));
                    ix_noninf = find(~isinf(normFluxOverNormProtein));
                    ix_intersect_Nan_inf=intersect(ix_nonNan,ix_noninf,'stable');
                    if(~isempty(ix_intersect_Nan_inf))
                        [alpha_fit,~,~,exitflag_t,~,~] = lsqlin(metData_mean_log(ix_intersect_Nan_inf,:),normFluxOverNormProtein(ix_intersect_Nan_inf)',[],[],[],[],repmat(0,l,1),repmat(4,l,1));
                    else
                        alpha_fit=NaN(l,1);
                        exitflag_t=NaN;
                    end
                    % save exit flag of alpha fit. 1 == good
                    exitflag = [exitflag;exitflag_t];
                    %% perform regulation analysis
                    %regulation analysis 1: only flux vs protein                    
                    mdl_glc_prot = fitlm(normalizedFlux_glc_log',normalizedProtein_glc_log','Intercept',false);
                    mdl_glut_prot = fitlm(normalizedFlux_glu_log',normalizedProtein_glu_log','Intercept',false);
                                   
                    %regulation analysis 2: only flux vs substrates 
                    metaboliteContribution = metData_mean_log*alpha_fit;                    
                    mdl_glc_met = fitlm(normalizedFlux_glc_log',metaboliteContribution(1:8),'Intercept',false);
                    mdl_glut_met = fitlm(normalizedFlux_glu_log',metaboliteContribution(9:end),'Intercept',false);
                    
                    %regulation analysis 3: combined flux vs (protein + substrates )                 
                    mdl_glc_prot_met = fitlm(normalizedFlux_glc_log',normalizedProtein_glc_log'+metaboliteContribution(1:8),'Intercept',false);
                    mdl_glut_prot_met = fitlm(normalizedFlux_glu_log',normalizedProtein_glu_log'+metaboliteContribution(9:end),'Intercept',false);
                    
                    % get regulation coefficients and SE
                    regCoeff_tmp=[real(mdl_glc_prot.Coefficients.Estimate),real(mdl_glut_prot.Coefficients.Estimate),real(mdl_glc_met.Coefficients.Estimate),real(mdl_glut_met.Coefficients.Estimate),real(mdl_glc_prot_met.Coefficients.Estimate),real(mdl_glut_prot_met.Coefficients.Estimate)];
                    regCoeff_SE_tmp=[real(mdl_glc_prot.Coefficients.SE),real(mdl_glut_prot.Coefficients.SE),real(mdl_glc_met.Coefficients.SE),real(mdl_glut_met.Coefficients.SE),real(mdl_glc_prot_met.Coefficients.SE),real(mdl_glut_prot_met.Coefficients.SE)];
                    
                    %save the connection between flux index, protein index,
                    %and regulatory coefficients:
                    %order regulatory coefficients: 1.protein_glc
                    %2.protein_glut 3.substrate_glc 4.substrate_glut
                    
                    indx_flux_protein_met = [indx_flux_protein_met;[i,ia(j)]];
                    regulatory_coefficient = [regulatory_coefficient;regCoeff_tmp];
                    regulatory_coefficient_SE = [regulatory_coefficient_SE;regCoeff_SE_tmp];
                    

                    %save the connection between flux index, protein index,
                    %metabolite index, and, fitting parameter alpha
                    indx_flux_protein_met_alpha_tmp = [repmat(i,l,1),repmat(ia(j),l,1),metabolites_intersect_ix,alpha_fit];
                    indx_flux_protein_met_alpha = [indx_flux_protein_met_alpha;indx_flux_protein_met_alpha_tmp];
                    
                    flux_summary = [flux_summary;[normalizedFlux_glc_log normalizedFlux_glu_log]];
                    protein_summary = [protein_summary;[normalizedProtein_glc_log normalizedProtein_glu_log]];
                    metabolite_summary = [metabolite_summary;metaboliteContribution'];
                end
            end 
        end
    end
end
regulatory_coefficient(find(regulatory_coefficient(:,1:2)==0)) = NaN;

%% save output of regulation analysis
data.regulationAnalysis.enzymeSaturation.exitflag = exitflag(:,1);
data.regulationAnalysis.enzymeSaturation.indices = indx_flux_protein_met;
data.regulationAnalysis.enzymeSaturation.connectionFluxProteinMetAlpha = indx_flux_protein_met_alpha;
data.regulationAnalysis.enzymeSaturation.regulatory_coefficient = regulatory_coefficient;
data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE = regulatory_coefficient_SE;
% save associated data
data.regulationAnalysis.enzymeSaturation.protein_data = protein_summary;
data.regulationAnalysis.enzymeSaturation.flux_data = flux_summary;
data.regulationAnalysis.enzymeSaturation.metabolite_data = metabolite_summary;
% save associated protein sector
data.regulationAnalysis.enzymeSaturation.assignedSector = data.protein.SectorAssignment(indx_flux_protein_met(:,2));

%% if there is no protein regulation coefficient in either condition, set all other coefficients to NaN as well
ix1_glc = isnan(data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,1));
ix1_glut = isnan(data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(:,2));

data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(ix1_glc,[3,5]) = NaN;
data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(ix1_glut,[4,6]) = NaN;

data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(ix1_glc,[1,3,5]) = NaN;
data.regulationAnalysis.enzymeSaturation.regulatory_coefficient_SE(ix1_glut,[2,4,6]) = NaN;


%% Next: calculate mean regulation coefficent for each reaction

uniqueReactions = unique(data.regulationAnalysis.enzymeSaturation.indices(:,1),'stable');
ReactionsIsTransport = zeros(length(uniqueReactions),1); %don't consider transport reactions

for i=1:length(uniqueReactions)
   ix1 = find(data.regulationAnalysis.enzymeSaturation.indices(:,1)==uniqueReactions(i,1));
   regCoef_tmp = data.regulationAnalysis.enzymeSaturation.regulatory_coefficient(ix1,:);
   for j=1:6
       ix_nonNan=find(~isnan(regCoef_tmp(:,j)));
       if(~isempty(ix_nonNan))
           regCoef(i,j) = mean(regCoef_tmp(ix_nonNan,j));
       else
           regCoef(i,j) = NaN;
       end
   end
   
   %moreover: ignore all transport reactions
   reactionName = data.flux.model.rxnNames{uniqueReactions(i,1),1};
   k = strfind(reactionName, 'transport');
   if(~isempty(k))
       ReactionsIsTransport(i,1) = 1;
   end
end

ix_noTransport = find(ReactionsIsTransport==0);

data.regulationAnalysis.enzymeSaturation.uniqueReactions_regCoeff = regCoef(ix_noTransport,:);
data.regulationAnalysis.enzymeSaturation.uniqueReactions_ix = uniqueReactions(ix_noTransport);


end