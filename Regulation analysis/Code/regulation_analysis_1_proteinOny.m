%% Karl Kochanowski

% regulation analysis: relate flux and protein abundance changes

% 1: perform regulation analysis as described in Kochanowski et al, MSB
% (2021)

function data = regulation_analysis_1_proteinOny(data)

%% underlying linear model
modelFun = @(p,x) p(1).*x;

%% get all relevant data
mueFlux = data.flux.fluxValues(7,:);
mueFlux_interpolated = data.flux.mue_interpolated;
limType_flux = data.flux.limType;
ix_glc_flux = find(strcmpi(limType_flux,'catabolic'));
ix_glu_flux = find(strcmpi(limType_flux,'anabolic'));

%get fluxes/proteins at exactly the same growth rate as the flux data
tmp_glc = abs(repmat(mueFlux_interpolated(ix_glc_flux)',1,8)-repmat(mueFlux(1:8),length(ix_glc_flux),1));
[~,ind_glc_overlap] = min(tmp_glc,[],1);

tmp_glu = abs(repmat(mueFlux_interpolated(ix_glu_flux)',1,8)-repmat(mueFlux(9:end),length(ix_glu_flux),1));

[~,ind_glu_overlap] = min(tmp_glu,[],1);
ind_glu_overlap = ind_glu_overlap+length(ix_glu_flux);

%% perform regulation analysis for all individual protein-reaction pairs, store results
indx_flux_protein = [];
regulatory_coefficient = [];
regulatory_coefficient_SE = [];
flux_summary = [];
protein_summary = [];

for i = 1:size(data.flux.fluxValues(:,1))
    
    fluxData = data.flux.fluxValues(i,:);
    fluxData_interpolated = data.flux.fluxValues_interpolated(i,:);
    
    %flag
    flag_flux_measured = data.flux.non_zero_flux(i,1);
    
    %first selection: only consider non-zero reactions, and ignore biomass
    %vector
    if (flag_flux_measured==1)&&(i~=7)
        
        %second selection: only consider reactions for which we also
        %measured at least one protein
        genesReaction = data.flux.model.genes(find(data.flux.model.rxnGeneMat(i,:)));
        
        [genes_intersect,ia,~] = intersect(data.protein.bnr,genesReaction);
        
        %check whether any of the proteins were actually quantified across
        %at least one of the limitations (flag > 0)
        flag_protein_quantified = data.protein.non_zero_protein_present(ia);
        
        ia_original = ia;
        ia = ia(find(flag_protein_quantified));
        
        %if all criteria are fulfilled: perform regulatory analysis
        %(separately for each limitation)
        
        if(length(ia)>0)
            proteinData = data.protein.relativeConc_interpolated(ia,[ind_glc_overlap ind_glu_overlap]);
            proteinNames = data.protein.gene(ia,1);
            
            %if flux direction negative in all cases: flip it
            if(max(fluxData)<0)
                fluxData_interpolated=fluxData_interpolated.*-1;
            end
            
            %normalize fluxes (by first data point)
            fluxData_norm_glc = fluxData_interpolated(ind_glc_overlap)./fluxData_interpolated(ind_glc_overlap(1,1));
            fluxData_log_norm_glc = log(abs(fluxData_norm_glc));
            
            fluxData_norm_glut = fluxData_interpolated(ind_glu_overlap)./fluxData_interpolated(ind_glu_overlap(1,1));
            fluxData_log_norm_glut = log(abs(fluxData_norm_glut));
            
            for j=1:length(proteinNames)
                %normalize protein data
                protein_log_norm_glc = log(proteinData(j,1:8)./proteinData(j,1));
                protein_log_norm_glut = log(proteinData(j,9:end)./proteinData(j,9));

                %calculate linear regression to get slopes               
                mdl_glc = fitlm(fluxData_log_norm_glc',protein_log_norm_glc','Intercept',false);
                mdl_glut = fitlm(fluxData_log_norm_glut',protein_log_norm_glut','Intercept',false);                
                
                indx_flux_protein = [indx_flux_protein;[i,ia(j)]];
                flux_summary = [flux_summary;[fluxData_log_norm_glc fluxData_log_norm_glut]];
                protein_summary = [protein_summary;[protein_log_norm_glc protein_log_norm_glut]];
                
                regulatory_coefficient = [regulatory_coefficient;[real(mdl_glc.Coefficients.Estimate) real(mdl_glut.Coefficients.Estimate)]];
                regulatory_coefficient_SE = [regulatory_coefficient_SE;[real(mdl_glc.Coefficients.SE) real(mdl_glut.Coefficients.SE)]];               
                
            end
        end        
    end   
end

regulatory_coefficient(find(regulatory_coefficient==0)) = NaN;

%% in some cases, several (identical) fluxes are connected to the same protein
%approach to find these examples: look at cases where the same protein id pops up in a row
protein_idx_change = diff(indx_flux_protein(:,2));
flag_change_protein = [1;find(protein_idx_change~=0)+1];

%% save regulation analysis results for all individual protein-reaction pairs
data.regulationAnalysis.fluxVSprotein.indx_flux_protein = indx_flux_protein(flag_change_protein,:);
data.regulationAnalysis.fluxVSprotein.regulatory_coefficient = regulatory_coefficient(flag_change_protein,:);
data.regulationAnalysis.fluxVSprotein.regulatory_coefficient_SE = regulatory_coefficient_SE(flag_change_protein,:);
data.regulationAnalysis.fluxVSprotein.assignedSector = data.protein.SectorAssignment(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(:,2));
data.regulationAnalysis.fluxVSprotein.flux_data = flux_summary(flag_change_protein,:);
data.regulationAnalysis.fluxVSprotein.protein_data = real(protein_summary(flag_change_protein,:));

%% calculate average regulation analysis coefficients for each reaction
% Extract unique reactions
uniqueReactions = unique(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(:,1),'stable');
% flag to exclude transport reactions
ReactionsIsTransport = zeros(length(uniqueReactions),1);
% information on 
uniqueReaction_Sector = {};
for i = 1:length(uniqueReactions)
   ix1 = find(data.regulationAnalysis.fluxVSprotein.indx_flux_protein(:,1) == uniqueReactions(i,1));
   regCoef_glc_tmp = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient(ix1,1);
   ix_nonNan_glc=find(~isnan(regCoef_glc_tmp));
   if(~isempty(ix_nonNan_glc))
       regCoef_glc(i,1) = mean(regCoef_glc_tmp(ix_nonNan_glc));
   else
       regCoef_glc(i,1) = NaN;
   end
   
   regCoef_glut_tmp = data.regulationAnalysis.fluxVSprotein.regulatory_coefficient(ix1,2);
   ix_nonNan_glut = find(~isnan(regCoef_glut_tmp));
   if(~isempty(ix_nonNan_glut))
       regCoef_glut(i,1) = mean(regCoef_glut_tmp(ix_nonNan_glut));
   else
       regCoef_glut(i,1) = NaN;
   end
   
   %moreover: ignore all transport reactions
   reactionName = data.flux.model.rxnNames{uniqueReactions(i,1),1};
   k = strfind(reactionName, 'transport');
   if(~isempty(k))
       ReactionsIsTransport(i,1)=1;
   end
   
   %moreover: assign sector to each unique reaction
   sectorAssignments = data.regulationAnalysis.fluxVSprotein.assignedSector(ix1,1);
   
   unique_sectors = unique(sectorAssignments);
   if(length(unique_sectors)>1)
       n_sectorAssignment = zeros(length(unique_sectors));
      for j=1:length(unique_sectors)
          n_sectorAssignment(j,1) = length(find(strcmpi(sectorAssignments,unique_sectors(j))));
      end
      [~,ix_group_sectorAssignment] = max(n_sectorAssignment);
      if(length(ix_group_sectorAssignment)>1)
          group_sectorAssignment = {'X'};
      else
          group_sectorAssignment = unique_sectors(ix_group_sectorAssignment);
      end
   else
       group_sectorAssignment = unique_sectors;
   end
   uniqueReaction_Sector(i,1) = group_sectorAssignment;
end

%do not consider transport reactions
ix_noTransport = find(ReactionsIsTransport==0);

%% save extracted mean regulation coefficients for unique reactions
data.regulationAnalysis.fluxVSprotein.uniqueReactions_ix = uniqueReactions(ix_noTransport);
data.regulationAnalysis.fluxVSprotein.uniqueReactions_regCoeff = [regCoef_glc(ix_noTransport,1) regCoef_glut(ix_noTransport,1)];
data.regulationAnalysis.fluxVSprotein.uniqueReactions_SectorAssignment = uniqueReaction_Sector(ix_noTransport,1);


end