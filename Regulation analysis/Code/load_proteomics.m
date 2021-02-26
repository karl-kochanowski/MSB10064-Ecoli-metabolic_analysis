%% Karl Kochanowski

% 1: load proteomics data from Hui 2015, MSB
% 2: interpolate proteomics data to enable matching the data to the
% flux/metabolomics data

% Note: only consider proteins that were detected in all titration steps
% in a given limitation

function data = load_proteomics()

%% load proteome
[~,gene,~] = xlsread('Hui et al 2015 proteomics',1,'A6:A1058');

%catabolic lim
[mu1,~,~] = xlsread('Hui et al 2015 proteomics',1,'B5:F5');
[proteinRel1,~,~] = xlsread('Hui et al 2015 proteomics',1,'B6:F1058');

%anabolic lim
[mu2,~,~] = xlsread('Hui et al 2015 proteomics',1,'J5:N5');
[proteinRel2,~,~] = xlsread('Hui et al 2015 proteomics',1,'J6:N1058');

l1 = {'catabolic'};
l2 = {'anabolic'};

protein.mue = horzcat(mu1(end:-1:1),mu2(end:-1:1));
protein.gene = gene;
protein.relativeConc = horzcat(proteinRel1(:,end:-1:1),proteinRel2(:,end:-1:1));
protein.limType = [repmat(l1,1,5),repmat(l2,1,5)];
data.protein = protein;

% %load absolute proteomics data for 
% [absproteome_data,~,~]=xlsread('2015-07-18 NCM-MG abs abundance conversion.xlsx',1,'F6:F848');
% [~,absproteome_gene,~]=xlsread('2015-07-18 NCM-MG abs abundance conversion.xlsx',1,'A6:A848');
% 
% absproteome_data(find(isnan(absproteome_data)))=0;
% absproteome_gene=absproteome_gene(find(absproteome_data~=0),1);
% absproteome_data=absproteome_data(find(absproteome_data~=0),1);
% 
% [g1,ia,ib]=intersect(protein.gene,absproteome_gene);
% 
% data.protein.absQuantAvailable=zeros(length(protein.gene),1);
% data.protein.absQuantAvailable(ia,1)=1;
% data.protein.absQuantFactor=ones(length(protein.gene),1);
% data.protein.absQuantFactor(ia,1)=absproteome_data(ib,1);

%% interpolate proteome (to be able to relate the data to the flux/metabolome data)
mueProt = data.protein.mue;
% growth rate range for interpolation
x_range = [0.2:0.001:1.05];

data.protein.mueInterpolated = [x_range,x_range];
data.protein.limType = [repmat({'catabolic'},1,length(x_range)),repmat({'anabolic'},1,length(x_range))];

% there are five conditions (i.e. gradual titration steps) per limitation
% for each limitation: only consider proteins that were detected in all
% five conditions
for i = 1:size(data.protein.relativeConc,1)
    proteinData = data.protein.relativeConc(i,:);
        
    ixIsNotNanGlc = ~isnan(proteinData(1:5));
    ixIsNotNanGlc2 = find(ixIsNotNanGlc);
    if(length(ixIsNotNanGlc2)>4)
        prot_interpol_glc = interp1(mueProt(ixIsNotNanGlc),proteinData(ixIsNotNanGlc),x_range,'linear','extrap');
        data.protein.non_zero_protein(i,1) = 1;
    else
        prot_interpol_glc = 0.*x_range;
        data.protein.non_zero_protein(i,1) = 0;
    end
    ixIsNotNanGlut = ~isnan(proteinData(6:10)); 
    ixIsNotNanGlut2 = find(ixIsNotNanGlut)+5;
    if(length(ixIsNotNanGlut2)>4)
        prot_interpol_glut = interp1(mueProt(ixIsNotNanGlut2),proteinData(ixIsNotNanGlut2),x_range,'linear','extrap');
        data.protein.non_zero_protein(i,2) = 1;
    else
        prot_interpol_glut = 0.*x_range;
        data.protein.non_zero_protein(i,2) = 0;
    end

    data.protein.relativeConc_interpolated(i,:) = [prot_interpol_glc,prot_interpol_glut];
end
data.protein.non_zero_protein_present = sum(data.protein.non_zero_protein,2);

%% load gene identifiers (FBA model: b-numbers, proteomcs: gene names)

[~,geneIDmapping,~]=xlsread('Ecoli-mapping-identifiers.xlsx',1,'C2:D4500');
for i=1:length(data.protein.gene)
   ix=find(strcmpi(geneIDmapping(:,1),data.protein.gene(i,1)));
   if(~isempty(ix))
       data.protein.bnr(i,1)=geneIDmapping(ix,2);
   else
       data.protein.bnr(i,1)={''};
   end
       
end

%% load sector assignment
[~,data.protein.SectorAssignment,~]=xlsread('Hui et al 2015 proteomics.xlsx',1,'Y6:Y1058');

end