%% 2021-01-03 Karl Kochanowski

% plot results of FVA

function plot_FVA_results(res)

% get flux data
flux = res.FBA.flux;
% get min and max FVA
FVAmin = res.FBA.minFVA;
FVAmax = res.FBA.maxFVA;

ix_nonzero = find(flux ~= 0);
ix_zero = find(flux == 0);

ratio_min = FVAmin./flux;
ratio_max = FVAmax./flux;

ratio_min(ix_zero) = NaN;
ratio_max(ix_zero) = NaN;

%% plot FVA Version 1: x-axis = median flux, y-axis: median FVAmin and FVAmax
% only consider non-zero fluxes
figure('Name','flux versus FVA V1')
line([-60 60],[-60 60],'Color','k','LineStyle','--');
hold on;
line([-60 60],[0 0],'Color','k');
line([0 0],[-60 60],'Color','k');
box on;
for i = 1:size(flux,1)
    if(sum(flux(i,:))~=0)
       ix1 = find(flux(i,:)~=0);
       if(length(ix1) > 0)
           
       flux_median = median(flux(i,ix1));
       FVA_min_median = median(FVAmin(i,ix1));
       FVA_max_median = median(FVAmax(i,ix1));
       
       FVA_min_median(find(FVA_min_median<-60)) = -60;
       FVA_max_median(find(FVA_max_median>60)) = 60;
                  
       plot(flux_median,flux_median,'ok');
       line([flux_median flux_median],[FVA_min_median FVA_max_median],'Color','r');
       end
    end
    
end
axis([-60 60 -60 60]);
set(gca,'FontSize',12,'XTick',[-60:20:60],'YTick',[-60:20:60]);
xlabel('median flux [mmol/h/gCDW]');
ylabel('FVA range [mmol/h/gCDW]');


% %% plot FVA Version 2: log-scale x-axis = median flux, y-axis: median FVAmin and FVAmax
% % only consider non-zero fluxes
% figure('Name','flux versus FVA V2')
% % line([0.0001 60],[-60 60],'Color','k','LineStyle','--');
% % hold on;
% % line([-60 60],[0 0],'Color','k');
% % line([0 0],[-60 60],'Color','k');
% % box on;
% for i = 1:size(flux,1)
%     if(sum(flux(i,:))~=0)
%        ix1 = find(flux(i,:)~=0);
%        if(length(ix1) == 16)
%            
%        flux_median = median(flux(i,ix1));
%        FVA_min_median = median(FVAmin(i,ix1));
%        FVA_max_median = median(FVAmax(i,ix1));
%        
%        if(flux_median > 0)
%            FVA_min_median(find(FVA_min_median<=0)) = 0.0001;
%            FVA_max_median(find(FVA_max_median>60)) = 60;
%        else
%            FVA_min_median(find(FVA_min_median==0)) = 0.0001;
%            FVA_max_median(find(FVA_max_median==0)) = 0.0001;
%        end
%                   
%        loglog(abs(flux_median),abs(flux_median),'ok');
%        hold on;
%        line([flux_median flux_median],[FVA_min_median FVA_max_median],'Color','r');
%        end
%     end
%     
% end
% axis([0.0001 60 0.0001 60]);
% 
% 



end