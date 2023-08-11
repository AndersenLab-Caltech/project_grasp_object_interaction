function [] = plot_error_matrix( Error_Matrix, Labels, Title )
% figure(); 
 
% 
% imagesc(1-err_bin)
% colorbar
% 
% xticks(phase_change); 
% xticklabels([data.task.phaseNames]);
% yticks(phase_change);
% yticklabels([data.task.phaseNames]);
% title('Classification accuracy when classifying per bin'); 

figure('units','normalized','outerposition',[0 0 0.65 0.65]);
colormap('jet')
g = imagesc(Error_Matrix); 
colorbar;  
%set(gca,'clim',[0 0.3]);
set(gca, 'xtick', 1:length(Error_Matrix)) 
set(gca, 'ytick', 1:length(Error_Matrix)) 

%xlabel('Class Labels');

%Labels
%Name = image2class_simple(Labels)
xticklabels([Labels]);
yticklabels([Labels]);
xlabel('True Class')
ylabel('Predicted Class')
 set(gca,'fontsize',13)

title(Title,'FontSize' , 15);

end

