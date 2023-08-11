function [ Error_matrix, Labels ] = misclassification( real_label,classified_label )

%Function which returns an error matrix which shows the classification accuracy. 
%The returned matrix is of size [n,n], where n is the number of classes present
%in the dataset. The diagonal of the matrix indicicates the percentage of 
%correctly classified data

temp = [real_label, classified_label];

% unique finds the number of uniques values of an array -> find out how
% many different labels we have to compute the matrix 
Error_matrix = zeros(length(unique(temp))); 
labels = unique(temp);
%labels
for i = 1:size(Error_matrix,1)
    
    for j = 1:size(Error_matrix,2)
        
        if ismember(labels(j), unique(real_label))
            Error_matrix(i,j) = (length(temp(temp(:,1) == labels(j) & temp(:,2) == labels(i))))/(nnz(temp(:,1) ==labels(j)));
        else 
            Error_matrix(i,j) = 0;
        
        end
    end 
end 

% figure(); imagesc(Error_matrix); colorbar; set(gca,'clim',[0 1]);
% xlabel('Class Labels');
  Labels= num2cell(labels);
% for i = 1:length(Labels_tmp)
%     Labels{2*i -1} = Labels_tmp{i};
%     Labels{2*i} = 0;
% end 
Labels = cellfun(@(x)num2str(x),Labels,'un',0);
%Labels
%xticklabels({'0', '3', '5', '6','8', '9'});
% xticklabels(Labels);
% yticklabels(Labels);
% 
% %set(gca, 'Xtick', Labels);
% %set(gca, 'Ytick', Labels);
% title('Error Matrix');

end 
