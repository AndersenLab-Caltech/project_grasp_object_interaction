function [] = print_text(Figure_labels)

disp('-------------------------------------------------------------------------------------------------------------')
disp('Analyzing following data :')
disp('-------------------------------------------------------------------------------------------------------------')
disp(['Subject           : ' Figure_labels.subject])
disp(['Thresholding      : ' Figure_labels.threshold])
disp(['Cue type          : ' Figure_labels.cue_type])
disp(['Unit region       : ' Figure_labels.unit_region])
disp(['Trial Type        : ' Figure_labels.Data_type])
LogicalStr = {'false', 'true'};
disp(['Using tuned units : ' LogicalStr{Figure_labels.tuned_units+1}])
disp(['Saving Figures    : ' LogicalStr{Figure_labels.save_figures+1}])
disp('-------------------------------------------------------------------------------------------------------------')
disp(['Grasps used : ' Figure_labels.Grasp_Names_c])


end

