function [] = plot_pie_chart(numbers_cue, Labels_pie, Figure_labels,num_tuned_tot, varargin )
    %I REALYY SHOULD COMMENT MY FUNCTIONS - DUH
    
    
    
    
    [varargin,x_label_text] = Utilities.argkeyval('xlabel_text',varargin, {' '});
    Utilities.argempty(varargin);
    %figure()

    FontSize = 35; 
    idx_remove = find(numbers_cue == 0,1);
    if ~isempty(idx_remove)
        numbers_cue(idx_remove) = []; 
        Labels_pie(idx_remove) = []; 
    end 
    
    f =figure('units','normalized','outerposition',[0 0 1 1]);
    h = pie(numbers_cue);
    T = h(strcmpi(get(h, 'Type'), 'text')); 
     for i =1:length(T)
        % T(i).String = strcat(Labels_pie(i),  T(i).String); 
        %decomment if want text and percentage inside the pie chart. it too
        %small right now, so I don't do it. 
        
 %        T(i).String = {Labels_pie(i),  T(i).String}; 
         T(i).FontSize = FontSize; 
         T(i).FontWeight = 'bold';
     end

%    T = strcat(Labels_pie, T); 
    P = cell2mat(get(T, 'Position')); 
    set(T, {'Position'}, num2cell(P*0.6,2)); 
   % text(P(:,2), P(:,2), Labels_pie(:), 'FontSize', 14);
    
    dim = [.45 0.82 .1 .1];  %  annotation('textbox',dim,'String',x_label_text,'FitBoxToText','on');
    annotation('textbox',dim,'String',x_label_text,'FitBoxToText','on' , 'LineStyle', 'none', 'FontSize', FontSize, 'FontWeight', 'bold');
 
    %title([ Figure_labels.subject ' -  ' Figure_labels.unit_region ' Region (#Units = ' num2str(num_tuned_tot)  ...
    %     , ', #Grasps = ' num2str(length(Figure_labels.Grasp_Names_c  )) ')'], 'FontSize', FontSize);
   
    title([ Figure_labels.subject ' - ' Figure_labels.cue_type ' - Region: ' Figure_labels.unit_region ' - #Units = ' num2str(num_tuned_tot)  ...
         ], 'FontSize', FontSize);
  
    colors_here = utile.get_color_rgb_codes(Labels_pie);
    colormap(cell2mat(colors_here'));
    legend(Labels_pie','FontSize', FontSize-5, 'location' ,'best');
    %gridLegend(h, 4, Labels_pie)
    %utile.columnlegend(1, Labels_pie','FontSize', FontSize);    
    
    %fig_title = (['Subject ' Figure_labels.subject '- Unit tuning for ' Figure_labels.unit_region ' Region #Units = ' num2str(num_tuned_tot) '.jpg']); 
    %fig_title2 = ([Figure_labels.subject '- Unit tuning for ' Figure_labels.unit_region ' Region #Units = ' num2str(num_tuned_tot) '.fig']); 

    fig_title = (['Subject ' Figure_labels.subject '- Unit gerneral population tuning for ' Figure_labels.unit_region ' Region #grasps = ' num2str(length(Figure_labels.Grasp_Names_c)) '.jpg']); 
    fig_title2 = ([Figure_labels.subject '- Unit general population tuning for ' Figure_labels.unit_region ' Region #grasps = ' num2str(length(Figure_labels.Grasp_Names_c)) '.fig']); 
    fig_title_svg = ([Figure_labels.subject '- Action Unit general population tuning for ' Figure_labels.unit_region ' Region #grasps = ' num2str(length(Figure_labels.Grasp_Names_c)) ' ' Figure_labels.Data_type '.svg']); 

    if Figure_labels.save_figures
        % saveFigure(f,fig_title);
        % savefig(f,fig_title2);
         saveas(f, fig_title_svg); 
         disp(['Figure for Subject ' Figure_labels.subject '- Unit tuning for ' Figure_labels.unit_region ' Region.jpg' ' saved']);
    end 

end

