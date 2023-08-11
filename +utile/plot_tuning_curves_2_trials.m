function [] = plot_tuning_curves_2_trials(groups,time_phase_labels,Feature_number,unit_number_c, Figure_labels, varargin)

    [varargin,names] = Utilities.argkeyval('Group_names',varargin, Figure_labels.Grasp_Names_c');
    [varargin,figure_title_supp] = Utilities.argkeyval('Figure_title',varargin,'');
    [varargin,save_figure_name] = Utilities.argkeyval('save_figure_name',varargin,'');
    [varargin,separated_per_condition] = Utilities.argkeyval('separated_per_condition',varargin, false);
    [varargin,max_firing_rate] = Utilities.argkeyval('y_axis_max_firing_rate',varargin, false);

    Utilities.argempty(varargin);

    if separated_per_condition 
        condition = size(groups{1},2)/Figure_labels.fixed_trial_number;
        %warning('Bad coding right now but will fix later, assumes that there are 8 trials per condition'); 
        
    else
        condition = 1; 
    end 
    f =   figure('units','normalized','outerposition',[0 0 1 1]);
    FontSize = 20; 
    color = utile.get_colors_per_grasp(names);
    special_color = utile.get_color_rgb_codes(names); 
    
%     if max_firing_rate
%         max_FR = max(max(cell2mat(groups)));
%     else
%         max_FR = max(max(cell2mat(groups))) - max(max(cell2mat(groups)))/2; 
%     end 

    for con = 1:condition
        
        if separated_per_condition 
            subplot(round(condition/2), round(condition/2-1), con)
            groups_temp =  cellfun(@(x) x(:,((con-1)*8 +1):con*8), groups, 'UniformOutput', false); 
            FontSize = 10; 

        else
            groups_temp = groups; 
        end 
        
        
        for j =1:size(groups,1)

            data_tmp = permute(groups_temp{j}, [2 1]); 
            %compute confidence interval with 1000 x bootstrapping
            ci = bootci(1000, {@mean,data_tmp});
            Mean_FR = squeeze(mean(data_tmp));

            %to correctly plot confidence interval on the figure substract mean
            %FR
            err_ci(1,:) = ci(2,:) - Mean_FR; 
            err_ci(2,:) = Mean_FR - ci(1,:);  

            timer_period = 0.05; 
            time_idx = ([0:length(Mean_FR)-1]*timer_period)'; 

            %Plot Error bar with bootstrapping confidence interval
            ER =utile.shadedErrorBar(time_idx,Mean_FR, err_ci,'lineprops',{color{j},'markerfacecolor',color{j}});
            ER.mainLine.Color = special_color{j};
            ER.patch.FaceColor = special_color{j};
            ER.edge(1).Color = special_color{j};
            ER.edge(2).Color = special_color{j};

            hold on
            h{j} = plot(time_idx,Mean_FR,color{j}, 'Linewidth', 3);
            h{j}.Color = special_color{j};
            hold on

        end 

        xlabel('time(s)', 'FontSize', FontSize); 
        ylabel('Firing rate (Hz)','FontSize', FontSize); 
        %ylim([0, max_FR]);
        
        time_phase_labels_double = [time_phase_labels; time_phase_labels+4];
        phase_time_idx = arrayfun(@(x) find(time_phase_labels_double == x,1), unique(time_phase_labels_double));
        %Plot lines indicating phase changes (ITI, Cue, Delay, Action)
        for  i = 1:length(phase_time_idx) 
             yL = get(gca,'YLim');
             phase_change = phase_time_idx(i);
             l(i) = line([time_idx(phase_change), time_idx(phase_change)],yL,'Color','k','LineStyle','--','Linewidth', 2);        
        end

        time_points = (time_idx(phase_time_idx));
        ylim_val=get(gca,'YLim');
       
        xlim([0, time_idx(end)]);
        %only plot legend and phase names once
        if con == condition
            legend([l(1), h{:}],[{' Phase Start'},names{:}], 'Location', 'best', 'FontSize', FontSize);
            %plot vertical text (names) of phases with small offset
            hold on 
            ht = text(time_points+0.3, repmat(ylim_val(2)-ylim_val(2)/4, [1 length(time_points)]), {'ITI', 'Image Cue', 'Delay' 'Action', 'ITI 2', 'Image Cue 2', 'Delay 2' 'Action 2', },'FontSize',FontSize+5);
            set(ht,'Rotation',90)
        end 
        
        set(gca,'fontsize',FontSize)
        %title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' ' Figure_labels.Data_type  ' - ' Figure_labels.cue_type figure_title_supp],'FontSize',FontSize+5)
        title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' Condition ' num2str(con)],'FontSize',FontSize+5)

        %fig_title = ([' Unit '  num2str(Feature_number) 'Legend  tuning curve.jpg']); 
        %fig_title_svg = (['Unit '  num2str(Feature_number) ' ' Figure_labels.Data_type '.svg']); 
    
    end 
    
    
    if Figure_labels.save_figures
      saveas(f, save_figure_name); 
      disp(['Figure for Unit '  num2str(Feature_number) ' saved']);    
    end 
end

