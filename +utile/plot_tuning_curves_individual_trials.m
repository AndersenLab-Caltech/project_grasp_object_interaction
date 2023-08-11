function [] = plot_tuning_curves_individual_trials(groups,time_phase_labels,Feature_number,unit_number_c, Figure_labels, varargin)

    %Plotting the firing rate of individual trials (instead of plotting the mean with
    %bootstrapping 95% confidence interval)

    [varargin,names] = Utilities.argkeyval('Group_names',varargin, Figure_labels.Grasp_Names_c');
    [varargin,figure_title_supp] = Utilities.argkeyval('Figure_title',varargin,'');
    [varargin,save_figure_name] = Utilities.argkeyval('save_figure_name',varargin,'');
    [varargin,separated_per_condition] = Utilities.argkeyval('separated_per_condition',varargin, false);
    [varargin,two_trials] = Utilities.argkeyval('two_trials',varargin, false);

    Utilities.argempty(varargin);

    if separated_per_condition 
        condition = size(groups{1},2)/Figure_labels.fixed_trial_number;
        %warning('Bad coding right now but will fix later, assumes that there are 8 trials per condition'); 
        
    else
        condition = 1; 
    end
    
    max_FR = max(max(cell2mat(groups)));

    %f = figure(); 
    f =   figure('units','normalized','outerposition',[0 0 1 1]);
    FontSize = 20; 
    color = utile.get_colors_per_grasp(names);
    special_color = utile.get_color_rgb_codes(names); 
    for con = 1:condition
        
        if separated_per_condition 
            subplot(round(condition/2), 2, con)
            groups_temp =  cellfun(@(x) x(:,((con-1)*8 +1):con*8), groups, 'UniformOutput', false); 
            FontSize = 10; 

        else
            groups_temp = groups; 
        end 
        
        for j =1:size(groups_temp,1)

           data_tmp = permute(groups_temp{j}, [2 1]); 
           Mean_FR = squeeze(mean(data_tmp));

           timer_period = 0.05; 
           time_idx = ([0:length(Mean_FR)-1]*timer_period)'; 

          %Plot firing rates for each trial
          for trial = 1:size(data_tmp,1)
              hold on
              p{j} = plot(time_idx,data_tmp(trial,:), color{j},'LineStyle','-', 'Linewidth', 0.8);
              p{j}.Color = special_color{j};
          end

          hold on
          h{j} = plot(time_idx,Mean_FR,color{j}, 'Linewidth', 2);
          h{j}.Color = special_color{j};
          hold on

        end 
   
        xlabel('time(s)', 'FontSize', FontSize); 
        ylabel('Firing rate (Hz)','FontSize', 55); 
        ylim([0, max_FR])

        if two_trials
            time_phase_labels_double = [time_phase_labels; time_phase_labels+4];
            phase_time_idx = arrayfun(@(x) find(time_phase_labels_double == x,1), unique(time_phase_labels_double));
            Phase_names = {'ITI', 'Image Cue', 'Delay' 'Action','ITI 2', 'Image Cue 2', 'Delay 2' 'Action 2'};
        else
            phase_time_idx = arrayfun(@(x) find(time_phase_labels == x,1), unique(time_phase_labels));
            Phase_names = {'ITI', 'Image Cue', 'Delay' 'Action'};
        end 
        

        for  i = 1:length(phase_time_idx) %(data.phase_time(2) - data.phase_time(1))+1
             yL = get(gca,'YLim');
             phase_change = phase_time_idx(i);
             l(i) = line([time_idx(phase_change), time_idx(phase_change)],yL,'Color','k','LineStyle','--','Linewidth', 2);
        end 
        
        time_points = (time_idx(phase_time_idx)); 
        get_ylim=get(gca,'YLim');
        xlim([0, time_idx(end)]);
        
        if con == condition
            legend([l(1), h{:}],[{' Phase Start'},names{:}], 'Location', 'best', 'FontSize', FontSize );
            %plot vertical names of phases
            hold on 
            ht = text(time_points+0.3, repmat(get_ylim(2)-get_ylim(2)/4, [1 length(time_points)]), Phase_names,'FontSize',FontSize +5);
            set(ht,'Rotation',90)
        end
 
        set(gca,'fontsize',FontSize)
        if condition == 1
            title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' All '],'FontSize',FontSize+5)
        else
            title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' Condition ' Figure_labels.Grasp_Names_c{con}],'FontSize',FontSize+5)
        end         
    end 
    
    %fig_title = ([' Unit '  num2str(Feature_number) 'Legend  tuning curve.jpg']); 
    %fig_title_svg = (['Unit '  num2str(Feature_number) ' ' Figure_labels.Data_type '.svg']); 
    if Figure_labels.save_figures
     saveas(f, save_figure_name); 
     disp(['Figure for Unit '  num2str(Feature_number) ' saved']);    
    end 
end

