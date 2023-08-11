function [] = plot_tuning_curves(groups,time_phase_labels,Feature_number,unit_number_c, Figure_labels, varargin)

    [varargin,names] = Utilities.argkeyval('Group_names',varargin, Figure_labels.Grasp_Names_c');
    [varargin,figure_title_supp] = Utilities.argkeyval('Figure_title',varargin,'');
    [varargin,save_figure_name] = Utilities.argkeyval('save_figure_name',varargin,'');
    [varargin,separated_per_condition] = Utilities.argkeyval('separated_per_condition',varargin, false);
    [varargin,Audio_file] = Utilities.argkeyval('Audio_file',varargin, []);
    [varargin,flagSvgFig] = Utilities.argkeyval('SvgFig',varargin, false);

    Utilities.argempty(varargin);

    if separated_per_condition 
        condition = size(groups{1},2)/Figure_labels.fixed_trial_number;
        %warning('Bad coding right now but will fix later, assumes that there are 8 trials per condition'); 
        
    else
        condition = 1; 
    end 
    %f =   figure('units','normalized','outerposition',[0 0 0.4 0.4]);
        f =   figure('units','normalized','outerposition',[0 0 1 1]);
        %f =   figure('units','normalized','outerposition',[0.2 0 0.3 1]);

    FontSize = 20; 
    special_color = utile.get_color_rgb_codes(names); 
    
    max_FR = max(max(cell2mat(groups)));% - max(max(cell2mat(groups)))/4; 

    for con = 1:condition
        %con
        if separated_per_condition 
            %subplot(round(condition/2), round(condition/2-1), con)
            subplot(round(condition/2), 2, con)
            try
                groups_temp =  cellfun(@(x) x(:,((con-1)*Figure_labels.fixed_trial_number +1):con*Figure_labels.fixed_trial_number), groups, 'UniformOutput', false); 
            catch
                error('Adapt the Group names!');
            end
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
            %figure();
            %ER =utile.shadedErrorBar(time_idx,Mean_FR, err_ci);
            %figure();
            ER =utile.shadedErrorBar(time_idx,Mean_FR, err_ci);
            %ER =utile.shadedErrorBar(time_idx,Mean_FR, err_ci,'lineprops',{color{j},'markerfacecolor',color{j}});
            ER.mainLine.Color = special_color{j};
            ER.patch.FaceColor = special_color{j};
            ER.edge(1).Color = special_color{j};
            ER.edge(2).Color = special_color{j};

            hold on
            h{j} = plot(time_idx,Mean_FR, 'Linewidth', 3);
            %h{j} = plot(time_idx,Mean_FR,color{j}, 'Linewidth', 3);
            h{j}.Color = special_color{j};
            hold on

        end 

        xlabel('time(s)', 'FontSize', FontSize -4); 
        ylabel('Firing rate (Hz)','FontSize', FontSize -4); 
        %ylim([0, max_FR/2])
        
        
        %Plot audio file 
        if ~isempty(Audio_file)
            hold on 
            yyaxis right
            a = plot(time_idx, Audio_file{con},'k', 'Linewidth',3);
            ylabel('Audio signal')
            ylim([0,1])
        end

        phase_time_idx = arrayfun(@(x) find(time_phase_labels == x,1), unique(time_phase_labels));
        %Plot lines indicating phase changes (ITI, Cue, Delay, Action)
        PhaseNames = {'ITI', 'Cue', 'Delay1', 'ImaginedSpeech', 'Delay2', 'Speech'};
        for  i = 1:length(phase_time_idx) 
             yL = get(gca,'YLim');
             phase_change = phase_time_idx(i);
             %l(i) = line([time_idx(phase_change), time_idx(phase_change)],yL,'Color','k','LineStyle','--','Linewidth', 2);    
             l(i) = xline(time_idx(phase_change), '--k', PhaseNames{i}, 'linewidth', 2, 'Fontsize', FontSize); 
             hold on
        end

        time_points = (time_idx(phase_time_idx));
        ylim_val=get(gca,'YLim');
       
        xlim([0, time_idx(end)]);
        %only plot legend and phase names once
        if con == condition
            if  isempty(Audio_file)
                legend([l(1), h{:}],[{' Phase Start'},names{:}], 'Location', 'best', 'FontSize', FontSize-4, 'Orientation', 'horizontal');
            else
                %legend([l(1),a, h{:}],[{' Phase Start','Mean Audio file' },names{:}], 'Location', 'best', 'FontSize', FontSize-3,'Orientation', 'horizontal');

            end 
            %plot vertical text (names) of phases with small offset
            hold on 
            %ht = text(time_points+0.3, repmat(ylim_val(2)-ylim_val(2)/4, [1 length(time_points)]), {'ITI', 'Image Cue', 'Delay' 'Action'},'FontSize',FontSize);
            %set(ht,'Rotation',90)
            
            
        end 
        
        set(gca,'fontsize',FontSize)
        %title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' ' Figure_labels.Data_type  ' - ' Figure_labels.cue_type figure_title_supp],'FontSize',FontSize+5)
        if condition == 1
            %title([' Unit ' num2str(Feature_number) ' ' Figure_labels.unit_region ' All '],'FontSize',FontSize+5)
            title([' Unit ' num2str(Feature_number) ' ' Figure_labels.unit_region ' All '],'FontSize',FontSize+5)

        else
            %title([' Unit ' num2str(Feature_number) ' ' Figure_labels.Unit_nsp{unit_number_c} ' Condition ' Figure_labels.Grasp_Names_c{con}],'FontSize',FontSize+5)
            title([' Unit ' num2str(Feature_number) ' ' Figure_labels.unit_region ' Condition ' num2str(con)],'FontSize',FontSize+5)

        end 

    end 
    
    
    if Figure_labels.save_figures
        saveas(f, [save_figure_name '.jpg']); 

      if flagSvgFig
        saveas(f, [save_figure_name '.fig'])
        saveas(f, [save_figure_name '.svg'])

      end 
      disp(['Figure for Unit '  num2str(Feature_number) ' saved']);    
    end 
end

