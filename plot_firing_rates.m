% plot firing rates of individual neurons
clc
clear all
close all 

% saveFolder = 'C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data';
% 
spike_sorting_type = 'unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
%taskName = 'GraspObject_Varied_Size'; % varied object/aperature sizes 
%taskName = 'GraspObject_5050'; % 50% Go, 50% No-Go task
taskName = 'GraspObject_Combined'; % all grasp/object combinations task
subject_id = 's3';
session_date = {'20250424'}; % 0830, 0921, 0929, 1005, 1030
% 
% DataName = ['Table_' subject_id '_' taskName spike_sorting_type '.mat'];
% Data = load(fullfile(saveFolder,DataName));

%uniqueSessionDays = unique(Data.session_date);

% for n_session = 1:num(numel(uniqueSessionDays))
% 
%     ind = ismember(Data.session_date,uniqueSessionDays{n_session});
%     DataSession = Data(ind,:);
% 
% end 

%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s2_20230720_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s2_20230725_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s2_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s2\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s2_20230824_unsorted_aligned_thr_-4.5_GraspObject');

% 4S_Action Data (indiv session)
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject_4S_Action\unsorted_aligned_thr_-4.5\s3_20230921_unsorted_aligned_thr_-4.5_GraspObject_4S_Action');

% Shuffled images data
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\Table_s3_GraspObject_Shuffled_unsorted_aligned_thr_-4.5');

% LOAD DATA
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName '_' spike_sorting_type]);
Data = Data.Go_data;

Data = Data(strcmp(Data.session_date, session_date), :); % pull desired session

% % add Aperature Size column
% sizeKeywords = ['Small', 'Medium', 'Large'];
% Data.Aperature_Size = cell(height(Data),1);
% % Loop through each label and extract the size information
% for i = 1:height(Data)
%     % Use regular expression to find the size keyword after the last underscore
%     tokens = regexp(Data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');
% 
%     if ~isempty(tokens)
%         % tokens is a cell array; extract the size keyword from it
%         Data.Aperature_Size{i} = tokens{1}{1};
%     end
% end

% remove faulty data
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20231016'};
elseif strcmp(subject_id, 's3')
    error_session = {};
end 

if ~isempty(error_session)
    condition = cellfun(@(x) strcmp(x, error_session), Data.session_date);
    Data = Data(~condition,:);
end

brainAreas = Data.frPerChannel{7}; % 6 for original task, 7 for Ripple?
phase_time_idx = Data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);
%uniqueAperatureSize = unique(Data.Aperature_Size);

%%

keyboard

%%
% analyzing fr for each grasp separated by modality
for n_brain = 3%:5 %:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    allTrials = cell(1,numChannels);

    for n_channel = 1:numChannels
        figure('units','normalized','outerposition',[0 0 0.5 1])
        sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            
            go_ind = cell2mat(Data.TrialCue) == 1;
            grasp_go_idx = logical(grasp_ind .* go_ind);
            
            fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
            CueType_name = Data.TrialType(grasp_go_idx);
            fr_sep_cue_type_mean = cell2mat(cellfun(@(x) mean(fr_grasp(:,ismember(CueType_name, x)), 2), uniqueCueTypes, 'UniformOutput', false)');
            fr_sep_cue_type_trial = cellfun(@(x) fr_grasp(:,ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

            % max_fr = cellfun(max(fr_sep_cue_type_trial);
            % ylim([0 max_fr]);
            
            % UNCOMMENT HERE
            subplot(numel(uniqueGraspTypes), 1, n_grasp);
            hold on;
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
            end
            
            %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
            err_bar = {};
            color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};
            if strcmp(taskName, 'GraspObject_Combined')
                color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};
            end

            for n_cueType = 1:numel(uniqueCueTypes) % analyzing modalities (H, HO, O)
                dataTmp = fr_sep_cue_type_trial{n_cueType}';

                %allTrials{n_channel} = dataTmp; % finding artifact
                
                % figure;
                % for i = 1:12
                %     subplot(4, 3, i); % 4 rows, 3 columns layout
                %     plot(1:174, dataTmp(i, :));
                %     title(['Channel - ' num2str(n_channel) 'Trial - ' num2str(i)]);
                %     xlabel('Timebins');
                %     ylabel('Value');
                % end


                ci = bootci(1000, {@mean,dataTmp});
                Mean_FR = squeeze(mean(dataTmp));
    
                %to correctly plot confidence interval on the figure substract mean
                %FR
                err_ci(1,:) = ci(2,:) - Mean_FR; 
                err_ci(2,:) = Mean_FR - ci(1,:); 
                
                % N = size(dataTmp, 1);
                % sem = std(dataTmp) / sqrt(N);  % standard error of the mean
                % CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
                % yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

                % ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci,'lineprops',color_info{n_cueType},'transparent',true);
                % UNCOMMENT FROM HERE
                ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci);

                hold on
                err_bar{n_cueType} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_cueType},'LineWidth',2);

                ER.mainLine.Color = color_info{n_cueType};
                ER.patch.FaceColor = color_info{n_cueType};
                ER.edge(1).Color = color_info{n_cueType};
                ER.edge(2).Color = color_info{n_cueType};

            end 
            title([uniqueGraspTypes{n_grasp} ' Grasp']);
            %legend(plotHandles, uniqueCueTypes);
            legend([err_bar{:}], uniqueCueTypes','Interpreter', 'none');
             set(gca, 'FontSize', 12)

            hold off;
        end
        xlabel('Time');
        ylabel('Average Firing Rate');
        set(gca, 'FontSize', 12);
        
    end
end

%% Finding artifact
% Convert cell array to a 3D matrix (12x174x63)
data3D = cat(3, allTrials{:});

% Compute the mean along the third dimension
meanTrials = mean(data3D, 3);

figure;
for i = 1:12
    subplot(4, 3, i); % 4 rows, 3 columns layout
    plot(1:174, meanTrials(i, :));
    title(['Trial - ' num2str(i)]);
    xlabel('Timebins');
    ylabel('Value');
end
%%
% analyzing fr for each grasp separated by aperature size
for n_brain = 1%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        figure('units','normalized','outerposition',[0 0 0.5 1])
        sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            
            go_ind = cell2mat(Data.TrialCue) == 1;
            grasp_go_idx = logical(grasp_ind .* go_ind);
            
            fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
            Size_name = Data.Aperature_Size(grasp_go_idx);
            fr_sep_size_type_mean = cell2mat(cellfun(@(x) mean(fr_grasp(:,ismember(Size_name, x)), 2), uniqueAperatureSize, 'UniformOutput', false)');
            fr_sep_size_type_trial = cellfun(@(x) fr_grasp(:,ismember(Size_name, x)), uniqueAperatureSize, 'UniformOutput', false);

            % max_fr = cellfun(max(fr_sep_cue_type_trial);
            % ylim([0 max_fr]);

            subplot(numel(uniqueGraspTypes), 1, n_grasp);
            hold on;
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
            end
            
            %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
            err_bar = {};
            color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

            for n_SizeType = 1:numel(uniqueAperatureSize) % analyzing sizes (S, M, L)
                dataTmp = fr_sep_size_type_trial{n_SizeType}';

                ci = bootci(1000, {@mean,dataTmp});
                Mean_FR = squeeze(mean(dataTmp));
    
                %to correctly plot confidence interval on the figure substract mean
                %FR
                err_ci(1,:) = ci(2,:) - Mean_FR; 
                err_ci(2,:) = Mean_FR - ci(1,:); 
                
                % N = size(dataTmp, 1);
                % sem = std(dataTmp) / sqrt(N);  % standard error of the mean
                % CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
                % yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

                % ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci,'lineprops',color_info{n_cueType},'transparent',true);
                ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci);

                hold on
                err_bar{n_SizeType} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_SizeType},'LineWidth',2);

                ER.mainLine.Color = color_info{n_SizeType};
                ER.patch.FaceColor = color_info{n_SizeType};
                ER.edge(1).Color = color_info{n_SizeType};
                ER.edge(2).Color = color_info{n_SizeType};

            end 
            title([uniqueGraspTypes{n_grasp} ' Grasp']);
            %legend(plotHandles, uniqueCueTypes);
            legend([err_bar{:}], uniqueAperatureSize','Interpreter', 'none');
             set(gca, 'FontSize', 12)

            hold off;
        end
        xlabel('Time');
        ylabel('Average Firing Rate');
        set(gca, 'FontSize', 12);
        
    end
end

%% SfN plots - all grasps on 1 plot

% Analyzing fr for each grasp, ignoring cue modality and plotting on a single plot
for n_brain = 1%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 2%1:numChannels
        figure('units','normalized','outerposition',[0 0 0.15 0.2])
        %sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)],'FontWeight','bold');
        
        hold on; % Hold the plot for multiple grasps

        % Colors for grasp type
        distinctColors = {[0.2, 0.13, 0.53], [0.067, 0.467, 0.2], [0.53, 0.8, 0.93], [0.53, 0.13, 0.33]}; % Purple, Green, Light Blue, Dark Pink
        plotHandles = []; % Initialize array to store plot handles for the legend
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            
            % Get indices for the current grasp type, ignoring cue type
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            go_ind = cell2mat(Data.TrialCue) == 1; % Assuming you want to consider only "Go" trials
            grasp_go_idx = logical(grasp_ind .* go_ind); % Combined index for grasp type and "Go" trials
            
            % Extract firing rate data for this grasp, ignoring cue modality
            fr_grasp = squeeze(frData(n_channel, :, grasp_go_idx)); 
            
            % Calculate the mean firing rate across all trials for the current grasp
            Mean_FR = mean(fr_grasp, 2)';
            
            % Compute the confidence intervals using bootstrapping
            ci = bootci(1000, {@mean, fr_grasp'});
            err_ci(1, :) = ci(2, :) - Mean_FR; 
            err_ci(2, :) = Mean_FR - ci(1, :); 
            
            % Plotting each grasp type on the same plot with distinct colors
            ER = utile.shadedErrorBar(1:length(Mean_FR), Mean_FR, err_ci, 'lineprops', {'Color', distinctColors{n_grasp}, 'LineWidth', 2});
            ER.mainLine.Color = distinctColors{n_grasp}; % Set color for the main line
            ER.patch.FaceColor = distinctColors{n_grasp}; % Match shaded area color to line color
            ER.edge(1).LineStyle = 'none'; % remove edge line
            ER.edge(2).LineStyle = 'none';
            % ER.edge(1).Color = distinctColors{n_grasp}; % Edge colors
            % ER.edge(2).Color = distinctColors{n_grasp};
            
            % Store the handle for the main line for the legend
            plotHandles(n_grasp) = ER.mainLine;
            legendEntries{n_grasp} = [uniqueGraspTypes{n_grasp}]; % Legend entry for each grasp type
        end
        
        % Add vertical lines for phase changes
        for n_phase = 1:numPhases
            xline(phase_changes(n_phase), 'k--', 'LineWidth', 1.5); %phaseNames{n_phase},
        end
        
        % Set labels and legend
        %xlabel('Timebins (50 ms)');
        xlim([30 134]) % 174 (# timebins) + 5 (buffer)
        xticks([0]); % corresponds to phases and end of trial
        %ylabel('Average Firing Rate (Hz)');
        ylim([8 32]);
        yticks([10 20 30]);
        %legend(plotHandles, {'L', 'MW', 'PP', 'S3F'}, 'Interpreter', 'none', 'Location', 'best');
        set(gca, 'FontSize', 12);
        
        hold off; % Release the plot hold
    end
end


%% Testing for sig grasp type differences (grasp-specific units)
p_values = []; % Initialize array to store p-values for each channel

cue_phase_bins = 42:81; % Define the timebins corresponding to the Cue phase

for n_brain = 5 % Specify the brain area (e.g., 4 = AIP)  
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        % Prepare to store data and labels for statistical testing
        all_data = [];
        group_labels = [];
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            % Get indices for the current grasp type, ignoring cue type
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            go_ind = cell2mat(Data.TrialCue) == 1; % Consider only "Go" trials
            grasp_go_idx = logical(grasp_ind .* go_ind); % Combined index for grasp type and "Go" trials
            
            % Extract firing rate data for this grasp, focusing on Cue phase
            fr_grasp = squeeze(frData(n_channel, cue_phase_bins, grasp_go_idx)); 
            
            % Average across time within the Cue phase
            mean_fr_cue_phase = mean(fr_grasp, 1); % Mean across time (cue phase)
            
            % Concatenate mean firing rates and labels for each grasp
            all_data = [all_data; mean_fr_cue_phase']; % Trials x 1
            group_labels = [group_labels; repmat(n_grasp, numel(mean_fr_cue_phase), 1)]; % Label for each grasp type
        end
        
        % Perform one-way ANOVA to test for differences across grasp types
        [p, ~, stats] = anova1(all_data, group_labels, 'off'); % 'off' to suppress plot
        p_values = [p_values; p]; % Store p-value for this channel
    end
end

% Display results
GS_significant_channels = find(p_values < 0.05); % Channels with significant grasp differences
disp(['GS Significant channels: ' num2str(GS_significant_channels')]);
percent_significant = (numel(GS_significant_channels)/numChannels)*100;
disp(['% Significant channels in ' brainAreas{n_brain} ' for Grasp Types: ' num2str(percent_significant)]);

%% SfN - separating out modality for one grasp

% analyzing fr for each grasp separated by modality
for n_brain = 1%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 22%1:numChannels
        % figure('units','normalized','outerposition',[0 0 0.25 0.3])
        % sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
        
        for n_grasp = 4%1:numel(uniqueGraspTypes)
            
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            
            go_ind = cell2mat(Data.TrialCue) == 1;
            grasp_go_idx = logical(grasp_ind .* go_ind);
            
            fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
            CueType_name = Data.TrialType(grasp_go_idx);
            fr_sep_cue_type_mean = cell2mat(cellfun(@(x) mean(fr_grasp(:,ismember(CueType_name, x)), 2), uniqueCueTypes, 'UniformOutput', false)');
            fr_sep_cue_type_trial = cellfun(@(x) fr_grasp(:,ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

            % max_fr = cellfun(max(fr_sep_cue_type_trial);
            % ylim([0 max_fr]);

            % Create a new figure for each grasp and cue modality
            figure('units', 'normalized', 'outerposition', [0 0 0.15 0.2])
            %sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel) ' - ' uniqueGraspTypes{n_grasp}], 'FontWeight','bold');
            hold on; % Hold the plot to add multiple cue types
            
            %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
            err_bar = {};
            color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

            for n_cueType = 1:numel(uniqueCueTypes) % analyzing modalities (HO, H, O)
                dataTmp = fr_sep_cue_type_trial{n_cueType}';

                ci = bootci(1000, {@mean,dataTmp});
                Mean_FR = squeeze(mean(dataTmp));
    
                %to correctly plot confidence interval on the figure substract mean
                %FR
                err_ci(1,:) = ci(2,:) - Mean_FR; 
                err_ci(2,:) = Mean_FR - ci(1,:); 
                
                % N = size(dataTmp, 1);
                % sem = std(dataTmp) / sqrt(N);  % standard error of the mean
                % CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
                % yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

                % ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci,'lineprops',color_info{n_cueType},'transparent',true);
                ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci);

                hold on;
                err_bar{n_cueType} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_cueType},'LineWidth',2);

                ER.mainLine.Color = color_info{n_cueType};
                ER.patch.FaceColor = color_info{n_cueType};
                ER.edge(1).LineStyle = 'none'; % remove edge line
                ER.edge(2).LineStyle = 'none';

            end 

            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', 'LineWidth', 1.5);
            end

             % Add labels and legend
            %xlabel('Timebins (50 ms)');
            xlim([30 134]) % 174 (# timebins) + 5 (buffer)
            xticks([0]); % corresponds to phases and end of trial
            %ylabel('Average Firing Rate (Hz)');
            ylim([10 38]);
            yticks([15 25 35]);
            %legend([err_bar{:}], {'H','H+O','O'}, 'Interpreter', 'none', 'Location', 'best');
            set(gca, 'FontSize', 12);

            hold off; % Release the plot hold for this figure
        end
    end
end

%% SfN plots - ignore separation of grasps and just look at separation of cue

% Analyzing fr for all grasps combined, separated by cue modality
for n_brain = 3%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        figure('units', 'normalized', 'outerposition', [0 0 0.15 0.2])
        %sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel) ' - All Grasps '], 'FontWeight', 'bold');
        hold on;
        
        % Combine activity for all grasps
        all_grasp_idx = cell2mat(Data.TrialCue) == 1; % "Go" trials only
        
        % Extract firing rate data for all grasps
        fr_combined = squeeze(frData(n_channel, :, all_grasp_idx)); 
        CueType_name = Data.TrialType(all_grasp_idx);
        
        % Separate by cue modality
        fr_sep_cue_type_trial = cellfun(@(x) fr_combined(:, ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

        % Define colors for cue types
        color_info = {[.1176 .5333 .8980], [0.8471 0.1059 0.3765], [1 0.7569 0.0275]};
        err_bar = {}; % Initialize error bar handles for legend

        for n_cueType = 4%1:2%numel(uniqueCueTypes) % Loop through cue modalities (e.g., HO, H, O)
            dataTmp = fr_sep_cue_type_trial{n_cueType}';
            
            % Compute confidence intervals using bootstrapping
            ci = bootci(1000, {@mean, dataTmp});
            Mean_FR = squeeze(mean(dataTmp));
            
            % Calculate error bars
            err_ci(1, :) = ci(2, :) - Mean_FR; 
            err_ci(2, :) = Mean_FR - ci(1, :); 
            
            % Plotting each cue modality on the same graph
            ER = utile.shadedErrorBar(1:length(dataTmp), mean(dataTmp), err_ci, 'lineprops', {'Color', color_info{n_cueType}, 'LineWidth', 2});
            hold on;
            err_bar{n_cueType} = plot(1:length(dataTmp), mean(dataTmp), 'Color', color_info{n_cueType}, 'LineWidth', 2);
            
            % Set color properties for the shaded error bars
            ER.mainLine.Color = color_info{n_cueType};
            ER.patch.FaceColor = color_info{n_cueType};
            ER.edge(1).LineStyle = 'none'; % remove edge line
            ER.edge(2).LineStyle = 'none';
        end

        % Add vertical lines for phase changes
        for n_phase = 1:numPhases
            xline(phase_changes(n_phase), 'k--', 'LineWidth', 1.5);
        end

        % Add labels and legend
        %xlabel('Timebins (50 ms)');
        xlim([30 134]) % 174 (# timebins) + 5 (buffer)
        xticks([0]); % corresponds to phases and end of trial
        %ylabel('Average Firing Rate (Hz)');
        ylim([8 32]);
        yticks([10 20 30]);
        legend([err_bar{:}], {'G','G+O'}, 'Interpreter', 'none', 'Location', 'best'); %uniqueCueTypes
        set(gca, 'FontSize', 12);

        hold off; % Release the plot hold for this figure
    end
end

%% figuring out what is wrong with Object only for Combined
% too many channels, plotting all channel's FR on one plot

for n_brain = 3 % Modify as needed
    
    frData = Data.frPerChannel{n_brain}; % Extract firing rate data
    numChannels = size(frData, 1);
    channelsPerFigure = 16; % 4x4 layout
    numFigures = ceil(numChannels / channelsPerFigure); % Number of figures required
    
    for figIdx = 1:numFigures
        figure('units', 'normalized', 'outerposition', [0 0 1 1]) % Fullscreen figure
        
        % Loop over 16 channels per figure
        for n_channel_sub = 1:channelsPerFigure
            n_channel = (figIdx - 1) * channelsPerFigure + n_channel_sub;
            
            % Ensure we don't exceed the total number of channels
            if n_channel > numChannels
                break;
            end
            
            subplot(4, 4, n_channel_sub);
            hold on;
            
            % Get all "Go" trials
            all_grasp_idx = cell2mat(Data.TrialCue) == 1; 
            
            % Extract firing rate data for all trials
            fr_combined = squeeze(frData(n_channel, :, all_grasp_idx)); 
            
            % Separate trials by cue modality
            CueType_name = Data.TrialType(all_grasp_idx);
            fr_sep_cue_type_trial = cellfun(@(x) fr_combined(:, ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

            % Define colors for cue types
            color_info = {[.1176 .5333 .8980], [0.8471 0.1059 0.3765], [1 0.7569 0.0275]};

            for n_cueType = 3 % Modify to plot specific cue modalities if needed
                dataTmp = fr_sep_cue_type_trial{n_cueType}';

                % Plot each trial separately
                for trial_idx = 1:size(dataTmp, 1)
                    plot(1:size(dataTmp, 2), dataTmp(trial_idx, :), 'LineWidth', 0.5);
                end
            end

            % Add vertical lines for phase changes
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', 'LineWidth', 1);
            end

            % Formatting
            xlim([0 174]);
            %ylim([8 32]);
            title(['Ch ' num2str(n_channel)]);
            set(gca, 'FontSize', 10);
            hold off;
        end
    end
end

%% all chan's FR per trial

frData = Data.frPerChannel{1}; 
numChannels = size(frData, 1);
%numTrials = size(frData,3); % for entire dataset
uniqueTrialTypes = unique(Data.TrialType);

trialsPerFigure = 18; % 3x6 grid
% numFigures = ceil(numTrials / trialsPerFigure); % Number of figures
% required for complete dataset

for n_type = 1:length(uniqueTrialTypes)
    trialType = uniqueTrialTypes{n_type}; % Current trial type
    
    % Find indices of trials belonging to this type
    trialTypeIdx = strcmp(Data.TrialType, trialType);
    numTrials = sum(trialTypeIdx);
    
    % Extract firing rate data for this trial type
    fr_combined = frData(:, :, trialTypeIdx); % (channels x timebins x selected trials)
    
    numFigures = ceil(numTrials / trialsPerFigure); % Number of figures required

    for figIdx = 1:numFigures
        figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
        sgtitle(['Trial Type: ' trialType], 'FontSize', 14, 'FontWeight', 'bold');
       
        
        % Loop through 18 trials per figure
        for trialSubIdx = 1:trialsPerFigure
            trialIdx = (figIdx - 1) * trialsPerFigure + trialSubIdx;
            
            % Ensure we don't exceed the total number of trials
            if trialIdx > numTrials
                break;
            end
            
            subplot(3, 6, trialSubIdx);
            hold on;
            
            % Plot each channel's firing rate for this trial
            for n_channel = 1:numChannels
                plot(1:size(fr_combined, 2), squeeze(fr_combined(n_channel, :, trialIdx)), 'LineWidth', 1);
            end
            
            % Add vertical lines for phase changes
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', 'LineWidth', 1.5);
            end
    
            % Formatting
            xlim([0 174]); 
            %ylim([8 32]); 
            title(['Trial ' num2str(trialIdx)]);
            set(gca, 'FontSize', 10);
            hold off;
        end
    end
end

%% testing for sig modality differences (modality-specific units)
% Analyzing fr for all grasps combined, separated by cue modality
p_values = []; % Initialize array to store p-values for each channel

cue_phase_bins = 42:81; % Define the timebins corresponding to the Cue phase

for n_brain = 5 % 1 = SMG, 4 = AIP, 5 = M1  
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        % Combine activity for all grasps
        all_grasp_idx = cell2mat(Data.TrialCue) == 1; % "Go" trials only
        
        % Extract firing rate data for all grasps
        fr_combined = squeeze(frData(n_channel, :, all_grasp_idx)); 
        CueType_name = Data.TrialType(all_grasp_idx);
        
        % Separate by cue modality
        fr_sep_cue_type_trial = cellfun(@(x) fr_combined(cue_phase_bins, ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);
        
        % Prepare data for statistical testing
        all_data = [];
        group_labels = [];
        
        for n_cueType = 1:numel(uniqueCueTypes) % Loop through cue modalities (e.g., HO, H, O)
            dataTmp = fr_sep_cue_type_trial{n_cueType}'; % Transpose for trials x timebins
            
            % Average across the cue phase
            mean_fr_cue_phase = mean(dataTmp, 2); % Average across timebins
            
            all_data = [all_data; mean_fr_cue_phase]; % Concatenate mean firing rates across trials
            group_labels = [group_labels; repmat(n_cueType, size(dataTmp, 1), 1)]; % Label data by cue type
        end
        
        % Perform one-way ANOVA to test for differences across cue modalities
        [p, ~, stats] = anova1(all_data, group_labels, 'off'); % 'off' to suppress plot
        p_values = [p_values; p]; % Store p-value for this channel
    end
end

% Display results
MS_significant_channels = find(p_values < 0.05); % Channels with significant modality differences
disp(['MS Significant channels:' num2str(MS_significant_channels')]);
percent_significant = (numel(MS_significant_channels)/numChannels)*100;
disp(['% Significant channels ' brainAreas{n_brain} ' for Modality : ' num2str(percent_significant)]);

%% determining overlap between MS and GS units within sessions
% Find common channels between modality-specific and grasp-specific
[overlapping_channels, idx_in_MS, idx_in_GS] = intersect(MS_significant_channels, GS_significant_channels);

% Count the number of overlapping channels
num_overlapping_channels = numel(overlapping_channels);

% determine percentage of overlap for each specific group
MS_GS_per_overlap_MS_pop = (num_overlapping_channels/numel(MS_significant_channels))*100;
MS_GS_per_overlap_GS_pop = (num_overlapping_channels/numel(GS_significant_channels))*100;

% Display results
disp(['% channels encoding both modality and grasp in ' brainAreas{n_brain} ' for MS pop: ' num2str(MS_GS_per_overlap_MS_pop)]);
disp(['% channels encoding both modality and grasp in ' brainAreas{n_brain} ' for GS pop: ' num2str(MS_GS_per_overlap_GS_pop)]);
%% mean across entire trial 
% I don't really like this because it is easy to lose things when averaging
% activity across entire trial. Averaging across Cue seems more accurate to
% the claim.

% Analyzing fr for all grasps combined, separated by cue modality
p_values = []; % Initialize array to store p-values for each channel

for n_brain = 1 % 1 = SMG, 4 = AIP, 5 = M1  
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        % Combine activity for all grasps
        all_grasp_idx = cell2mat(Data.TrialCue) == 1; % "Go" trials only
        
        % Extract firing rate data for all grasps
        fr_combined = squeeze(frData(n_channel, :, all_grasp_idx)); 
        CueType_name = Data.TrialType(all_grasp_idx);
        
        % Separate by cue modality
        fr_sep_cue_type_trial = cellfun(@(x) fr_combined(:, ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);
        
        % Prepare data for statistical testing
        all_data = [];
        group_labels = [];
        
        for n_cueType = 1:numel(uniqueCueTypes) % Loop through cue modalities (e.g., HO, H, O)
            dataTmp = fr_sep_cue_type_trial{n_cueType}'; % Transpose for trials x timebins
            all_data = [all_data; mean(dataTmp, 1)]; % Concatenate mean firing rates across trials
            group_labels = [group_labels; repmat(n_cueType, size(dataTmp, 1), 1)]; % Label data by cue type
        end
        
        % Perform one-way ANOVA to test for differences across cue modalities
        [p, ~, stats] = anova1(all_data, group_labels, 'off'); % 'off' to suppress plot
        p_values = [p_values; p]; % Store p-value for this channel
        
    end
end

% Display results
significant_channels = find(p_values < 0.05); % Channels with significant modality differences
percent_significant = (numel(significant_channels)/numChannels)*100;
disp(['% Significant channels ' brainAreas(n_brain) ' : ' num2str(percent_significant)]);

%% significance totals

%% SMG - nearly equivalent prop resp to grasp and modality
% modality-specifc during trial (responding differently to at least 1 modality during entire trial)
MS_per_sig_SMG_t = [73.0159 76.1905 50 53.2258 51.6129];
MS_ave_SMG_t = mean(MS_per_sig_SMG_t); % 60.8090
MS_std_SMG_t = std(MS_per_sig_SMG_t); % 12.6936

% modality-specifc during Cue (responding differently to at least 1 modality during Cue phase)
MS_per_sig_SMG_c = [57.1429 69.8413 50 54.8387 56.4516];
MS_ave_SMG_c = mean(MS_per_sig_SMG_c); % 57.6549
MS_std_SMG_c = std(MS_per_sig_SMG_c); % 7.3612

% grasp-specific during Cue (responding differently to at least 1 grasp during the
% Cue phase)
GS_per_sig_SMG = [61.9048 60.3175 43.5484 58.0645 58.0645];
GS_ave_SMG = mean(GS_per_sig_SMG); % 56.3799
GS_std_SMG = std(GS_per_sig_SMG); % 7.3545

% overlap between MS and GS units within sessions for MS pop
MG_overlap_SMG_MS_pop = [69.4444 59.0909 48.3871 50 60];
MG_ave_SMG_MS_pop = mean(MG_overlap_SMG_MS_pop); % 57.3845
MG_std_SMG_MS_pop = std(MG_overlap_SMG_MS_pop); % 8.5246

% overlap between MS and GS units within sessions for GS pop
MG_overlap_SMG_GS_pop = [64.1026 68.4211 55.5556 47.2222 58.3333];
MG_ave_SMG_GS_pop = mean(MG_overlap_SMG_GS_pop); % 58.7270
MG_std_SMG_GS_pop = std(MG_overlap_SMG_GS_pop); % 8.1463

% similar amt of overlap among the channels, with about 60% encoding both
% grasp and modality (~40 +/- 8% indep pop for each (modality- and
% grasp-specific in SMG) 

%% AIP - we get much higher prop responding to modality than specific grasps
% modality-specific during trial
MS_per_sig_AIP_t = [24.4898 73.4694 72.5 48.7805 62.7907 46.1538];
MS_ave_AIP_t = mean(MS_per_sig_AIP_t); % 60.7389
MS_std_AIP_t = std(MS_per_sig_AIP_t); % 12.8484

% modality-specifc during Cue
MS_per_sig_AIP_c = [63.4921 63.2653 52.5 46.3415 34.8837 43.5897];
MS_ave_AIP_c = mean(MS_per_sig_AIP_c); % 48.1160
MS_std_AIP_c = std(MS_per_sig_AIP_c); % 10.5765

% grasp-specific during Cue
GS_per_sig_AIP = [24.4898 35 26.8293 30.2326 48.7179];
GS_ave_AIP = mean(GS_per_sig_AIP); % 33.0539
GS_std_AIP = std(GS_per_sig_AIP); % 9.6073

% overlap between MS and GS units within sessions for MS pop
MG_overlap_AIP_MS_pop = [32.2581 42.8571 15.7895 33.3333 52.9412];
MG_ave_AIP_MS_pop = mean(MG_overlap_AIP_MS_pop); % 35.4358
MG_std_AIP_MS_pop = std(MG_overlap_AIP_MS_pop); % 13.8023

% overlap between MS and GS units within sessions for GS pop
MG_overlap_AIP_GS_pop = [83.333 64.2857 27.2727 38.4615 47.3684];
MG_ave_AIP_GS_pop = mean(MG_overlap_AIP_GS_pop); % 52.1443
MG_std_AIP_GS_pop = std(MG_overlap_AIP_GS_pop); % 22.0725

% most of the grasp-specific units also encode modality, but the reverse is
% not true. Modality-specific units engage another population (~65 +/- 14%)
% in AIP

%% M1 - nearly equiv prop resp to grasp and modality
% modality-specific; this might be articially inflated bc M1 responds much
% less to Object overall
MS_per_sig_M1_t = [85.7143 77.4194 75 63.4921 70.3125];
MS_ave_M1_t = mean(MS_per_sig_M1_t); % 74.3877
MS_std_M1_t = std(MS_per_sig_M1_t); % 8.2641

% modality-specifc during Cue
MS_per_sig_M1_c = [69.8413 75.8065 70.3125 50.7937 73.4375];
MS_ave_M1_c = mean(MS_per_sig_M1_c); % 68.0383
MS_std_M1_c = std(MS_per_sig_M1_c); % 9.9410

% grasp-specific during Cue
GS_per_sig_M1 = [63.4921 72.5806 60.9375 65.0794 76.5625];
GS_ave_M1 = mean(GS_per_sig_M1); % 67.7304
GS_std_M1 = std(GS_per_sig_M1); % 6.5701

% overlap between MS and GS units within sessions for MS pop
MG_overlap_M1_MS_pop = [61.3636 74.4681 62.2222 71.875 74.4681];
MG_ave_M1_MS_pop = mean(MG_overlap_M1_MS_pop); % 68.8794
MG_std_M1_MS_pop = std(MG_overlap_M1_MS_pop); % 6.5621

% overlap between MS and GS units within sessions for GS pop
MG_overlap_M1_GS_pop = [67.5 77.7778 71.7949 56.0976 71.4286];
MG_ave_M1_GS_pop = mean(MG_overlap_M1_GS_pop); % 68.9198
MG_std_M1_GS_pop = std(MG_overlap_M1_GS_pop); % 8.0537

% similar amt of overlap among the channels, with about 70% encoding both
% grasp and modality (~30 +/- 7% indep pop for each (modality- and
% grasp-specific in M1)

%% plotting overlap data

brainRegions = {'SMG', 'AIP', 'M1'};
percentages_grasp = [GS_ave_SMG, GS_ave_AIP, GS_ave_M1];
percentages_modality = [MS_ave_SMG_c, MS_ave_AIP_c, MS_ave_M1_c];
std_dev_grasp = [GS_std_SMG, GS_std_AIP, GS_std_M1]; 
std_dev_modality = [MS_std_SMG_c, MS_std_AIP_c, MS_std_M1_c]; 

% Overlap data
percent_overlap_grasp = [MG_ave_SMG_GS_pop, MG_ave_AIP_GS_pop, MG_ave_M1_GS_pop];
percent_overlap_modality = [MG_ave_SMG_MS_pop, MG_ave_AIP_MS_pop, MG_ave_M1_MS_pop];
std_dev_overlap_grasp = [MG_std_SMG_GS_pop, MG_std_AIP_GS_pop, MG_std_M1_GS_pop]; 
std_dev_overlap_modality = [MG_std_SMG_MS_pop, MG_std_AIP_MS_pop, MG_std_M1_MS_pop];

% Create the figure and subplots
figure;

% Top subplot for grasp-specific and modality-specific percentages with standard deviations
subplot(2, 1, 1);

% X-axis positions for bars
x = 1:length(brainRegions);
barWidth = 0.35; % Width of bars

% Bars for grasp-specific and modality-specific percentages
hold on;
b1 = bar(x - barWidth/2, percentages_grasp, barWidth, 'FaceColor', [0, 0.3529, 0.7098], 'EdgeColor', 'k'); % Grasp-specific
b2 = bar(x + barWidth/2, percentages_modality, barWidth, 'FaceColor', [0.8327, 0.1961, 0.1255], 'EdgeColor', 'k'); % Modality-specific

% Add error bars for standard deviations
errorbar(x - barWidth/2, percentages_grasp, std_dev_grasp, 'k.', 'LineWidth', 1.5);
errorbar(x + barWidth/2, percentages_modality, std_dev_modality, 'k.', 'LineWidth', 1.5);

% Customize the subplot
set(gca, 'XTick', x, 'XTickLabel', brainRegions);
xlabel('Brain Region');
xlim([min(x) - 0.5, max(x) + 0.5]);
ylabel('Percentage (%)');
ylim([0 100])
yticks([0 50 100])
title('Percentage of Total Units');
legend({'Grasp-specific', 'Modality-specific'}, 'Location', 'best');
hold off;

% Bottom subplot for overlap percentages with standard deviations
subplot(2, 1, 2);

% X-axis positions for overlap bars
barWidth = 0.35; % Width of bars

% Bars for overlap percentages
hold on;
b3 = bar(x - barWidth/2, percent_overlap_grasp, barWidth, 'FaceColor', [0, 0.3529, 0.7098], 'EdgeColor', 'k'); % Overlap in grasp-specific
b4 = bar(x + barWidth/2, percent_overlap_modality, barWidth, 'FaceColor', [0.8327, 0.1961, 0.1255], 'EdgeColor', 'k'); % Overlap in modality-specific

% Add error bars for standard deviations
errorbar(x - barWidth/2, percent_overlap_grasp, std_dev_overlap_grasp, 'k.', 'LineWidth', 1.5);
errorbar(x + barWidth/2, percent_overlap_modality, std_dev_overlap_modality, 'k.', 'LineWidth', 1.5);

% Customize the subplot
set(gca, 'XTick', x, 'XTickLabel', brainRegions);
xlabel('Brain Region');
xlim([min(x) - 0.5, max(x) + 0.5]);
ylabel('Percentage (%)');
ylim([0 100])
yticks([0 50 100])
title('Percentage Overlap between Grasp and Modality');
legend({'Overlap in Grasp-specific', 'Overlap in Modality-specific'}, 'Location', 'best');
hold off;
