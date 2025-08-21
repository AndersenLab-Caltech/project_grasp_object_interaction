% plot firing rates of individual neurons
clc
clear all
close all 

% saveFolder = 'C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data';
% 
spike_sorting_type = 'unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
%taskName = 'GraspObject_Varied_Size'; % varied object/aperture sizes 
%taskName = 'GraspObject_5050'; % 50% Go, 50% No-Go task
taskName = 'GraspObject_Combined'; % all grasp/object combinations task
subject_id = 's2';
session_date = {'20250728'}; % 0830, 0921, 0929, 1005, 1030 

% LOAD DATA
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName '_' spike_sorting_type]);
Data = Data.Go_data;

%Data = Data(strcmp(Data.session_date, session_date), :); % pull desired session, comment out when calculating specific units across sessions

% add Aperture Size column
sizeKeywords = ['Small', 'Medium', 'Large'];
Data.Aperture_Size = cell(height(Data),1);
% Loop through each label and extract the size information
for i = 1:height(Data)
    % Use regular expression to find the size keyword after the last underscore
    tokens = regexp(Data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');

    if ~isempty(tokens)
        % tokens is a cell array; extract the size keyword from it
        Data.Aperture_Size{i} = tokens{1}{1};
    end
end

if strcmp(taskName, 'GraspObject_Combined')
    Data.TrialType(strcmp(Data.TrialType, 'Unknown')) = {'Combined'}; % adds in Combined as Trial type
    % add in column with Object Type for Combined trials and original trial types (H, HO, O with Associated)
    % Loop through each label and extract the object information
    for i = 1:height(Data)
        % Use regular expression to find the size keyword after the last underscore
        tokens = regexp(Data.LabelNames{i}, '_(deck|block|rod|ball)$', 'tokens');
        
        if ~isempty(tokens)
            % tokens is a cell array; extract the size keyword from it
            Data.ObjectType{i} = tokens{1}{1};
        else 
            Data.ObjectType{i} = 'Associated';
        end
    end
    color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
end

% remove faulty data
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20231016'};
elseif strcmp(subject_id, 's3')
    error_session = {'20250212'};
elseif strcmp(subject_id, 's4')
    error_session = {'20240613'};
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
sessions_all = unique(Data.session_date);
numSessions = numel(sessions_all);

uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);
if strcmp(taskName, 'GraspObject_Varied_Size')
    uniqueApertureSize = unique(Data.Aperture_Size);
end
if strcmp(taskName, 'GraspObject_Combined')
    uniqueObjectType = unique(Data.ObjectType);
    % remove "Associated" => usual grasp with usual object
    uniqueObjectType = uniqueObjectType(2:5);
end

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
%% analyzing fr for each grasp separated by aperture size/object type
if strcmp(taskName, 'GraspObject_Varied_Size')
    for n_brain = 4%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
        
        frData = Data.frPerChannel{n_brain};
        numChannels = size(frData, 1);
        
        for n_channel = 1:numChannels
            figure('units','normalized','outerposition',[0 0 .25 .35]);%[0 0 0.5 1])
            sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
            
            %for n_grasp = 1:numel(uniqueGraspTypes)
                
                %grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
                
                go_ind = cell2mat(Data.TrialCue) == 1;
                %grasp_go_idx = logical(grasp_ind .* go_ind);
                
                %fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
                fr_go = squeeze(frData(n_channel,:,go_ind));
                Size_name = Data.Aperture_Size(go_ind);%grasp_go_idx);
                fr_sep_size_type_mean = cell2mat(cellfun(@(x) mean(fr_go(:,ismember(Size_name, x)), 2), uniqueApertureSize, 'UniformOutput', false)');
                fr_sep_size_type_trial = cellfun(@(x) fr_go(:,ismember(Size_name, x)), uniqueApertureSize, 'UniformOutput', false);
    
                % max_fr = cellfun(max(fr_sep_cue_type_trial);
                % ylim([0 max_fr]);
    
                %subplot(numel(uniqueGraspTypes), 1, n_grasp);
                hold on;
                for n_phase = 1:numPhases
                    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
                end
                
                %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
                err_bar = {};
                %color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569.0275]}; %conditions
                color_info = {[.4688 .3672 .9375],[.9922 .3789 0],[.9961 .6875 0]}; %sizes
    
                for n_SizeType = [1,3] %1:numel(uniqueApertureSize) % analyzing sizes (S, M, L)
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
    
                    ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci);
    
                    hold on
                    err_bar{n_SizeType} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_SizeType},'LineWidth',2);
    
                    ER.mainLine.Color = color_info{n_SizeType};
                    ER.patch.FaceColor = color_info{n_SizeType};
                    ER.edge(1).Color = 'none';
                    ER.edge(2).Color = 'none';
    
                end 
                % title([uniqueGraspTypes{n_grasp} ' Grasp']);
                %legend(plotHandles, uniqueCueTypes);
                % legend([err_bar{:}], uniqueApertureSize','Interpreter', 'none');
                 set(gca, 'FontSize', 12)
    
                hold off;
            %end
            % legend([err_bar{:}], uniqueApertureSize','Interpreter', 'none');
            legend([err_bar{:}], [{'Large'} {'Small'}],'Interpreter', 'none');
            xlabel('Timebins (50 ms)');
            xlim([0 180]);
            ylabel('Average Firing Rate');
            set(gca, 'FontSize', 12);
            
        end
    end
elseif strcmp(taskName, 'GraspObject_Combined')
    for n_brain = 1%:length(brainAreas) % 1:5 for AN, 1:3 for FG, [1, 3:6] for GB
        
        frData = Data.frPerChannel{n_brain};
        numChannels = size(frData, 1);
        
        for n_channel = 1:numChannels
            figure('units','normalized','outerposition',[0 0 .17 .25]);%[0 0 0.5 1])
            sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
            
            %for n_grasp = 1:numel(uniqueGraspTypes)
                
                %grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
                
                go_ind = cell2mat(Data.TrialCue) == 1;
                %grasp_go_idx = logical(grasp_ind .* go_ind);
                
                %fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
                fr_go = squeeze(frData(n_channel,:,go_ind));
                Object_name = Data.ObjectType(go_ind);%grasp_go_idx);
                fr_sep_object_type_mean = cell2mat(cellfun(@(x) mean(fr_go(:,ismember(Object_name, x)), 2), uniqueObjectType, 'UniformOutput', false)');
                fr_sep_object_type_trial = cellfun(@(x) fr_go(:,ismember(Object_name, x)), uniqueObjectType, 'UniformOutput', false);
    
                % max_fr = cellfun(max(fr_sep_cue_type_trial);
                % ylim([0 max_fr]);
    
                %subplot(numel(uniqueGraspTypes), 1, n_grasp);
                hold on;
                for n_phase = 1:numPhases
                    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
                end
                
                %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
                err_bar = {};
                %color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569.0275]}; %conditions
                color_info = {[.3906 .5586 .9961],[.4688 .3672 .9375],[.8594 .1484 .4961],[.9922 .3789 0]}; % objects (including Assoc.) [.9961 .6875 0],
    
                for n_type = 1:numel(uniqueObjectType) % analyzing sizes (S, M, L)
                    dataTmp = fr_sep_object_type_trial{n_type}';
    
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
    
                    ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),err_ci);
    
                    hold on
                    err_bar{n_type} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_type},'LineWidth',2);
    
                    ER.mainLine.Color = color_info{n_type};
                    ER.patch.FaceColor = color_info{n_type};
                    ER.edge(1).Color = 'none';
                    ER.edge(2).Color = 'none';
    
                end 
                % title([uniqueGraspTypes{n_grasp} ' Grasp']);
                %legend(plotHandles, uniqueCueTypes);
                % legend([err_bar{:}], uniqueApertureSize','Interpreter', 'none');
                 set(gca, 'FontSize', 12)
    
                hold off;
            %end
            % legend([err_bar{:}], uniqueApertureSize','Interpreter', 'none');
            legend([err_bar{:}], uniqueObjectType,'Interpreter', 'none');
            xlabel('Timebins (50 ms)');
            xlim([0 180]);
            ylabel('Average Firing Rate');
            set(gca, 'FontSize', 12);
            
        end
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

phase_bins = 42:81; % Define the timebins corresponding to the Cue phase

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
            fr_grasp = squeeze(frData(n_channel, phase_bins, grasp_go_idx)); 
            
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



%% GS and MS Parameters
phase_bins = 94:174; % cue = 42:82, action = 94:174
alpha = 0.05;

% Preallocate storage
GS_percent_per_session = cell(numel(brainAreas), 1);
MS_percent_per_session = cell(numel(brainAreas), 1);
overlap_MS_pop = cell(numel(brainAreas), 1);
overlap_GS_pop = cell(numel(brainAreas), 1);

for n_session = 1:numSessions
    session_date = sessions_all{n_session};
    SessionData = Data(strcmp(Data.session_date, session_date), :);
    fprintf('\n===== Session %d =====\n', n_session)

    for n_brain = 1%[1,4,5]
        frData = SessionData.frPerChannel{n_brain};
        numChannels = size(frData, 1);

        % --- Grasp-specific analysis ---
        p_values_grasp = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            all_data = [];
            group_labels = [];

            for n_grasp = 1:numel(uniqueGraspTypes)
                grasp_ind = ismember(SessionData.GraspType, uniqueGraspTypes{n_grasp});
                go_ind = cell2mat(SessionData.TrialCue) == 1;
                grasp_go_idx = grasp_ind & go_ind;

                fr_grasp = squeeze(frData(n_channel, phase_bins, grasp_go_idx));
                mean_fr = mean(fr_grasp, 1);

                all_data = [all_data; mean_fr'];
                group_labels = [group_labels; repmat(n_grasp, numel(mean_fr), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_grasp(n_channel) = p;
        end

        GS_sig_idx = find(p_values_grasp < alpha);
        GS_percent = (numel(GS_sig_idx) / numChannels) * 100;
        GS_percent_per_session{n_brain}(n_session) = GS_percent;

        % --- Modality-specific analysis ---
        p_values_modality = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            go_idx = cell2mat(SessionData.TrialCue) == 1;
            fr_combined = squeeze(frData(n_channel, :, go_idx));
            cue_names = SessionData.TrialType(go_idx);

            fr_by_modality = cellfun(@(x) fr_combined(phase_bins, ismember(cue_names, x)), uniqueCueTypes, 'UniformOutput', false);

            all_data = [];
            group_labels = [];

            for n_cue = 1:numel(uniqueCueTypes)
                fr_tmp = fr_by_modality{n_cue}';
                mean_fr = mean(fr_tmp, 2);

                all_data = [all_data; mean_fr];
                group_labels = [group_labels; repmat(n_cue, size(fr_tmp, 1), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_modality(n_channel) = p;
        end

        MS_sig_idx = find(p_values_modality < alpha);
        MS_percent = (numel(MS_sig_idx) / numChannels) * 100;
        MS_percent_per_session{n_brain}(n_session) = MS_percent;

        % --- Overlap analysis ---
        [overlap_channels, ~, ~] = intersect(MS_sig_idx, GS_sig_idx);
        overlap_pct_MS = (numel(overlap_channels) / numel(MS_sig_idx)) * 100;
        overlap_pct_GS = (numel(overlap_channels) / numel(GS_sig_idx)) * 100;

        overlap_MS_pop{n_brain}(n_session) = overlap_pct_MS;
        overlap_GS_pop{n_brain}(n_session) = overlap_pct_GS;

        % --- Output per region and session ---
        fprintf('Region: %s | Grasp: %.2f%% | Modality: %.2f%% | Overlap (MS): %.2f%% | Overlap (GS): %.2f%%\n', ...
            brainAreas{n_brain}, GS_percent, MS_percent, overlap_pct_MS, overlap_pct_GS);
    end
end

% Summary statistics
for n_brain = 1%:numel(brainAreas)
    fprintf('\n--- %s ---\n', brainAreas{n_brain});
    fprintf('Grasp: %.2f ± %.2f %%\n', mean(GS_percent_per_session{n_brain}, 'omitnan'), std(GS_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Modality: %.2f ± %.2f %%\n', mean(MS_percent_per_session{n_brain}, 'omitnan'), std(MS_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Overlap (MS pop): %.2f ± %.2f %%\n', mean(overlap_MS_pop{n_brain}, 'omitnan'), std(overlap_MS_pop{n_brain}, 'omitnan'));
    fprintf('Overlap (GS pop): %.2f ± %.2f %%\n', mean(overlap_GS_pop{n_brain}, 'omitnan'), std(overlap_GS_pop{n_brain}, 'omitnan'));
end

% plotting overlap data
% === Colors ===
graspColor         = [0, 0.3529, 0.7098];        % Blue
graspOverlapColor  = [0.0, 0.2, 0.45];           % Darker blue
modalityColor      = [0.8327, 0.1961, 0.1255];   % Red
modalityOverlapColor = [0.5, 0.1, 0.1];          % Darker red

% === Brain regions ===
brainRegions = {'SMG'};%, 'AIP', 'M1'};
regionIdx = [1];%,4,5];
x = 1:length(regionIdx);
barWidth = 0.35;

% === Compute values ===
total_GS = NaN(1, length(x));
total_MS = NaN(1, length(x));
overlap_abs_GS = NaN(1, length(x));
overlap_abs_MS = NaN(1, length(x));
GS_std = NaN(1, length(x));
MS_std = NaN(1, length(x));

for i = 1:length(regionIdx)
    r = regionIdx(i);
    
    gs_vals = GS_percent_per_session{r};
    ms_vals = MS_percent_per_session{r};
    ol_gs_vals = overlap_GS_pop{r};
    ol_ms_vals = overlap_MS_pop{r};

    total_GS(i) = mean(gs_vals, 'omitnan');
    total_MS(i) = mean(ms_vals, 'omitnan');
    GS_std(i) = std(gs_vals, 'omitnan');
    MS_std(i) = std(ms_vals, 'omitnan');

    overlap_abs_GS(i) = total_GS(i) * mean(ol_gs_vals, 'omitnan') / 100;
    overlap_abs_MS(i) = total_MS(i) * mean(ol_ms_vals, 'omitnan') / 100;
end

% === Plot ===
figure('Color','w','Position',[100, 100, 800, 450]); 
hold on;

for i = 1:length(x)
    % ----------------------
    % Grasp-specific bar (left bar)
    % ----------------------
    xG = x(i) - barWidth/2;
    totalG = total_GS(i);
    overlapG = overlap_abs_GS(i);
    gsOnly = totalG - overlapG;

    % Bottom (Overlap: Grasp + Modality)
    rectangle('Position', [xG, 0, barWidth, overlapG], ...
        'FaceColor', graspOverlapColor, 'EdgeColor', 'k');

    % Top (Grasp-only)
    rectangle('Position', [xG, overlapG, barWidth, gsOnly], ...
        'FaceColor', graspColor, 'EdgeColor', 'k');

    % Error bar centered on grasp bar
    errorbar(xG + barWidth/2, totalG, GS_std(i), 'k.', 'LineWidth', 1.5);

    % ----------------------
    % Modality-specific bar (right bar)
    % ----------------------
    xM = x(i) + barWidth/2;
    totalM = total_MS(i);
    overlapM = overlap_abs_MS(i);
    msOnly = totalM - overlapM;

    % Bottom (Overlap: Modality + Grasp)
    rectangle('Position', [xM, 0, barWidth, overlapM], ...
        'FaceColor', modalityOverlapColor, 'EdgeColor', 'k');

    % Top (Modality-only)
    rectangle('Position', [xM, overlapM, barWidth, msOnly], ...
        'FaceColor', modalityColor, 'EdgeColor', 'k');

    % Error bar centered on modality bar
    errorbar(xM + barWidth/2, totalM, MS_std(i), 'k.', 'LineWidth', 1.5);
end

set(gca, 'XTick', x, 'XTickLabel', brainRegions, 'FontSize', 12);
ylabel('Percentage of Units (%)');
ylim([0 100]);
yticks(0:25:100);
title('Grasp- and Modality-specific Units during Action Phase');

% Legend
h1 = patch(NaN, NaN, graspColor);
h2 = patch(NaN, NaN, graspOverlapColor);
h3 = patch(NaN, NaN, modalityColor);
h4 = patch(NaN, NaN, modalityOverlapColor);
legend([h2 h1 h4 h3], {'Grasp+Modality', 'Grasp-only', 'Modality+Grasp', 'Modality-only'}, ...
    'Location', 'northeastoutside');

%% GS, Size Parameters 
phase_bins = 94:174; % cue = 42:82, action = 94:174
alpha = 0.05;

% Preallocate storage
GS_percent_per_session = cell(numel(brainAreas), 1);
MS_percent_per_session = cell(numel(brainAreas), 1);
size_percent_per_session = cell(numel(brainAreas), 1);
overlap_MS_pop = cell(numel(brainAreas), 1);
overlap_GS_pop = cell(numel(brainAreas), 1);
overlap_size_pop = cell(numel(brainAreas), 1);

for n_session = 1:numSessions
    session_date = sessions_all{n_session};
    SessionData = Data(strcmp(Data.session_date, session_date), :);
    fprintf('\n===== Session %d =====\n', n_session)

    for n_brain = 1% [1,4,5] 
        frData = SessionData.frPerChannel{n_brain};
        numChannels = size(frData, 1);

        % --- Grasp-specific analysis ---
        p_values_grasp = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            all_data = [];
            group_labels = [];

            for n_grasp = 1:numel(uniqueGraspTypes)
                grasp_ind = ismember(SessionData.GraspType, uniqueGraspTypes{n_grasp});
                go_ind = cell2mat(SessionData.TrialCue) == 1;
                grasp_go_idx = grasp_ind & go_ind;

                fr_grasp = squeeze(frData(n_channel, phase_bins, grasp_go_idx));
                mean_fr = mean(fr_grasp, 1);

                all_data = [all_data; mean_fr'];
                group_labels = [group_labels; repmat(n_grasp, numel(mean_fr), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_grasp(n_channel) = p;
        end

        GS_sig_idx = find(p_values_grasp < alpha);
        GS_percent = (numel(GS_sig_idx) / numChannels) * 100;
        GS_percent_per_session{n_brain}(n_session) = GS_percent;

        % --- Size-specific analysis ---
        p_values_size = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            go_idx = cell2mat(SessionData.TrialCue) == 1;
            fr_combined = squeeze(frData(n_channel, :, go_idx));
            size_names = SessionData.Aperture_Size(go_idx);

            fr_by_size = cellfun(@(x) fr_combined(phase_bins, ismember(size_names, x)), uniqueApertureSize, 'UniformOutput', false);

            all_data = [];
            group_labels = [];

            for n_size = 1:numel(uniqueApertureSize)
                fr_tmp = fr_by_size{n_size}';
                mean_fr = mean(fr_tmp, 2);

                all_data = [all_data; mean_fr];
                group_labels = [group_labels; repmat(n_size, size(fr_tmp, 1), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_size(n_channel) = p;
        end

        size_sig_idx = find(p_values_size < alpha);
        size_percent = (numel(size_sig_idx) / numChannels) * 100;
        size_percent_per_session{n_brain}(n_session) = size_percent;

        % --- Overlap analysis ---
        [overlap_channels, ~, ~] = intersect(size_sig_idx, GS_sig_idx);
        overlap_pct_size = (numel(overlap_channels) / numel(size_sig_idx)) * 100;
        overlap_pct_GS = (numel(overlap_channels) / numel(GS_sig_idx)) * 100;

        overlap_size_pop{n_brain}(n_session) = overlap_pct_size;
        overlap_GS_pop{n_brain}(n_session) = overlap_pct_GS;

        % --- Output per region and session ---
        fprintf('Region: %s | Grasp: %.2f%% | Size: %.2f%% | Overlap (Size): %.2f%% | Overlap (GS): %.2f%%\n', ...
            brainAreas{n_brain}, GS_percent, size_percent, overlap_pct_size, overlap_pct_GS);
    end
end

% Summary statistics
for n_brain = 1%:numel(brainAreas)
    fprintf('\n--- %s ---\n', brainAreas{n_brain});
    fprintf('Grasp: %.2f ± %.2f %%\n', mean(GS_percent_per_session{n_brain}, 'omitnan'), std(GS_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Modality: %.2f ± %.2f %%\n', mean(size_percent_per_session{n_brain}, 'omitnan'), std(size_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Overlap (Size pop): %.2f ± %.2f %%\n', mean(overlap_size_pop{n_brain}, 'omitnan'), std(overlap_size_pop{n_brain}, 'omitnan'));
    fprintf('Overlap (GS pop): %.2f ± %.2f %%\n', mean(overlap_GS_pop{n_brain}, 'omitnan'), std(overlap_GS_pop{n_brain}, 'omitnan'));
end

% plotting overlap data
% === Colors ===
graspColor         = [0, 0.3529, 0.7098];        % Blue
graspOverlapColor  = [0.0, 0.2, 0.45];           % Darker blue
modalityColor      = [0.8327, 0.1961, 0.1255];   % Red
modalityOverlapColor = [0.5, 0.1, 0.1];          % Darker red
sizeColor            = [0.2039, 0.6902, 0.2902]; % Bright green 
sizeOverlapColor     = [0.1059, 0.3451, 0.1451]; % Darker green 


% === Brain regions ===
brainRegions = {'SMG'}; %, 'AIP', 'M1'};
regionIdx = [1];%, 4, 5];
x = 1:length(regionIdx);
barWidth = 0.35;

% === Compute values ===
total_GS = NaN(1, length(x));
total_size = NaN(1, length(x));
overlap_abs_GS = NaN(1, length(x));
overlap_abs_size = NaN(1, length(x));
GS_std = NaN(1, length(x));
size_std = NaN(1, length(x));

for i = 1:length(regionIdx)
    r = regionIdx(i);
    
    gs_vals = GS_percent_per_session{r};
    size_vals = size_percent_per_session{r};
    ol_gs_vals = overlap_GS_pop{r};
    ol_size_vals = overlap_size_pop{r};

    total_GS(i) = mean(gs_vals, 'omitnan');
    total_size(i) = mean(size_vals, 'omitnan');
    GS_std(i) = std(gs_vals, 'omitnan');
    size_std(i) = std(size_vals, 'omitnan');

    overlap_abs_GS(i) = total_GS(i) * mean(ol_gs_vals, 'omitnan') / 100;
    overlap_abs_size(i) = total_size(i) * mean(ol_size_vals, 'omitnan') / 100;
end

% === Plot ===
figure('Color','w','Position',[100, 100, 800, 450]); 
hold on;

for i = 1:length(x)
    % ----------------------
    % Grasp-specific bar (left bar)
    % ----------------------
    xG = x(i) - barWidth/2;
    totalG = total_GS(i);
    overlapG = overlap_abs_GS(i);
    gsOnly = totalG - overlapG;

    % Bottom (Overlap: Grasp + Modality)
    rectangle('Position', [xG, 0, barWidth, overlapG], ...
        'FaceColor', graspOverlapColor, 'EdgeColor', 'k');

    % Top (Grasp-only)
    rectangle('Position', [xG, overlapG, barWidth, gsOnly], ...
        'FaceColor', graspColor, 'EdgeColor', 'k');

    % Error bar centered on grasp bar
    errorbar(xG + barWidth/2, totalG, GS_std(i), 'k.', 'LineWidth', 1.5);

    % ----------------------
    % Size-specific bar (right bar)
    % ----------------------
    xS = x(i) + barWidth/2;
    totalS = total_size(i);
    overlapS = overlap_abs_size(i);
    sizeOnly = totalS - overlapS;

    % Bottom (Overlap: Size + Grasp)
    rectangle('Position', [xS, 0, barWidth, overlapS], ...
        'FaceColor', sizeOverlapColor, 'EdgeColor', 'k');

    % Top (Modality-only)
    rectangle('Position', [xS, overlapS, barWidth, sizeOnly], ...
        'FaceColor', sizeColor, 'EdgeColor', 'k');

    % Error bar centered on modality bar
    errorbar(xS + barWidth/2, totalS, size_std(i), 'k.', 'LineWidth', 1.5);
end

set(gca, 'XTick', x, 'XTickLabel', brainRegions, 'FontSize', 12);
ylabel('Percentage of Units (%)');
ylim([0 100]);
yticks(0:25:100);
title('Grasp- and Size-specific Units during Action Phase');

% Legend
h1 = patch(NaN, NaN, graspColor);
h2 = patch(NaN, NaN, graspOverlapColor);
h3 = patch(NaN, NaN, sizeColor);
h4 = patch(NaN, NaN, sizeOverlapColor);
legend([h2 h1 h4 h3], {'Grasp+Size', 'Grasp-only', 'Size+Grasp', 'Size-only'}, ...
    'Location', 'northeastoutside');
%% GS, Object Parameters 
phase = 'Cue';
%phase = 'Action';
if strcmp(phase, 'Cue')
    phase_bins = 42:82; 
elseif strcmp(phase, 'Action')
    phase_bins = 94:174; 
else 
    keyboard
    error([phase ' not added to pipeline, must add timebins'])
end
alpha = 0.05;

% Preallocate storage
GS_percent_per_session = cell(numel(brainAreas), 1);
MS_percent_per_session = cell(numel(brainAreas), 1);
object_percent_per_session = cell(numel(brainAreas), 1);
overlap_MS_pop = cell(numel(brainAreas), 1);
overlap_GS_pop = cell(numel(brainAreas), 1);
overlap_object_pop = cell(numel(brainAreas), 1);

for n_session = 1:numSessions
    session_date = sessions_all{n_session};
    SessionData = Data(strcmp(Data.session_date, session_date), :);
    fprintf('\n===== Session %d =====\n', n_session)

    for n_brain = [1,4,5] 
        frData = SessionData.frPerChannel{n_brain};
        numChannels = size(frData, 1);

        % --- Grasp-specific analysis ---
        p_values_grasp = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            all_data = [];
            group_labels = [];

            for n_grasp = 1:numel(uniqueGraspTypes)
                grasp_ind = ismember(SessionData.GraspType, uniqueGraspTypes{n_grasp});
                go_ind = cell2mat(SessionData.TrialCue) == 1;
                grasp_go_idx = grasp_ind & go_ind;

                fr_grasp = squeeze(frData(n_channel, phase_bins, grasp_go_idx));
                mean_fr = mean(fr_grasp, 1);

                all_data = [all_data; mean_fr'];
                group_labels = [group_labels; repmat(n_grasp, numel(mean_fr), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_grasp(n_channel) = p;
        end

        GS_sig_idx = find(p_values_grasp < alpha);
        GS_percent = (numel(GS_sig_idx) / numChannels) * 100;
        GS_percent_per_session{n_brain}(n_session) = GS_percent;

        % --- Object-specific analysis ---
        p_values_object = NaN(numChannels, 1);

        for n_channel = 1:numChannels
            go_idx = cell2mat(SessionData.TrialCue) == 1;
            fr_combined = squeeze(frData(n_channel, :, go_idx));
            object_labels = {'deck','block','rod','ball'};
            object_names = SessionData.ObjectType(go_idx);

            fr_by_object = cellfun(@(x) fr_combined(phase_bins, ismember(object_names, x)), uniqueObjectType, 'UniformOutput', false);

            all_data = [];
            group_labels = [];

            for n_object = 1:numel(uniqueObjectType)
                fr_tmp = fr_by_object{n_object}';
                mean_fr = mean(fr_tmp, 2);

                all_data = [all_data; mean_fr];
                group_labels = [group_labels; repmat(n_object, size(fr_tmp, 1), 1)];
            end

            [p, ~] = anova1(all_data, group_labels, 'off');
            p_values_object(n_channel) = p;
        end

        object_sig_idx = find(p_values_object < alpha);
        object_percent = (numel(object_sig_idx) / numChannels) * 100;
        object_percent_per_session{n_brain}(n_session) = object_percent;

        % --- Overlap analysis ---
        [overlap_channels, ~, ~] = intersect(object_sig_idx, GS_sig_idx);
        overlap_pct_object = (numel(overlap_channels) / numel(object_sig_idx)) * 100;
        overlap_pct_GS = (numel(overlap_channels) / numel(GS_sig_idx)) * 100;

        overlap_object_pop{n_brain}(n_session) = overlap_pct_object;
        overlap_GS_pop{n_brain}(n_session) = overlap_pct_GS;

        % --- Output per region and session ---
        fprintf('Region: %s | Grasp: %.2f%% | Object: %.2f%% | Overlap (Object): %.2f%% | Overlap (GS): %.2f%%\n', ...
            brainAreas{n_brain}, GS_percent, object_percent, overlap_pct_object, overlap_pct_GS);
    end
end

% Summary statistics
for n_brain = 1:numel(brainAreas)
    fprintf('\n--- %s ---\n', brainAreas{n_brain});
    fprintf('Grasp: %.2f ± %.2f %%\n', mean(GS_percent_per_session{n_brain}, 'omitnan'), std(GS_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Modality: %.2f ± %.2f %%\n', mean(object_percent_per_session{n_brain}, 'omitnan'), std(object_percent_per_session{n_brain}, 'omitnan'));
    fprintf('Overlap (Size pop): %.2f ± %.2f %%\n', mean(overlap_object_pop{n_brain}, 'omitnan'), std(overlap_object_pop{n_brain}, 'omitnan'));
    fprintf('Overlap (GS pop): %.2f ± %.2f %%\n', mean(overlap_GS_pop{n_brain}, 'omitnan'), std(overlap_GS_pop{n_brain}, 'omitnan'));
end

% plotting overlap data
% === Colors ===
graspColor         = [0, 0.3529, 0.7098];        % Blue
graspOverlapColor  = [0.0, 0.2, 0.45];           % Darker blue
modalityColor      = [0.8327, 0.1961, 0.1255];   % Red
modalityOverlapColor = [0.5, 0.1, 0.1];          % Darker red
sizeColor            = [0.2039, 0.6902, 0.2902]; % Bright green 
sizeOverlapColor     = [0.1059, 0.3451, 0.1451]; % Darker green 
objectColor          = [0.9608, 0.7333, 0.1176]; % Bright yellow
objectOverlapColor   = [0.5294, 0.3961, 0.0627]; % Darker yellow 


% === Brain regions ===
brainRegions = {'SMG', 'AIP', 'M1'};
regionIdx = [1, 4, 5];
x = 1:length(regionIdx);
barWidth = 0.35;

% === Compute values ===
total_GS = NaN(1, length(x));
total_object = NaN(1, length(x));
overlap_abs_GS = NaN(1, length(x));
overlap_abs_object = NaN(1, length(x));
GS_std = NaN(1, length(x));
object_std = NaN(1, length(x));

for i = 1:length(regionIdx)
    r = regionIdx(i);
    
    gs_vals = GS_percent_per_session{r};
    object_vals = object_percent_per_session{r};
    ol_gs_vals = overlap_GS_pop{r};
    ol_object_vals = overlap_object_pop{r};

    total_GS(i) = mean(gs_vals, 'omitnan');
    total_object(i) = mean(object_vals, 'omitnan');
    GS_std(i) = std(gs_vals, 'omitnan');
    object_std(i) = std(object_vals, 'omitnan');

    overlap_abs_GS(i) = total_GS(i) * mean(ol_gs_vals, 'omitnan') / 100;
    overlap_abs_object(i) = total_object(i) * mean(ol_object_vals, 'omitnan') / 100;
end

% === Plot ===
figure('Color','w','Position',[100, 100, 800, 450]); 
hold on;

for i = 1:length(x)
    % ----------------------
    % Grasp-specific bar (left bar)
    % ----------------------
    xG = x(i) - barWidth/2;
    totalG = total_GS(i);
    overlapG = overlap_abs_GS(i);
    gsOnly = totalG - overlapG;

    % Bottom (Overlap: Grasp + Modality)
    rectangle('Position', [xG, 0, barWidth, overlapG], ...
        'FaceColor', graspOverlapColor, 'EdgeColor', 'k');

    % Top (Grasp-only)
    rectangle('Position', [xG, overlapG, barWidth, gsOnly], ...
        'FaceColor', graspColor, 'EdgeColor', 'k');

    % Error bar centered on grasp bar
    errorbar(xG + barWidth/2, totalG, GS_std(i), 'k.', 'LineWidth', 1.5);

    % ----------------------
    % Object-specific bar (right bar)
    % ----------------------
    xObj = x(i) + barWidth/2;
    totalObj = total_object(i);
    overlapObj = overlap_abs_object(i);
    objectOnly = totalObj - overlapObj;

    % Bottom (Overlap: Object + Grasp)
    rectangle('Position', [xObj, 0, barWidth, overlapObj], ...
        'FaceColor', objectOverlapColor, 'EdgeColor', 'k');

    % Top (Modality-only)
    rectangle('Position', [xObj, overlapObj, barWidth, objectOnly], ...
        'FaceColor', objectColor, 'EdgeColor', 'k');

    % Error bar centered on modality bar
    errorbar(xObj + barWidth/2, totalObj, object_std(i), 'k.', 'LineWidth', 1.5);
end

set(gca, 'XTick', x, 'XTickLabel', brainRegions, 'FontSize', 12);
ylabel('Percentage of Units (%)');
ylim([0 100]);
yticks(0:25:100);
title(['Grasp- and Object-specific Units during ' phase]);

% Legend
h1 = patch(NaN, NaN, graspColor);
h2 = patch(NaN, NaN, graspOverlapColor);
h3 = patch(NaN, NaN, objectColor);
h4 = patch(NaN, NaN, objectOverlapColor);
legend([h2 h1 h4 h3], {'Grasp+Object', 'Grasp-only', 'Object+Grasp', 'Object-only'}, ...
    'Location', 'northeastoutside');