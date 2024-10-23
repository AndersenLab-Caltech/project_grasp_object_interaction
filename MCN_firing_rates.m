% Plot firing rates of individual neurons

% Clear command window, workspace, and close all figures
clc
clear all
close all 

% Define spike sorting type, task name, subject ID, and session date
spike_sorting_type = 'unsorted_aligned_thr_-4.5';
taskName = 'GraspObject_4S_Action';
subject_id = 's3';
session_date = {'20231030'}; 

% Load the data file based on subject ID, task name, and spike sorting type
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName '_' spike_sorting_type]);
Data = Data.Go_data;

% Filter the data to keep only the desired session
Data = Data(strcmp(Data.session_date, session_date), :); 

% Extract brain area information and phase indices from the data
brainAreas = Data.frPerChannel{6};
phase_time_idx = Data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

% Identify unique grasp types and cue types from the data
uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);

% Begin analysis of firing rates for each brain area
for n_brain = 1 % Looking only at SMG
    
    % Extract firing rate data for the current brain area
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    % Iterate through each channel in the brain area
    for n_channel = 1:numChannels
        figure('units','normalized','outerposition',[0 0 0.5 1])
        sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
        
        % Iterate through each unique grasp type
        for n_grasp = 1:numel(uniqueGraspTypes)
            
            % Find indices for the current grasp type
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            
            % Find indices for GO cue trials
            go_ind = cell2mat(Data.TrialCue) == 1;
            
            % Combine grasp type and GO cue indices
            grasp_go_idx = logical(grasp_ind .* go_ind);
            
            % Extract firing rate data for the current grasp type and channel
            fr_grasp = squeeze(frData(n_channel,:,grasp_go_idx)); 
            
            % Get the cue type names for the selected trials
            CueType_name = Data.TrialType(grasp_go_idx);

            % Calculate the mean firing rate for each cue type

            % Step 1: Initialize a cell array to store mean firing rates for each cue type
            mean_firing_rate_per_cue_type = cell(size(uniqueCueTypes));

            % Step 2: Iterate through each unique cue type
            for n_cueType = 1:numel(uniqueCueTypes)
                
                % Step 3: Find the indices for the current cue type
                cue_type_ind = ismember(CueType_name, uniqueCueTypes{n_cueType});

                % Step 4: Extract firing rates corresponding to the current cue type
                fr_for_current_cue_type = fr_grasp(:, cue_type_ind);

                % Step 5: Compute the mean firing rate across trials for this cue type
                mean_firing_rate_per_cue_type{n_cueType} = mean(fr_for_current_cue_type, 2);
            end
            
            % Convert cell array to a matrix for easier handling
            fr_sep_cue_type_mean = cell2mat(mean_firing_rate_per_cue_type');

            % Extract firing rate data separated by cue type for later use
            fr_sep_cue_type_trial = cell(size(uniqueCueTypes));
            for n_cueType = 1:numel(uniqueCueTypes)
                fr_sep_cue_type_trial{n_cueType} = fr_grasp(:, ismember(CueType_name, uniqueCueTypes{n_cueType}));
            end

            % Plot firing rates for each phase of the trial
            subplot(numel(uniqueGraspTypes), 1, n_grasp);
            hold on;
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
            end
            
            % Initialize variables for plotting error bars
            err_bar = {};
            color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

            % Plot firing rates separated by cue type
            for n_cueType = 1:numel(uniqueCueTypes)
                % Extract data for current cue type
                dataTmp = fr_sep_cue_type_trial{n_cueType}';

                % Calculate confidence intervals using bootstrapping
                ci = bootci(1000, {@mean,dataTmp});
                Mean_FR = squeeze(mean(dataTmp));
                err_ci(1,:) = ci(2,:) - Mean_FR; 
                err_ci(2,:) = Mean_FR - ci(1,:); 

                % Plot the firing rates with error bars
                ER = utile.shadedErrorBar(1:length(dataTmp), mean(dataTmp), err_ci);
                hold on
                err_bar{n_cueType} = plot(1:length(dataTmp), mean(dataTmp), 'Color', color_info{n_cueType}, 'LineWidth', 2);

                % Set plot colors for the shaded error bars
                ER.mainLine.Color = color_info{n_cueType};
                ER.patch.FaceColor = color_info{n_cueType};
                ER.edge(1).Color = color_info{n_cueType};
                ER.edge(2).Color = color_info{n_cueType};
            end 

            % Set plot title and legend
            title([uniqueGraspTypes{n_grasp} ' Grasp']);
            legend([err_bar{:}], uniqueCueTypes, 'Interpreter', 'none');
            set(gca, 'FontSize', 12)
            hold off;
        end
        
        % Label the x and y axes
        xlabel('Time');
        ylabel('Average Firing Rate');
        set(gca, 'FontSize', 12);
    end
end
