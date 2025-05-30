% plot firing rates of individual neurons per size
clc
clear all
close all %- closes all the figures

% saveFolder = 'C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data';
% 
spike_sorting_type = 'unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
taskName = 'GraspObject_Varied_Size'; % varied object/aperature sizes 
subject_id = 's3';
session_date = {'20240521'}; 
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

% add Aperature Size column
sizeKeywords = ['Small', 'Medium', 'Large'];
Data.Aperature_Size = cell(height(Data),1);
% Loop through each label and extract the size information
for i = 1:height(Data)
    % Use regular expression to find the size keyword after the last underscore
    tokens = regexp(Data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');
    
    if ~isempty(tokens)
        % tokens is a cell array; extract the size keyword from it
        Data.Aperature_Size{i} = tokens{1}{1};
    end
end

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

brainAreas = Data.frPerChannel{7};
phase_time_idx = Data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);
uniqueAperatureSize = unique(Data.Aperature_Size);

%%
keyboard

%%
% analyzing fr for each grasp separated by modality
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
            CueType_name = Data.TrialType(grasp_go_idx);
            fr_sep_cue_type_mean = cell2mat(cellfun(@(x) mean(fr_grasp(:,ismember(CueType_name, x)), 2), uniqueCueTypes, 'UniformOutput', false)');
            fr_sep_cue_type_trial = cellfun(@(x) fr_grasp(:,ismember(CueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

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