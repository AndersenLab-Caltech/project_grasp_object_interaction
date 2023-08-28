% plot firing rates of individual neurons
clc
clear all
%close all - closes all the figures

% saveFolder = 'C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data';
% 
spike_sorting_type = 'unsorted_aligned_thr_-4.5';
taskName = 'GraspObject';
subject_id = 's2';
session_dates = {'20230824'}; 
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

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');
Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s2\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s2_20230824_unsorted_aligned_thr_-4.5_GraspObject');
Data = Data.Go_data;

brainAreas = Data.frPerChannel{6};
phase_time_idx = Data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
blub = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(blub) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);

for n_brain = 1:3 %:length(brainAreas) 1:5 for AN, 1:3 for FG
    
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
            cueType_name = Data.TrialType(grasp_go_idx);
            fr_sep_cue_type_mean = cell2mat(cellfun(@(x) mean(fr_grasp(:,ismember(cueType_name, x)), 2), uniqueCueTypes, 'UniformOutput', false)');
            fr_sep_cue_type_trial = cellfun(@(x) fr_grasp(:,ismember(cueType_name, x)), uniqueCueTypes, 'UniformOutput', false);

            subplot(numel(uniqueGraspTypes), 1, n_grasp);
            hold on;
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
            end
            
            %chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
            err_bar = {};
            for n_cueType = 1:numel(uniqueCueTypes)
                dataTmp = fr_sep_cue_type_trial{n_cueType}';
                N = size(dataTmp, 1);
                sem = std(dataTmp) / sqrt(N);  % standard error of the mean
                CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
                yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
           



                %%figure(); 
                %subplot(2,1,1)
               % plot(1:length(dataTmp), mean(dataTmp), 'LineWidth', 2);
                %hold on 
               % plot(1:length(dataTmp), mean(dataTmp) +yCI95 , 'LineWidth', 2);
                %subplot(2,1,2)
                errorbar(mean(dataTmp), yCI95);
                hold on 
                err_bar{n_cueType} = plot(1:length(dataTmp), mean(dataTmp), 'LineWidth',2);
            end 

            %ylim([0 max(chan_fr)]);
            title([uniqueGraspTypes{n_grasp} ' Grasp']);
            hold off;
        end
        xlabel('Time');
        ylabel('Average Firing Rate');
        legend([err_bar{:}],uniqueCueTypes, 'Location', 'Best', 'Interpreter', 'none');
    end
end


%% CI example

% x = 1:100;                                          % Create Independent Variable
% y = randn(50,100);                                  % Create Dependent Variable ‘Experiments’ Data
% 
% N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
% yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
% ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
% 
% CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
% yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% 
% figure
% plot(x, yMean)                                      % Plot Mean Of All Experiments
% hold on
% plot(x, yCI95+yMean)                                % Plot 95% Confidence Intervals Of All Experiments
% hold off
% grid