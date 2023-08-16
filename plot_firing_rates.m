% plot firing rates of individual neurons
clc
clear all
%close all - closes all the figures

% saveFolder = 'C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data';
% 
% spike_sorting_type = '_unsorted_aligned_thr_-4.5';
% taskName = 'GraspObject';
% subject_id = 's3';
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


Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');


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

for n_brain = 1%:length(brainAreas)
    
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    
    for n_channel = 1:numChannels
        figure;
        sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            
            grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
            
            go_ind = cell2mat(Data.TrialCue) == 1;
            grasp_go_idx = logical(grasp_ind .* go_ind);
            
            fr_grasp = frData(:,:,grasp_go_idx);
            cueType_name = Data.TrialType(grasp_go_idx);
            fr_sep_cue_type_mean = cellfun(@(x) mean(fr_grasp(:,:,ismember(cueType_name, x)), 3), uniqueCueTypes, 'UniformOutput', false);
            
            subplot(numel(uniqueGraspTypes), 1, n_grasp);
            hold on;
            for n_phase = 1:numPhases
                xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
            end
            
            chan_fr = cell2mat(cellfun(@(x) x(n_channel,:), fr_sep_cue_type_mean, 'UniformOutput',false));
            plot(chan_fr', 'LineWidth', 2);
            
            title([uniqueGraspTypes{n_grasp} ' Grasp']);
            hold off;
        end
        xlabel('Time');
        ylabel('Average Firing Rate');
        legend(uniqueCueTypes, 'Location', 'Best', 'Interpreter', 'none');
    end
end
