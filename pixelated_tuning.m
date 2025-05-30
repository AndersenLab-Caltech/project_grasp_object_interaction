clc
clear all
close all

subject_id = 's2';
unit_region = 'SMG';

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
taskName = 'GraspObject_Shuffled'; % shuffled images
flag_4S = true; % true = updated 4S action phase; false = original 2S action phase
flag_shuffled = true; % true = shuffled images task

if ~flag_4S
    TaskCue = 'GraspObject';
    min_timebin_length = 134; % NOT VALID FOR 20230831    
elseif ~flag_shuffled
    TaskCue = 'GraspObject_4S_Action';
    min_timebin_length = 174;
else
    TaskCue = 'GraspObject_Shuffled';
    min_timebin_length = 174; 
end 

% load data
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' TaskCue spike_sorting_type]);
Go_data = Data.Go_data;

% remove faulty sessions, if any
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20231016'};
elseif strcmp(subject_id, 's3')
    error_session = {'20231207','20231212'}; % tester sessions with varying # of trials
end 

if ~isempty(error_session)
    %condition = cellfun(@(x) strcmp(x, error_session), Go_data.session_date);
    condition = ismember(Go_data.session_date, error_session);
    Go_data = Go_data(~condition,:);
end

flag_Shuffled = false; % true = analyze pixelated images

flagRegressionTuning = true;

if flagRegressionTuning
    analysis_type = 'LinearRegression';
else
    analysis_type = 'KruskalWallis';
end

flagBinPerBin = true;
multipleComparePhase = true;
flagTunedChannels = true;
flagSaveData = true;

% add Shuffled column
% find indices of shuffled
Shuffled_idx = find(contains(Go_data.LabelNames, 'Shuffled'));
% create logical vector
shuffled = false(size(Go_data.LabelNames));
shuffled(Shuffled_idx) = true;
% add column to Data table
Go_data.Shuffled = shuffled;

% chose cue type:
% taskCuesAll = {'Regular Images', 'Pixelated Images'}; 
taskCuesAll = {'Hand', 'Hand-Object', 'Object'};
sessions_all = unique(Go_data.session_date);
numSessions = numel(sessions_all);
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

numUnitsPerSession = zeros(numSessions,1);
%% Analysis
for n_session = 1:numSessions

    disp(['Classification session ' sessions_all{n_session} ]);  

    %find idx of current session day
    idxThisSession = ismember(Go_data.session_date, sessions_all(n_session));
        

    %extract data from selected brain area

     if strcmp('SMG', unit_region)
            SessionData = Go_data.SMG_Go(idxThisSession,:);
        elseif strcmp('PMV', unit_region) 
            SessionData = Go_data.PMV_Go(idxThisSession,:);
        elseif strcmp('S1', unit_region)
            SessionData = Go_data.S1X_Go(idxThisSession,:);
         elseif strcmp('M1', unit_region) 
            SessionData = Go_data.M1_Go(idxThisSession,:);
        elseif strcmp('AIP', unit_region)
            SessionData = Go_data.AIP_Go(idxThisSession,:);
        else
            error([unit_region ' does not exist '])
     end

     % skip session days that are empty - relevant for S1 session 20230810
     if isempty(SessionData{1})
        continue
     end
    
    %labels 
    sessionLabels = Go_data.GoLabels(idxThisSession,:);
    
    %trialType
    trialTypeSession = Go_data.TrialType(idxThisSession,:);

    %get idx for shuffled or unshuffled trials
    ShuffledImage_idx = Go_data.Shuffled(idxThisSession,:);

    timePhaseLabels = Go_data.time_phase_labels(idxThisSession);


    if flag_Shuffled
        SessionData = SessionData(ShuffledImage_idx);
        sessionLabels = sessionLabels(ShuffledImage_idx);
        timePhaseLabels = timePhaseLabels(ShuffledImage_idx);
        trialTypeSession = trialTypeSession(ShuffledImage_idx);
        taskImageType = {'Pixelated Images'};
    else
        SessionData = SessionData(~ShuffledImage_idx);
        sessionLabels = sessionLabels(~ShuffledImage_idx);
        timePhaseLabels = timePhaseLabels(~ShuffledImage_idx);
        trialTypeSession = trialTypeSession(~ShuffledImage_idx);
        taskImageType = 'Regular Images';
    end
     
    %seperate data according to cue modality => Shuffled vs not

    unTrialType = unique(Go_data.TrialType);
    unType = unique(Go_data.Shuffled);
    numUnitsPerSession(n_session) = size(SessionData{1},2);
    % loop through cue modalities 
    for n_type = 1:numel(unTrialType) % unType for Shuffled vs Regular; unTrialType for H,H+O,H

        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
        %imageTypeIdx = ismember(imageType, unType(n_type)); % if
        %separating data by only shuffled v reg

         if flagTunedChannels
            %Compute index of units that are tuned
            if flagRegressionTuning
            
                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin,numTunedChannelsPerCategory,~,~,p_per_phase] ...
                     = classification.getRegressionTunedChannels_paper(SessionData(trialTypeIdx),sessionLabels(trialTypeIdx), ...
                 timePhaseLabels(trialTypeIdx), 'multcompare', multipleComparePhase, 'BinperBinTuning', flagBinPerBin);

                condToTest = arrayfun(@(x) preproc.image2class_simple(x),  unique(sessionLabels), 'UniformOutput', false);

                if nnz(sumBin) ~= 0
                    figure();
                    plot(sumBin);
                end

                tuned_channels_per_type{n_type,n_session} = numTunedChannelsPerCategory;
            else

                tuned_channels_per_graps{n_type,n_session} = [];
                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin]= classification.getTunedChannels(SessionData(trialTypeIdx),sessionLabels(trialTypeIdx), ...
                timePhaseLabels(trialTypeIdx), 'multcompare', multipleComparePhase,'removeITItuning', 'false', 'BinperBinTuning', flagBinPerBin);
                sumBin = sumBin';

                % tuned_channels_per_type{n_type,n_session} = [];
                % [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin]= classification.getTunedChannels(SessionData(imageTypeIdx),sessionLabels(imageTypeIdx), ...
                % timePhaseLabels(imageTypeIdx), 'multcompare', multipleComparePhase,'removeITItuning', 'false', 'BinperBinTuning', flagBinPerBin);
                % sumBin = sumBin';
            end

                if nnz(sumBin) > 0
                sum_bin_all{n_type, n_session } = sumBin;

            else
                sum_bin_all{n_type, n_session } = [];

                end

            tuned_channels_per_phase{n_type,n_session} = sumPhase;
            tuned_channels_per_phase_vector{n_type,n_session} = tunedChannelsPhase;
         
         end
    end 
       
end 

shufLabel = ["Regular", "Pixelated"];
shufLabel = shufLabel(flag_Shuffled + 1);

% Create the filename using the brain region and analysis type
filename = "tuned_channels_" + TaskCue + '_' + unit_region + "_" + analysis_type + "_" + shufLabel + ".mat";

directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
full_path = fullfile(directory, filename);

% Save the relevant variables with the dynamic filename
save(full_path, 'sum_bin_all', 'tuned_channels_per_phase', 'tuned_channels_per_phase_vector','numUnitsPerSession');

%% load analyzed data
shufLabel = ["Regular", "Pixelated"];
shufLabel = shufLabel(flag_Shuffled + 1);

directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
analysis_type = 'LinearRegression'; % 'LinearRegression' or 'KW'
filename = "tuned_channels_" + TaskCue + '_' + unit_region + "_" + analysis_type + "_" + shufLabel + ".mat";
full_path = fullfile(directory, filename);
load(full_path);
%% PLOTS
if flag_Shuffled
    taskImageType = {'Pixelated Images'};
else
    taskImageType = 'Regular Images';
end

% bar graph w/ CIs
phase_yCI95 = [];
phase_tuned_mean = [];

sessionToInclude = 1:numel(numUnitsPerSession);

figure('units','normalized','outerposition',[0 0 .3 0.45]);
for n_type = 1:numel(taskCuesAll)
    dataTmp = cell2mat(tuned_channels_per_phase(n_type,sessionToInclude)')*100;
    percentage_tuned = dataTmp./numUnitsPerSession(sessionToInclude);
    yCI95tmp = utile.calculate_CI(percentage_tuned);
    phase_yCI95(n_type,:) = yCI95tmp(2,:);
    phase_tuned_mean(n_type,:) = mean(percentage_tuned,1);
    % figure()
    subplot(1,numel(taskCuesAll),n_type)

    hold on
    bar(phase_tuned_mean(n_type,:),'FaceColor',color_info{n_type});

    hold on
    errorbar(phase_tuned_mean(n_type,:),phase_yCI95(n_type,:),'Color','k');

    title(taskCuesAll(n_type));
    xticks(1:numel(phaseNames));
    xticklabels(phaseNames);
    xtickangle(45);
    ylabel('% of Tuned Units');
    ylim([0 70]);
    sgtitle(['Tuned Units in ' unit_region ' -- ' taskImageType])
    set(gca, 'FontSize', 12);
end

% line plot w/o CIs bc little data rn => update when we have more data
% for n_type = 1:numel(unTrialType) % unType for Shuffled v reg; unTrialType for modality
%     tunedUnitsPerTypeBin(n_type,:)  = sum(cell2mat(sum_bin_all(n_type,:)),2);
% 
% end 
% 
% figure('units','normalized','outerposition',[0 0 0.3 0.45]);
% plot((tunedUnitsPerTypeBin'./sum(numUnitsPerSession))*100,'LineWidth',2);
% hold on
% for n_phase = 1:numPhases
%     xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5,'FontSize',12);
% end
% title(['Tuned Units Throughout Trial in ' unit_region ' -- ' taskImageType]);
% xlabel('Time Bins (50 ms)');
% xlim([0 (min_timebin_length + 5)])
% ylabel('% of Total Units');
% ylim([0 100]);
% legend(taskCuesAll, 'Location', 'Best','FontSize',12);
% set(gca, 'FontSize', 12);
% hold off
%% plotting only Cue phase
if flag_Shuffled
    taskImageType = {'Pixelated Images'};
else
    taskImageType = 'Regular Images';
end

% Extract Cue phase index
cueIdx = find(strcmp(phaseNames, 'Cue')); % Adjust if needed

sessionToInclude = 1:numel(numUnitsPerSession);
%sessionToInclude = sessionToInclude(2:end); % take out first session (maybe bc novel, barely any tuning)
numConditions = numel(taskCuesAll);

% Preallocate arrays
phase_yCI95 = zeros(numConditions, 1);
phase_tuned_mean = zeros(numConditions, 1);

% Collect data for Cue phase
figure('units','normalized','outerposition',[0 0 .13 0.25]); 
hold on

for n_type = 1:numel(taskCuesAll)

    dataTmp = cell2mat(tuned_channels_per_phase(n_type, sessionToInclude)') * 100;
    percentage_tuned = dataTmp ./ numUnitsPerSession(sessionToInclude);
    
    % Compute mean and CI
    yCI95tmp = utile.calculate_CI(percentage_tuned(:, cueIdx)); 
    phase_yCI95(n_type) = yCI95tmp(2);
    phase_tuned_mean(n_type) = mean(percentage_tuned(:, cueIdx), 1);

    % Bar plot for mean
    bar(n_type, phase_tuned_mean(n_type), 'FaceColor', color_info{n_type}, 'BarWidth', 0.6);
    
    % Error bars for CI
    errorbar(n_type, phase_tuned_mean(n_type), phase_yCI95(n_type), 'k', 'LineWidth', 1.5);
    
    % Scatter plot of individual session data points
    scatter(repmat(n_type, numel(sessionToInclude), 1), percentage_tuned(:, cueIdx), 50, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Formatting
xticks(1:numConditions);
xticklabels([]);%taskCuesAll);
xtickangle(45);
%ylabel('% of Tuned Units');
ylim([0 70]);
yticks([0 35 70]);
title(unit_region); %' -- ' phaseNames{2} ' -- ' taskImageType])
set(gca, 'FontSize', 12);

hold off

%% for line plot w/ 95% CI
per_bin_yCI95 = [];
per_bin_tuned_mean = [];
%sessionToInclude = setdiff(1:numSessions,1);

% code for empty/missing session data
rowsToKeep = numUnitsPerSession ~= 0;
numUnitsPerSession = numUnitsPerSession(rowsToKeep);
sessionToInclude = 1:numel(numUnitsPerSession);
colsToKeep = true(1,numSessions);
for n_session = 1:numSessions
    if all(cellfun('isempty',sum_bin_all(:,n_session)))
        colsToKeep(n_session) = false;
    end
end
sum_bin_all = sum_bin_all(:,colsToKeep);

figure('units','normalized','outerposition',[0 0 0.3 0.5]);
err_bar = {};
for n_type = 1:numel(taskCuesAll)
    dataTmp = cell2mat(sum_bin_all(n_type,sessionToInclude))*100;
    percentage_tuned = dataTmp./(numUnitsPerSession(sessionToInclude)');
    yCI95tmp = utile.calculate_CI(percentage_tuned');
    per_bin_yCI95(n_type,:) = yCI95tmp(2,:);
    per_bin_tuned_mean(n_type,:) = mean(percentage_tuned,2);

    hold on
    err_bar{n_type} = plot(1:length(dataTmp),per_bin_tuned_mean(n_type,:),'Color', color_info{n_type},'LineWidth',2);

    ER = utile.shadedErrorBar(1:length(dataTmp),per_bin_tuned_mean(n_type,:),per_bin_yCI95(n_type,:));
    ER.mainLine.Color = color_info{n_type};
    ER.patch.FaceColor = color_info{n_type};
    ER.edge(1).Color = color_info{n_type};
    ER.edge(2).Color = color_info{n_type};
end

for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5,'FontSize',12);
end

title(append('Tuned Units Throughout Trial in ', unit_region,' - ',taskImageType));
xlabel('Time Bins (50 ms)');
xlim([0 (min_timebin_length + 5)]) % 5 chosen as a buffer
ylabel('% of Total Units');
ylim([0 70]);
legend([err_bar{:}], taskCuesAll,'Location', 'Best','Interpreter', 'none','FontSize',12);
set(gca, 'FontSize', 12);