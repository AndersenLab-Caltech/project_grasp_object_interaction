clc
clear all
close all 

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
taskName = 'GraspObject_Varied_Size'; % varied object/aperture sizes
%taskName = 'GraspObject_GB_Images'; % for GB
%taskName = 'GraspObject_Combined'; % for Combined task
subject_id = 's2';

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName spike_sorting_type]);

% keyboard

Go_data = Data.Go_data;

% remove faulty sessions, if any
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20231016'};
elseif strcmp(subject_id, 's3')
    error_session = {'20250212'};
elseif strcmp(subject_id, 's4')
    error_session = {'20240613'};
end 

if ~isempty(error_session)
    condition = cellfun(@(x) strcmp(x, error_session), Go_data.session_date);
    Go_data = Go_data(~condition,:);
end

if strcmp(taskName, 'GraspObject_Varied_Size')
    % add Aperture Size column
    sizeKeywords = ['Small', 'Medium', 'Large'];
    Go_data.Aperture_Size = cell(height(Go_data),1);
    % Loop through each label and extract the size information
    for i = 1:height(Go_data)
        % Use regular expression to find the size keyword after the last underscore
        tokens = regexp(Go_data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');
        
        if ~isempty(tokens)
            % tokens is a cell array; extract the size keyword from it
            Go_data.Aperture_Size{i} = tokens{1}{1};
        end
    end
end
if strcmp(taskName, 'GraspObject_Combined')
    Go_data.TrialType(strcmp(Go_data.TrialType, 'Unknown')) = {'Combined'}; % adds in Combined as Trial type
    % add in column with Object Type for Combined trials and original trial types (H, HO, O with Associated)
    % Loop through each label and extract the object information
    for i = 1:height(Go_data)
        % Use regular expression to find the size keyword after the last underscore
        tokens = regexp(Go_data.LabelNames{i}, '_(deck|block|rod|ball)$', 'tokens');
        
        if ~isempty(tokens)
            % tokens is a cell array; extract the size keyword from it
            Go_data.ObjectType{i} = tokens{1}{1};
        else 
            Go_data.ObjectType{i} = 'Associated';
        end
    end
    %color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
end
brainAreas = Go_data.frPerChannel{6}; % 7 for Ripple, 6 for Blackrock
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phaseTimeTmp = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
phaseNames = {'ITI', 'Planning', 'Delay', 'Action'};
% color_info = {[0.3359, 0.7031, 0.9101],[0.8984 0.6211 0],[0.8320 0.3672 0],[0.7969, 0.4726, 0.6523],[0, 0.6171, 0.4492]}; % SMG, PMV, S1, AIP, M1
color_info = {[0.3359, 0.7031, 0.9101],[0.8984 0.6211 0],[0.8320 0.3672 0],[0.7969, 0.4726, 0.6523],[0, 0.6171, 0.4492],[.9961 .6875 0]}; % SMG, PMV, S1, AIP, M1, dlPFC
sessions_all = unique(Go_data.session_date);
numSessions = numel(sessions_all);
uniqueCueTypes = {'Hand','Hand-Object','Object'};
if strcmp(taskName, 'GraspObject_Combined')
    uniqueCueTypes = {'Combined','Hand','Hand-Object','Object'};
end

flagGoTrials = true; 

n_regions = length(brainAreas);
num_timebins = 174;

keyboard

%% analysis and plot for all brain regions together
% Initialize storage
all_errTest_timebin_all_regions = NaN(numSessions, num_timebins, n_regions);
first_sig_idx_all = NaN(1, n_regions);
first_sig_perc_all = NaN(1, n_regions);
peak_Cue_idx_all = NaN(1, n_regions);
peak_Cue_perc_all = NaN(1, n_regions);
peak_perc_all = NaN(1, n_regions);
peak_idx_all = NaN(1, n_regions);

% Create figure
figure('units','normalized','outerposition',[0 0 0.6 0.5]);
hold on;
plot_handles = gobjects(n_regions, 1);

for n_region = [1,3,4,5,6] %1:n_regions % GB: [1,3,4,5,6]
    unit_region = brainAreas{n_region};
    disp(['Processing brain region: ' unit_region]);

    all_errTest_timebin = NaN(numSessions, num_timebins);  % Per region

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
            elseif strcmp('dlPFC', unit_region)
                SessionData = Go_data.dlPFC_Go(idxThisSession,:);
            else
                error([unit_region ' does not exist '])
         end
    
         % skip session days that are empty - relevant for S1 session 20230810
         if isempty(SessionData{1})
            continue
         end
    
         % Z-scoring => calc the mean of the FR of each unit across all trials per timebin, not per phase)
        num_trials = length(SessionData);
        [num_timebins, num_units] = size(SessionData{1});
        % reconfigure matrix to store all trials
        all_data = NaN(num_timebins,num_units,num_trials); % 174 x 63 x 142
        for t = 1:num_trials
            all_data(:,:,t) = SessionData{t}; % (timebins x units) for each trial
        end
    
        % Compute the mean and SD across trials
        mean_fr = mean(all_data, 3); % Result is (timebins x units)
        std_fr = std(all_data, 0, 3); % (timebins x units)
        std_fr(std_fr == 0) = 1; % Avoid division by zero by setting std_fr to 1 where it's zero
    
        % Z-score normalization: (X - mean) / std
        z_scored_fr = (all_data - mean_fr) ./ std_fr; % (timebins x units x trials)
    
        % Initialize the cell array
        z_scored_data = cell(num_trials, 1);
    
        % Fill the cell array with z-scored data
        for t = 1:num_trials
            z_scored_data{t} = z_scored_fr(:,:,t); % Extract each trial's matrix
        end
    
        %labels 
        sessionLabels = Go_data.GoLabels(idxThisSession,:);
        
        %trialType
        trialTypeSession = Go_data.TrialType(idxThisSession,:);
    
        % grasp labels
        graspTypeSession = Go_data.GraspType(idxThisSession,:);

        % ApertureSize
        %apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);

        % object labels
        objectTypeSession = Go_data.ObjectType(idxThisSession,:);
      
        %get idx for Go or NoGo trials
        GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
        time_phase_labels = Go_data.time_phase_labels(idxThisSession);
        
    
        if flagGoTrials
            %SessionData = SessionData(GoNoGoidx);
            SessionData = z_scored_data(GoNoGoidx);
            sessionLabels = sessionLabels(GoNoGoidx);
            time_phase_labels = time_phase_labels(GoNoGoidx);
            trialTypeSession = trialTypeSession(GoNoGoidx);
            graspTypeSession = graspTypeSession(GoNoGoidx);
            %apertureSizeSession = apertureSizeSession(GoNoGoidx);
            objectTypeSession = objectTypeSession(GoNoGoidx);
        else
            SessionData = SessionData(~GoNoGoidx);
            sessionLabels = sessionLabels(~GoNoGoidx);
            time_phase_labels = time_phase_labels(~GoNoGoidx);
            trialTypeSession = trialTypeSession(~GoNoGoidx);
        end
         
        % seperate data according to cue condition 
        unTrialType = unique(Go_data.TrialType);
    
        % separate data according to grasp
        unGraspType = unique(Go_data.GraspType);

        % separate data according to object 
        unObjectType = unique(Go_data.ObjectType); % includes "Associated" currently

        % % only keeping Small and Large
        % % Define which sizes to keep
        % sizesToKeep = {'Small', 'Large'};
        % 
        % % Find indices of trials belonging to 'small' or 'large' sessions
        % SLsizeIdx = ismember(apertureSizeSession, sizesToKeep);
        % 
        % % Extract trials that correspond to 'small' or 'large' sessions
        % SessionData = SessionData(SLsizeIdx);
        % sessionLabels = sessionLabels(SLsizeIdx);
        % time_phase_labels = time_phase_labels(SLsizeIdx);
        % trialTypeSession = trialTypeSession(SLsizeIdx);
        % graspTypeSession = graspTypeSession(SLsizeIdx);
        % apertureSizeSession = apertureSizeSession(SLsizeIdx);
        % %objectTypeSession = objectTypeSession(SLsizeIdx);
        % 
        % % seperate data according to size 
        % unAperture = unique(apertureSizeSession);
    
        % code for CMs of grasps and modalities separately
        sessionLabels_modality = trialTypeSession;
        sessionLabels_grasp = graspTypeSession;
        %sessionLabels_size = apertureSizeSession;
        sessionLabels_object = objectTypeSession;
    
        % Convert modality labels ('Hand', 'HandObject', 'Object') to numerical values
        modality_labels = {'Hand', 'Hand_Object', 'Object'}; %'Combined'
        grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
        %size_labels = {'Small', 'Medium', 'Large'};
        object_labels = {'deck','block','rod','ball'};
    
        sessionLabels_modality_num = zeros(size(sessionLabels_modality));  % Initialize numerical labels
        sessionLabels_grasp_num = zeros(size(sessionLabels_grasp));
        %sessionLabels_size_num = zeros(size(sessionLabels_size));
        sessionLabels_object_num = zeros(size(sessionLabels_object));
    
        % Loop through labels and assign numerical values
        for i = 1:length(modality_labels)
            sessionLabels_modality_num(strcmp(sessionLabels_modality, modality_labels{i})) = i;
        end  
        for i = 1:length(grasp_labels)
            sessionLabels_grasp_num(strcmp(sessionLabels_grasp, grasp_labels{i})) = i;
        end 
        % for i = 1:length(size_labels)
        %     sessionLabels_size_num(strcmp(sessionLabels_size, size_labels{i})) = i;
        % end
        for i = 1:length(object_labels)
            sessionLabels_object_num(strcmp(sessionLabels_object, object_labels{i})) = i;
        end 
    
        combinedTrialIdx = ismember(objectTypeSession, object_labels); % removes "Associated" object type trials

        errTest_timebin = NaN(num_timebins, 1);
    
        % loop through timebins 
        for n_bin = 1:num_timebins
       
            % Extract data for this timebin across all trials
            data_per_timebin = cell2mat(cellfun(@(x) x(n_bin, :), SessionData(combinedTrialIdx), 'UniformOutput', false)); %(combinedTrialIdx)
            
            % Run classification
            [~, errTestTmp, ~, ~, ~] = classification.LDA_classification_rep(data_per_timebin, sessionLabels_grasp_num(combinedTrialIdx), ...
                'flagErrorMatrix', false, 'PCA_variance', 95, 'flagLeaveOneOut', true, 'flagRandomPerm', false); %(combinedTrialIdx)
        
            % Store classification accuracy
            errTest_timebin(n_bin) = (1 - mean(errTestTmp)) * 100;
    
            set(gca, 'FontSize', 12);
        end    
        
        % Store the classification accuracy for each session
        all_errTest_timebin(n_session, :) = errTest_timebin;
    
    end
    % Store all sessions for current region
    all_errTest_timebin_all_regions(:, :, n_region) = all_errTest_timebin;

    % Compute stats
    CI95 = utile.calculate_CI(all_errTest_timebin); 
    mean_acc = mean(all_errTest_timebin, 1, 'omitnan');
    ci_upper = CI95(2, :);
    ci_lower = CI95(1, :);
    chance = 1 / numel(unGraspType) * 100;

    first_sig_idx = find(mean_acc(phase_changes(2):end) + ci_lower(phase_changes(2):end) > chance, 1, 'first') + 42;
    if isempty(first_sig_idx)
        first_sig_idx = 1;
    end
    first_sig_perc = mean_acc(first_sig_idx);

    [~, cue_peak_rel_idx] = max(mean_acc(phase_changes(2):phase_changes(3)));
    peak_Cue_idx = phase_changes(2) + cue_peak_rel_idx - 1;
    peak_Cue_perc = mean_acc(peak_Cue_idx);

    [peak_perc, peak_idx] = max(mean_acc);

    % Store region summary stats
    first_sig_idx_all(n_region) = first_sig_idx;
    first_sig_perc_all(n_region) = first_sig_perc;
    peak_Cue_idx_all(n_region) = peak_Cue_idx;
    peak_Cue_perc_all(n_region) = peak_Cue_perc;
    peak_perc_all(n_region) = peak_perc;
    peak_idx_all(n_region) = peak_idx;

    % Plotting
    ER = utile.shadedErrorBar(1:num_timebins, mean_acc, ci_upper, 'lineprops', '-b');
    ER.mainLine.Color = color_info{n_region};
    ER.mainLine.LineWidth = 2;
    ER.patch.FaceColor = color_info{n_region};
    ER.edge(1).LineStyle = 'none';
    ER.edge(2).LineStyle = 'none';
    plot_handles(n_region) = ER.mainLine;
    
end

% Finalize plot
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5, 'FontSize', 12);
end
xlim([0 179]);
xtickangle(45);
xlabel('Timebin');
ylabel('Classification Accuracy [%]');
title('Object Classification Accuracy Over Time');
yline(chance, '--r', 'LineWidth', 1.5);
ylim([0 100]);
yticks([0 25 50 75 100]);
legend(plot_handles, brainAreas{[1,3,4,5,6]}, 'Location', 'best');
set(gca, 'FontSize', 12);
hold off;

%% TEST analysis for timing of onset and peak
% Initialize storage
all_errTest_timebin_all_regions = NaN(numSessions, num_timebins, n_regions);
all_errTest_timebin_shuf_all_regions = NaN(numSessions, num_timebins, n_regions);
first_sig_timebin_Cue_all = NaN(1, n_regions);
first_sig_perc_all = NaN(1, n_regions);
peak_Cue_timebin_all = NaN(1, n_regions);
peak_Cue_perc_all = NaN(1, n_regions);
peak_Action_timebin_all = NaN(1, n_regions);
peak_Action_perc_all = NaN(1, n_regions);

for n_region = 1 %1:n_regions % GB: [1,3,4,5,6]
    unit_region = brainAreas{n_region};
    disp(['Processing brain region: ' unit_region]);

    all_errTest_timebin = NaN(numSessions, num_timebins);  % Per region
    all_errTest_timebin_shuf = NaN(numSessions, num_timebins);

    for n_session = 1%:numSessions
    
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
            elseif strcmp('dlPFC', unit_region)
                SessionData = Go_data.dlPFC_Go(idxThisSession,:);
            else
                error([unit_region ' does not exist '])
         end
    
         % skip session days that are empty - relevant for S1 session 20230810
         if isempty(SessionData{1})
            continue
         end
    
         % Z-scoring => calc the mean of the FR of each unit across all trials per timebin, not per phase)
        num_trials = length(SessionData);
        [num_timebins, num_units] = size(SessionData{1});
        % reconfigure matrix to store all trials
        all_data = NaN(num_timebins,num_units,num_trials); % 174 x 63 x 142
        for t = 1:num_trials
            all_data(:,:,t) = SessionData{t}; % (timebins x units) for each trial
        end
    
        % Compute the mean and SD across trials
        mean_fr = mean(all_data, 3); % Result is (timebins x units)
        std_fr = std(all_data, 0, 3); % (timebins x units)
        std_fr(std_fr == 0) = 1; % Avoid division by zero by setting std_fr to 1 where it's zero
    
        % Z-score normalization: (X - mean) / std
        z_scored_fr = (all_data - mean_fr) ./ std_fr; % (timebins x units x trials)
    
        % Initialize the cell array
        z_scored_data = cell(num_trials, 1);
    
        % Fill the cell array with z-scored data
        for t = 1:num_trials
            z_scored_data{t} = z_scored_fr(:,:,t); % Extract each trial's matrix
        end
    
        %labels 
        sessionLabels = Go_data.GoLabels(idxThisSession,:);
        
        %trialType
        trialTypeSession = Go_data.TrialType(idxThisSession,:);
    
        % grasp labels
        graspTypeSession = Go_data.GraspType(idxThisSession,:);
      
        %get idx for Go or NoGo trials
        GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
        time_phase_labels = Go_data.time_phase_labels(idxThisSession);
        
    
        if flagGoTrials
            %SessionData = SessionData(GoNoGoidx);
            SessionData = z_scored_data(GoNoGoidx);
            sessionLabels = sessionLabels(GoNoGoidx);
            time_phase_labels = time_phase_labels(GoNoGoidx);
            trialTypeSession = trialTypeSession(GoNoGoidx);
            graspTypeSession = graspTypeSession(GoNoGoidx);
            %apertureSizeSession = apertureSizeSession(GoNoGoidx);
            %objectTypeSession = objectTypeSession(GoNoGoidx);
        else
            SessionData = SessionData(~GoNoGoidx);
            sessionLabels = sessionLabels(~GoNoGoidx);
            time_phase_labels = time_phase_labels(~GoNoGoidx);
            trialTypeSession = trialTypeSession(~GoNoGoidx);
        end
         
        % seperate data according to cue condition 
        unTrialType = unique(Go_data.TrialType);
    
        % separate data according to grasp
        unGraspType = unique(Go_data.GraspType);
    
        % code for CMs of grasps and modalities separately
        sessionLabels_modality = trialTypeSession;
        sessionLabels_grasp = graspTypeSession;
    
        % Convert modality labels ('Hand', 'HandObject', 'Object') to numerical values
        modality_labels = {'Hand', 'Hand_Object', 'Object'}; %'Combined'
        grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
    
        sessionLabels_modality_num = zeros(size(sessionLabels_modality));  % Initialize numerical labels
        sessionLabels_grasp_num = zeros(size(sessionLabels_grasp));
    
        % Loop through labels and assign numerical values
        for i = 1:length(modality_labels)
            sessionLabels_modality_num(strcmp(sessionLabels_modality, modality_labels{i})) = i;
        end  
        for i = 1:length(grasp_labels)
            sessionLabels_grasp_num(strcmp(sessionLabels_grasp, grasp_labels{i})) = i;
        end 

        errTest_timebin = NaN(num_timebins, 1);
    
        % loop through Cue 
        for n_bin = 42:82 % 42:82 % cue 94:174 % action 
       
            % Extract data for this timebin across all trials
            data_per_cue_timebins = cell2mat(cellfun(@(x) x(n_bin, :), SessionData, 'UniformOutput', false)); 
            
            % Run classification
            [~, errTestTmp, ~, ~, ~] = classification.LDA_classification_rep(data_per_cue_timebins, sessionLabels_grasp_num, ...
                'flagErrorMatrix', false, 'PCA_variance', 95, 'flagLeaveOneOut', true, 'flagRandomPerm', false); % regular classification
        
            % Store classification accuracy
            errTest_timebin(n_bin) = (1 - mean(errTestTmp)) * 100;

            % Run classification on shuffled data
            [~, errTestTmp_shuf, ~, ~, ~] = classification.LDA_classification_rep(data_per_cue_timebins, sessionLabels_grasp_num, ...
                'flagErrorMatrix', false, 'PCA_variance', 95, 'flagLeaveOneOut', true, 'flagRandomPerm', true); % shuffled classification
        
            % Store shuffled classification accuracy
            errTest_timebin_shuf(n_bin) = (1 - mean(errTestTmp_shuf)) * 100;
    
        end    
        
        % Store the classification accuracy for each session
        all_errTest_timebin(n_session, :) = errTest_timebin;

        % Store the shuffled classification accuracy for each session
        all_errTest_timebin_shuf(n_session, :) = errTest_timebin_shuf;
    
    end
    % Store all sessions for current region
    all_errTest_timebin_all_regions(:, :, n_region) = all_errTest_timebin;

    % Store all shuffled sessions for current region
    all_errTest_timebin_shuf_all_regions(:, :, n_region) = all_errTest_timebin_shuf;

    % Compute stats
    CI95 = utile.calculate_CI(all_errTest_timebin); 
    mean_acc = mean(all_errTest_timebin, 1, 'omitnan');
    ci_upper = CI95(2, :);
    ci_lower = CI95(1, :);

    % Compute shuffled stats
    shuf_CI95 = utile.calculate_CI(all_errTest_timebin_shuf); 
    shuf_mean_acc = mean(all_errTest_timebin_shuf, 1, 'omitnan');
    shuf_ci_upper = shuf_CI95(2, :);
    shuf_ci_lower = shuf_CI95(1, :);
    chance = shuf_mean_acc(phase_changes(2):phase_changes(3)) + shuf_ci_upper(phase_changes(2):phase_changes(3));

    first_sig_timebin_Cue = find(mean_acc(phase_changes(2):phase_changes(3)) + ci_lower(phase_changes(2):phase_changes(3)) > chance, 1, 'first') + 42;
    if isempty(first_sig_timebin_Cue)
        first_sig_timebin_Cue = 1;
    end
    first_sig_perc = mean_acc(first_sig_timebin_Cue);

    [~, cue_peak_rel_idx] = max(mean_acc(phase_changes(2):phase_changes(3)));
    peak_Cue_timebin = phase_changes(2) + cue_peak_rel_idx - 1;
    peak_Cue_perc = mean_acc(peak_Cue_timebin);

    % Store region summary stats
    first_sig_timebin_Cue_all(n_region) = first_sig_timebin_Cue;
    first_sig_perc_all(n_region) = first_sig_perc;
    peak_Cue_timebin_all(n_region) = peak_Cue_timebin;
    peak_Cue_perc_all(n_region) = peak_Cue_perc;
    
end

%% Save everything to one file
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
%filename = "decoded_objects_per_timebin_" + taskName + "_ALL_REGIONS_LDA_" + goLabel + "_z_scored.mat";
filename = "decoded_grasps_per_timebin_" + taskName + "_ALL_REGIONS_Shuffled_LDA_" + goLabel + "_z_scored.mat"; % shuffled data
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
save(fullfile(directory, filename), ...
    'all_errTest_timebin_all_regions', ...
    'first_sig_idx_all', 'first_sig_perc_all', ...
    'peak_Cue_idx_all', 'peak_Cue_perc_all', ...
    'peak_idx_all','peak_perc_all',...
    'brainAreas');

keyboard
%% analysis and plot for individual regions 
%Go_data = Data.Go_data;

% remove faulty sessions, if any
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20250618'}; % remove this date bc not enough trials per grasp for object decoding per grasp
elseif strcmp(subject_id, 's3')
    error_session = {'20250212'};
elseif strcmp(subject_id, 's4')
    error_session = {'20240613'};
end

if ~isempty(error_session)
    condition = cellfun(@(x) strcmp(x, error_session), Go_data.session_date);
    Go_data = Go_data(~condition,:);
end

unit_region = 'dlPFC';
%brainAreas = Go_data.frPerChannel{6};
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phaseTimeTmp = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
%color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Cue Conditions
%color_info = {[0.2, 0.13, 0.53], [0.067, 0.467, 0.2], [0.53, 0.8, 0.93], [0.53, 0.13, 0.33]}; % grasps: Purple, Green, Light Blue, Dark Pink
color_info = {[.3906 .5586 .9961],[.4688 .3672 .9375],[.8594 .1484 .4961],[.9922 .3789 0]}; % objects

sessions_all = unique(Go_data.session_date);
numSessions = numel(sessions_all);
uniqueCueTypes = {'Hand','Hand-Object','Object'};
if strcmp(taskName, 'GraspObject_Combined')
    uniqueCueTypes = {'Combined','Hand','Hand-Object','Object'};
    %color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
end
if strcmp(taskName, 'GraspObject_Varied_Size')
    % add Aperture Size column
    sizeKeywords = ['Small', 'Medium', 'Large'];
    Go_data.Aperture_Size = cell(height(Go_data),1);
    % Loop through each label and extract the size information
    for i = 1:height(Go_data)
        % Use regular expression to find the size keyword after the last underscore
        tokens = regexp(Go_data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');
        
        if ~isempty(tokens)
            % tokens is a cell array; extract the size keyword from it
            Go_data.Aperture_Size{i} = tokens{1}{1};
        end
    end
end

flagGoTrials = true; 

taskSizesAll = {'Small', 'Medium', 'Large'};
graspTypesAll = {'PalmarPinch', 'MediumWrap', 'Sphere3Finger', 'Lateral'};

% Initialize a matrix to store classification accuracy for all sessions
%all_errTest_timebin = NaN(numSessions, 174);
%all_errTest_timebin = NaN(numSessions, numel(uniqueCueTypes), 174); % for sep by Cue Modality
all_errTest_timebin = NaN(numSessions, numel(graspTypesAll), 174); % for sep by Grasp

% Initialize a cell array to store confusion matrices for each phase
%confMatAllSessions = cell(numSessions,numPhases,1); 
%figure(); % for CM

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
        elseif strcmp('dlPFC', unit_region)
            SessionData = Go_data.dlPFC_Go(idxThisSession,:);
        else
            error([unit_region ' does not exist '])
     end

     % skip session days that are empty - relevant for S1 session 20230810
     if isempty(SessionData{1})
        continue
     end

    % Z-scoring => calc the mean of the FR of each unit across all trials per timebin, not per phase)
    num_trials = length(SessionData);
    [num_timebins, num_units] = size(SessionData{1});
    % reconfigure matrix to store all trials
    all_data = NaN(num_timebins,num_units,num_trials); % 174 x 63 x 142
    for t = 1:num_trials
        all_data(:,:,t) = SessionData{t}; % (timebins x units) for each trial
    end

    % Compute the mean and SD across trials
    mean_fr = mean(all_data, 3); % Result is (timebins x units)
    std_fr = std(all_data, 0, 3); % (timebins x units)
    std_fr(std_fr == 0) = 1; % Avoid division by zero by setting std_fr to 1 where it's zero

    % Z-score normalization: (X - mean) / std
    z_scored_fr = (all_data - mean_fr) ./ std_fr; % (timebins x units x trials)

    % Initialize the cell array
    z_scored_data = cell(num_trials, 1);

    % Fill the cell array with z-scored data
    for t = 1:num_trials
        z_scored_data{t} = z_scored_fr(:,:,t); % Extract each trial's matrix
    end

    %labels 
    sessionLabels = Go_data.GoLabels(idxThisSession,:);
    
    %trialType
    trialTypeSession = Go_data.TrialType(idxThisSession,:);

    % grasp labels
    graspTypeSession = Go_data.GraspType(idxThisSession,:);

    %AperatureSize
    %apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);

    % object labels
    objectTypeSession = Go_data.ObjectType(idxThisSession,:);
  
    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    time_phase_labels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        %SessionData = SessionData(GoNoGoidx);
        SessionData = z_scored_data(GoNoGoidx);
        sessionLabels = sessionLabels(GoNoGoidx);
        time_phase_labels = time_phase_labels(GoNoGoidx);
        trialTypeSession = trialTypeSession(GoNoGoidx);
        graspTypeSession = graspTypeSession(GoNoGoidx);
        %apertureSizeSession = apertureSizeSession(GoNoGoidx);
        objectTypeSession = objectTypeSession(GoNoGoidx);
    else
        SessionData = SessionData(~GoNoGoidx);
        sessionLabels = sessionLabels(~GoNoGoidx);
        time_phase_labels = time_phase_labels(~GoNoGoidx);
        trialTypeSession = trialTypeSession(~GoNoGoidx);
    end
     
    % seperate data according to cue condition 
    unTrialType = unique(Go_data.TrialType);

    % separate data according to grasp
    unGraspType = unique(Go_data.GraspType);

    % separate data according to object
    unObjectType = unique(Go_data.ObjectType);

    % seperate data according to size 
    %unAperture = unique(Go_data.Aperture_Size);

    % % only keeping Small and Large
    % % Define which sizes to keep
    % sizesToKeep = {'Small', 'Large'};
    % 
    % % Find indices of trials belonging to 'small' or 'large' sessions
    % SLsizeIdx = ismember(apertureSizeSession, sizesToKeep);
    % 
    % % Extract trials that correspond to 'small' or 'large' sessions
    % SessionData = SessionData(SLsizeIdx);
    % sessionLabels = sessionLabels(SLsizeIdx);
    % time_phase_labels = time_phase_labels(SLsizeIdx);
    % trialTypeSession = trialTypeSession(SLsizeIdx);
    % graspTypeSession = graspTypeSession(SLsizeIdx);
    % apertureSizeSession = apertureSizeSession(SLsizeIdx);
    % 
    % % seperate data according to size 
    % unAperture = unique(apertureSizeSession);

    % code for CMs of grasps and modalities separately
    sessionLabels_modality = trialTypeSession;
    sessionLabels_grasp = graspTypeSession;
    %sessionLabels_size = apertureSizeSession;
    sessionLabels_object = objectTypeSession;

    % Convert labels to numerical values
    modality_labels = {'Hand', 'Hand_Object', 'Object'};
    grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
    %size_labels = {'Small', 'Large'};
    object_labels = {'ball','block','deck','rod'};

    sessionLabels_modality_num = zeros(size(sessionLabels_modality));  % Initialize numerical labels
    sessionLabels_grasp_num = zeros(size(sessionLabels_grasp));
    %sessionLabels_size_num = zeros(size(sessionLabels_size));
    sessionLabels_object_num = zeros(size(sessionLabels_object));

    % Loop through labels and assign numerical values
    for i = 1:length(modality_labels)
        sessionLabels_modality_num(strcmp(sessionLabels_modality, modality_labels{i})) = i;
    end  
    for i = 1:length(grasp_labels)
        sessionLabels_grasp_num(strcmp(sessionLabels_grasp, grasp_labels{i})) = i;
    end 
    % for i = 1:length(size_labels)
    %     sessionLabels_size_num(strcmp(sessionLabels_size, size_labels{i})) = i;
    % end
    for i = 1:length(object_labels)
        sessionLabels_object_num(strcmp(sessionLabels_object, object_labels{i})) = i;
    end 

    combinedTrialIdx = ismember(objectTypeSession, object_labels); % removes "Associated" object type trials
    SessionData = SessionData(combinedTrialIdx); % get data for combined trials only
    sessionLabels_object_num = sessionLabels_object_num(combinedTrialIdx);
    sessionLabels_grasp_num = sessionLabels_grasp_num(combinedTrialIdx);

    errTest_timebin = NaN(num_timebins, 1);
    
    for n_type = 1:numel(object_labels) 

        % find idx of trial type 
        %trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
        %trialGraspIdx = ismember(graspTypeSession(combinedTrialIdx), unGraspType(n_type));
        trialObjectIdx = ismember(objectTypeSession(combinedTrialIdx), object_labels(n_type));

        % loop through timebins 
        for n_bin = 1:num_timebins
       
            % Extract data for this timebin across all trials
            %data_per_timebin_per_grasp = cell2mat(cellfun(@(x) x(n_bin, :), SessionData(trialGraspIdx), 'UniformOutput', false));
            %data_per_timebin_per_cue = cell2mat(cellfun(@(x) x(n_bin, :), SessionData(trialTypeIdx), 'UniformOutput', false));
            data_per_timebin_per_object = cell2mat(cellfun(@(x) x(n_bin, :), SessionData(trialObjectIdx), 'UniformOutput', false));

            % Run classification
            [~, errTestTmp, ~, ~, ~] = classification.LDA_classification_rep(data_per_timebin_per_object, sessionLabels_grasp_num(trialObjectIdx), ...
                'flagErrorMatrix', false, 'PCA_variance', 95, 'flagLeaveOneOut', true);
        
            % Store classification accuracy
            errTest_timebin(n_bin) = (1 - mean(errTestTmp)) * 100;
    
            set(gca, 'FontSize', 12);
        end    
    
        % Store the classification accuracy for each session
        all_errTest_timebin(n_session, n_type, :) = errTest_timebin;
    end
end

% === Plot classification accuracy over time for each Trial/Grasp Type ===
figure('units','normalized','outerposition',[0 0 0.28 0.23]); % laptop: [0 0 0.46 0.37]
hold on;

% first_sig_idx = zeros(1, numel(unTrialType));
% first_sig_perc = zeros(1, numel(unTrialType));
% peak_Cue_idx = zeros(1, numel(unTrialType));
% peak_Cue_perc = zeros(1, numel(unTrialType));
% peak_perc_all = NaN(1, numel(unTrialType));
% peak_idx_all = NaN(1, numel(unTrialType));

first_sig_idx = zeros(1, numel(unGraspType));
first_sig_perc = zeros(1, numel(unGraspType));
peak_Cue_idx = zeros(1, numel(unGraspType));
peak_Cue_perc = zeros(1, numel(unGraspType));
peak_perc_all = NaN(1, numel(unGraspType));
peak_idx_all = NaN(1, numel(unGraspType));

CI95 = utile.calculate_CI(all_errTest_timebin); % (2 x numTypes x timebins)
mean_acc = squeeze(mean(all_errTest_timebin, 1, 'omitnan')); % (numTypes x timebins)
ci_upper = squeeze(CI95(2,:,:)); % (numTypes x timebins)
ci_lower = squeeze(CI95(1,:,:)); % (numTypes x timebins)

chance = 1 / numel(unGraspType) * 100;

legend_handles = gobjects(numel(unGraspType), 1); % store line handles
legend_labels = cell(numel(unGraspType), 1);      % store labels

for n_type = 1:numel(unGraspType)
    % Find first significant bin (Cue phase or later)
    sig_search_idx = phase_changes(2):num_timebins;
    sig_idx = find(mean_acc(n_type, sig_search_idx) + ci_lower(n_type, sig_search_idx) > chance, 1, 'first');
    if ~isempty(sig_idx)
        first_sig_idx(n_type) = sig_search_idx(sig_idx);
        first_sig_perc(n_type) = mean_acc(n_type, first_sig_idx(n_type));
    end

    % Find peak in Cue phase
    cue_range = phase_changes(2):phase_changes(3)-1;
    [~, rel_peak_idx] = max(mean_acc(n_type, cue_range));
    peak_Cue_idx(n_type) = cue_range(1) + rel_peak_idx - 1;
    peak_Cue_perc(n_type) = mean_acc(n_type, peak_Cue_idx(n_type));

    [peak_perc, peak_idx] = max(mean_acc(n_type,:));

    peak_perc_all(n_type) = peak_perc;
    peak_idx_all(n_type) = peak_idx;

    % Plot shaded error
    ER = utile.shadedErrorBar(1:num_timebins, mean_acc(n_type,:), ...
                              ci_upper(n_type,:), ...
                              'lineprops', '-','transparent', true);
    ER.mainLine.Color = color_info{n_type};
    ER.mainLine.LineWidth = 2;
    ER.patch.FaceColor = color_info{n_type};
    ER.edge(1).LineStyle = 'none'; 
    ER.edge(2).LineStyle = 'none';

    % Store for legend
    legend_handles(n_type) = ER.mainLine;
    %legend_labels{n_type} = uniqueCueTypes{n_type};
    %legend_labels{n_type} = unGraspType{n_type};
    legend_labels{n_type} = object_labels{n_type};

    % % Plot first significant timepoint
    % plot(first_sig_idx(n_type), first_sig_perc(n_type), 'v', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{n_type}, 'MarkerEdgeColor', 'k');
    % 
    % % Plot peak timepoint
    % plot(peak_Cue_idx(n_type), peak_Cue_perc(n_type), 'o', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{n_type}, 'MarkerEdgeColor', 'k');
end

% Add phase lines
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5,'FontSize',12);
end

xlim([0 num_timebins + 1]);
ylim([0 100]);
yline(chance,'--r','LineWidth',1.5);
yticks([0 25 50 75 100]);
xlabel('Timebin');
ylabel('Classification Accuracy [%]');
title(['Grasp Classification Over Time - ' unit_region]);
legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;


% save variables

goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);

% Create the filename using the brain region and analysis type
filename = "decoded_grasps_per_timebin_per_object_" + taskName + '_' + unit_region + "_LDA_" + goLabel + "_z_scored.mat"; % when z-scoring

directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
full_path = fullfile(directory, filename);

% Save the relevant variables with the dynamic filename
save(full_path, 'all_errTest_timebin','first_sig_idx','first_sig_perc','peak_Cue_idx','peak_Cue_perc','peak_idx_all','peak_perc_all');

%% load Data
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
%filename = "decoded_grasp_per_timebin_per_cue_" + taskName + '_' + unit_region + "_LDA_" + goLabel + "_z_scored.mat"; % when z-scoring
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
%filename = "decoded_objects_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

%% load multiple so can plot variable decoding together
% load Data GRASP
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% GRASP %
G_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
G_first_sig_idx_all = first_sig_idx_all;
G_first_sig_perc_all = first_sig_perc_all;
G_peak_Cue_idx_all = peak_Cue_idx_all;
G_peak_Cue_perc_all = peak_Cue_perc_all;
%G_peak_perc_all = peak_perc_all; % not in original pipeline
%G_peak_idx_all = peak_idx_all;

%% load multiple so can plot variable decoding together
% load Data OBJECT/SIZE
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_sizeExtremes_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; % SIZE
%filename = "decoded_objects_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; % OBJECT
full_path = fullfile(directory, filename);
load(full_path);

% OBJECT/SIZE EXTREMES %
X_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
X_first_sig_idx_all = first_sig_idx_all;
X_first_sig_perc_all = first_sig_perc_all;
X_peak_Cue_idx_all = peak_Cue_idx_all;
X_peak_Cue_perc_all = peak_Cue_perc_all;
%X_peak_perc_all = peak_perc_all; % not in original pipeline
%X_peak_idx_all = peak_idx_all;

%% load SHUFFLED to replace chance line
% load Data SHUFFLED GRASP
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% SHUFFLED GRASP %
SG_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
SG_first_sig_idx_all = first_sig_idx_all;
SG_first_sig_perc_all = first_sig_perc_all;
SG_peak_Cue_idx_all = peak_Cue_idx_all;
SG_peak_Cue_perc_all = peak_Cue_perc_all;
SG_peak_perc_all = peak_perc_all;
SG_peak_idx_all = peak_idx_all;

%% load SHUFFLED to replace chance line
% load Data SHUFFLED OBJECTS/SIZE
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
%filename = "decoded_objects_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
filename = "decoded_sizeExtremes_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% SHUFFLED OBJECT/SIZE %
SX_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
SX_first_sig_idx_all = first_sig_idx_all;
SX_first_sig_perc_all = first_sig_perc_all;
SX_peak_Cue_idx_all = peak_Cue_idx_all;
SX_peak_Cue_perc_all = peak_Cue_perc_all;
SX_peak_perc_all = peak_perc_all;
SX_peak_idx_all = peak_idx_all;

%% plot grasp and X classification together
% Determine number of items to be plotted together (decoded items, ie. grasp, object, size, shuffled)
n_decoded = 4;
%chance = 1 / 4 * 100; 
numPhases = numel(phase_changes);

% Create figure
figure('units','normalized','outerposition',[0 0 0.08 0.163]); %.137 or .08
hold on;
plot_handles = gobjects(n_decoded, 1);
%color_info = {[0.3359, 0.7031, 0.9101],[0.8984 0.6211 0],[0.8320 0.3672 0],[0.7969, 0.4726, 0.6523],[0, 0.6171, 0.4492],[.9961 .6875 0]}; % SMG, PMV, S1, AIP, M1, dlPFC
% colors for grasp only
%color_info = {[0, 0.3529, 0.7098],[0 0 0],[0.75 0.75 0.75]};
% colors for object
%color_info = {[0, 0.3529, 0.7098],[0.9608, 0.7333, 0.1176],[0.75 0.75 0.75]}; % blue for grasp, yellow for object (pull from FR), gray for shuffled
% colors for size
color_info = {[0, 0.3529, 0.7098],[0.2039, 0.6902, 0.2902],[0.75 0.75 0.75]}; % blue for grasp, green for size (pull from FR), gray for shuffled

for n_region = 4 %1:n_regions
    % PLOT GRASP
    G_mean_acc = mean(G_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    G_CI95 = utile.calculate_CI(G_all_errTest_timebin_all_regions(:,:,n_region));
    G_ci_upper = G_CI95(2, :);
    G_ci_lower = G_CI95(1, :);

    G_ER = utile.shadedErrorBar(1:num_timebins, G_mean_acc, G_ci_upper, 'lineprops', '-b');
    G_ER.mainLine.Color = color_info{1};
    G_ER.mainLine.LineWidth = 2;
    G_ER.patch.FaceColor = color_info{1};
    G_ER.edge(1).LineStyle = 'none';
    G_ER.edge(2).LineStyle = 'none';
    plot_handles(1) = G_ER.mainLine;

    hold on;

    % PLOT SHUFFLED GRASP
    SG_mean_acc = mean(SG_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    SG_CI95 = utile.calculate_CI(SG_all_errTest_timebin_all_regions(:,:,n_region));
    SG_ci_upper = SG_CI95(2, :);
    SG_ci_lower = SG_CI95(1, :);

    SG_ER = utile.shadedErrorBar(1:num_timebins, SG_mean_acc, SG_ci_upper, 'lineprops', '-.');
    SG_ER.mainLine.Color = color_info{3};
    SG_ER.mainLine.LineWidth = 2;
    SG_ER.patch.FaceColor = color_info{3};
    SG_ER.edge(1).LineStyle = 'none';
    SG_ER.edge(2).LineStyle = 'none';
    plot_handles(2) = SG_ER.mainLine;

    hold on;

    % % PLOT X
    % X_mean_acc = mean(X_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    % X_CI95 = utile.calculate_CI(X_all_errTest_timebin_all_regions(:,:,n_region));
    % X_ci_upper = X_CI95(2, :);
    % X_ci_lower = X_CI95(1, :);
    % 
    % X_ER = utile.shadedErrorBar(1:num_timebins, X_mean_acc, X_ci_upper, 'lineprops', '-b');
    % X_ER.mainLine.Color = color_info{2};
    % X_ER.mainLine.LineWidth = 2;
    % X_ER.patch.FaceColor = color_info{2};
    % X_ER.edge(1).LineStyle = 'none';
    % X_ER.edge(2).LineStyle = 'none';
    % plot_handles(3) = X_ER.mainLine;
    % 
    % hold on;
    % 
    % % PLOT SHUFFLED X
    % SX_mean_acc = mean(SX_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    % SX_CI95 = utile.calculate_CI(SX_all_errTest_timebin_all_regions(:,:,n_region));
    % SX_ci_upper = SX_CI95(2, :);
    % SX_ci_lower = SX_CI95(1, :);
    % 
    % SX_ER = utile.shadedErrorBar(1:num_timebins, SX_mean_acc, SX_ci_upper, 'lineprops', '--');
    % SX_ER.mainLine.Color = color_info{3};
    % SX_ER.mainLine.LineWidth = 2;
    % SX_ER.patch.FaceColor = color_info{3};
    % SX_ER.edge(1).LineStyle = 'none';
    % SX_ER.edge(2).LineStyle = 'none';
    % plot_handles(4) = SX_ER.mainLine;


    % % Plot first significant timepoint for GRASP
    % plot(G_first_sig_idx_all(n_region), G_first_sig_perc_all(n_region), 'v', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{1}, 'MarkerEdgeColor', 'k');
    % 
    % % Plot peak timepoint for GRASP
    % plot(G_peak_Cue_idx_all(n_region), G_peak_Cue_perc_all(n_region), 'o', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{1}, 'MarkerEdgeColor', 'k');
    % 
    % % Plot first significant timepoint for X
    % plot(X_first_sig_idx_all(n_region), X_first_sig_perc_all(n_region), 'v', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{2}, 'MarkerEdgeColor', 'k');
    % 
    % % Plot peak timepoint for X
    % plot(X_peak_Cue_idx_all(n_region), X_peak_Cue_perc_all(n_region), 'o', ...
    %     'MarkerSize', 6, 'MarkerFaceColor', color_info{2}, 'MarkerEdgeColor', 'k');
end

% Finalize plot
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', 'LineWidth', 1.5, 'FontSize', 12); %phaseNames{n_phase}
end
%xlim([0 num_timebins]);
xtickangle();
%xlabel('Time (s)');
xlim([30 134])
xticks([1 42 83 124]); %165
xticklabels([0]); %8
%ylabel('Classification Accuracy [%]');
%title([brainAreas(n_region)]);
%yline(chance, '--r', 'LineWidth', 1.5);
ylim([15 65]); %40 80 ; 15 80
yticks([25 65]); %50 80  ; 25 80
%legend(plot_handles, [{'Grasp'},{'Shuffled Grasp'},{'Object'},{'Shuffled Object'}], 'Location', 'best'); %[{'Grasp'},{'Shuffled Grasp'},{'Object'},{'Shuffled Object'}]
set(gca, 'FontSize', 10);
hold off;

%% plot all regions once loaded

% Determine number of regions
n_regions = 3;
chance = 1 / 4 * 100; 
numPhases = numel(phase_changes);

% Create figure
figure('units','normalized','outerposition',[0 0 0.15 0.18]);
hold on;
plot_handles = gobjects(n_regions, 1);
color_info = {[0.3359, 0.7031, 0.9101],[0.8984 0.6211 0],[0.8320 0.3672 0],[0.7969, 0.4726, 0.6523],[0, 0.6171, 0.4492],[.9961 .6875 0]}; % SMG, PMV, S1, AIP, M1, dlPFC

for n_region = [3,5]%:n_regions % GB: [1,3,4,5] %1:n_regions
    mean_acc = mean(all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    CI95 = utile.calculate_CI(all_errTest_timebin_all_regions(:,:,n_region));
    ci_upper = CI95(2, :);
    ci_lower = CI95(1, :);

    ER = utile.shadedErrorBar(1:num_timebins, mean_acc, ci_upper, 'lineprops', '-b');
    ER.mainLine.Color = color_info{n_region};
    ER.mainLine.LineWidth = 2;
    ER.patch.FaceColor = color_info{n_region};
    ER.edge(1).LineStyle = 'none';
    ER.edge(2).LineStyle = 'none';
    plot_handles(n_region) = ER.mainLine;

    % Plot first significant timepoint
    plot(first_sig_idx_all(n_region), first_sig_perc_all(n_region), 'v', ...
        'MarkerSize', 6, 'MarkerFaceColor', color_info{n_region}, 'MarkerEdgeColor', 'k');

    % Plot peak timepoint
    plot(peak_Cue_idx_all(n_region), peak_Cue_perc_all(n_region), 'o', ...
        'MarkerSize', 6, 'MarkerFaceColor', color_info{n_region}, 'MarkerEdgeColor', 'k');
end

% Finalize plot
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5, 'FontSize', 12);
end
xlim([0 num_timebins]);
xtickangle(45);
xticklabels([]);
%xlabel('Timebin');
%ylabel('Classification Accuracy [%]');
%title('Grasp Classification Accuracy Over Time');
yline(chance, '--r', 'LineWidth', 1.5);
ylim([0 100]);
yticks([0 50 100]);
legend(plot_handles([3,5]), brainAreas([3,5]), 'Location', 'best'); % GB: ([1,3:5])
set(gca, 'FontSize', 12);
hold off;

%% plot single brain area

% figure('units','normalized','outerposition',[0 0 0.35 0.3]); % for desktop [0 0 0.35 0.3], for laptop [0 0 0.4 0.4]
% CI95 = utile.calculate_CI(all_errTest_timebin); 
% mean_acc = mean(all_errTest_timebin, 1, 'omitnan');       % (1 x numTimebins)
% ci_upper = CI95(2, :); % upper bound of the 95% CI
% ci_lower = CI95(1, :); % lower bound of the 95% CI
% chance = 1/(numel(unGraspType))*100;
% 
% first_sig_idx = find(mean_acc(phase_changes(2):end) + ci_lower(phase_changes(2):end) > chance, 1, 'first'); % limited to Cue phase and beyond, accounts for lower CI bound for significance
% first_sig_idx = first_sig_idx + 42; % to account for limiting to Cue phase
% first_sig_perc = mean_acc(first_sig_idx);
% 
% [~, cue_peak_rel_idx] = max(mean_acc(phase_changes(2):phase_changes(3))); % relative index
% peak_Cue_idx = phase_changes(2) + cue_peak_rel_idx - 1; % convert to absolute index in mean_acc
% peak_Cue_perc = mean_acc(peak_Cue_idx);
% 
% ER = utile.shadedErrorBar(1:num_timebins, mean_acc, ci_upper, 'lineprops', '-b');
% ER.mainLine.Color = color_info{1};
% ER.patch.FaceColor = color_info{1};
% ER.edge(1).LineStyle = 'none'; 
% ER.edge(2).LineStyle = 'none';
% hold on;
% % Plot first significant timepoint
% h1 = plot(first_sig_idx, first_sig_perc, 'v', ...
%     'MarkerSize', 5, 'MarkerFaceColor', color_info{1}, 'MarkerEdgeColor', 'k');
% 
% % Plot peak in Cue phase
% h2 = plot(peak_Cue_idx, peak_Cue_perc, 'o', ...
%     'MarkerSize', 5, 'MarkerFaceColor', color_info{1}, 'MarkerEdgeColor', 'k');
% 
% % Add text labels
% text(first_sig_idx + 2, first_sig_perc, num2str(first_sig_idx), ...
%     'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
% text(peak_Cue_idx + 2, peak_Cue_perc, num2str(peak_Cue_idx), ...
%     'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
% 
% for n_phase = 1:numPhases
%     xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5,'FontSize',12);
% end
% xlim([0 179]);
% xtickangle(45);
% xlabel('Timebin');
% ylabel('Classification Accuracy [%]');
% title(['Grasp Classification Accuracy Over Time - ' unit_region]); 
% yline(chance,'--r','LineWidth',1.5);
% ylim([0 100]);
% yticks([0 25 50 75 100]);
% legend([ER.mainLine, h1, h2], {'Mean Accuracy', 'First Sig. Timebin', 'Peak Timebin'}, ...
%        'Location', 'southeast');
% hold off;
% set(gca, 'FontSize', 12);

%% onset timing per region
% Initialize storage
first_sig_timebin_Cue_all = NaN(numSessions, n_regions); % numSessions
first_sig_perc_Cue_all = NaN(numSessions, n_regions);
peak_Cue_timebin_all = NaN(numSessions, n_regions);
peak_Cue_perc_all = NaN(numSessions, n_regions);
first_sig_timebin_Action_all = NaN(numSessions, n_regions); % numSessions
first_sig_perc_Action_all = NaN(numSessions, n_regions);
peak_Action_timebin_all = NaN(numSessions, n_regions);
peak_Action_perc_all = NaN(numSessions, n_regions);

%loop through each region
for n_region = 1:numel(brainAreas)
    % Compute stats
    % % GRASP % %
    regionData = G_all_errTest_timebin_all_regions(:,:,n_region);
    shuf_regionData = SG_all_errTest_timebin_all_regions(:,:,n_region);
    % % SIZE/SHAPE % %
    % regionData = X_all_errTest_timebin_all_regions(:,:,n_region);
    % shuf_regionData = SX_all_errTest_timebin_all_regions(:,:,n_region);

    first_sig_timebin_Cue_sessions = NaN(1, numSessions); %numSessions
    first_sig_perc_Cue_session = NaN(1, numSessions);
    peak_Cue_timebin_sessions = NaN(1, numSessions);
    peak_Cue_perc_sessions = NaN(1, numSessions);
    first_sig_timebin_Action_sessions = NaN(1, numSessions); %numSessions
    first_sig_perc_Action_session = NaN(1, numSessions);
    peak_Action_timebin_sessions = NaN(1, numSessions);
    peak_Action_perc_sessions = NaN(1, numSessions);

    for n_session = 1:numSessions
        CI95 = utile.calculate_CI(regionData); 
        ci_upper = CI95(2, :);
        ci_lower = CI95(1, :);
        perSessionData = regionData(n_session,:);
        
        % Compute shuffled stats
        shuf_CI95 = utile.calculate_CI(shuf_regionData); 
        shuf_ci_upper = shuf_CI95(2, :);
        shuf_ci_lower = shuf_CI95(1, :);
        shuf_perSessionData = shuf_regionData(n_session,:);

        % % CUE % %

        chance_Cue = shuf_perSessionData(phase_changes(2):phase_changes(3)) + shuf_ci_upper(phase_changes(2):phase_changes(3));
        
        first_sig_timebin_Cue = find(perSessionData(phase_changes(2):phase_changes(3)) + ci_lower(phase_changes(2):phase_changes(3)) > chance_Cue, 1, 'first') + 42;
        if isempty(first_sig_timebin_Cue)
            first_sig_timebin_Cue = 1;
        end
        first_sig_perc_Cue = perSessionData(first_sig_timebin_Cue);
        
        [~, cue_peak_rel_idx] = max(perSessionData(phase_changes(2):phase_changes(3)));
        peak_Cue_timebin = phase_changes(2) + cue_peak_rel_idx - 1;
        peak_Cue_perc = perSessionData(peak_Cue_timebin);
        
        % Store region summary stats
        first_sig_timebin_Cue_sessions(n_session) = first_sig_timebin_Cue;
        first_sig_perc_Cue_session(n_session) = first_sig_perc_Cue;
        peak_Cue_timebin_sessions(n_session) = peak_Cue_timebin;
        peak_Cue_perc_sessions(n_session) = peak_Cue_perc;

        % % ACTION % %
        chance_Action = shuf_perSessionData(phase_changes(4):num_timebins) + shuf_ci_upper(phase_changes(4):num_timebins);
        
        first_sig_timebin_Action = find(perSessionData(phase_changes(4):num_timebins) + ci_lower(phase_changes(4):num_timebins) > chance_Action, 1, 'first') + 94;
        if isempty(first_sig_timebin_Action)
            first_sig_timebin_Action = 1;
        end
        first_sig_perc_Action = perSessionData(first_sig_timebin_Action);
        
        [~, action_peak_rel_idx] = max(perSessionData(phase_changes(4):num_timebins));
        peak_Action_timebin = phase_changes(4) + action_peak_rel_idx - 1;
        peak_Action_perc = perSessionData(peak_Action_timebin);
        
        % Store region summary stats
        first_sig_timebin_Action_sessions(n_session) = first_sig_timebin_Action;
        first_sig_perc_Action_session(n_session) = first_sig_perc_Action;
        peak_Action_timebin_sessions(n_session) = peak_Action_timebin;
        peak_Action_perc_sessions(n_session) = peak_Action_perc;


    end

    % Store region summary stats - CUE
    first_sig_timebin_Cue_all(:,n_region) = first_sig_timebin_Cue_sessions;
    first_sig_perc_Cue_all(:,n_region) = first_sig_perc_Cue_session;
    peak_Cue_timebin_all(:,n_region) = peak_Cue_timebin_sessions;
    peak_Cue_perc_all(:,n_region) = peak_Cue_perc_sessions;
    
    % Store region summary stats - ACTION
    first_sig_timebin_Action_all(:,n_region) = first_sig_timebin_Action_sessions;
    first_sig_perc_Action_all(:,n_region) = first_sig_perc_Action_session;
    peak_Action_timebin_all(:,n_region) = peak_Action_timebin_sessions;
    peak_Action_perc_all(:,n_region) = peak_Action_perc_sessions;

end

%% box plots
figure('units','normalized','outerposition',[0 0 0.5 0.75])
%color_info = [0.3359, 0.7031, 0.9101;0.8984 0.6211 0;0.8320 0.3672 0;0.7969, 0.4726, 0.6523;0, 0.6171, 0.4492;.9961 .6875 0]; % SMG, PMV, S1, AIP, M1, dlPFC

% % AN % %
% groupOrder = ["AIP","SMG","M1","PMV","S1"];
% colorOrder = [0.7969, 0.4726, 0.6523;0.3359, 0.7031, 0.9101;0, 0.6171, 0.4492;0.8984 0.6211 0;0.8320 0.3672 0]; % AIP,SMG,M1,PMV,S1
% % GB % %
groupOrder = ["AIP","SMG","M1","S1","dlPFC","PMV"];
colorOrder = [0.7969, 0.4726, 0.6523;0.3359, 0.7031, 0.9101;0, 0.6171, 0.4492;0.8320 0.3672 0; .9961 .6875 0; 0.8984 0.6211 0]; % AIP,SMG,M1,S1,dlPFC
% % FG % %
% groupOrder = ["SMG","PMV","S1"];
% colorOrder = [0.3359, 0.7031, 0.9101;0.8984 0.6211 0;0.8320 0.3672 0]; % SMG,PMV,S1
%sgtitle('Planning Stats');
[~, reorderedIdx] = ismember(brainAreas, groupOrder);
% Planning onset timing
subplot(2,2,1);
h = boxplot(first_sig_timebin_Cue_all, brainAreas,"Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = first_sig_timebin_Cue_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Planning Onset Response');
set(h,{'linew'},{1.5})
xlim([42 82]);
xticks([42 52 62 72 82]);
xticklabels([2 2.5 3 3.5 4]);
xlabel('Time (s)');
ylabel('Brain Region');

% Planning onset percentage
subplot(2,2,2);
h = boxplot(first_sig_perc_Cue_all, brainAreas,"Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = first_sig_perc_Cue_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Planning Onset Percentage');
set(h,{'linew'},{1.5})
xlim([0 100]);
xticks([0 20 40 60 80 100]);
xlabel('Percentage (%)');
ylabel('Brain Region');

% Planning peak timing
subplot(2,2,3);
h = boxplot(peak_Cue_timebin_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = peak_Cue_timebin_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Planning Peak Response');
set(h,{'linew'},{1.5})
xlim([42 82]);
xticks([42 52 62 72 82]);
xticklabels([2 2.5 3 3.5 4]);
xlabel('Time (s)');
ylabel('Brain Region');

% Planning peak percentage
subplot(2,2,4);
h = boxplot(peak_Cue_perc_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = peak_Cue_perc_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Planning Peak Percentage');
set(h,{'linew'},{1.5})
xlim([0 100]);
xticks([0 20 40 60 80 100]);
xlabel('Percentage (%)');
ylabel('Brain Region');

% Action plots
figure('units','normalized','outerposition',[0 0 0.5 0.75])
%sgtitle('Planning Stats');
% Action onset timing
subplot(2,2,1);
h = boxplot(first_sig_timebin_Action_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = first_sig_timebin_Action_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Action Onset Response');
set(h,{'linew'},{1.5})
xlim([94 174]);
xticks([94 104 114 124 134 144 154 164 174]);
xticklabels([4.5 5 5.5 6 6.5 7 7.5 8 8.5]);
xlabel('Time (s)');
ylabel('Brain Region');

% Action onset percentage
subplot(2,2,2);
h = boxplot(first_sig_perc_Action_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder);%, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = first_sig_perc_Action_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Action Onset Percentage');
set(h,{'linew'},{1.5})
xlim([0 100]);
xticks([0 20 40 60 80 100]);
xlabel('Percentage (%)');
ylabel('Brain Region');

% Action peak timing
subplot(2,2,3);
h = boxplot(peak_Action_timebin_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = peak_Action_timebin_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Action Peak Response');
set(h,{'linew'},{1.5})
xlim([94 174]);
xticks([94 104 114 124 134 144 154 164 174]);
xticklabels([4.5 5 5.5 6 6.5 7 7.5 8 8.5]);
xlabel('Time (s)');
ylabel('Brain Region');

% Action peak percentage
subplot(2,2,4);
h = boxplot(peak_Action_perc_all, brainAreas, "Orientation", "horizontal", "BoxStyle","outline","Colors",colorOrder, "GroupOrder",groupOrder); %, 'Notch', 'on', 'Labels', {'SMG','PMv','S1','AIP','M1'})
hold on
% Scatter plot of individual session data points
for n_region = 1:numel(brainAreas)
    x_vals = peak_Action_perc_all(:, n_region);   % data values along *x-axis*
    y_vals = repmat(reorderedIdx(n_region), length(x_vals), 1); % region index along *y-axis*
    
    scatter(x_vals, y_vals, 10, 'k', 'o', ...
            'MarkerFaceAlpha', 0.2);
end
title('Action Peak Percentage');
set(h,{'linew'},{1.5})
xlim([0 100]);
xticks([0 20 40 60 80 100]);
xlabel('Percentage (%)');
ylabel('Brain Region');

%% load multiple subjects to plot together
% load GRASP data for S2
subject_id = 's2';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S2 GRASP %
G2_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
G2_first_sig_idx_all = first_sig_idx_all;
G2_first_sig_perc_all = first_sig_perc_all;
G2_peak_Cue_idx_all = peak_Cue_idx_all;
G2_peak_Cue_perc_all = peak_Cue_perc_all;
%G_peak_perc_all = peak_perc_all; % not in original pipeline
%G_peak_idx_all = peak_idx_all;

% load GRASP data for S3
subject_id = 's3';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S3 GRASP %
G3_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
G3_first_sig_idx_all = first_sig_idx_all;
G3_first_sig_perc_all = first_sig_perc_all;
G3_peak_Cue_idx_all = peak_Cue_idx_all;
G3_peak_Cue_perc_all = peak_Cue_perc_all;
%G_peak_perc_all = peak_perc_all; % not in original pipeline
%G_peak_idx_all = peak_idx_all;

% load GRASP data for S4
subject_id = 's4';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S4 GRASP %
G4_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
G4_first_sig_idx_all = first_sig_idx_all;
G4_first_sig_perc_all = first_sig_perc_all;
G4_peak_Cue_idx_all = peak_Cue_idx_all;
G4_peak_Cue_perc_all = peak_Cue_perc_all;
%G_peak_perc_all = peak_perc_all; % not in original pipeline
%G_peak_idx_all = peak_idx_all;

% % load Data OBJECT/SIZE
% goLabel = ["NoGo", "Go"];
% goLabel = goLabel(flagGoTrials + 1);
% directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
% filename = "decoded_sizeExtremes_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; % SIZE
% %filename = "decoded_objects_per_timebin_" + taskName + '_ALL_REGIONS' + "_LDA_" + goLabel + "_z_scored.mat"; % OBJECT
% full_path = fullfile(directory, filename);
% load(full_path);
% 
% % OBJECT/SIZE EXTREMES %
% X_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
% X_first_sig_idx_all = first_sig_idx_all;
% X_first_sig_perc_all = first_sig_perc_all;
% X_peak_Cue_idx_all = peak_Cue_idx_all;
% X_peak_Cue_perc_all = peak_Cue_perc_all;
% %X_peak_perc_all = peak_perc_all; % not in original pipeline
% %X_peak_idx_all = peak_idx_all;

% load SHUFFLED GRASP data for S2
subject_id = 's2';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S2 SHUFFLED GRASP %
SG2_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
SG2_first_sig_idx_all = first_sig_idx_all;
SG2_first_sig_perc_all = first_sig_perc_all;
SG2_peak_Cue_idx_all = peak_Cue_idx_all;
SG2_peak_Cue_perc_all = peak_Cue_perc_all;
SG2_peak_perc_all = peak_perc_all;
SG2_peak_idx_all = peak_idx_all;

% load SHUFFLED GRASP data for S3
subject_id = 's3';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S3 SHUFFLED GRASP %
SG3_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
SG3_first_sig_idx_all = first_sig_idx_all;
SG3_first_sig_perc_all = first_sig_perc_all;
SG3_peak_Cue_idx_all = peak_Cue_idx_all;
SG3_peak_Cue_perc_all = peak_Cue_perc_all;
SG3_peak_perc_all = peak_perc_all;
SG3_peak_idx_all = peak_idx_all;

% load SHUFFLED GRASP data for S4
subject_id = 's4';
goLabel = ["NoGo", "Go"];
goLabel = goLabel(flagGoTrials + 1);
directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
filename = "decoded_grasps_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
full_path = fullfile(directory, filename);
load(full_path);

% S2 SHUFFLED GRASP %
SG4_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
SG4_first_sig_idx_all = first_sig_idx_all;
SG4_first_sig_perc_all = first_sig_perc_all;
SG4_peak_Cue_idx_all = peak_Cue_idx_all;
SG4_peak_Cue_perc_all = peak_Cue_perc_all;
SG4_peak_perc_all = peak_perc_all;
SG4_peak_idx_all = peak_idx_all;

% % load SHUFFLED to replace chance line
% % load Data SHUFFLED OBJECTS/SIZE
% goLabel = ["NoGo", "Go"];
% goLabel = goLabel(flagGoTrials + 1);
% directory = ['C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction\analyzedData\' subject_id];
% %filename = "decoded_objects_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
% filename = "decoded_sizeExtremes_per_timebin_" + taskName + '_ALL_REGIONS_Shuffled' + "_LDA_" + goLabel + "_z_scored.mat"; 
% full_path = fullfile(directory, filename);
% load(full_path);
% 
% % SHUFFLED OBJECT/SIZE %
% SX_all_errTest_timebin_all_regions = all_errTest_timebin_all_regions;
% SX_first_sig_idx_all = first_sig_idx_all;
% SX_first_sig_perc_all = first_sig_perc_all;
% SX_peak_Cue_idx_all = peak_Cue_idx_all;
% SX_peak_Cue_perc_all = peak_Cue_perc_all;
% SX_peak_perc_all = peak_perc_all;
% SX_peak_idx_all = peak_idx_all;

%% plot altogether
% Determine number of items to be plotted together (decoded items, ie. grasp, object, size, shuffled)
n_decoded = 4;
%chance = 1 / 4 * 100; 
numPhases = numel(phase_changes);
n_subjects = 3;

% Create figure
figure('units','normalized','outerposition',[0 0 0.195 0.22]); %.137 or .08
hold on;
plot_handles = gobjects(n_decoded, 1);
%color_info = {[0.3359, 0.7031, 0.9101],[0.8984 0.6211 0],[0.8320 0.3672 0],[0.7969, 0.4726, 0.6523],[0, 0.6171, 0.4492],[.9961 .6875 0]}; % SMG, PMV, S1, AIP, M1, dlPFC
% colors for grasp only
color_info = {[0, 0.3529, 0.7098],[0 0 0],[0.75 0.75 0.75]};
% colors for object
%color_info = {[0, 0.3529, 0.7098],[0.9608, 0.7333, 0.1176],[0.75 0.75 0.75]}; % blue for grasp, yellow for object (pull from FR), gray for shuffled
% colors for size
%color_info = {[0, 0.3529, 0.7098],[0.2039, 0.6902, 0.2902],[0.75 0.75 0.75]}; % blue for grasp, green for size (pull from FR), gray for shuffled

for n_region = 5 %1:n_regions
    % PLOT S2 GRASP
    G2_mean_acc = mean(G2_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    G2_CI95 = utile.calculate_CI(G2_all_errTest_timebin_all_regions(:,:,n_region));
    G2_ci_upper = G2_CI95(2, :);
    G2_ci_lower = G2_CI95(1, :);

    G2_ER = utile.shadedErrorBar(1:num_timebins, G2_mean_acc, G2_ci_upper, 'lineprops', '-.');
    G2_ER.mainLine.Color = color_info{1};
    G2_ER.mainLine.LineWidth = 2;
    G2_ER.patch.FaceColor = color_info{1};
    G2_ER.edge(1).LineStyle = 'none';
    G2_ER.edge(2).LineStyle = 'none';
    plot_handles(1) = G2_ER.mainLine;

    hold on;

    % PLOT S2 SHUFFLED GRASP
    SG2_mean_acc = mean(SG2_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    SG2_CI95 = utile.calculate_CI(SG2_all_errTest_timebin_all_regions(:,:,n_region));
    SG2_ci_upper = SG2_CI95(2, :);
    SG2_ci_lower = SG2_CI95(1, :);

    SG2_ER = utile.shadedErrorBar(1:num_timebins, SG2_mean_acc, SG2_ci_upper, 'lineprops', '-.');
    SG2_ER.mainLine.Color = color_info{3};
    SG2_ER.mainLine.LineWidth = 2;
    SG2_ER.patch.FaceColor = color_info{3};
    SG2_ER.edge(1).LineStyle = 'none';
    SG2_ER.edge(2).LineStyle = 'none';
    plot_handles(2) = SG2_ER.mainLine;

    hold on;

    % PLOT S3 GRASP
    G3_mean_acc = mean(G3_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    G3_CI95 = utile.calculate_CI(G3_all_errTest_timebin_all_regions(:,:,n_region));
    G3_ci_upper = G3_CI95(2, :);
    G3_ci_lower = G3_CI95(1, :);

    G3_ER = utile.shadedErrorBar(1:num_timebins, G3_mean_acc, G3_ci_upper, 'lineprops', '-b');
    G3_ER.mainLine.Color = color_info{1};
    G3_ER.mainLine.LineWidth = 2;
    G3_ER.patch.FaceColor = color_info{1};
    G3_ER.edge(1).LineStyle = 'none';
    G3_ER.edge(2).LineStyle = 'none';
    plot_handles(3) = G3_ER.mainLine;

    hold on;

    % PLOT S3 SHUFFLED GRASP
    SG3_mean_acc = mean(SG3_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    SG3_CI95 = utile.calculate_CI(SG3_all_errTest_timebin_all_regions(:,:,n_region));
    SG3_ci_upper = SG3_CI95(2, :);
    SG3_ci_lower = SG3_CI95(1, :);

    SG3_ER = utile.shadedErrorBar(1:num_timebins, SG3_mean_acc, SG3_ci_upper, 'lineprops', '-b');
    SG3_ER.mainLine.Color = color_info{3};
    SG3_ER.mainLine.LineWidth = 2;
    SG3_ER.patch.FaceColor = color_info{3};
    SG3_ER.edge(1).LineStyle = 'none';
    SG3_ER.edge(2).LineStyle = 'none';
    plot_handles(4) = SG3_ER.mainLine;

    hold on;

    % PLOT S4 GRASP
    G4_mean_acc = mean(G4_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    G4_CI95 = utile.calculate_CI(G4_all_errTest_timebin_all_regions(:,:,n_region));
    G4_ci_upper = G4_CI95(2, :);
    G4_ci_lower = G4_CI95(1, :);

    G4_ER = utile.shadedErrorBar(1:num_timebins, G4_mean_acc, G4_ci_upper, 'lineprops', '.');
    G4_ER.mainLine.Color = color_info{1};
    G4_ER.mainLine.LineWidth = 2;
    G4_ER.patch.FaceColor = color_info{1};
    G4_ER.edge(1).LineStyle = 'none';
    G4_ER.edge(2).LineStyle = 'none';
    plot_handles(5) = G4_ER.mainLine;

    hold on;

    % PLOT S4 SHUFFLED GRASP
    SG4_mean_acc = mean(SG4_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    SG4_CI95 = utile.calculate_CI(SG4_all_errTest_timebin_all_regions(:,:,n_region));
    SG4_ci_upper = SG4_CI95(2, :);
    SG4_ci_lower = SG4_CI95(1, :);

    SG4_ER = utile.shadedErrorBar(1:num_timebins, SG4_mean_acc, SG4_ci_upper, 'lineprops', '.');
    SG4_ER.mainLine.Color = color_info{3};
    SG4_ER.mainLine.LineWidth = 2;
    SG4_ER.patch.FaceColor = color_info{3};
    SG4_ER.edge(1).LineStyle = 'none';
    SG4_ER.edge(2).LineStyle = 'none';
    plot_handles(6) = SG4_ER.mainLine;

    hold on;

    % % PLOT X
    % X_mean_acc = mean(X_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    % X_CI95 = utile.calculate_CI(X_all_errTest_timebin_all_regions(:,:,n_region));
    % X_ci_upper = X_CI95(2, :);
    % X_ci_lower = X_CI95(1, :);
    % 
    % X_ER = utile.shadedErrorBar(1:num_timebins, X_mean_acc, X_ci_upper, 'lineprops', '-b');
    % X_ER.mainLine.Color = color_info{2};
    % X_ER.mainLine.LineWidth = 2;
    % X_ER.patch.FaceColor = color_info{2};
    % X_ER.edge(1).LineStyle = 'none';
    % X_ER.edge(2).LineStyle = 'none';
    % plot_handles(3) = X_ER.mainLine;
    % 
    % hold on;
    % 
    % % PLOT SHUFFLED X
    % SX_mean_acc = mean(SX_all_errTest_timebin_all_regions(:,:,n_region), 1, 'omitnan');
    % SX_CI95 = utile.calculate_CI(SX_all_errTest_timebin_all_regions(:,:,n_region));
    % SX_ci_upper = SX_CI95(2, :);
    % SX_ci_lower = SX_CI95(1, :);
    % 
    % SX_ER = utile.shadedErrorBar(1:num_timebins, SX_mean_acc, SX_ci_upper, 'lineprops', '--');
    % SX_ER.mainLine.Color = color_info{3};
    % SX_ER.mainLine.LineWidth = 2;
    % SX_ER.patch.FaceColor = color_info{3};
    % SX_ER.edge(1).LineStyle = 'none';
    % SX_ER.edge(2).LineStyle = 'none';
    % plot_handles(4) = SX_ER.mainLine;

end

% Finalize plot
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase},'LineWidth', 1.5, 'FontSize', 10); %phaseNames{n_phase}
end
xlim([0 174]);
xtickangle();
xlabel('Time (s)');
xticks([0 42 83 124 165]); %165
xticklabels([0 2 4 6 8]); %8
ylabel('Classification Accuracy [%]');
title('Grasp Classification');
%yline(chance, '--r', 'LineWidth', 1.5);
ylim([0 100]); %40 80 ; 15 80
yticks([0 25 50 75 100]); %50 80  ; 25 80
%legend(plot_handles, [{'S2'},{'S2 Shuffled'},{'S3'},{'S3 Shuffled'},{'S4'},{'S4 Shuffled'},], 'Location', 'best'); %[{'Grasp'},{'Shuffled Grasp'},{'Object'},{'Shuffled Object'}]
set(gca, 'FontSize', 10);
hold off;