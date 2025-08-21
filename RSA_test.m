clc
clear all
close all

subject_id = 's3';
unit_region = 'AIP';

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
flag_4S = true; % true = updated 4S action phase; false = original 2S action phase
flag_shuffled = false; % true = shuffled images task
flag_varied_sizes = false; % true for varied sizes task
flag_GB_images = false; % true for task using images of GB's own hands and real objects
flag_5050 = false; % true for 50% Go 50% NoGo trials
flag_combined = true; % true for combinations task

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

if flag_varied_sizes
    TaskCue = 'GraspObject_Varied_Size';
    min_timebin_length = 174;
end
if flag_GB_images
    TaskCue = 'GraspObject_GB_Images';
    min_timebin_length = 174;
end

if flag_5050
    TaskCue = 'GraspObject_5050';
    min_timebin_length = 174;
end

if flag_combined
    TaskCue = 'GraspObject_Combined';
    min_timebin_length = 174;
end

% Task Variables
% 4S data
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' TaskCue spike_sorting_type]);
Go_data = Data.Go_data;
color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

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

% add Aperture Size column
if strcmp(TaskCue, 'GraspObject_Varied_Size')
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

if strcmp(TaskCue, 'GraspObject_Combined')
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
    color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
end

flagGoTrials = true; % false = No-Go

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

%chose cue type:
taskCuesAll = {'Hand', 'Hand-Object', 'Object'};
if flag_combined
    taskCuesAll = {'Combined','Hand', 'Hand-Object', 'Object'};
end
sessions_all = unique(Go_data.session_date);
numSessions = numel(sessions_all);
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phase_bin_ranges = {
    1:phase_changes(2)-1;                   % ITI
    phase_changes(2):phase_changes(3)-1;    % Cue
    phase_changes(3):phase_changes(4)-1;    % Delay
    phase_changes(4):numel(phase_time_idx)  % Action
};
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

numUnitsPerSession = zeros(numSessions,1);
%%
% Example: assume 64 trials, 4 grasps × 4 objects
n_grasps = 4;
n_objects = 4;
n_conditions = n_grasps * n_objects; % 16

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

    %labels 
    sessionLabels = Go_data.GoLabels(idxThisSession,:);
    
    %trialType
    trialTypeSession = Go_data.TrialType(idxThisSession,:);

    % grasp labels
    graspTypeSession = Go_data.GraspType(idxThisSession,:);

    %ApertureSize
    %apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);

    % object labels
    objectTypeSession = Go_data.ObjectType(idxThisSession,:);

    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    timePhaseLabels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        SessionData = SessionData(GoNoGoidx);
        sessionLabels = sessionLabels(GoNoGoidx);
        timePhaseLabels = timePhaseLabels(GoNoGoidx);
        trialTypeSession = trialTypeSession(GoNoGoidx);
        graspTypeSession = graspTypeSession(GoNoGoidx);
        %apertureSizeSession = apertureSizeSession(GoNoGoidx);
        objectTypeSession = objectTypeSession(GoNoGoidx);
    else
        SessionData = z_scored_data(~GoNoGoidx);
        sessionLabels = sessionLabels(~GoNoGoidx);
        timePhaseLabels = timePhaseLabels(~GoNoGoidx);
        trialTypeSession = trialTypeSession(~GoNoGoidx);
    end
     
    % seperate data according to cue modality 
    unTrialType = unique(Go_data.TrialType);

    numUnitsPerSession(n_session) = size(SessionData{1},2);

    sessionLabels_modality = trialTypeSession;
    sessionLabels_grasp = graspTypeSession;
    %sessionLabels_size = apertureSizeSession;
    sessionLabels_object = objectTypeSession;

    % Convert word labels (ie.,'Hand', 'HandObject', 'Object') to numerical values
    modality_labels = {'Combined','Hand', 'Hand_Object', 'Object'}; % 'Combined'
    grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
    %size_labels = {'Small', 'Medium', 'Large'};
    object_labels = {'block','rod','deck','ball'}; % numbering mirrors associated grasp

    % Initialize numerical labels
    sessionLabels_modality_num = zeros(size(sessionLabels_modality)); 
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

    % create condition ID for each trial
    condition_ids = (sessionLabels_grasp_num - 1) * n_objects + sessionLabels_object_num; % labeled 1-16 for each trial

    % loop through cue modalities 
    for n_type = 1 % Combined dataset only

        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
        sessionLabels_object_num = sessionLabels_object_num(trialTypeIdx); % TEST for object tuning, idx so only Combined dataset used
        sessionLabels_grasp_num = sessionLabels_grasp_num(trialTypeIdx);
        SessionData = SessionData(trialTypeIdx); % TEST keep data only for Combined trials
        timePhaseLabels = timePhaseLabels(trialTypeIdx);
        condition_ids = condition_ids(trialTypeIdx);
        n_trials = numel(SessionData);

        for n_phase = [2,4]
            phase_window = phase_bin_ranges{n_phase};
            n_neurons = size(SessionData{1},2);
            avg_FR = zeros(n_trials, n_neurons);
            for i = 1:n_trials
                avg_FR(i, :) = mean(SessionData{i}(phase_window, :), 1);  % mean over timebins
            end
        end
    end
end
%% Bootstrapped RSA
n_bootstraps = 1000;
min_trials_per_condition = min(histcounts(condition_ids, 0.5:1:(n_conditions+0.5)));

all_dissimilarity = zeros(n_conditions, n_conditions, n_bootstraps);

for b = 1:n_bootstraps
    sampled_vectors = zeros(n_conditions, n_neurons);
    
    for c = 1:n_conditions
        trials_c = find(condition_ids == c);
        sampled_idx = randsample(trials_c, min_trials_per_condition, false);
        sampled_vectors(c, :) = mean(avg_FR(sampled_idx, :), 1);
    end
    
    dissim = 1 - corr(sampled_vectors');  % 1 - Pearson correlation
    dissim(1:n_conditions+1:end) = 0;
    all_dissimilarity(:,:,b) = dissim;
end

mean_dissimilarity = mean(all_dissimilarity, 3);  % average RDM
% Force symmetry
mean_dissimilarity = (mean_dissimilarity + mean_dissimilarity') / 2;

%% visualize RDM
figure;
imagesc(mean_dissimilarity);
axis square;
colorbar;
title('Mean Representational Dissimilarity Matrix (Bootstrapped)');
xlabel('Condition'); ylabel('Condition');

%% MDS visualization
[Y, ~] = mdscale(mean_dissimilarity, 2);

figure; hold on;
for i = 1:n_conditions
    g = mod(i-1, n_grasps) + 1;
    o = floor((i-1) / n_grasps) + 1;
    label = sprintf('%s-%s', grasp_labels{g}, object_labels{o});
    scatter(Y(i,1), Y(i,2), 100, 'filled');
    text(Y(i,1), Y(i,2), label, ...
        'VerticalAlignment','bottom','HorizontalAlignment','right');
end
xlabel('MDS Dim 1'); ylabel('MDS Dim 2');
title('MDS Embedding with Grasp-Object Labels');
axis equal;

%% color coding
object_per_condition = mod((1:n_conditions) - 1, n_grasps) + 1;  % [1 x 16]
grasp_per_condition = floor(((1:n_conditions) - 1) / n_grasps) + 1;  % [1 x 16]

grasp_colors = {[0.2, 0.13, 0.53], [0.067, 0.467, 0.2], [0.53, 0.8, 0.93], [0.53, 0.13, 0.33]}; % Purple, Green, Light Blue, Dark Pink
%object_colors = {[.9961 .6875 0],[.3906 .5586 .9961],[.4688 .3672 .9375],[.8594 .1484 .4961],[.9922 .3789 0]}; % objects (including Assoc.)

figure('units','normalized','outerposition',[0 0 0.25 0.33]); hold on;
for i = 1:n_conditions
    g = object_per_condition(i);
    scatter(Y(i,1), Y(i,2), 100, 'filled', 'MarkerFaceColor', grasp_colors{g});
end
xlabel('MDS Dim 1'); ylabel('MDS Dim 2');
title([sessions_all{n_session} ' - MDS: Grasp Grouping - ' phaseNames{n_phase}]);
legend(grasp_labels, 'Location', 'bestoutside'); 
axis equal;

% plot by object type
figure('units','normalized','outerposition',[0 0 0.25 0.33]); hold on;
for i = 1:n_conditions
    o = grasp_per_condition(i);
    scatter(Y(i,1), Y(i,2), 100, 'filled', 'MarkerFaceColor', grasp_colors{o});
end
% plot one point for each object to use in legend
object_handles = gobjects(n_objects,1);
for o = 1:n_objects
    object_handles(o) = scatter(nan, nan, 100, 'filled', 'MarkerFaceColor', grasp_colors{o});
end

% Then create the legend
legend(object_handles, object_labels, 'Location', 'bestoutside');

xlabel('MDS Dim 1');
ylabel('MDS Dim 2');
title([sessions_all{n_session} ' - MDS: Object Grouping - ' phaseNames{n_phase}]);
axis equal;

%% plotting all sessions
n_grasps = 4;
n_objects = 4;
n_conditions = n_grasps * n_objects; % 16
grasp_colors = {[0.2, 0.13, 0.53], [0.067, 0.467, 0.2], [0.53, 0.8, 0.93], [0.53, 0.13, 0.33]}; % Purple, Green, Light Blue, Dark Pink
object_per_condition = mod((1:n_conditions) - 1, n_grasps) + 1;  % [1 x 16]
grasp_per_condition = floor(((1:n_conditions) - 1) / n_grasps) + 1;  % [1 x 16]
n_bootstraps = 1000;


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

    %labels 
    sessionLabels = Go_data.GoLabels(idxThisSession,:);
    
    %trialType
    trialTypeSession = Go_data.TrialType(idxThisSession,:);

    % grasp labels
    graspTypeSession = Go_data.GraspType(idxThisSession,:);

    %ApertureSize
    %apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);

    % object labels
    objectTypeSession = Go_data.ObjectType(idxThisSession,:);

    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    timePhaseLabels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        SessionData = SessionData(GoNoGoidx);
        sessionLabels = sessionLabels(GoNoGoidx);
        timePhaseLabels = timePhaseLabels(GoNoGoidx);
        trialTypeSession = trialTypeSession(GoNoGoidx);
        graspTypeSession = graspTypeSession(GoNoGoidx);
        %apertureSizeSession = apertureSizeSession(GoNoGoidx);
        objectTypeSession = objectTypeSession(GoNoGoidx);
    else
        SessionData = z_scored_data(~GoNoGoidx);
        sessionLabels = sessionLabels(~GoNoGoidx);
        timePhaseLabels = timePhaseLabels(~GoNoGoidx);
        trialTypeSession = trialTypeSession(~GoNoGoidx);
    end
     
    % seperate data according to cue modality 
    unTrialType = unique(Go_data.TrialType);

    numUnitsPerSession(n_session) = size(SessionData{1},2);

    sessionLabels_modality = trialTypeSession;
    sessionLabels_grasp = graspTypeSession;
    %sessionLabels_size = apertureSizeSession;
    sessionLabels_object = objectTypeSession;

    % Convert word labels (ie.,'Hand', 'HandObject', 'Object') to numerical values
    modality_labels = {'Combined','Hand', 'Hand_Object', 'Object'}; % 'Combined'
    grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
    %size_labels = {'Small', 'Medium', 'Large'};
    object_labels = {'block','rod','deck','ball'}; % numbering mirrors associated grasp

    % Initialize numerical labels
    sessionLabels_modality_num = zeros(size(sessionLabels_modality)); 
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

    % create condition ID for each trial
    condition_ids = (sessionLabels_grasp_num - 1) * n_objects + sessionLabels_object_num; % labeled 1-16 for each trial

    % loop through cue modalities 
    for n_type = 1 % Combined dataset only

        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
        sessionLabels_object_num = sessionLabels_object_num(trialTypeIdx); % TEST for object tuning, idx so only Combined dataset used
        sessionLabels_grasp_num = sessionLabels_grasp_num(trialTypeIdx);
        SessionData = SessionData(trialTypeIdx); % TEST keep data only for Combined trials
        timePhaseLabels = timePhaseLabels(trialTypeIdx);
        condition_ids = condition_ids(trialTypeIdx);
        n_trials = numel(SessionData);
        figure('Name', sessions_all{n_session}, 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.6 0.6]);
        n_axes = 4;  % e.g., Cue-Grasp, Action-Grasp, Cue-Object, Action-Object
        t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

        
        for n_phase = [2,4]
            phase_window = phase_bin_ranges{n_phase};            
            n_neurons = size(SessionData{1},2);
            avg_FR = zeros(n_trials, n_neurons);
            for i = 1:n_trials
                avg_FR(i, :) = mean(SessionData{i}(phase_window, :), 1);  % mean over timebins
            end
        
            % --- RSA with Bootstrapping ---
            min_trials_per_condition = min(histcounts(condition_ids, 0.5:1:(n_conditions+0.5)));
            all_dissimilarity = zeros(n_conditions, n_conditions, n_bootstraps);
        
            for b = 1:n_bootstraps
                sampled_vectors = zeros(n_conditions, n_neurons);
                for c = 1:n_conditions
                    trials_c = find(condition_ids == c);
                    sampled_idx = randsample(trials_c, min_trials_per_condition, false);
                    sampled_vectors(c, :) = mean(avg_FR(sampled_idx, :), 1);
                end
                dissim = 1 - corr(sampled_vectors');  % 1 - Pearson correlation
                dissim(1:n_conditions+1:end) = 0;
                all_dissimilarity(:,:,b) = dissim;
            end
        
            mean_dissimilarity = mean(all_dissimilarity, 3);
            mean_dissimilarity = (mean_dissimilarity + mean_dissimilarity') / 2;
            [Y, ~] = mdscale(mean_dissimilarity, 3);
        
            % plotting
            % Row index: 1 for Cue (phase 2), 2 for Action (phase 4)
            row_idx = find([2, 4] == n_phase);
            % --- MDS by Grasp ---
            nexttile(t, (row_idx - 1) * 2 + 1);  % Left column
            hold on;
            
            for i = 1:n_conditions
                g = grasp_per_condition(i);
                scatter(Y(i,1), Y(i,2), 100, 'filled', 'MarkerFaceColor', grasp_colors{g});
            end
            grasp_handles = gobjects(n_grasps,1);
            for g = 1:n_grasps
                grasp_handles(g) = scatter(nan, nan, 100, 'filled', 'MarkerFaceColor', grasp_colors{g});
            end
            legend(grasp_handles, grasp_labels, 'Location', 'eastoutside');
            title([sessions_all{n_session} ' - ' phaseNames{n_phase} ' (Grasp)']);
            axis equal; xlabel('Dim 1'); ylabel('Dim 2');
        
            % --- MDS by Object ---
            nexttile(t, (row_idx - 1) * 2 + 2);  % Right column
            hold on;
            for i = 1:n_conditions
                o = object_per_condition(i);
                scatter(Y(i,1), Y(i,2), 100, 'filled', 'MarkerFaceColor', grasp_colors{o});
            end
            legend(object_labels, 'Location', 'eastoutside');
            title([sessions_all{n_session} ' - ' phaseNames{n_phase} ' (Object)']);
            axis equal; xlabel('Dim 1'); ylabel('Dim 2');
        end
    end
end
