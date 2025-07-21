clc
clear all
close all

subject_id = 's3';
unit_region = 'SMG';

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
    
    % Z-scoring => calc the mean of the FR of each unit across all trials per timebin, not per phase)
    num_trials = length(SessionData);
    [num_timebins, num_units] = size(SessionData{1});
    % reconfigure matrix to store all trials
    all_data = NaN(num_timebins,num_units,num_trials); 
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

    %ApertureSize
    %apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);

    % object labels
    objectTypeSession = Go_data.ObjectType(idxThisSession,:);

    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    timePhaseLabels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        SessionData = SessionData(GoNoGoidx);
        %SessionData = z_scored_data(GoNoGoidx);
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

    % loop through cue modalities 
    for n_type = 1 % Combined dataset only

        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
        sessionLabels_object_num = sessionLabels_object_num(trialTypeIdx); % TEST for object tuning, idx so only Combined dataset used
        sessionLabels_grasp_num = sessionLabels_grasp_num(trialTypeIdx);
        SessionData = SessionData(trialTypeIdx); % TEST keep data only for Combined trials
        timePhaseLabels = timePhaseLabels(trialTypeIdx);
    end
end

%% loop through phases
for n_phase = 1:numPhases
    phase_window = phase_bin_ranges{n_phase}; % example: 500–1000ms (50ms bins)
    avg_FR = squeeze(mean(FR(phase_window, :, :), 1)); % [n_trials x n_neurons]
end

