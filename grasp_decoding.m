clc
clear all
close all 

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
taskName = 'GraspObject_Varied_Size'; % varied object/aperture sizes
%taskName = 'GraspObject_GB_Images'; % for GB
%taskName = 'GraspObject_Combined'; % for Combined task
subject_id = 's3';

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName spike_sorting_type]);

%% Analysis

if ~strcmp(taskName, 'GraspObject_Varied_Size')
    Go_data = Data.Go_data;

    color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

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
        color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
    end
    
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
    
    unit_region = 'M1';
    brainAreas = Go_data.frPerChannel{7};
    phase_time_idx = Go_data.time_phase_labels{1,1};
    numPhases = numel(unique(phase_time_idx));
    phaseTimeTmp = diff(phase_time_idx);
    phase_changes(1) = 1;
    phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
    phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
    %color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Conditions
    
    %define brain region. 
    % loop to session days 
        % separate data into go and no go
            % seperate data according to cue modality 
                %seperate according to task phases
                    %decode 
    
    
    
    sessions_all = unique(Go_data.session_date);
    numSessions = numel(sessions_all);
    uniqueCueTypes = {'Hand','Hand-Object','Object'};
    if strcmp(taskName, 'GraspObject_Combined')
        uniqueCueTypes = {'Combined','Hand','Hand-Object','Object'};
    end
    
    flagGoTrials = true; %if true, extract Go trials, if false, extract NoGo trials

    % Initialize a matrix to store classification accuracy for all sessions
    all_errTest = NaN(numSessions, numPhases, numel(uniqueCueTypes)); % for sep by Cue Modality
    %all_errTest = NaN(numSessions, numPhases, numel(uniqueGraspTypes)); % for sep by Grasp Types
    
    % Initialize a cell array to store confusion matrices for each phase
    confMatAllSessions = cell(numSessions,numPhases,numel(uniqueCueTypes)); %cell(numSessions,numPhases); 1 bc only running for Cue rn
    figure(); 
    
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
    
        % loop through cue conditions
        for n_type = 1:numel(unTrialType) 
            
            % find idx of trial type 
            trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
    
            % code for CMs of grasps and modalities separately
            sessionLabels_modality = trialTypeSession;
            sessionLabels_grasp = graspTypeSession;
            sessionLabels_object = objectTypeSession;

            % Convert modality labels ('Hand', 'HandObject', 'Object') to numerical values
            modality_labels = {'Hand', 'Hand_Object', 'Object'}; % 'Combined',
            grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
            object_labels = {'block','rod','deck','ball'};

            sessionLabels_modality_num = zeros(size(sessionLabels_modality));  % Initialize numerical labels
            sessionLabels_grasp_num = zeros(size(sessionLabels_grasp));
            sessionLabels_object_num = zeros(size(sessionLabels_object));

            % Loop through labels and assign numerical values
            for i = 1:length(modality_labels)
                sessionLabels_modality_num(strcmp(sessionLabels_modality, modality_labels{i})) = i;
            end  
            for i = 1:length(grasp_labels)
                sessionLabels_grasp_num(strcmp(sessionLabels_grasp, grasp_labels{i})) = i;
            end 
            for i = 1:length(object_labels)
                sessionLabels_object_num(strcmp(sessionLabels_object, object_labels{i})) = i;
            end 
    
            % loop through task phases 
            for n_phase = 1:numPhases
           
                data_per_phase_per_cue = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData(trialTypeIdx),time_phase_labels(trialTypeIdx), 'UniformOutput', false));
                %data_per_phase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false)); %1:128 = H and HO
    
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % original
                [errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels_grasp_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % testing with my method to see if we still get same decoding
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels_object_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true);
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(data_per_phase,sessionLabels_grasp_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                errTest(n_phase,n_type) =  (1-mean(errTestTmp))*100; % sub 1 for n_type when separating by condition first
                %title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase}]) % ' - ' unTrialType{n_type}])
                

                % Store the confusion matrix for this phase and session
                confMatAllSessions{n_session, n_phase, n_type} = confMat; % sub 1 for n_type when separating by condition first
    
                %data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false));
                
                % % confusion matrices
                % if n_phase ~= 1
                %     %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true); % original
                %     [errTrain, errTestTmp,ConfMat] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels_grasp_num,'flagErrorMatrix', true, 'PCA_variance', 95,'flagLeaveOneOut', true); % CM for grasps/modalities separately (can switch out sessionLabels_grasp/modality)
                %     % error matrices are generating properly but are not
                %     % saving values into ConfMat, rather the individual
                %     % classifications, so look into how to get the values
                %     % found when flagging ErrorMatrix
                %     % Store the confusion matrix
                %     confMatAllSessions{n_session, n_phase} = ConfMat;
                %     title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                % end

                % confusion matrices
                % if n_phase ~= 1
                %     %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true); % original
                %     [errTrain, errTestTmp,trueLabels,predictedLabels,cm] = classification.LDA_classification_rep(z_data_per_phase_per_all,sessionLabels_grasp_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %     % Store the confusion matrix
                %     confMatAllSessions{n_session, n_phase} = cm;
                %     title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                % end
                
                set(gca, 'FontSize', 12);
            end   
        end 
        
        % Store the classification accuracy for each session
        all_errTest(n_session, :, :) = errTest; % add ,:) aka 3rd dimension when separating by CONDITION
    
        subplot(4,2,n_session) % code using var names instead of hard 
    
      %  figure();
        b = bar(errTest);
        colors = cell2mat(color_info');
        for n_color = 1:length(b)
            b(n_color).FaceColor = colors(n_color,:);
        end
        
        %yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100)
        yline(1/numel(unique(unGraspType))*100)
        ylim([0 100])
        xticklabels(phaseNames)
        
        % if n_session ==3
        %     legend(unTrialType, 'Interpreter','none')
        % end 
        ylabel('Classification percentage [%]')
        title(sessions_all{n_session})
        %SessionData = 
    
    end
    legend(unTrialType, 'Interpreter','none')
    sgtitle(unit_region)

else 
    Go_data = Data.Go_data;

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
    
    % remove faulty sessions, if any
    error_session = {};
    if strcmp(subject_id, 's2')
        error_session = {'20231016'};
    elseif strcmp(subject_id, 's3')
        error_session = {};
    end 
    
    if ~isempty(error_session)
        condition = cellfun(@(x) strcmp(x, error_session), Go_data.session_date);
        Go_data = Go_data(~condition,:);
    end
    
    unit_region = 'SMG';
    brainAreas = Go_data.frPerChannel{7};
    phase_time_idx = Go_data.time_phase_labels{1,1};
    numPhases = numel(unique(phase_time_idx));
    phaseTimeTmp = diff(phase_time_idx);
    phase_changes(1) = 1;
    phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
    phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
    
    uniqueGraspTypes = unique(Go_data.GraspType);
    uniqueCueTypes = unique(Go_data.TrialType);
    uniqueApertureSize = unique(Go_data.Aperture_Size); % comment out when not varied sizes
    
    %define brain region. 
    % loop to session days 
        % separate data into go and no go
            % seperate data according to cue modality 
                % separate according to grasp
                    %seperate according to task phases
                        %decode 
    
    
    
    sessions_all = unique(Go_data.session_date);
    numSessions = numel(sessions_all);
    
    
    flagGoTrials = true; %if true, extract Go trials, if false, extract NoGo trials
    
    % Initialize a matrix to store classification accuracy for all sessions
    %all_errTest = NaN(numSessions, numPhases, numel(uniqueCueTypes)); % for sep by Cue Modality
    %all_errTest = NaN(numSessions, numPhases, numel(uniqueGraspTypes)); % for sep by Grasp Types
    %all_errTest = NaN(numSessions, numPhases, numel(uniqueApertureSize)); % for sep by Aperture Size
    all_errTest = NaN(numSessions, numPhases); % unseparated
    
    % Initialize a cell array to store confusion matrices for each phase
    %confMatAllSessions = cell(numSessions,numPhases,numel(uniqueCueTypes)); %grasp/cueType ; cell(numSessions,numPhases); 1 bc only running for Cue rn
    
    % Original code
    figure(); 
    
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
    
        % grasp labels
        graspTypeSession = Go_data.GraspType(idxThisSession,:);
        
        %trialType
        trialTypeSession = Go_data.TrialType(idxThisSession,:);
    
        %ApertureSize
        apertureSizeSession = Go_data.Aperture_Size(idxThisSession,:);
    
        %get idx for Go or NoGo trials
        GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
        time_phase_labels = Go_data.time_phase_labels(idxThisSession);
        
    
        if flagGoTrials
            %SessionData = SessionData(GoNoGoidx);
            z_scored_data = z_scored_data(GoNoGoidx);
            sessionLabels = sessionLabels(GoNoGoidx);
            time_phase_labels = time_phase_labels(GoNoGoidx);
            trialTypeSession = trialTypeSession(GoNoGoidx);
            graspTypeSession = graspTypeSession(GoNoGoidx);
            apertureSizeSession = apertureSizeSession(GoNoGoidx);
        else
            SessionData = SessionData(~GoNoGoidx);
            sessionLabels = sessionLabels(~GoNoGoidx);
            time_phase_labels = time_phase_labels(~GoNoGoidx);
            trialTypeSession = trialTypeSession(~GoNoGoidx);
            apertureSizeSession = apertureSizeSession(~GoNoGoidx);
        end
         
        %seperate data according to cue modality 
    
        unTrialType = unique(Go_data.TrialType);
    
        %seperate data according to size 
    
        unAperture = unique(Go_data.Aperture_Size);

        % separate data according to grasp

        unGraspType = unique(Go_data.GraspType);

        % only keeping Small and Large
        % Define which sizes to keep
        sizesToKeep = {'Small', 'Large'};

        % Find indices of trials belonging to 'small' or 'large' sessions
        SLsizeIdx = ismember(apertureSizeSession, sizesToKeep);

        % Extract trials that correspond to 'small' or 'large' sessions
        z_scored_data = z_scored_data(SLsizeIdx);
        sessionLabels = sessionLabels(SLsizeIdx);
        time_phase_labels = time_phase_labels(SLsizeIdx);
        trialTypeSession = trialTypeSession(SLsizeIdx);
        graspTypeSession = graspTypeSession(SLsizeIdx);
        apertureSizeSession = apertureSizeSession(SLsizeIdx);
        uniqueApertureSize = unique(apertureSizeSession);

        % loop through cue modalities/sizes/grasps
        %for n_type = 1:numel(unGraspType) %1:numel(unGraspType) %1:numel(unTrialType) %n_size = 1:numel(unAperture)       %
            
            % find idx of trial type 
            %trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
    
            % find idx of size type 
            %trialSizeIdx = ismember(apertureSizeSession, unAperture(n_size));

            % find idx of grasp type 
            %trialGraspIdx = ismember(graspTypeSession, unGraspType(n_type));
            
            sessionLabels_modality = trialTypeSession;
            sessionLabels_grasp = graspTypeSession;
            sessionLabels_size = apertureSizeSession;
             
            % Convert modality labels ('Hand', 'HandObject', 'Object') to numerical values
            modality_labels = {'Hand', 'Hand_Object', 'Object'};
            grasp_labels = {'Lateral', 'MediumWrap', 'PalmarPinch', 'Sphere3Finger'};
            size_labels = {'Small', 'Medium', 'Large'};

            sessionLabels_modality_num = zeros(size(sessionLabels_modality));  % Initialize numerical labels
            sessionLabels_grasp_num = zeros(size(sessionLabels_grasp));
            sessionLabels_size_num = zeros(size(sessionLabels_size));

            % Loop through labels and assign numerical values
            for i = 1:length(modality_labels)
                sessionLabels_modality_num(strcmp(sessionLabels_modality, modality_labels{i})) = i;
            end  
            for i = 1:length(grasp_labels)
                sessionLabels_grasp_num(strcmp(sessionLabels_grasp, grasp_labels{i})) = i;
            end 
            for i = 1:length(size_labels)
                sessionLabels_size_num(strcmp(sessionLabels_size, size_labels{i})) = i;
            end
    
            % loop through task phases 
            for n_phase = 1:numPhases
           
                % data_per_phase_per_cue = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData(trialTypeIdx),time_phase_labels(trialTypeIdx),'UniformOutput', false)); % original (sep by cue data)
                %data_per_phase_per_size = cell2mat(arrayfun(@(x,y)
                    %mean(x{1,1}(y{:}==
                    %n_phase,:),1),SessionData(trialSizeIdx),time_phase_labels(trialSizeIdx),
                    %'UniformOutput', false)); % sep by size data
                %data_per_phase_per_grasp = cell2mat(arrayfun(@(x,y)
                    %mean(x{1,1}(y{:}==
                    %n_phase,:),1),SessionData(trialGraspIdx),time_phase_labels(trialGraspIdx),
                    %'UniformOutput', false)); % sep by grasp data
                data_per_phase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),z_scored_data,time_phase_labels,'UniformOutput', false)); % all data

                %z_data_per_phase_per_condition = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),z_scored_data(trialTypeIdx),time_phase_labels(trialTypeIdx),'UniformOutput', false)); % sep z-scored data by condition per phase
                %z_data_per_phase_per_grasp = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),z_scored_data(trialGraspIdx),time_phase_labels(trialGraspIdx),'UniformOutput', false)); % sep by grasp data

                % [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % original
                % [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels_size_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying size from modality
                
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(z_data_per_phase_per_condition,sessionLabels_grasp_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying grasp from modality using z-scored data
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(z_data_per_phase_per_condition,sessionLabels_size_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying size from modality using z-scored data
                %[errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(z_data_per_phase_per_grasp,sessionLabels_size_num(trialGraspIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying size from grasp using z-scored data
                %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_grasp,sessionLabels_size_num(trialGraspIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true);
                [errTrain, errTestTmp,~,~,confMat] = classification.LDA_classification_rep(data_per_phase,sessionLabels_size_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                errTest(n_phase,1) =  (1-mean(errTestTmp))*100; % sub 1 for : when separating by condition
    
                % Store the confusion matrix for this phase and session
                confMatAllSessions{n_session, n_phase, 1} = confMat; % sub 1 for : when separating by condition
                %title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ' - ' unTrialType{n_type} ]) %grasp/trialType; works for TrialType but not GraspType?
    
                %data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels,'UniformOutput', false)); % original
                %z_data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),z_scored_data,time_phase_labels, 'UniformOutput', false));
                % [errTrain, errTestTmp] = classification.LDA_classification_rep(z_data_per_phase_per_all,sessionLabels_size_num,'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying size using z-scored data
                % errTest(n_phase) =  (1-mean(errTestTmp))*100;

                % confusion matrices
                % if n_phase ~= 1
                %     %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true); % original
                %     [errTrain, errTestTmp,trueLabels,predictedLabels,cm] = classification.LDA_classification_rep(z_data_per_phase_per_all,sessionLabels_grasp_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %     % Store the confusion matrix
                %     confMatAllSessions{n_session, n_phase} = cm;
                %     title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                % end
    
    
                % % confusion matrices for F30 (Object decoding in one session)
                % object_labels = sessionLabels(129:192,:); 
                % object_data = data_per_phase_per_all(129:192,:);
                % if n_phase ~= 1
                %     [errTrain, errTestTmp] = classification.LDA_classification_rep(object_data,object_labels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                % 
                %     title([unit_region ' - ' phaseNames{n_phase}])
                % end
                
                set(gca, 'FontSize', 12);
    
              
    
            end   
        %end 
        
        % Store the classification accuracy for each session
        all_errTest(n_session, :, 1) = errTest; % sub 1 for : when separating by condition
    
        subplot(4,2,n_session) % code using var names instead of hard 
    
      %  figure();
        bar(errTest)
        %yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100) % when using assigned session labels
        yline((1/numel(uniqueApertureSize))*100) % when classifying size
        %yline((1/numel(uniqueGraspTypes))*100);
        ylim([0 100])
        xticklabels(phaseNames)
    
        if n_session ==3
            %legend(unTrialType, 'Interpreter','none') % sep by cue
            legend(unGraspType, 'Interpreter','none') % sep by grasp
        end 
        ylabel('Classification percentage [%]')
        title(sessions_all{n_session})
        %SessionData = 
    
    end
    
    sgtitle(unit_region)

end

keyboard
%% averaging CMs

for n_phase = 2:numPhases
    % Initialize a figure for the current phase to display multiple confusion matrices (one for each grasp type)
    figure('units','normalized','outerposition',[0 0 .95 0.4]); 
    
    %for n_type = 1:numel(unTrialType) % grasp/trialType
        % Initialize an empty matrix for the summed confusion matrix for the current phase and grasp type
        confMatSize = size(confMatAllSessions{1, n_phase});%, n_type}); % Get size of the first confusion matrix for the phase and type
        confMatSum = zeros(confMatSize);  % Initialize the sum matrix with zeros

        % Sum the confusion matrices across sessions for the current phase and grasp type
        for n_session = 1:numSessions
            cmx = confMatAllSessions{n_session, n_phase};%, n_type};  
            if ~isempty(cmx)
                confMatSum = confMatSum + cmx;
            end
        end
        
        % Calculate the average confusion matrix by dividing the sum by the number of sessions
        confMatAvg = confMatSum / numSessions;
        
        % Normalize each row (class-wise) by dividing by the sum of the row
        confMatAvgNorm = round((confMatAvg ./ sum(confMatAvg, 2)) * 100);
        
        % Create subplot for the current grasp type
        %subplot(1, numel(unTrialType), n_type); % numel(unGrasp/TrialType)
        
        % Display the normalized average confusion matrix for this grasp type
        cmx = confusionchart(confMatAvgNorm);
        cmx.DiagonalColor = '#0072BD';
        cmx.OffDiagonalColor = '#0072BD';
        cmx.FontSize = 20;
        cmx.Title = sprintf('Normalized Average Classification\n%s - %s - %s', unit_region, phaseNames{n_phase});%, unTrialType{n_type}); %grasp/trialType
    %end
end
%% poster
for n_phase = 2%:numPhases
    % Initialize an empty matrix for the summed confusion matrix for the current phase
    confMatSize = size(confMatAllSessions{1, n_phase}); % Get size of the first confusion matrix for the phase
    confMatSum = zeros(confMatSize);  % Initialize the sum matrix with zeros

    % Sum the confusion matrices across sessions for the current phase
    for n_session = 1:numSessions
        cmx = confMatAllSessions{n_session, n_phase};  
        if ~isempty(cmx)
            confMatSum = confMatSum + cmx;
        end
    end

    % Calculate the average confusion matrix by dividing the sum by the number of sessions
    confMatAvg = confMatSum / numSessions;

    % Normalize each row (class-wise) by dividing by the sum of the row
    confMatAvgNorm = round((confMatAvg ./ sum(confMatAvg, 2)) * 100);

    % Display the normalized average confusion matrix
    figure('Position', [100 100 350 350]);
    cmx = confusionchart(confMatAvgNorm);
    cmx.DiagonalColor = '#0072BD';
    cmx.OffDiagonalColor = '#0072BD';
    cmx.FontSize = 25;
    %cmx.Title = sprintf('Normalized Average Classification Across Sessions\n%s - %s', unit_region, phaseNames{n_phase});
end

%%
keyboard

%% ANOVA
[numSessions, numPhases, numModalities] = size(all_errTest);

% Loop through each phase
for n_phase = 1:numPhases
    % Extract data for the current phase
    phaseData = squeeze(all_errTest(:, n_phase, :)); % Shape: [numSessions, numModalities]

    % Reshape data for ANOVA
    data = [];
    modality = [];

    % Loop through each modality
    for n_modality = 1:numModalities
        % Append the classification percentage
        data = [data; phaseData(:, n_modality)];
        % Append modality identifier
        modality = [modality; repmat(n_modality, numSessions, 1)];
    end

    % Convert modality to categorical for ANOVA
    modality = categorical(modality);

    % Create a table with the reshaped data
    dataTable = table(data, modality);

    % Perform ANOVA for the current phase
    [p, tbl, stats] = anova1(dataTable.data, dataTable.modality, 'off');

    % Display ANOVA results
    disp(['Phase ' num2str(n_phase) ':']);
    disp(['p-value = ' num2str(p)]);

    % Perform post-hoc tests if ANOVA is significant
    if p < 0.05
        [c, m, h, gnames] = multcompare(stats, 'CType', 'tukey-kramer');
        disp('Post-hoc test results:');
        disp(c);
    end
end

%% if RM ANOVA is required:
% % ANOVA
% [numSessions, numPhases, numModalities] = size(all_errTest);
% 
% for n_phase = 1:numPhases
%     % Extract data for the current phase
%     phaseData = squeeze(all_errTest(:, n_phase, :)); % Shape: [numSessions, numModalities]
% 
%     % Convert data to table format
%     modality = repmat((1:numModalities)', numSessions, 1);
%     data = phaseData(:);
% 
%     % Create a table with the reshaped data
%     dataTable = table(data, categorical(modality), 'VariableNames', {'Data', 'Modality'});
% 
%     % Specify the within-subject factor (the modalities)
%     WithinDesign = table({'1'; '2'; '3'}, 'VariableNames', {'Modality'});
% 
%     % Fit a repeated measures model
%     rm = fitrm(dataTable, 'Data ~ Modality', 'WithinDesign', WithinDesign);
% 
%     % Perform repeated measures ANOVA
%     ranovatbl = ranova(rm);
% 
%     % Display the ANOVA table
%     disp(['Phase ' num2str(n_phase) ':']);
%     disp(ranovatbl);
% 
%     % Perform post-hoc tests if ANOVA is significant
%     if ranovatbl.pValue(1) < 0.05
%         [c, m, h, gnames] = multcompare(ranova(rm), 'CType', 'tukey-kramer'); %multcompare doesn't work with rm??
%         disp('Post-hoc test results:');
%         disp(c);
%     end
% end

%% plot all classification percentages together with SD
% Calculate the mean and standard deviation across sessions
mean_errTest = squeeze(mean(all_errTest, 1, 'omitnan')); % 2D matrix: [numPhases, numTrialTypes]
std_errTest = squeeze(std(all_errTest, 0, 1, 'omitnan')); % Same dimensions

% Plot the average performance with error bars
figure('units','normalized','outerposition',[0 0 .2 0.25]); % .2 .25 when on desktop, .35 0.4 when on laptop
hold on;

% Create a bar plot with phases on the x-axis
bar_handle = bar(mean_errTest);
%color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % for Conditions
color_info = {[.3632 .2266 .6055],[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]}; % Combinations task (purple at beginning)
%color_info = {[0.2, 0.13, 0.53], [0.067, 0.467, 0.2], [0.53, 0.8, 0.93], [0.53, 0.13, 0.33]}; % grasps: Purple, Green, Light Blue, Dark Pink
colors = cell2mat(color_info');
for n_color = 1:length(bar_handle)
    bar_handle(n_color).FaceColor = colors(n_color,:);
end

% Add error bars
numPhases = size(mean_errTest, 1);
for n_type = 1:numel(unTrialType) %(unTrialType) unGraspType
    % Get the X positions for the current group of bars
    xPositions = bar_handle(n_type).XEndPoints; 
    % Plot the error bars
    errorbar(xPositions, mean_errTest(:, n_type), std_errTest(:, n_type), 'k.', 'LineWidth', 1.5); 
end
% xPositions = bar_handle.XEndPoints; % these two lines are for unseparated data (ie. the entire dataset or Combined dataset)
% errorbar(xPositions, mean_errTest(:, :), std_errTest(:, :), 'k.', 'LineWidth', 1.5);

% Customize plot appearance
chance = 1/(numel(unGraspType))*100; %grasp or object
%chance = 1/(numel(uniqueApertureSize))*100;
yline(chance,'--k','LineWidth',1.5);
ylim([0 100]);
xticks(1:numel(phaseNames));%numPhases);
xticklabels(phaseNames);
ylabel('Classification percentage [%]');
yticks([0 50 100]);
legend(unTrialType, 'Interpreter', 'none');
%legend(unGraspType, 'Interpreter', 'none');
title(unit_region);
set(gca, 'FontSize', 12);
hold off;

%% plot average confusion matrix per phase
% Calculate the average confusion matrix across sessions for each phase
for n_phase = 1:numPhases
    confMatAllSessions{n_phase} = confMatAllSessions{n_phase} / numSessions;
end

for n_phase = 2:numPhases
    figure;
    confusionchart(confMatAllSessions{n_phase}, unique(sessionLabels), 'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
    title(['Average Confusion Matrix - ' unit_region ' - ' phaseNames{n_phase}]);
    set(gca, 'FontSize', 12);
end


%% F30 grasp plot
% figure('Position',[500 500 225 300]);
% region_cue = errTest(2,1:2); % (phase, conditions)
% h = bar(region_cue,'FaceColor','flat');
% h.CData(1,:) = [0.1176, 0.5333, 0.8980];
% h.CData(2,:) = [0.8471, 0.1059, 0.3765];
% %h.CData(3,:) = [1 .7569 .0275];
% yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100,'LineStyle','--','LineWidth',2)
% ylim([0 100])
% xticklabels({'G','G+O'})
% xlim([0.5, 2.5]);
% %ylabel('Classification percentage [%]')
% yticks([0 50 100])
% title([unit_region ' - Cue'])
% set(gca, 'FontSize', 15);
%% F30 cue condition plot
% figure('Position',[500 500 225 300]);
% region_cue = [77, 88, 91]; % value during Cue
% h = bar(region_cue,'FaceColor','flat');
% h.CData(1,:) = [0.7969, 0.4726, 0.6523];
% h.CData(2,:) = [0.3359, 0.7031, 0.9101];
% h.CData(3,:) = [0, 0.6171, 0.4492];
% yline(50,'LineStyle','--','LineWidth',2)
% ylim([0 100])
% xticklabels({'AIP','SMG','M1'})
% xlim([0.5, 3.5]);
% ylabel('Classification percentage [%]')
% yticks([0 50 100])
% % title(unit_region)
% set(gca, 'FontSize', 15);