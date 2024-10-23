clc
clear all
close all 

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
%taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
taskName = 'GraspObject_Varied_Size'; % varied object/aperature sizes
%taskName = 'GraspObject_GB_Images'; % for GB
subject_id = 's3';

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName spike_sorting_type]);

%%

if ~strcmp(taskName, 'GraspObject_Varied_Size')
    Go_data = Data.Go_data;
    
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
    brainAreas = Go_data.frPerChannel{6};
    phase_time_idx = Go_data.time_phase_labels{1,1};
    numPhases = numel(unique(phase_time_idx));
    phaseTimeTmp = diff(phase_time_idx);
    phase_changes(1) = 1;
    phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
    phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
    
    %define brain region. 
    % loop to session days 
        % separate data into go and no go
            % seperate data according to cue modality 
                %seperate according to task phases
                    %decode 
    
    
    
    sessions_all = unique(Go_data.session_date);
    numSessions = numel(sessions_all);
    
    flagGoTrials = true; %if true, extract Go trials, if false, extract NoGo trials
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
            SessionData = SessionData(GoNoGoidx);
            sessionLabels = sessionLabels(GoNoGoidx);
            time_phase_labels = time_phase_labels(GoNoGoidx);
            trialTypeSession = trialTypeSession(GoNoGoidx);
            graspTypeSession = graspTypeSession(GoNoGoidx);
        else
            SessionData = SessionData(~GoNoGoidx);
            sessionLabels = sessionLabels(~GoNoGoidx);
            time_phase_labels = time_phase_labels(~GoNoGoidx);
            trialTypeSession = trialTypeSession(~GoNoGoidx);
        end
         
        %seperate data according to cue modality 
    
        unTrialType = unique(Go_data.TrialType);
    
        % loop through cue modalities/sizes
        for n_type = 1:numel(unTrialType) 
            
            % find idx of trial type 
            trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
    
            % code for CMs of grasps and modalities separately
            sessionLabels_modality = trialTypeSession;
            sessionLabels_grasp = graspTypeSession;

            % Convert modality labels ('Hand', 'HandObject', 'Object') to numerical values
            modality_labels = {'Hand', 'Hand_Object', 'Object'};
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
    
            % loop through task phases 
            for n_phase = 1:numPhases
           
                data_per_phase_per_cue = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData(trialTypeIdx),time_phase_labels(trialTypeIdx), 'UniformOutput', false));
                %data_per_phase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false));
    
                % [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % original
                [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels_grasp_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % testing with my method to see if we still get same decoding
                %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase,sessionLabels_modality_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                errTest(n_phase,n_type) =  (1-mean(errTestTmp))*100;
    
                data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false));
                
                % confusion matrices
                % if n_phase ~= 1
                %     [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true); % original
                %     %[errTrain, errTestTmp,~,~,cm] =
                %     %classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels_grasp_num,'flagErrorMatrix', true, 'PCA_variance', 95,'flagLeaveOneOut', true); % CM for grasps/modalities separately (can switch out sessionLabels_grasp/modality)
                % 
                %     % Store the confusion matrix
                %     %confMatAllSessions{n_session, n_phase} = cm;
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
        end 
        
        % Store the classification accuracy for each session
        %all_errTest(n_session, :, :) = errTest;
    
        subplot(4,2,n_session) % code using var names instead of hard 
    
      %  figure();
        bar(errTest)
        yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100)
        ylim([0 100])
        xticklabels(phaseNames)
    
        if n_session ==3
            legend(unTrialType, 'Interpreter','none')
        end 
        ylabel('Classification percentage [%]')
        title(sessions_all{n_session})
        %SessionData = 
    
    end
    
    sgtitle(unit_region)

else 
    Go_data = Data.Go_data;

    % add Aperature Size column
    sizeKeywords = ['Small', 'Medium', 'Large'];
    Go_data.Aperature_Size = cell(height(Go_data),1);
    % Loop through each label and extract the size information
    for i = 1:height(Go_data)
        % Use regular expression to find the size keyword after the last underscore
        tokens = regexp(Go_data.LabelNames{i}, '_(Small|Medium|Large)$', 'tokens');
        
        if ~isempty(tokens)
            % tokens is a cell array; extract the size keyword from it
            Go_data.Aperature_Size{i} = tokens{1}{1};
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
    
    unit_region = 'M1';
    brainAreas = Go_data.frPerChannel{6};
    phase_time_idx = Go_data.time_phase_labels{1,1};
    numPhases = numel(unique(phase_time_idx));
    phaseTimeTmp = diff(phase_time_idx);
    phase_changes(1) = 1;
    phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
    phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
    
    uniqueGraspTypes = unique(Go_data.GraspType);
    uniqueCueTypes = unique(Go_data.TrialType);
    uniqueAperatureSize = unique(Go_data.Aperature_Size); % comment out when not varied sizes
    
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
    all_errTest = NaN(numSessions, numPhases, numel(uniqueCueTypes)); % for sep by Cue Modality
    %all_errTest = NaN(numSessions, numPhases, numel(uniqueGraspTypes)); % for sep by Grasp Types
    
    % Initialize a cell array to store confusion matrices for each phase
    confMatAllSessions = cell(numSessions,1); %cell(numSessions,numPhases); 1 bc only running for Cue rn
    
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
        
        %labels 
        sessionLabels = Go_data.GoLabels(idxThisSession,:);
    
        % grasp labels
        graspTypeSession = Go_data.GraspType(idxThisSession,:);
        
        %trialType
        trialTypeSession = Go_data.TrialType(idxThisSession,:);
    
        %AperatureSize
        aperatureSizeSession = Go_data.Aperature_Size(idxThisSession,:);
    
        %get idx for Go or NoGo trials
        GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
        time_phase_labels = Go_data.time_phase_labels(idxThisSession);
        
    
        if flagGoTrials
            SessionData = SessionData(GoNoGoidx);
            sessionLabels = sessionLabels(GoNoGoidx);
            time_phase_labels = time_phase_labels(GoNoGoidx);
            trialTypeSession = trialTypeSession(GoNoGoidx);
            graspTypeSession = graspTypeSession(GoNoGoidx);
            aperatureSizeSession = aperatureSizeSession(GoNoGoidx);
        else
            SessionData = SessionData(~GoNoGoidx);
            sessionLabels = sessionLabels(~GoNoGoidx);
            time_phase_labels = time_phase_labels(~GoNoGoidx);
            trialTypeSession = trialTypeSession(~GoNoGoidx);
            aperatureSizeSession = aperatureSizeSession(~GoNoGoidx);
        end
         
        %seperate data according to cue modality 
    
        unTrialType = unique(Go_data.TrialType);
    
        %seperate data according to cue modality 
    
        unAperature = unique(Go_data.Aperature_Size);

        % separate data according to grasp

        unGraspType = unique(Go_data.GraspType);
    
        % loop through cue modalities/sizes/grasps
        for n_type = 1:numel(unTrialType) %1:numel(unGraspType) %1:numel(unTrialType) %n_size = 1:numel(unAperature)       %
            
            % find idx of trial type 
            trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
    
            % find idx of size type 
            %trialSizeIdx = ismember(aperatureSizeSession, unAperature(n_size));

            % find idx of grasp type 
            %trialGraspIdx = ismember(graspTypeSession, unGraspType(n_type));
            
            sessionLabels_modality = trialTypeSession;
            sessionLabels_grasp = graspTypeSession;
            sessionLabels_size = aperatureSizeSession;
             
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
           
                data_per_phase_per_cue = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData(trialTypeIdx),time_phase_labels(trialTypeIdx),'UniformOutput', false)); % original (sep by cue data)
                %data_per_phase_per_size = cell2mat(arrayfun(@(x,y)
                    %mean(x{1,1}(y{:}==
                    %n_phase,:),1),SessionData(trialSizeIdx),time_phase_labels(trialSizeIdx),
                    %'UniformOutput', false)); % sep by size data
                %data_per_phase_per_grasp = cell2mat(arrayfun(@(x,y)
                    %mean(x{1,1}(y{:}==
                    %n_phase,:),1),SessionData(trialGraspIdx),time_phase_labels(trialGraspIdx),
                    %'UniformOutput', false)); % sep by grasp data
                %data_per_phase = cell2mat(arrayfun(@(x,y)
                    %mean(x{1,1}(y{:}==
                    %n_phase,:),1),SessionData,time_phase_labels,
                    %'UniformOutput', false)); % all data
    
                % [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % original
                [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels_size_num(trialTypeIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true); % classifying grasp from modality
                %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_grasp,sessionLabels_size_num(trialGraspIdx),'flagErrorMatrix', false, 'PCA_variance', 95,'flagLeaveOneOut', true);
                %[errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase,sessionLabels_modality_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
                errTest(n_phase,n_type) =  (1-mean(errTestTmp))*100;
    
                % Store the confusion matrix for this phase and session
                %confMatAllSessions{n_phase, n_session} = confMat;
    
                data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false));
                % confusion matrices
                % if n_phase ~= 1
                %     [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true); % original
                %     %[errTrain, errTestTmp,~,~,cm] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels_grasp_num, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                %     % Store the confusion matrix
                %     %confMatAllSessions{n_session, n_phase} = cm;
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
        end 
        
        % Store the classification accuracy for each session
        all_errTest(n_session, :, :) = errTest;
    
        subplot(4,2,n_session) % code using var names instead of hard 
    
      %  figure();
        bar(errTest)
        %yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100) % when using assigned session labels
        yline((1/numel(uniqueAperatureSize))*100) % when classifying size
        %yline((1/numel(uniqueGraspTypes))*100);
        ylim([0 100])
        xticklabels(phaseNames)
    
        if n_session ==3
            legend(unTrialType, 'Interpreter','none') % sep by cue
            %legend(unGraspType, 'Interpreter','none') % sep by grasp
        end 
        ylabel('Classification percentage [%]')
        title(sessions_all{n_session})
        %SessionData = 
    
    end
    
    sgtitle(unit_region)

end

%% averaging CMs
% % Initialize an empty matrix for the summed confusion matrix
% confMatSize = size(confMatAllSessions{1, n_phase}); % Get size of the first confusion matrix for phase 2
% confMatSum = zeros(confMatSize);  % Initialize the sum matrix with zeros
% 
% % Sum the confusion matrices across sessions
% for n_session = 1:numSessions
%     cm = confMatAllSessions{n_session, n_phase};  
%     if ~isempty(cm)
%         confMatSum = confMatSum + cm;
%     end
% end
% 
% % Calculate the average confusion matrix by dividing the sum by the number of sessions
% confMatAvg = confMatSum / numSessions;
% 
% % Normalize each row (class-wise) by dividing by the sum of the row
% confMatAvgNorm = round((confMatAvg ./ sum(confMatAvg, 2))*100);
% 
% % Display the normalized average confusion matrix
% figure('Position',[100 100 350 350]);
% cm = confusionchart(confMatAvgNorm);
% cm.DiagonalColor = '#0072BD';
% cm.OffDiagonalColor = '#0072BD';
% cm.FontSize = 20;
% %title(['Normalized Average Classification Across Sessions - ' unit_region ' - ' phaseNames{n_phase}]);


%% ANOVA
% [numSessions, numPhases, numModalities] = size(all_errTest);
% 
% % Loop through each phase
% for n_phase = 1:numPhases
%     % Extract data for the current phase
%     phaseData = squeeze(all_errTest(:, n_phase, :)); % Shape: [numSessions, numModalities]
% 
%     % Reshape data for ANOVA
%     data = [];
%     modality = [];
% 
%     % Loop through each modality
%     for n_modality = 1:numModalities
%         % Append the classification percentage
%         data = [data; phaseData(:, n_modality)];
%         % Append modality identifier
%         modality = [modality; repmat(n_modality, numSessions, 1)];
%     end
% 
%     % Convert modality to categorical for ANOVA
%     modality = categorical(modality);
% 
%     % Create a table with the reshaped data
%     dataTable = table(data, modality);
% 
%     % Perform ANOVA for the current phase
%     [p, tbl, stats] = anova1(dataTable.data, dataTable.modality, 'off');
% 
%     % Display ANOVA results
%     disp(['Phase ' num2str(n_phase) ':']);
%     disp(['p-value = ' num2str(p)]);
% 
%     % Perform post-hoc tests if ANOVA is significant
%     if p < 0.05
%         [c, m, h, gnames] = multcompare(stats, 'CType', 'tukey-kramer');
%         disp('Post-hoc test results:');
%         disp(c);
%     end
% end

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
%%

% plot all classification percentages together with SD
% Calculate the mean and standard deviation across sessions
mean_errTest = squeeze(mean(all_errTest, 1, 'omitnan')); % 2D matrix: [numPhases, numTrialTypes]
std_errTest = squeeze(std(all_errTest, 0, 1, 'omitnan')); % Same dimensions

% Plot the average performance with error bars
figure('units','normalized','outerposition',[0 0 .2 0.25]);
hold on;

% Create a bar plot with phases on the x-axis
bar_handle = bar(mean_errTest); 

% Add error bars
numPhases = size(mean_errTest, 1);
for n_type = 1:numel(unTrialType) %(unGraspType) %(unTrialType)
    % Get the X positions for the current group of bars
    xPositions = bar_handle(n_type).XEndPoints; 
    % Plot the error bars
    errorbar(xPositions, mean_errTest(:, n_type), std_errTest(:, n_type), 'k.', 'LineWidth', 1.5); 
end

% Customize plot appearance
%chance = 1/(numel(uniqueGraspTypes))*100;
chance = 1/(numel(uniqueAperatureSize))*100;
yline(chance,'--k','LineWidth',1.5);
ylim([0 100]);
xticks(1:numPhases);
xticklabels(phaseNames);
ylabel('Classification percentage [%]');
yticks([0 50 100]);
%legend(unTrialType, 'Interpreter', 'none');
title(unit_region);
set(gca, 'FontSize', 12);
hold off;

% plot average confusion matrix per phase
% Calculate the average confusion matrix across sessions for each phase
for n_phase = 1:numPhases
    confMatAllSessions{n_phase} = confMatAllSessions{n_phase} / numSessions;
end

for n_phase = 1:numPhases
    figure;
    confusionchart(confMatAllSessions{n_phase}, unique(sessionLabels), 'RowSummary', 'row-normalized', 'ColumnSummary', 'column-normalized');
    title(['Average Confusion Matrix - ' unit_region ' - ' phaseNames{n_phase}]);
    set(gca, 'FontSize', 12);
end


%% F30 plot
% figure('Position',[500 500 400 300]);
% region_cue = errTest(4,1:2);
% h = bar(region_cue,'FaceColor','flat');
% h.CData(1,:) = [0.1176, 0.5333, 0.8980];
% h.CData(2,:) = [0.8471, 0.1059, 0.3765];
% h.CData(3,:) = [1 .7569 .0275];
% yline(1/numel(unique(sessionLabels(trialTypeIdx)))*100,'LineStyle','--','LineWidth',2)
% ylim([0 100])
% xticklabels({'G','G+O'})
% xlim([0.5, 2.5]);
% ylabel('Classification percentage [%]')
% yticks([0 50 100])
% title(unit_region)
% set(gca, 'FontSize', 15);
