%% GLMs for MCN

%% start with simpliest thing first (unit spiking per trial)

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
len_timebin = 0.05; % FRs are in 50 ms timebins

% Identify unique grasp types and cue types from the data
uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);
graspNumbers = 1:4;  % Corresponding numbers

% Determine outcomes for each unit
for n_brain = 1 % Looking only at SMG
    
    % Extract firing rate data for the current brain area
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    numTrials = size(frData, 3);

    % Initialize a cell array to store spike rate for each channel
    spike_all_data = cell(numChannels, 1); 

    % Initialize a cell array to store outcomes for each channel
    outcome = zeros(numChannels, numTrials); 
    
    % Iterate through each channel in the brain area
    for n_channel = 1:numChannels 
            
        % Extract firing rate data for the channel
        fr_data = squeeze(frData(n_channel,:,:));

        % Store the firing rate data in the cell array
        %fr_all_results{n_channel} = fr_data; 
    
        % Convert FR to spikes
        spike_all_data{n_channel} = fr_data * len_timebin; 

        % Sum spikes across trial
        outcome(n_channel,:) = sum(spike_all_data{n_channel,1},1);
        
    end  

end

% create Design Matrix
% designMatrix should be a 4 (grasps) x 216 (trials) matrix of 0s and 1s
% where 1s denote which grasp is occurring during that trial

designMatrix = zeros(numel(uniqueGraspTypes),numTrials);

for n_grasp = 1:(numel(uniqueGraspTypes))
    % Find indices for the current grasp type
    grasp_ind = ismember(Data.GraspType, uniqueGraspTypes{n_grasp});
        
    % Find indices for GO cue trials (exclude NOGO bc they should not
        % be included in the decoder/model of activity
    go_ind = cell2mat(Data.TrialCue) == 1;

    % Combine grasp type and GO cue indices
    grasp_go_idx = logical(grasp_ind .* go_ind);

    % Set the corresponding elements in the design matrix to 1
    designMatrix(n_grasp, grasp_go_idx) = 1;
end

% fit each unit to GLM and decode 

% Preallocate cell arrays to store the results
b_values = cell(numChannels, 1);
% The coefficients in b_values represent the strength and direction of 
    % the relationship between the predictors (in your design matrix) 
    % and the outcome variable. For a Poisson model, these coefficients
    % indicate how much the log of the expected count changes with a 
    % one-unit change in the predictor.
dev_values = cell(numChannels, 1);
stats_values = cell(numChannels, 1);
% You can check which coefficients are statistically significant by 
    % examining the p-values in stats.p. Small p-values (typically < 0.05) 
    % indicate that the corresponding predictor is significantly 
    % associated with the outcome variable.
    % Use stats.se to calculate confidence intervals for the coefficients, 
    % which can help you assess the precision of your estimates.
predicted_grasps = zeros(numChannels,numTrials);

for n_channel = 1:numChannels
    % Fit the GLM for the current channel
    [b, dev, stats] = glmfit(designMatrix', outcome(n_channel,:)', 'poisson');
    
    % Store the results in the preallocated cell arrays
    b_values{n_channel} = b;
    dev_values{n_channel} = dev;
    stats_values{n_channel} = stats;

    % Decode grasp types 
    for n_trial = 1:numTrials
        % Compute the likelihood for each grasp type for each trial
        likelihoods = zeros(numel(uniqueGraspTypes), 1);
        
        for n_grasp = 1:numel(uniqueGraspTypes)
            % Calculate the linear predictor for the current grasp type
            linear_predictor = b_values{n_channel}(1) + b_values{n_channel}(n_grasp+1);
            
            % Compute the Poisson likelihood (log-likelihood to ensure only positive values)
            lambda = exp(linear_predictor);
            likelihoods(n_grasp) = log(poisspdf(round(outcome(n_channel, n_trial)), lambda));
            
        end
        
        % Normalize to get posterior probabilities (uniform priors)
        posterior_probs = likelihoods / sum(likelihoods);
        
        % Assign the grasp type with the highest posterior probability
        [~, predicted_grasp] = max(posterior_probs);
        predicted_grasps(n_channel, n_trial) = predicted_grasp;
    end
end

% Accuracy of each channel
% Initialize actual grasp array for later
actual_grasps = zeros(1, numTrials); 

% Convert grasp names to numbers
for i = 1:numTrials
    actual_grasps(i) = graspNumbers(strcmp(Data.GraspType{i}, uniqueGraspTypes));
end

% Initialize a vector to store accuracy for each channel
accuracy = zeros(numChannels, 1);

% Loop over each channel to compute the accuracy
for n_channel = 1:numChannels
    % Compare predicted grasps with actual grasps
    correct_predictions = predicted_grasps(n_channel, :) == actual_grasps;
    
    % Calculate the percentage of correct predictions
    accuracy(n_channel) = sum(correct_predictions) / numTrials * 100;
end

min_acc = min(accuracy); % 2.7778
max_acc = max(accuracy); % 25
accuracy_range = [min_acc; max_acc];

%% trying to visualize

for n_channel = 1:numChannels
    figure;
    errorbar(1:length(b_values{n_channel}), b_values{n_channel}, stats_values{n_channel}.se, 'o');
    title(['Channel ' num2str(n_channel) ' Coefficient Estimates']);
    xlabel('Predictors');
    ylabel('Coefficient Estimate');
    grid on;
end

% most coefficients are so small bc baseline firing is high so the other
% coefficients just look like they hang around 0, some dip!!

%% plotting residuals
% Residuals (difference between observed and predicted outcomes) can be 
    % examined to assess the model's goodness-of-fit. You can plot 
    % residuals to check for patterns, which might indicate model 
    % misspecification.

for n_channel = 1:numChannels
    residuals = stats_values{n_channel}.resid;
    fitted_values = exp(b_values{n_channel}(1) + designMatrix' * b_values{n_channel}(2:end));
    figure;
    scatter(fitted_values, residuals);
    title(['Channel ' num2str(n_channel) ' Residual Plot']);
    xlabel('Fitted Values');
    ylabel('Residuals');
    grid on;
end

% I know the model isn't great so probably not much to be gleaned from
% these but ask Uri what the intuition is for reading these plots.

%% Confusion Matrix

% Convert designMatrix to actual_grasps
% [~, actual_grasps] = max(designMatrix, [], 1);  % 1x216 vector

for n_channel = 1:numChannels
    figure;
    
    % Plot the confusion matrix
    confusionchart(actual_grasps, predicted_grasps(n_channel,:));
    title(['Confusion Matrix for Predicted Grasp Types - Channel ' num2str(n_channel)]);
end

% shows as number of trials
% as expected, not a stellar decoder

%% Confusion matrix with percentages

for n_channel = 1:numChannels
    % Compute the confusion matrix
    cm = confusionmat(actual_grasps, predicted_grasps(n_channel,:));
    
    % Normalize the confusion matrix to get percentages
    cm_percentage = cm / sum(cm(:)) * 100;  % Normalize by total number of trials and convert to percentage
    
    % Create a new figure for each confusion matrix
    figure;
    
    % Define grasp type labels
    grasp_labels = arrayfun(@num2str, 1:numel(uniqueGraspTypes), 'UniformOutput', false);
    
    % Plot the confusion matrix as a heatmap with percentages
    heatmap(grasp_labels, grasp_labels, cm_percentage, ...
            'Title', ['Confusion Matrix for Predicted Grasp Types - Channel ' num2str(n_channel)], ...
            'XLabel', 'Predicted', 'YLabel', 'Actual', ...
            'ColorbarVisible', 'on');

end

%% AIC vs BIC model comparison 
% You can compute the Akaike Information Criterion (AIC) and Bayesian 
    % Information Criterion (BIC) from the deviance to compare different 
    % models. Lower values indicate a better trade-off between model fit 
    % and complexity.

for n_channel = 1:numChannels
    aic = 2 * length(b_values{n_channel}) + dev_values{n_channel};
    bic = log(numTrials) * length(b_values{n_channel}) + dev_values{n_channel};
    fprintf('Channel %d: AIC = %.2f, BIC = %.2f\n', n_channel, aic, bic);
end

%% Prediction over trials

for n_channel = 1:numChannels
    predicted_spikes = exp(b_values{n_channel}(1) + designMatrix' * b_values{n_channel}(2:end));
    actual_spikes = outcome(n_channel, :);
    figure;
    plot(actual_spikes, 'o');
    hold on;
    plot(predicted_spikes, 'x');
    title(['Channel ' num2str(n_channel) ' Actual vs Predicted Spikes']);
    xlabel('Trial');
    ylabel('Spike Count');
    legend('Actual', 'Predicted');
    grid on;
    hold off;
end

%% scale up to include timebins in predictor and parameters

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
len_timebin = 0.05; % FRs are in 50 ms timebins

% Identify unique grasp types and cue types from the data
uniqueGraspTypes = unique(Data.GraspType);
uniqueCueTypes = unique(Data.TrialType);

% Determine outcomes for each unit
for n_brain = 1 % Looking only at SMG
    
    % Extract firing rate data for the current brain area
    frData = Data.frPerChannel{n_brain};
    numChannels = size(frData, 1);
    numTrials = size(frData, 3);

    % Initialize a cell array to store spike rate for each channel
    spike_all_data = cell(numChannels, 1); % replace first 1 with numChannels when scaling up

    % Initialize a cell array to store outcomes for each channel
    outcome = zeros(numChannels, (numTrials*numel(phase_time_idx))); 
    
    % Iterate through each channel in the brain area
    for n_channel = 1:numChannels 
            
        % Extract firing rate data for the channel
        fr_data = squeeze(frData(n_channel,:,:));
    
        % Convert FR to spikes
        spike_all_data{n_channel} = fr_data * len_timebin; % I think this 
            % is the term we use for scaling up to include time? I think I
            % need to concatenate each cell on top of each other, such that
            % the columns get moved under the first row

        % Concatenate spikes per trial per timebin
        outcome(n_channel,:) = reshape(spike_all_data{n_channel}, 1, []); % each column is a unit (y)
        
    end  

end

% create Design Matrix
% designMatrix should be a 178 (4 grasps + 174 timebins) x 37584 (216 trials x 174 timebins) 
% matrix of 0s and 1s where 1s denote which grasp is occurring and at what
% timebin => not sure how to index for this yet

% Define the number of timebins and initialize the design matrix
numTimebins = numel(phase_time_idx);  
numGrasps = numel(uniqueGraspTypes);  
designMatrix = zeros(numGrasps + numTimebins, numTrials * numTimebins);

% Loop over each trial to fill the design matrix
for n_trial = 1:numTrials
    % Find the grasp type for the current trial
    grasp_type = find(ismember(uniqueGraspTypes, Data.GraspType{n_trial}));
    
    % Set the grasp indicator (rows 1 to numGrasps) for the current trial
    start_idx = (n_trial - 1) * numTimebins + 1;
    end_idx = n_trial * numTimebins;
    designMatrix(grasp_type, start_idx:end_idx) = 1;
    
    % Set the timebin indicators (rows numGrasps+1 to numGrasps+numTimebins)
    for n_timebin = 1:numTimebins
        timebin_idx = numGrasps + n_timebin;
        designMatrix(timebin_idx, start_idx + n_timebin - 1) = 1;
    end
end

designMatrix = designMatrix';

%% Create new matrix to be able to decode across time

I_timebin = designMatrix(:,5:end);
newDM = [];

for n_grasp = 1:numGrasps
    newDM = [newDM, (designMatrix(:,n_grasp) * ones(1,numTimebins)) .* I_timebin];
end


%% Fit GLM and Decode Grasp Types

% Preallocate cell arrays to store the results
b_values = cell(numChannels, 1);
dev_values = cell(numChannels, 1);
stats_values = cell(numChannels, 1);
predicted_grasps = zeros(numChannels, numTrials);

for n_channel = 1:numChannels
    % Fit the GLM for the current channel using the full design matrix
    [b, dev, stats] = glmfit(newDM, outcome(n_channel,:)', 'poisson', 'constant','off'); 

    % I get a 62x1 matrix with each cell containing a 696x1 double, I think
    % I'm doing something wrong here? I feel like each unit should contain
    % the info across all timebins for all trials? but this might be b
    % values for each grasp across timebins (4 x 174 = 696)?
    % yes, what I get from my GLM is estimated activity of each unit during
    % each grasp for any trial
    
    % Store the results in the preallocated cell arrays
    b_values{n_channel} = b;
    dev_values{n_channel} = dev;
    stats_values{n_channel} = stats;

   
end

%% Decode

%plot(1:174,exp(b_values{1}(1:174)),1:174,exp(b_values{1}(175:(174+174))))
% Initialize array
DM = zeros(numTimebins,(numGrasps*numTimebins));
logL = zeros(numTimebins,numGrasps,numChannels,numTrials);

% Loop through each channel
for n_channel = 1:numChannels
    %Loop through each grasp
    for n_grasp = 1:numGrasps
        
        % Index into newDM to obtain all trial indicators for that specific
        % grasp
        timebin_start_idx = ((n_grasp - 1) * numTimebins) + 1;
        timebin_end_idx = ((n_grasp - 1) * numTimebins) + numTimebins;
        DM = newDM(:,(timebin_start_idx:timebin_end_idx));

        for n_trial = 1:numTrials
            % want to index into outcome with timebins per trial? which trial? I think I need to loop through trial?
            trial_start_idx = (n_trial - 1) * numTimebins + 1;
            trial_end_idx = n_trial * numTimebins;
            % obtain linear predictor and log likelihoods
            b_val = b_values{n_channel}(timebin_start_idx:timebin_end_idx)';
            lambda = glmval(b_val',eye(174),'log','constant','off');
            logL(:,n_grasp,n_channel,n_trial) = log(poisspdf(round(outcome(n_channel,trial_start_idx:trial_end_idx))',lambda));
        end
    end 
end

%% Visualizing each channel's contribution
% I want to go through channel by channel and look at plot((1:174,exp(b_values{1}(1:174)),1:174,exp(b_values{1}(175:(174+174))))
% this is an example looking at the b-values for neuron 1 across time for
% grasp 1 and 2 (can mult b_vals by 20 to get FR/s?) => I want to do this 
% for each neuron, all 4 grasps, on a subplot that also plots the decode 
% that each unit is doing => the two should make sense together

% so I just need to create loops through every channel and trial and I can
% subplot so they're right next to each other

for n_channel = 1 %1:numChannels
    figure()
    %sgtitle([brainAreas{n_brain} ' - Channel ' num2str(n_channel)]);
    % Define specific colors for each grasp
    grasp_colors = [0 0.4470 0.7410;   % Blue for Grasp 1
                    0.8500 0.3250 0.0980; % Red for Grasp 2
                    0.9290 0.6940 0.1250; % Yellow for Grasp 3
                    0.4940 0.1840 0.5560];% Purple for Grasp 4
    for n_grasp = 1:numGrasps

        % grasp indices
        grasp1_idx = 1:numTimebins;
        grasp2_idx = (numTimebins + 1):(numTimebins*2);
        grasp3_idx = (numTimebins*2 + 1):(numTimebins*3);
        grasp4_idx = (numTimebins*3 + 1):(numTimebins*4);

        % Time index
        time_idx = 1:numTimebins;
        
        % convert b_values to FR
        est_FR = exp(b_values{n_channel})*20; % mult by 20 bc 50ms timebins

        % Plot figure
        subplot(2,1,1);
        hold on;
        plot(time_idx, est_FR(grasp1_idx), 'Color', grasp_colors(1,:));
        plot(time_idx, est_FR(grasp2_idx), 'Color', grasp_colors(2,:));
        plot(time_idx, est_FR(grasp3_idx), 'Color', grasp_colors(3,:));
        plot(time_idx, est_FR(grasp4_idx), 'Color', grasp_colors(4,:));
        hold on;
        for n_phase = 1:numPhases
            xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
        end

        legend(uniqueGraspTypes);
        title(['Activity Across Trial - Channel ' num2str(n_channel)])
        xlabel("Timebins (50 ms)")
        ylabel('Firing Rate (spikes/sec)')

        for n_trial = 1%:numTrials
            % Calculate log likelihood probability of each grasp for each
            % channel and trial 
            lp1 =(sum(cumsum(logL(:,1,n_channel,n_trial),1),3));
            lp2 =(sum(cumsum(logL(:,2,n_channel,n_trial),1),3));
            lp3 =(sum(cumsum(logL(:,3,n_channel,n_trial),1),3));
            lp4 =(sum(cumsum(logL(:,4,n_channel,n_trial),1),3));

            % Compute the Constant using max so that at least one grasp
            % will prove the winning (bc values are so small)
            C = max([lp1 lp2 lp3 lp4]')';

            % Calculate probabilities for each grasp
            p1 = exp(lp1 - C)./(exp(lp1 - C) + exp(lp2-C) + exp(lp3-C) +exp(lp4-C));
            p2 = exp(lp2 - C)./(exp(lp1 - C) + exp(lp2-C) + exp(lp3-C) +exp(lp4-C));
            p3 = exp(lp3 - C)./(exp(lp1 - C) + exp(lp2-C) + exp(lp3-C) +exp(lp4-C));
            p4 = exp(lp4 - C)./(exp(lp1 - C) + exp(lp2-C) + exp(lp3-C) +exp(lp4-C));

            % Grasp prediction
            % Assign the grasp type with the highest probability
            [~, predicted_grasp] = max([p1(174) p2(174) p3(174) p4(174)]); % chose the probability 
                % at the last timebin so that it accounted for all the
                % information during the trial?
            predicted_grasps(n_channel, n_trial) = predicted_grasp;

            % Plot probabilities
            subplot(2,1,2)
            cla; % clesr previous subplot so not all plotted on same fig
            plot(time_idx,p1,time_idx,p2,time_idx,p3,time_idx,p4)

            legend(uniqueGraspTypes);
            title(['Decoding Across Trial ' num2str(n_trial)])
            xlabel("Timebins (50 ms)")
            ylabel('Probability')

            % pause;
        end
    end
end

%% Leave one out (LOO) cross validation

for n_trial = 1%:numTrials
    % Indices for the current trial
    trial_start_idx = (n_trial - 1) * numTimebins + 1;
    trial_end_idx = n_trial * numTimebins;

    % Split the data into training and testing
    train_idx = setdiff(1:(numTrials * numTimebins), trial_start_idx:trial_end_idx);
    test_idx = trial_start_idx:trial_end_idx;

    % Training data
    train_outcome = outcome(:, train_idx);
    train_designMatrix = newDM(train_idx, :);

    % Testing data
    test_outcome = outcome(:, test_idx);
    test_designMatrix = newDM(test_idx, :);

    % Fit the GLM using the training data
    for n_channel = 1:numChannels
        [b, ~, ~] = glmfit(train_designMatrix, train_outcome(n_channel, :)', 'poisson', 'constant', 'off');

        % Predict the grasp on the test trial
        lambda = glmval(b, test_designMatrix, 'log', 'constant', 'off');
        logL_test = log(poisspdf(round(test_outcome(n_channel, :)'), lambda));

        % Since logL_test corresponds to a single trial's timebins, the grasp probabilities
        % need to be calculated across the timebins without further indexing into grasps.
        logL_grasps = zeros(numGrasps, 1);
        for n_grasp = 1:numGrasps           
            % Sum the log likelihoods across all timebins for this grasp
            logL_grasps(n_grasp) = sum(logL_test);
        end

        % Determine which grasp has the highest probability
        [~, predicted_grasp] = max(logL_grasps);
        predicted_grasps_LOO(n_channel, n_trial) = predicted_grasp;

        % Store the actual grasp for comparison
        actual_grasps_LOO(n_channel, n_trial) = find(ismember(uniqueGraspTypes, Data.GraspType{n_trial}));
    end
end


%% Accuracy of each channel

graspNumbers = 1:4;  % Corresponding numbers

% Initialize the numeric array
actual_grasps = zeros(1, numTrials);  % 1x216 array

% Convert grasp names to numbers
for i = 1:numTrials
    actual_grasps(i) = graspNumbers(strcmp(Data.GraspType{i}, uniqueGraspTypes));
end

% Initialize a vector to store accuracy for each channel
accuracy = zeros(numChannels, 1);

% Loop over each channel to compute the accuracy
for n_channel = 1:numChannels
    % Compare predicted grasps with actual grasps
    correct_predictions = predicted_grasps_LOO(n_channel, :) == actual_grasps;
    
    % Calculate the percentage of correct predictions
    accuracy(n_channel) = sum(correct_predictions) / numTrials * 100;
end

min_acc = min(accuracy); % 27.3148
max_acc = max(accuracy); % 60.6481
mean_acc = mean(accuracy); % 43.6380 for all grasp types

%% Grasp Type accuracies => much higher accuracies when broken up by grasp type,
    % maybe indicating tuning to specfic grasps. 

% Initialize a matrix to store the accuracy for each channel and each grasp type
accuracy_per_grasp = zeros(numChannels, numel(uniqueGraspTypes));

% Loop through each channel
for n_channel = 1:numChannels
    % Loop through each grasp type
    for n_grasp = 1:numel(uniqueGraspTypes)
        % Find the indices of trials that correspond to the current grasp type
        grasp_trials = actual_grasps == n_grasp;
        
        % Compare the predicted grasps to actual grasps for these trials
        correct_predictions = predicted_grasps(n_channel, grasp_trials) == actual_grasps(grasp_trials);
        
        % Calculate the accuracy for the current grasp type and channel
        accuracy_per_grasp(n_channel, n_grasp) = sum(correct_predictions) / sum(grasp_trials) * 100;
    end
end

%% Confusion Matrix

% Convert designMatrix to actual_grasps
% [~, actual_grasps] = max(designMatrix, [], 1);  % 1x216 vector

for n_channel = 1:numChannels
    figure;
    
    % Plot the confusion matrix
    confusionchart(actual_grasps, predicted_grasps(n_channel,:));
    title(['Confusion Matrix for Predicted Grasp Types - Channel ' num2str(n_channel)]);
end

% some units are just picking the same grasp type every time. maybe this
% makes sense bc we didn't average over trials before fitting?



%% everything after this is WRONG; take average FR across grasps
% Initialize an array to store the final mean value for each channel and grasp type
mean_fr_grasp_results = zeros(1, numel(uniqueGraspTypes)); %replace 1 w/ numChannels when scaling up

% Iterate over each channel and grasp type
for n_channel = 3 %1:numChannels
    for n_grasp = 1:numel(uniqueGraspTypes)
        
        % Extract the 174x48 matrix for the current channel and grasp type
        fr_grasp = fr_grasp_results{1, n_grasp}; % substitute 1 for n_channel when scaling
        
        % Calculate the mean of each column (resulting in a 1x48 row vector)
        mean_cols = mean(fr_grasp, 1);
        
        % Calculate the mean across the row vector to get a single value
        mean_fr_grasp = mean(mean_cols);
        
        % Store the single mean value
        mean_fr_grasp_results(1, n_grasp) = mean_fr_grasp; % substitute 1 for n_channel when scaling
        
    end
end

%% fit GLM

% y (outcome) must be a numeric array so can't use uniqueGraspTypes
% assign each a numerical value (1 = L, 2 = MW, 3 = PP, 4 = S3F)

% y should be number of spikes so divide fr by time (50ms) => y is sum of
% spikes across each trial

outcome1 = (1:numel(uniqueGraspTypes))';

designMatrix1 = mean_fr_grasp_results';
[x1,dev1,stats1] = glmfit(designMatrix1,outcome1,'poisson'); %constant = off

%% add time back in

fr_grasp_results; % cell of 174 x 48 for each grasp

% take average FR across trials per time per grasp
% Initialize an array to store the final mean value for each channel and grasp type
mean_fr_time_grasp_results = cell(1, numel(uniqueGraspTypes)); %replace 1 w/ numChannels when scaling up

% Iterate over each channel and grasp type
for n_channel = 3 %1:numChannels
    for n_grasp = 1:numel(uniqueGraspTypes)
        
        % Extract the 174x48 matrix for the current channel and grasp type
        fr_grasp = fr_grasp_results{1, n_grasp}; % substitute 1 for n_channel when scaling
        
        % Calculate the mean of each column (resulting in a 174x1 col vector)
        mean_fr_time_grasp = mean(fr_grasp, 2);
        
        % Trial averaged activity per grasp across time
        mean_fr_time_grasp_results{1, n_grasp} = mean_fr_time_grasp; % substitute 1 for n_channel when scaling
        
    end
end

% I think I want to add in phases for each cell and then concatenate such
% that I have 696 (174 x 4) x 2 (avg FR at that timepoint and trial phase),
% this is not just adding in time parameter though, this is adding in time
% parameter. I think the phase parameter should be added separately?

%phase_time_idx; % phases across trial (1 = ITI, 2 = Cue, 3 = Delay, 4 = Action)

time_idx = (1:numel(phase_time_idx))';

% Initialize an array to store the final mean value for each channel and grasp type
mean_fr_time_results = cell(1, numel(uniqueGraspTypes)); %replace 1 w/ numChannels when scaling up

% Loop through each cell in the 1x4 cell array
for n_grasp = 1:numel(uniqueGraspTypes)
    % Extract the 174x1 double from the current cell
    fr_data = mean_fr_time_grasp_results{n_grasp};
    
    % Concatenate the 174x1 fr_data with the 174x1 time_idx to form a 174x2 matrix
    updated_data = [fr_data, time_idx];
    
    % Store the updated 174x2 matrix back into the cell
    mean_fr_time_results{n_grasp} = updated_data;
end

% Initialize an empty array to store the concatenated results
designMatrix2 = [];

% Loop through each cell in mean_fr_time_grasp_results
for n_grasp = 1:numel(uniqueGraspTypes)
    % Concatenate each 174x2 matrix from the cell array vertically
    designMatrix2 = [designMatrix2; mean_fr_time_results{n_grasp}];
end

% designMatrix is now trial averaged FR across each timebin, one grasp on
% top of the other (696) x 2 (FR, time bin labels)

% Create the outcome variable
outcome2 = [ones(numel(phase_time_idx), 1); 
           2 * ones(numel(phase_time_idx), 1); 
           3 * ones(numel(phase_time_idx), 1); 
           4 * ones(numel(phase_time_idx), 1)];

% fit glm
[x2,dev2,stats2] = glmfit(designMatrix2,outcome2,'poisson');

%% add in time phases

% would I have to take the average FR of the phase? (I'll try it both ways)
phase_time_idx; % phases across trial (1 = ITI, 2 = Cue, 3 = Delay, 4 = Action)

% Initialize an array to store the final mean value for each channel and grasp type
wholeTrial_fr_phase_results = cell(1, numel(uniqueGraspTypes)); %replace 1 w/ numChannels when scaling up

% Loop through each cell in the 1x4 cell array
for n_grasp = 1:numel(uniqueGraspTypes)
    % Extract the 174x1 double from the current cell
    fr_data = mean_fr_time_grasp_results{n_grasp};
    
    % Concatenate the 174x1 fr_data with the 174x1 time_idx to form a 174x2 matrix
    updated_data = [fr_data, phase_time_idx];
    
    % Store the updated 174x2 matrix back into the cell
    wholeTrial_fr_phase_results{n_grasp} = updated_data;
end

% Initialize an empty array to store the concatenated results
designMatrix3 = [];

% Loop through each cell in mean_fr_time_grasp_results
for n_grasp = 1:numel(uniqueGraspTypes)
    % Concatenate each 174x2 matrix from the cell array vertically
    designMatrix3 = [designMatrix3; wholeTrial_fr_phase_results{n_grasp}];
end

% designMatrix2 is now trial averaged FR across each timebin, one grasp on
% top of the other (696) x 2 (FR, phase labels)

% fit glm
[x3,dev3,stats3] = glmfit(designMatrix3,outcome2,'poisson');

% very similar dev2 and dev3

%% okay now take the average across entire phase bc maybe that is necessary?

% Initialize the result array (4x4 double if there are 4 phases and 4 grasp types)
mean_fr_phase_grasp_results = cell(1, 4); % Replace 1 with numChannels if you're using more channels

% Iterate over each grasp type
for n_grasp = 1:numel(uniqueGraspTypes)
    
    % Extract the 174x1 double for the current grasp type
    fr_time = mean_fr_time_grasp_results{1, n_grasp}; % substitute 1 for n_channel when scaling
    
    % Initialize a 4x1 vector to store the phase averages
    phase_averages = zeros(numel(unique(phase_time_idx)), 1);
    
    % Loop over each phase
    for n_phase = 1:numel(unique(phase_time_idx))
        % Find the indices corresponding to the current phase
        phase_idx = (phase_time_idx == n_phase);
        
        % Calculate the mean firing rate for the current phase
        phase_averages(n_phase) = mean(fr_time(phase_idx));
    end
    
    % Store the result back in the cell array
    mean_fr_phase_grasp_results{1, n_grasp} = phase_averages; % substitute 1 for n_channel when scaling
end

% Initialize an array to store the final mean value for each channel and grasp type
mean_fr_phase_results = cell(1, numel(uniqueGraspTypes)); %replace 1 w/ numChannels when scaling up

% Loop through each cell in the 1x4 cell array
for n_grasp = 1:numel(uniqueGraspTypes)
    % Extract the 4x1 double from the current cell
    mean_phase_data = mean_fr_phase_grasp_results{n_grasp};
    
    % Concatenate the 174x1 fr_data with the 174x1 time_idx to form a 174x2 matrix
    updated_data = [mean_phase_data, unique(phase_time_idx)];
    
    % Store the updated 174x2 matrix back into the cell
    mean_fr_phase_results{n_grasp} = updated_data;
end

% Initialize an empty array to store the concatenated results
designMatrix4 = [];

% Loop through each cell in mean_fr_time_grasp_results
for n_grasp = 1:numel(uniqueGraspTypes)
    % Concatenate each 4x2 matrix from the cell array vertically
    designMatrix4 = [designMatrix4; mean_fr_phase_results{n_grasp}];
end

% designMatrix4 is now trial averaged FR across each phase, one grasp on
% top of the other (16) x 2 (FR, phase labels)

% Create the outcome variable
outcome4 = [ones(numel(unique(phase_time_idx)), 1); 
           2 * ones(numel(unique(phase_time_idx)), 1); 
           3 * ones(numel(unique(phase_time_idx)), 1); 
           4 * ones(numel(unique(phase_time_idx)), 1)];

% fit glm
[x4,dev4,stats4] = glmfit(designMatrix4,outcome4,'poisson');

% much smaller dev4 but visually seems like this unit is not contributing
% much