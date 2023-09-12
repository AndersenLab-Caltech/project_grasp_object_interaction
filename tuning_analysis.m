clc
clear all
close all



subject_id = 's2';
unit_region = 'SMG';

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject';
flag_4S = true; % true = updated 4S action phase; false = original 2S action phase

if ~flag_4S
    TaskCue = 'GraspObject';
    min_timebin_length = 134; % NOT VALID FOR 20230831    
else
    TaskCue = 'GraspObject_4S_Action';
    min_timebin_length = 174; 
end 
% 4S data
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_GraspObject_4S_Action_unsorted_aligned_thr_-4.5']);

% original 2S data
%Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_GraspObject_unsorted_aligned_thr_-4.5']);

Go_data = Data.Go_data;

flagGoTrials = true;

flagRegressionTuning = false;

flagBinPerBin = true;
multipleComparePhase = true;
flagTunedChannels = true;
flagSaveData = true;

%chose cue type:
taskCuesAll = {'Hand', 'Hand-Object', 'Object'};
sessions_all = unique(Go_data.session_date);
numSessions = numel(sessions_all);
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phase_changes_idx = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phase_changes_idx) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};
color_info = {[.8471 .1059 .3765],[.1176 .5333 .8980],[1 .7569 .0275]};

numUnitsPerSession = zeros(numSessions,1);
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

    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    timePhaseLabels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        SessionData = SessionData(GoNoGoidx);
        sessionLabels = sessionLabels(GoNoGoidx);
        timePhaseLabels = timePhaseLabels(GoNoGoidx);
        trialTypeSession = trialTypeSession(GoNoGoidx);
    else
        SessionData = SessionData(~GoNoGoidx);
        sessionLabels = sessionLabels(~GoNoGoidx);
        timePhaseLabels = timePhaseLabels(~GoNoGoidx);
        trialTypeSession = trialTypeSession(~GoNoGoidx);
    end
     
    %seperate data according to cue modality 

    unTrialType = unique(Go_data.TrialType);
    numUnitsPerSession(n_session) = size(SessionData{1},2);
    % loop through cue modalities 
    for n_type = 1:numel(unTrialType)

        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));

         if flagTunedChannels
            %Compute index of units that are tuned
            if flagRegressionTuning
            
                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin,numTunedChannelsPerCategory,~,~,p_per_phase] ...
                     = classification.getRegressionTunedChannels_paper(SessionData(trialTypeIdx),sessionLabels(trialTypeIdx), ...
                 timePhaseLabels, 'multcompare', multipleComparePhase, 'BinperBinTuning', flagBinPerBin);

                condToTest = arrayfun(@(x) preproc.image2class_simple(x),  unique(sessionLabels), 'UniformOutput', false);

                if nnz(sumBin) ~= 0
                    figure();
                    plot(sumBin);
                end

                tuned_channels_per_graps{n_type,n_session} = numTunedChannelsPerCategory;
            else

                tuned_channels_per_graps{n_type,n_session} = [];
                [tunedCombinedChannels, tunedChannelsPhase, tunedChannelsBin, sumPhase, sumBin]= classification.getTunedChannels(SessionData(trialTypeIdx),sessionLabels(trialTypeIdx), ...
                timePhaseLabels(trialTypeIdx), 'multcompare', multipleComparePhase,'removeITItuning', 'false', 'BinperBinTuning', flagBinPerBin);
                sumBin = sumBin';
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

% for n_cueType = 1:numel(taskCuesAll)
%     dataTmp = fr_sep_cue_type_trial{n_cueType}';
%     % N = size(dataTmp, 1);
%     % sem = std(dataTmp) / sqrt(N);  % standard error of the mean
%     % CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
%     % yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
%     % 
%     data_std = std(dataTmp);
%     data_mean = mean(dataTmp);
% 
%     %%figure(); 
%     %subplot(2,1,1)
%    % plot(1:length(dataTmp), mean(dataTmp), 'LineWidth', 2);
%     %hold on 
%    % plot(1:length(dataTmp), mean(dataTmp) +yCI95 , 'LineWidth', 2);
%     %subplot(2,1,2)
%     % figure();
%     % errorbar(mean(dataTmp), yCI95);
%     % hold on 
%     % err_bar{n_cueType} = plot(1:length(dataTmp), mean(dataTmp), 'LineWidth',2);
%     % figure();
%     %yCI95Shaded = mean(dataTmp) + yCI95;
% 
%     hold on
%     err_bar{n_cueType} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_cueType},'LineWidth',2);
% 
%     ER.mainLine.Color = color_info{n_cueType};
%     ER.patch.FaceColor = color_info{n_cueType};
%     ER.edge(1).Color = color_info{n_cueType};
%     ER.edge(2).Color = color_info{n_cueType};
%       %  err_bar{n_cueType} = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),yCI95(1,:),'lineprops','-b','transparent',true);
% 
%     % Capture the plot handle for this data series
%     %plotHandles = [plotHandles, err_bar{n_cueType}.mainLine];
% 
%    % hold on;
%     %plot(mean(dataTmp), 'LineWidth',2,)
% end

% err_bar = [];
for n_type = 1:numel(unTrialType) 
    if numSessions ~= 1
        tunedUnitsPerType(n_type,:) = sum(cell2mat(tuned_channels_per_phase(n_type,:)'));

    else
        tunedUnitsPerType(n_type,:) = cell2mat(tuned_channels_per_phase(n_type,:)');

    end 
end

% N = size(tunedUnitsPerType, 2);
% sem = std(tunedUnitsPerType(:,1)) / sqrt(N);  % standard error of the mean
% CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
% yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% err_bar(n_type) = [mean(dataTmp)./sum(numUnitsPerSession) yCI95(2,:)./sum(numUnitsPerSession)];

figure(); 
bar(tunedUnitsPerType'./sum(numUnitsPerSession));
hold on;
title(['Tuned Units in ' unit_region]);
xticks(1:numel(phaseNames));
xticklabels(phaseNames);
ylabel('% of Tuned Units');
ylim([0 0.7]);
legend(taskCuesAll, 'Location', 'Best', 'Interpreter', 'none');
set(gca, 'FontSize', 12);
hold off


% err_bar = {};
for n_type = 1:numel(unTrialType)
    tunedUnitsPerTypeBin(n_type,:)  = sum(cell2mat(sum_bin_all(n_type,:)),2);
    % dataTmp = tunedUnitsPerTypeBin(n_type,:)./sum(numUnitsPerSession);
    % N = size(dataTmp, 2);
    % sem = std(dataTmp) / sqrt(N);  % standard error of the mean
    % CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
    % yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    % ER = utile.shadedErrorBar(1:length(dataTmp),mean(dataTmp),yCI95(:,1));
    % err_bar{n_type} = plot(1:length(dataTmp),mean(dataTmp),'Color', color_info{n_type},'LineWidth',2);
    % 
    % ER.mainLine.Color = color_info{n_type};
    % ER.patch.FaceColor = color_info{n_type};
    % ER.edge(1).Color = color_info{n_type};
    % ER.edge(2).Color = color_info{n_type};
end 

figure(); 
plot(tunedUnitsPerTypeBin'./sum(numUnitsPerSession));
hold on
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5);
end
title(['Tuned Units Throughout Trial in ' unit_region]);
xlabel('Time (ms)');
xlim([0 (min_timebin_length + 5)])
ylabel('% of Tuned Units');
ylim([0 0.7]);
legend(taskCuesAll, 'Location', 'Best');
set(gca, 'FontSize', 12);
hold off
