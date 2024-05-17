clc
clear all
close all

subject_id = 's4';
unit_region = 'M1';

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
flag_4S = true; % true = updated 4S action phase; false = original 2S action phase
flag_shuffled = false; % true = shuffled images task

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
%% Regular task
% 4S data
Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' TaskCue spike_sorting_type]);
Go_data = Data.Go_data;

% remove faulty sessions, if any
error_session = {};
if strcmp(subject_id, 's2')
    error_session = {'20231016'};
elseif strcmp(subject_id, 's3')
    error_session = {};
elseif strcmp(subject_id, 's4')
    error_session = {};
end 

if ~isempty(error_session)
    condition = cellfun(@(x) strcmp(x, error_session), Go_data.session_date);
    Go_data = Go_data(~condition,:);
end

flagGoTrials = true; % false = No-Go

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
color_info = {[.1176 .5333 .8980],[.8471 .1059 .3765],[1 .7569 .0275]};

numUnitsPerSession = zeros(numSessions,1);

% Initialize cell arrays to store results
hand_ho_overlap_units_all = cell(numSessions,1);
% hand_only_units_all = cell(numSessions,1);
% ho_only_units_h_all = cell(numSessions,1);
object_ho_overlap_units_all = cell(numSessions,1);
% object_only_units_all = cell(numSessions,1);
% ho_only_units_o_all = cell(numSessions,1);
object_hand_overlap_units_all = cell(numSessions,1);
% object_only_units_h_all = cell(numSessions,1);
% hand_only_units_o_all = cell(numSessions,1);
object_hand_ho_overlap_units_all = cell(numSessions,1);

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
                 timePhaseLabels(trialTypeIdx), 'multcompare', multipleComparePhase, 'BinperBinTuning', flagBinPerBin);

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
    
    % calculating tuning overlap
    % H-HO overlap
    hand_ho_overlap_vector = (tuned_channels_per_phase_vector{1,n_session} == 1) & (tuned_channels_per_phase_vector{2,n_session} == 1); % this tells me the overlap between hand and hand-object units
    hand_ho_overlap_units = sum(hand_ho_overlap_vector, 1);
    hand_ho_overlap_units_all{n_session} = hand_ho_overlap_units;
    
    % hand_only_units = tuned_channels_per_phase{1,n_session} - hand_ho_overlap_units;
    % ho_only_units_h = tuned_channels_per_phase{2,n_session} - hand_ho_overlap_units;
    % hand_only_units_all{n_session} = hand_only_units;
    % ho_only_units_h_all{n_session} = ho_only_units_h;
    
    % O-HO overlap
    object_ho_overlap_vector = (tuned_channels_per_phase_vector{3,n_session} == 1) & (tuned_channels_per_phase_vector{2,n_session} == 1); % this tells me the overlap between object and hand-object units
    object_ho_overlap_units = sum(object_ho_overlap_vector, 1);
    object_ho_overlap_units_all{n_session} = object_ho_overlap_units;
    
    % object_only_units = tuned_channels_per_phase{3,n_session} - object_ho_overlap_units;
    % ho_only_units_o = tuned_channels_per_phase{2,n_session} - object_ho_overlap_units;
    % object_only_units_all{n_session} = object_only_units;
    % ho_only_units_o_all{n_session} = ho_only_units_o;
    
    % O-H overlap
    object_hand_overlap_vector = (tuned_channels_per_phase_vector{3,n_session} == 1) & (tuned_channels_per_phase_vector{1,n_session} == 1); % this tells me the overlap between object and hand units
    object_hand_overlap_units = sum(object_hand_overlap_vector, 1);
    object_hand_overlap_units_all{n_session} = object_hand_overlap_units;

    % object_only_units_h = tuned_channels_per_phase{3,n_session} - object_hand_overlap_units;
    % hand_only_units_o = tuned_channels_per_phase{1,n_session} - object_hand_overlap_units;
    % object_only_units_h_all{n_session} = object_only_units_h;
    % hand_only_units_o_all{n_session} = hand_only_units_o;
    
    % all 3 modalities overlap
    object_hand_ho_overlap_vector = (tuned_channels_per_phase_vector{3,n_session} == 1) & (tuned_channels_per_phase_vector{1,n_session} == 1) & (tuned_channels_per_phase_vector{2,n_session} == 1); % this tells me the overlap between object and hand units
    object_hand_ho_overlap_units = sum(object_hand_ho_overlap_vector, 1);
    object_hand_ho_overlap_units_all{n_session} = object_hand_ho_overlap_units;
end 

hand_ho_overlap_units_all = sum(cell2mat(hand_ho_overlap_units_all'));
% hand_only_units_all = cell2mat(hand_only_units_all');
% ho_only_units_h_all = cell2mat(ho_only_units_h_all');
object_ho_overlap_units_all = sum(cell2mat(object_ho_overlap_units_all'));
% object_only_units_all = cell2mat(object_only_units_all');
% ho_only_units_o_all = cell2mat(ho_only_units_o_all');
object_hand_overlap_units_all = sum(cell2mat(object_hand_overlap_units_all'));
% object_only_units_h_all = cell2mat(object_only_units_h_all');
% hand_only_units_o_all = cell2mat(hand_only_units_o_all');
object_hand_ho_overlap_units_all = sum(cell2mat(object_hand_ho_overlap_units_all'));

% hand_ho_overlap_units_all = sum(hand_ho_overlap_units_all);
% % hand_only_units_all = sum(hand_only_units_all);
% % ho_only_units_h_all = sum(ho_only_units_h_all);
% object_ho_overlap_units_all = sum(object_ho_overlap_units_all);
% % object_only_units_all = sum(object_only_units_all);
% % ho_only_units_o_all = sum(ho_only_units_o_all);
% object_hand_overlap_units_all = sum(object_hand_overlap_units_all);
% % object_only_units_h_all = sum(object_only_units_h_all);
% % hand_only_units_o_all = sum(hand_only_units_o_all);
% object_hand_ho_overlap_units_all = sum(object_hand_ho_overlap_units_all);

hand_total_units = sum(cell2mat(tuned_channels_per_phase(1,:)'));
ho_total_units = sum(cell2mat(tuned_channels_per_phase(2,:)'));
object_total_units = sum(cell2mat(tuned_channels_per_phase(3,:)'));



%% tuned units overlapping

% tuned_channels_per_phase_vector; % 3 (modality) x 5 (sessions)
% % I can compare within sessions how much overlap there is and then average
% % the sessions together to get average overlap
% 
% % H-HO overlap
% hand_ho_overlap_vector = (tuned_channels_per_phase_vector{1,1} == 1) & (tuned_channels_per_phase_vector{2,1} == 1); % this tells me the overlap between hand and hand-object units
% % I can sum and then substract from the total to get the venn diagram
% hand_ho_overlap_units = sum(hand_ho_overlap_vector, 1);
% 
% tuned_channels_per_phase; % total units for each modality
% hand_only_units = tuned_channels_per_phase{1,1} - hand_ho_overlap_units;
% ho_only_units_h = tuned_channels_per_phase{2,1} - hand_ho_overlap_units;
% 
% % next find average by iterating through each session and then finding the
% % mean
% 
% % O-HO overlap
% object_ho_overlap_vector = (tuned_channels_per_phase_vector{3,1} == 1) & (tuned_channels_per_phase_vector{2,1} == 1); % this tells me the overlap between object and hand-object units
% % I can sum and then substract from the total to get the venn diagram
% object_ho_overlap_units = sum(object_ho_overlap_vector, 1);
% 
% tuned_channels_per_phase; % total units for each modality
% object_only_units = tuned_channels_per_phase{3,1} - object_ho_overlap_units;
% ho_only_units_o = tuned_channels_per_phase{2,1} - object_ho_overlap_units;
% 
% % next find average by iterating through each session and then finding the
% % mean
% 
% % O-H overlap
% object_hand_overlap_vector = (tuned_channels_per_phase_vector{3,1} == 1) & (tuned_channels_per_phase_vector{1,1} == 1); % this tells me the overlap between object and hand units
% % I can sum and then substract from the total to get the venn diagram
% object_hand_overlap_units = sum(object_hand_overlap_vector, 1);
% 
% tuned_channels_per_phase; % total units for each modality
% object_only_units_h = tuned_channels_per_phase{3,1} - object_hand_overlap_units;
% hand_only_units_o = tuned_channels_per_phase{1,1} - object_hand_overlap_units;
% 
% % next find average by iterating through each session and then finding the
% % mean
% 
% % all 3 modalities overlap
% object_hand_ho_overlap_vector = (tuned_channels_per_phase_vector{3,1} == 1) & (tuned_channels_per_phase_vector{1,1} == 1) & (tuned_channels_per_phase_vector{2,1} == 1); % this tells me the overlap between object and hand units
% % I can sum and then substract from the total to get the venn diagram
% object_hand_ho_overlap_units = sum(object_hand_ho_overlap_vector, 1);
% 
% tuned_channels_per_phase; % total units for each modality
% object_only_units_h = tuned_channels_per_phase{3,1} - object_hand_overlap_units;
% hand_only_units_o = tuned_channels_per_phase{1,1} - object_hand_overlap_units;

%% example (requires Statistics and Machine Learning Toolbox) => unsure if
% % can handle 3 inputs
% % Sample data (replace with your own data)
% set1 = randi([0, 1], 1, 100); % Binary data for set 1
% set2 = randi([0, 1], 1, 100); % Binary data for set 2
% 
% % Create a logical array for the Venn diagram
% venn_data = [sum(set1 & ~set2), sum(~set1 & set2), sum(set1 & set2)];
% 
% % Create a Venn diagram using vennplot
% figure;
% vennplot(venn_data, 'FaceColor', {'r', 'g', 'b'}, 'FaceAlpha', 0.5);
% 
% % Add labels
% vennlabel({'Set 1', 'Set 2'});
% 
% % Add a title
% title('Venn Diagram');
% 
% % Adjust the display
% axis equal;

%% bar plot of tuned units w/ 95% CIs (work on getting them all on same plot)
%save('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\Workspaces\LinearRegression\s2\s2_GraspObject_2S_unsorted_aligned_thr_-4.5_SMG_Example.mat','sum_bin_all')
%ExampleSMG = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\Workspaces\LinearRegression\s2\s2_GraspObject_2S_unsorted_aligned_thr_-4.5_SMG_Example.mat');
% upload linear regression analysis
%load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\Workspaces\LinearRegression\s2\s2_GraspObject_2S_unsorted_aligned_thr_-4.5_PMV.mat');

phase_yCI95 = [];
phase_tuned_mean = [];
%sessionToInclude = setdiff(1:numSessions,1);

% code for empty/missing session data
rowsToKeep = numUnitsPerSession ~= 0;
numUnitsPerSession = numUnitsPerSession(rowsToKeep);
sessionToInclude = 1:numel(numUnitsPerSession);
colsToKeep = true(1,numSessions);
for n_session = 1:numSessions
    if all(cellfun('isempty',tuned_channels_per_phase(:,n_session)))
        colsToKeep(n_session) = false;
    end
end
tuned_channels_per_phase = tuned_channels_per_phase(:,colsToKeep);

figure('units','normalized','outerposition',[0 0 .38 0.38]);
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
    ylabel('% of Total Units');
    ylim([0 70]);
    sgtitle(['Tuned Units in ' unit_region])
    set(gca, 'FontSize', 12);
end




%% bar plot w/o CIs

for n_type = 1:numel(unTrialType) 
    if numSessions ~= 1
        tunedUnitsPerType(n_type,:) = sum(cell2mat(tuned_channels_per_phase(n_type,:)'));

    else
        tunedUnitsPerType(n_type,:) = cell2mat(tuned_channels_per_phase(n_type,:)');

    end 
end

figure('units','normalized','outerposition',[0 0 0.2 0.35]);
bar((tunedUnitsPerType'./sum(numUnitsPerSession))*100);
%bar((((tunedUnitsPerType')*8)./sum(numUnitsPerSession))*100);
%bar(tunedUnitsPerType');
hold on;
title(['Tuned Units in ' unit_region]);
xticks(1:numel(phaseNames));
xticklabels(phaseNames);
ylabel('% of Total Units');
%ylabel('# of Tuned Units');
ylim([0 70]);
%ylim([0 50]);
legend(taskCuesAll, 'Location', 'Best', 'Interpreter', 'none','FontSize',12);
set(gca, 'FontSize', 12);
hold off

% phaseNames = {'Action'};
% taskCuesAll = {'G', 'G+O'};
% %tuned_channels_per_phase = [23 36; 38 47]; % pulling out cue and action of H & H+O, specific session for proposal
% %tuned_channels_per_phase = [50; 55]; % pulling out action of H & H+O, specific session for proposal
% 
% tunedUnitsPerType = tuned_channels_per_phase;
% figure('Position',[500 500 400 300]);
% h = bar((tunedUnitsPerType'./sum(numUnitsPerSession))*100, 'FaceColor','flat');
% h.CData(1,:) = [0.1176, 0.5333, 0.8980];
% h.CData(2,:) = [0.8471, 0.1059, 0.3765];
% %bar((((tunedUnitsPerType')*8)./sum(numUnitsPerSession))*100);
% %bar(tunedUnitsPerType');
% hold on;
% title(unit_region);
% xticks(1:numel(taskCuesAll));
% xticklabels(taskCuesAll);
% xlim([0.5, 2.5]);
% ylabel('% of Total Units');
% %ylabel('# of Tuned Units');
% ylim([0 100]);
% yticks([0 50 100]);
% %ylim([0 50]);
% %legend(taskCuesAll, 'Location', 'Best', 'Interpreter', 'none','FontSize',12);
% set(gca, 'FontSize', 12);
% hold off

%% line plot w/o CIs

for n_type = 1:numel(unTrialType)
    tunedUnitsPerTypeBin(n_type,:)  = sum(cell2mat(sum_bin_all(n_type,:)),2);

end 

figure('units','normalized','outerposition',[0 0 0.3 0.45]);
plot((tunedUnitsPerTypeBin'./sum(numUnitsPerSession))*100,'LineWidth',2);
%plot((((tunedUnitsPerTypeBin')*8)./sum(numUnitsPerSession))*100,'LineWidth',2);
%plot(tunedUnitsPerTypeBin','LineWidth',2);
hold on
for n_phase = 1:numPhases
    xline(phase_changes(n_phase), 'k--', phaseNames{n_phase}, 'LineWidth', 1.5,'FontSize',12);
end
title(['Tuned Units Throughout Trial in ' unit_region]);
xlabel('Time Bins (50 ms)');
xlim([0 (min_timebin_length + 5)])
%xticks([0 50 100 150]);
ylabel('% of Total Units');
%ylabel('# of Tuned Units');
ylim([0 70]);
%yticks([0 20 40 60]);
%ylim([0 50]);
legend(taskCuesAll, 'Location', 'Best','FontSize',12);
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

figure('units','normalized','outerposition',[0 0 0.3 0.45]);
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

title(['Tuned Units Throughout Trial in ' unit_region]);
xlabel('Time Bins (50 ms)');
xlim([0 (min_timebin_length + 5)]) % 5 chosen as a buffer
ylabel('% of Total Units');
ylim([0 70]);
legend([err_bar{:}], taskCuesAll,'Location', 'Best','Interpreter', 'none','FontSize',12);
set(gca, 'FontSize', 12);
%%
% blub = percentage_tuned.*numUnitsPerSession';
% blub2= mean(blub');
% figure(); plot(zscore(blub2)); hold on; plot(zscore(per_bin_tuned_mean(n_type,:)))