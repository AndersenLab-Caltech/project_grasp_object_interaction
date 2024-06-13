clc
clear all
close all % closes all the figures

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
taskName = 'GraspObject_4S_Action';
%taskName = 'GraspObject_Shuffled'; % shuffled images
%taskName = 'GraspObject_Varied_Size'; % varied object/aperature sizes 
subject_id = 's4';

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_' taskName spike_sorting_type]);

%%
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

unit_region = 'dlPFC';
brainAreas = Go_data.frPerChannel{6};
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phaseTimeTmp = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

uniqueGraspTypes = unique(Go_data.GraspType);
uniqueCueTypes = unique(Go_data.TrialType);
%uniqueAperatureSize = unique(Data.Aperature_Size); % comment out when not varied sizes

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

    %get idx for Go or NoGo trials
    GoNoGoidx =  logical(cell2mat(Go_data.TrialCue(idxThisSession,:)));
    time_phase_labels = Go_data.time_phase_labels(idxThisSession);
    

    if flagGoTrials
        SessionData = SessionData(GoNoGoidx);
        sessionLabels = sessionLabels(GoNoGoidx);
        time_phase_labels = time_phase_labels(GoNoGoidx);
        trialTypeSession = trialTypeSession(GoNoGoidx);
    else
        SessionData = SessionData(~GoNoGoidx);
        sessionLabels = sessionLabels(~GoNoGoidx);
        time_phase_labels = time_phase_labels(~GoNoGoidx);
        trialTypeSession = trialTypeSession(~GoNoGoidx);
    end
     
    %seperate data according to cue modality 

    unTrialType = unique(Go_data.TrialType);

    % loop through cue modalities 

    for n_type = 1:numel(unTrialType)
        
        % find idx of trial type 
        trialTypeIdx = ismember(trialTypeSession, unTrialType(n_type));
    
        % loop through task phases 
        for n_phase = 1:numPhases
            data_per_phase_per_cue = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData(trialTypeIdx),time_phase_labels(trialTypeIdx), 'UniformOutput', false));
    
            [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_cue,sessionLabels(trialTypeIdx), 'flagErrorMatrix', false, 'PCA_variance', 95, 'flagLeaveOneOut', true);
            errTest(n_phase,n_type) =  (1-mean(errTestTmp))*100;
           
            data_per_phase_per_all = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),SessionData,time_phase_labels, 'UniformOutput', false));
            % % confusion matrices
            % if n_phase ~= 1
            %     [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
            % 
            %     title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
            % end


            % % confusion matrices for F30
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
