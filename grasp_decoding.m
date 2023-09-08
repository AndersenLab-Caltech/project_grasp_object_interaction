clc
clear all
%close all - closes all the figures

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
taskName = 'GraspObject';
subject_id = 's2';

% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230803_unsorted_aligned_thr_-4.5_GraspObject');
%Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230724_unsorted_aligned_thr_-4.5_GraspObject');
% Data = load('C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\s3\Data\IndividualFiles\GraspObject\unsorted_aligned_thr_-4.5\s3_20230721_unsorted_aligned_thr_-4.5_GraspObject');

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_GraspObject_unsorted_aligned_thr_-4.5']);

%%
Go_data = Data.Go_data;
unit_region = 'PMV';
brainAreas = Go_data.frPerChannel{6};
phase_time_idx = Go_data.time_phase_labels{1,1};
numPhases = numel(unique(phase_time_idx));
phaseTimeTmp = diff(phase_time_idx);
phase_changes(1) = 1;
phase_changes(2:numPhases) = find(phaseTimeTmp) + 1;
phaseNames = {'ITI', 'Cue', 'Delay', 'Action'};

uniqueGraspTypes = unique(Go_data.GraspType);
uniqueCueTypes = unique(Go_data.TrialType);

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
            if n_phase ~= 1
                [errTrain, errTestTmp] = classification.LDA_classification_rep(data_per_phase_per_all,sessionLabels, 'flagErrorMatrix', true, 'PCA_variance', 95, 'flagLeaveOneOut', true);
                
                title([ sessions_all{n_session} ' - ' unit_region ' - ' phaseNames{n_phase} ])
            end

        end 

    end 

    subplot(3,2,n_session) % code using var names instead of hard 

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
