clc
clear all
%close all



subject_id = 's3';
unit_region = 'S1';

spike_sorting_type = '_unsorted_aligned_thr_-4.5';
taskName = 'GraspObject';

Data = load(['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\Table_' subject_id '_GraspObject_unsorted_aligned_thr_-4.5']);

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




for n_type = 1:numel(unTrialType)
    tunedUnitsPerType(n_type,:) = sum(cell2mat(tuned_channels_per_phase(n_type,:)'));
end 

figure(); bar(tunedUnitsPerType')



for n_type = 1:numel(unTrialType)
    tunedUnitsPerTypeBin(n_type,:)  = sum(cell2mat(sum_bin_all(n_type,:)),2);

end 

figure(); 
plot(tunedUnitsPerTypeBin'./sum(numUnitsPerSession));
legend(taskCuesAll)
