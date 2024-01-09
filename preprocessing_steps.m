% Cleaned up preprocessing file for Speech paper. For S2 and S3

clc 
clear all
close all

%subject_id = 's2';  % s2 or p3 or n1
subject_id = 's3';  % s2 or p3 or n1

subject = hst.Subject(subject_id);
flag_dPCA = false; 
flag_4S = true; 
flag_shuffled = true; % true for shuffled images

if strcmp(subject_id, 's2')
    %session_dates = {'20230831','20230907'};
    session_dates = {'20231201'};
elseif strcmp(subject_id, 's3')
    session_dates = {'20231207'};
else 
    error('unknown subject')
end 

save_data = true; 
flagRemoveTrials = true; 

spike_sorting_type = 'unsorted_aligned_thr_-4.5';
%spike_sorting_type = 'unsorted_aligned_noratefilt_4.5';
%spike_sorting_type = 'unsorted_aligned_noratefilt';
%spike_sorting_type = 'sorting'; % I did not rethreshold the session before spike sorting... idk if that will work? 


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



save_data_pathway = ['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\IndividualFiles\' TaskCue '\' spike_sorting_type];

%create the folder if it does not exist in the path yet
if ~exist(save_data_pathway)
    mkdir(save_data_pathway)
end


flagRemoveErrorTrials = true; %remove trials were errors occured
session_date_idx = 1:length(session_dates);

classification_result = [];
classification_result_com = [];
last_session_date = 0;

for n_session = session_date_idx 

    disp('-----------------------------------------------------------')
    disp(['Session Day        : ' session_dates{n_session}])
    disp(['Spike sorting type : ' spike_sorting_type])
    disp(['Cue type           : ' TaskCue])
    disp('-----------------------------------------------------------')
    disp(['Current save mode is : ' num2str(save_data)]); 

    session = hst.Session(session_dates{n_session}, subject);
    taskfiles = session.getTaskFiles('GraspObject');
    %pwd = 'C:\Users\macthurston\Documents\GitHub\project_grasp_object_interaction'; % location on this computer

    %Extract correct datasets from excel files 
    filename = [subject_id '_good_trials_' TaskCue '.xlsx'];
    
    [numb,txt,raw] = xlsread(fullfile(pwd,'ExcelFiles',filename));

    session_date = str2double(session_dates{n_session});
    good_datablocks = numb(session_date==numb(:,1), 2:end);
        
    good_datablocks(isnan(good_datablocks)) = [];
    
    if (length(good_datablocks) <1)
        error ('Add to excel sheet, No good dataset present, skip or check for problem')
    elseif good_datablocks == 0
        error('No good dataset present for this day, has been checked, skip')  
    elseif any(ismember(good_datablocks, 1000))
        warning('There is a problem with one of the datasets, ask Spencer? Removing it for now');
        good_datablocks(ismember(good_datablocks,1000)) = []; 
    else
        disp('Everything seems ok');
    end

    
    disp(['Good_datasets : ' num2str(good_datablocks)]); 
       
 
    %%
    %to verify length of time bins
    time_bin_size = zeros(1,length(good_datablocks)); 
    
    %Either sorting (sorted spikes), noisy (includes noise units), unsorted
    %(unsorted data), smoothed (halfkerner = 0 and causal = 0, how much the
    %neuronal data is smoothed). 

    individual_runs = cell(1,length(good_datablocks));
    featdef = cell(1,length(good_datablocks));
    
    %predeclare variables
    TrialNumber = [];
    LabelNames = {};
    TrialType = {};
    TrialCue = [];
    GoLabels = [];
    session_date = [];
    cueType = {};
    time_phase_labels = [];
    time_trial = [];

    for n_dataset = 1:length(good_datablocks) 
    
            %load task 
            task = hst.Task(taskfiles{good_datablocks(n_dataset)});  
            
            %idx of trials to include                                  
            if flagRemoveErrorTrials
                data_subset = setdiff(1:task.numTrials,  preproc.errorTrialsPerTaskfile(task, subject_id));                
            else 
                data_subset = 1:task.numTrials; 
            end 
            
            numTrials = length(data_subset);
            
            %extract and process neural data
            individual_runs{n_dataset} = preproc.get_neural_data_paper(task,...
                'spikes',spike_sorting_type,...
                'ratefilt',true,...
                'trials',data_subset,...
                'min_timebin_length', min_timebin_length);
                       
                       
            %save relevant variables for analysis 
            time_bin_size(n_dataset) = size(individual_runs{n_dataset}.fr_adapted,1);
                
            individual_runs{n_dataset}.task.hDebug= []; %remove saving debugger, created errors later when loading the file
            LabelNames_ind = {task.trialparams(:).action}';
            individual_runs{n_dataset}.Labels = LabelNames_ind;
            LabelNames = vertcat(LabelNames, LabelNames_ind(data_subset));
            
            GraspType = cellfun(@(x) strsplit(x, '_'), LabelNames,'UniformOutput',false);
            GraspType = cellfun(@(x) x{1,1}, GraspType,'UniformOutput',false);

            TrialCue_ind = {task.trialparams(:).cue}';
            individual_runs{n_dataset}.Cue = TrialCue_ind;
            TrialCue = vertcat(TrialCue, TrialCue_ind(data_subset));
            
            TrialType_ind = {task.trialparams(:).action}';
            for i = 1:numel(TrialType_ind)
                label = TrialType_ind{i};  % Extract the current trial label
                if contains(label, 'Hand') && ~contains(label, 'Object')
                    TrialType_ind{i} = 'Hand';
                elseif contains(label, 'Object') && ~contains(label, 'Hand')
                    TrialType_ind{i} = 'Object';
                elseif contains(label, 'Hand') && contains(label, 'Object')
                    TrialType_ind{i} = 'Hand_Object';
                else
                    % Handle cases where the label doesn't match any predefined pattern
                    TrialType_ind{i} = 'Unknown';  
                end
            end
            individual_runs{n_dataset}.TrialType = TrialType_ind;
            TrialType = vertcat(TrialType, TrialType_ind(data_subset)); 

            
            cueType_ind = {task.trialparams(:).cueType}';
            cueType = vertcat(cueType, cueType_ind(data_subset));
            session_date_val = {session_dates{n_session}};
         
            session_date_ind =  repmat(session_date_val,numTrials , 1);
            session_date = vertcat(session_date, session_date_ind);
            time_phase_labels_ind = repmat({individual_runs{n_dataset}.phase_labels'},numTrials , 1);
            time_phase_labels = vertcat(time_phase_labels, time_phase_labels_ind);
            time_trial_ind = repmat({individual_runs{n_dataset}.relt},numTrials , 1);
            time_trial = vertcat(time_trial, time_trial_ind);
            
    end
  
    %when combining datasets of different blocks, select features that are
    %present in each block. e.g. ratefilt might kick features our in block
    %1 but not block 2 -> allows combining datasets
    if (length(good_datablocks) >1)
        [combined_all, IA, IB] = intersect(individual_runs{1,1}.featdef(:,{'nsp', 'channel' ,'unit'}),individual_runs{1,2}.featdef(:,{ 'nsp','channel' ,'unit'}));
        c1 = individual_runs{1,1}.featdef(IA, 1:5); %KICK OUT DATASET
        c2 = individual_runs{1,2}.featdef(IB, 1:5); 
        disp('----------------------------------------------------------------');
        disp(['Original number of units for dataset 1: ' num2str(size(individual_runs{1,1}.featdef,1)) ' . Original number of units for dataset 2: ' num2str(size(individual_runs{1,2}.featdef,1))]);
        disp(['Number of kept units: ' num2str(size(c1,1))]);
    else %keep same format
        [combined_all, IA, IB] = intersect(individual_runs{1,1}.featdef(:,{'nsp', 'channel' ,'unit'}),individual_runs{1,1}.featdef(:,{ 'nsp', 'channel' ,'unit'}));
    end

    %if there are more than 2 dataset, repeat so that all combined dataset
    %include the same features
    if length(good_datablocks) > 2  
        for kk = 3:(length(good_datablocks))
            [combined_all, IA, IB] = intersect(combined_all, individual_runs{1,kk}.featdef(:,{ 'nsp', 'channel' ,'unit'})); 
            disp(['Original number of units for dataset' num2str(kk) ' : ' num2str(size(individual_runs{1,kk}.featdef,1))]);
            disp(['Number of kept units: ' num2str(size(combined_all,1))]);
         end   
    end 

    % Adapt the length of the bins to the shortest one for fr_adapted
    %%
    featdef_adapted = [];
    SMG_Go = [];
    PMV_Go = [];
    S1X_Go = [];
    AIP_Go = [];
    M1_Go  = [];
    dPca_data_tmp = [];

    for n_dataset = 1:length(good_datablocks) 
        
        featdef_ind = individual_runs{1,n_dataset}.featdef;
        feature_index_to_keep= ismember(featdef_ind(:,{ 'nsp', 'channel' ,'unit'}),combined_all); 
        %find index of features present in both datasets
        
        %adapt featureset and data according to features to keep
        featdef_ind_sub = featdef_ind(feature_index_to_keep,:);
        dataset_channel = featdef_ind_sub.dataset_channel;   
        channel = featdef_ind_sub.channel;   
        data_ind = individual_runs{1, n_dataset}.fr_adapted;
        data_ind = data_ind(:,feature_index_to_keep,:); 
        
        pca_tmp = permute(data_ind,[2,1,3]);
    
         dPca_data_tmp = cat(3,dPca_data_tmp, pca_tmp);
        
        %separate channels according to brain area
        if strcmp(subject_id, 's2')
            SMG_idx = channel <= 96 .* ismember(featdef_ind_sub.nsp_name, 'APX');
            PMV_idx = logical((channel > 96 & channel <= 224) .* ismember(featdef_ind_sub.nsp_name, 'APX'));
            S1_idx = channel <= 96   & ismember(featdef_ind_sub.nsp_name, 'S1X_S1');
            AIP_idx = dataset_channel < 0; %does not exist for s2
            M1_idx = dataset_channel < 0; %does not exist for s2

            %SMG_idx = dataset_channel <= 96 .* ismember(featdef_ind_sub.nsp_name, 'APX');
            %PMV_idx1 = (dataset_channel > 96 & dataset_channel <= 224) .* ismember(featdef_ind_sub.nsp_name, 'APX');
            %S1_idx2 = dataset_channel > 225  & ismember(featdef_ind_sub.nsp_name, 'S1X_S1');
            %AIP_idx = dataset_channel < 0; %does not exist for s2
            %M1_idx = dataset_channel < 0; %does not exist for s2

            % if nnz(PMV_idx) ~= nnz(PMV_idx1)
            %     keyboard
            %     %problem with separating channels into appropriate brian
            %     %area
            % end

        elseif strcmp(subject_id, 's3')
            %implement dataset channels for s3 for SMG, PMv and S1
            SMG_idx =  ismember(dataset_channel,[129:160, 225:256]);
            PMV_idx =  ismember(dataset_channel,[161:224]);
            S1_idx = ismember(dataset_channel,[(1:128) + 256]);
            AIP_idx =  ismember(dataset_channel, [1:32, 97:128]);
            M1_idx =  ismember(dataset_channel,[33:96]);
        else
            keyboard
            % add participant
        end
        
        if nnz(unique(time_bin_size)) > 1
            keyboard
        end
        
        Brain_idx = [SMG_idx, PMV_idx, S1_idx, AIP_idx, M1_idx];
        Brain_areas = {'SMG', 'PMV', 'S1', 'AIP', 'M1'};
        %separate dataset channels per brain area 
        SMG_Go_ind = arrayfun(@(x) data_ind(:,SMG_idx,x), 1:individual_runs{1, n_dataset}.numTrials,'UniformOutput', false)';
        PMV_Go_ind = arrayfun(@(x) data_ind(:,PMV_idx,x), 1:individual_runs{1, n_dataset}.numTrials,'UniformOutput', false)';
        S1_Go_ind = arrayfun(@(x) data_ind(:,S1_idx,x), 1:individual_runs{1, n_dataset}.numTrials,'UniformOutput', false)';
        AIP_Go_ind = arrayfun(@(x) data_ind(:,AIP_idx,x), 1:individual_runs{1, n_dataset}.numTrials,'UniformOutput', false)';
        M1_Go_ind = arrayfun(@(x) data_ind(:,M1_idx,x), 1:individual_runs{1, n_dataset}.numTrials,'UniformOutput', false)';
   
        %concatenate together
        SMG_Go = vertcat(SMG_Go, SMG_Go_ind);
        PMV_Go = vertcat(PMV_Go, PMV_Go_ind);
        S1X_Go = vertcat(S1X_Go, S1_Go_ind);
        AIP_Go = vertcat(AIP_Go, AIP_Go_ind);
        M1_Go  = vertcat(M1_Go, M1_Go_ind);
    end 
     
  
    TrialNumber = 1:length(SMG_Go); 
    GoLabels = cell(size(TrialNumber))';

    if strcmp(subject_id, 's2')
        Go_data = [array2table(TrialNumber') cell2table(LabelNames) cell2table(GraspType) cell2table(TrialType) array2table(TrialCue) ...
                    cell2table(cueType)  cell2table(SMG_Go) cell2table(PMV_Go) cell2table(S1X_Go) cell2table(GoLabels) ...
                    cell2table(session_date) cell2table(time_phase_labels) cell2table(time_trial)];

    elseif strcmp(subject_id, 's3')
        Go_data = [array2table(TrialNumber') cell2table(LabelNames) cell2table(GraspType) cell2table(TrialType) array2table(TrialCue) ...
                    cell2table(cueType)  cell2table(SMG_Go) cell2table(PMV_Go) cell2table(S1X_Go) cell2table(AIP_Go) cell2table(M1_Go) cell2table(GoLabels) ...
                    cell2table(session_date) cell2table(time_phase_labels) cell2table(time_trial)];
    end 
    
      


        if flag_dPCA 
            % pred data for dPCA
            LabelsInd = (preproc.image2class_simple(LabelNames))';
            numCueTypes = numel(unique(cueType));
            numClasses = numel(unique(LabelsInd));
            %separate fr per cue type
            frPerCueType =  arrayfun(@(x) dPca_data_tmp(:,:,ismember(cueType, x)), unique(cueType), 'UniformOutput', false);
            LabelsPerCueType = arrayfun(@(x) LabelsInd(ismember(cueType, x)), unique(cueType), 'UniformOutput', false);
                
            if strcmp(subject_id, 's2')
                trialPerWord = 8;
            elseif strcmp(subject_id, 's3')
                trialPerWord = 16;
            end 
            %format for dPCA analysis
            dPCAData = zeros(size(frPerCueType{1},1),numCueTypes, numClasses, size(frPerCueType{1},2), trialPerWord);

       for n_cue = 1:numCueTypes

            %firing rate per condition
            frPerCondition = arrayfun(@(x) frPerCueType{n_cue}(:,:,LabelsPerCueType{n_cue} == x), unique(LabelsInd), 'UniformOutput', false);

            %dPCA requires that the data has the same number of repetitions each
            %time -> if one dataset does not have enough repetitions, complete it by
            %the average of the other trials

            trSize = cell2mat(cellfun(@(x) size(x,3), frPerCondition, 'UniformOutput', false));
           
            %find groups with missing trials
            trAdapt = find(trSize ~= trialPerWord);

            for n_tr = 1:length(trAdapt)
                idxToChange = trAdapt(n_tr);
                %calculate average FR per condition
                meanToAdd = nanmean(frPerCondition{idxToChange},3);
                lenToChange = trialPerWord - trSize(idxToChange);
                %add average FR to get 16 trials per repetitions
                for toAdd = 1:lenToChange
                  frPerCondition{idxToChange}(:,:,trSize(idxToChange) + toAdd) = meanToAdd;
                end 
            end 

            for n_classes = 1:numClasses
                dPCAData(:,n_cue,n_classes,:,:) = frPerCondition{n_classes};
            end 

        end 

            dPCAPerArea = arrayfun(@(x) dPCAData(Brain_idx(:,x),:,:,:,:), 1:size(Brain_idx,2), 'UniformOutput', false);
            dPCA = cell(size(TrialNumber))';
            dPCA(1:5) = dPCAPerArea;
            dPCA{6} = Brain_areas;

             Go_data = [Go_data cell2table(dPCA)];
        end


    frPerChannelAll = arrayfun(@(x) dPca_data_tmp(Brain_idx(:,x),:,:), 1:size(Brain_idx,2), 'UniformOutput', false);
    frPerChannel = cell(size(TrialNumber))';
    frPerChannel(1:length(frPerChannelAll)) = frPerChannelAll;
    frPerChannel{length(frPerChannelAll)+1} = Brain_areas;
    
    Go_data = [Go_data cell2table(frPerChannel)];
 
   filename_save = [subject_id '_' session_dates{n_session} '_' spike_sorting_type '_' TaskCue '.mat'];
   
   disp(['Processed file '  filename_save])

   filename_save = fullfile(save_data_pathway,filename_save);
   
    
   save(filename_save, 'Go_data', 'individual_runs', '-v7.3'); 

   

end  

keyboard
 
%%
%combine data together
    
% TO DO: vertically concatenate tables together
% think about that last dataaset only has 156 on average timebins -> where
% exactly do we lose timebins, e.g. in whcih phase? could be all slightly
% misaligned to be honest.. but it is what it is. I also should average
% time phases for each session individually for classification purposes at
% least (does not really work for bin per bin regression analysis). 
 

% Combine datasets into one
%subject_id = 's3';  % s2 or p3 or n1
%spike_sorting_type = 'unsorted_aligned_thr_-4.5';
%spike_sorting_type = 'unsorted_aligned_noratefilt_4.5';
%TaskCue = 'GraspObject_4S_Action';
%spike_sorting_type = 'sorting_aligned_thr_-4.5';
%spike_sorting_type = 'sorting_aligned_noratefilt_4.5';

%spike_sorting_type = 'unsorted_aligned_noratefilt';
%spike_sorting_type = 'unsorted_aligned';

%save_data_pathway = ['D:\Users\Sarah\Documents\Saved_Data\InternalSpeechPaper\' subject_id '\Data\IndividualFiles\' spike_sorting_type];
%save_data_pathway = ['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subject_id '\Data\IndividualFiles\' TaskCue '\' spike_sorting_type];
save_data_pathway =  ['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data\IndividualFiles\' TaskCue '\' spike_sorting_type];

datafiles = dir([save_data_pathway '\*.mat']);

for n_session = 1:length(datafiles)
    
    Data = load(fullfile(datafiles(n_session).folder, datafiles(n_session).name));
    if n_session == 1
        Go_data = Data.Go_data;       
    else
            Go_data = vertcat(Go_data, Data.Go_data);
    end 
    
end 

Go_data.TrialNumber = (1:size(Go_data,1))';
%transform name of condition into code
Go_data.GoLabels = (preproc.image2class_simple(Go_data.LabelNames))';


filename_save = ['Table_' subject_id '_' TaskCue '_' spike_sorting_type   '.mat'];

save_data_pathway =  ['C:\Users\macthurston\OneDrive - Kaiser Permanente\CaltechData\GraspObject_project\' subject_id '\Data'];

filename_save = fullfile(save_data_pathway,filename_save);

save(filename_save, 'Go_data', '-v7.3'); 
disp([ 'saved ' filename_save])
