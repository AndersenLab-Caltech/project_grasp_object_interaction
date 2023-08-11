% Cleaned up preprocessing file for Speech paper. For S2 and S3

clc 
clear all
%close all

subject_id = 's2';  % s2 or p3 or n1
subject_id = 's3';  % s2 or p3 or n1

subject = hst.Subject(subject_id);

flag8TrialBlocks = false; 
flag16TrialBlocks = true; 
flag_dPCA = true; 

if strcmp(subject_id, 's2')
    session_dates = {'20210712', '20210722','20210729','20210923','20210930','20211011','20211018','20211027','20211103', '20220323'};
    flag8TrialBlocks = false;
    flag16TrialBlocks = false;
elseif strcmp(subject_id, 's3')
    session_dates = {'20230106','20230706','20230712','20230713'};
    
    if flag8TrialBlocks
        session_dates = {'20230106','20230106','20230706','20230706','20230706','20230706','20230712','20230712','20230712',...
            '20230713','20230713','20230713','20230713'};
    end
    
    if flag16TrialBlocks
        session_dates = {'20230106','20230706','20230706','20230712','20230713','20230713'};
    end


    %for 20230706, only the 1st run does not have S1 data. I could rewrite
    %the code to include S1 data for the 3 other runs. keep in mind!
else 
    error('unknown subject')
end 

save_data = true; 
flagRemoveTrials = true; 


%spike_sorting_type = 'sorting_aligned_thr_-4.5_noSmoothing';
%spike_sorting_type = 'sorting_aligned_thr_-4.5_smoothed';
%spike_sorting_type = 'sorting_aligned_noratefilt_4.5';
%spike_sorting_type = 'unsorted_aligned_thr_-4.0';
%spike_sorting_type = 'sorting_aligned_thr_-4.5';
%spike_sorting_type = 'sorting_aligned_noratefilt_4.5';

%Might explain why they looked so different over different session days. It
%could be that thresholding it will improve the accuracy that I have? idk

%spike_sorting_type = 'unsorted_aligned_noratefilt_4.5';
%spike_sorting_type = 'unsorted_aligned_noratefilt';
%spike_sorting_type = 'sorting'; % I did not rethreshold the session before spike sorting... idk if that will work? 

if strcmp(subject_id, 's2')
    spike_sorting_type = 'sorting_aligned_thr_-4.5';
elseif strcmp(subject_id, 's3')
    spike_sorting_type = 'sorting';
    spike_sorting_type = 'sorting_noratefilt';

end 

%spike_sorting_type = 'unsorted_aligned_noratefilt';


%GraspCue = 'SpeechTraining';

if strcmp(subject_id, 's2')
    TaskCue = 'Speech';
elseif strcmp(subject_id, 's3')
    TaskCue = 'Speech_WrittenCue';
    
    if flag8TrialBlocks
        TaskCue = 'Speech_WrittenCue_8TrialBlocks';
    end
    
    if flag16TrialBlocks
        TaskCue = 'Speech_WrittenCue_16TrialBlocks';
    end
else
    keyboard
end


save_data_pathway = ['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subject_id '\Data\IndividualFiles\' TaskCue '\' spike_sorting_type];
audio_data_pathway = ['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subject_id '\Data\AudioFiles'];

if ~exist(save_data_pathway)
    mkdir(save_data_pathway)
end

min_timebin_length = 157;

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
    taskfiles = session.getTaskFiles('Speech');

    %Extract correct datasets from excel files 
    filename = [subject_id '_good_trials_' TaskCue '.xlsx'];
    
    [numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech_paper\ExcelFiles\' filename]); 
   
    session_date = str2double(session_dates{n_session});
    good_datasets = numb(session_date==numb(:,1), 2:end);
    
    if last_session_date ~= session_date
        dayIdx = 1;
    else
        dayIdx = dayIdx+1;
    end
    dayIdx
    last_session_date = session_date;
    good_datasets = good_datasets(dayIdx,:);
    
    good_datasets(isnan(good_datasets)) = [];
    
    if (length(good_datasets) <1)
        error ('Add to excel sheet, No good dataset present, skip or check for problem')
    elseif good_datasets == 0
        error('No good dataset present for this day, has been checked, skip')  
    elseif any(ismember(good_datasets, 1000))
        warning('There is a problem with one of the datasets, ask Spencer? Removing it for now');
        good_datasets(ismember(good_datasets,1000)) = []; 
    
    elseif good_datasets == 100
        disp(['No auditory dataset for session ' num2str(session_date)])       
        continue
    else
        disp('Everything seems ok');
    end

    
    disp(['Good_datasets : ' num2str(good_datasets)]); 
       
 
    %%
    %to verify length of time bins
    time_bin_size = zeros(1,length(good_datasets)); 
    
    %Either sorting (sorted spikes), noisy (includes noise units), unsorted
    %(unsorted data), smoothed (halfkerner = 0 and causal = 0, how much the
    %neuronal data is smoothed). 

    individual_runs = cell(1,length(good_datasets));
    featdef = cell(1,length(good_datasets));
    
    %predeclare variables
    TrialNumber = [];
    LabelNames = {};
    GoLabels = [];
    session_date = [];
    run_labels = [];
    cueType = {};
    time_phase_labels = [];
    time_trial = [];
    
    %predeclare variables for audio data
    
    Audio_raw = {};
    Audio_envelope = {};
    Audio_time = {};
   
    blubData = [];
    blubLabels = [];

    for n_dataset = 1:length(good_datasets) 
    
            %load task 
            task = hst.Task(taskfiles{good_datasets(n_dataset)});  
            
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
                       
            %to be removed - verification
            neuralData = individual_runs{n_dataset}.fr_adapted;
            neuralDataInternal = neuralData(individual_runs{n_dataset}.phase_labels == 4,:,:);
            neuralDataInternalSMG = squeeze(mean(neuralDataInternal(:,individual_runs{n_dataset}.featdef.nsp  == 1,:),1));
            Labels = vertcat(task.trialparams.Image_number);            
            [a,errTest] = classification.classification_simple(neuralDataInternalSMG',Labels(data_subset),'PCA_percentage', 95, 'LeaveOneOut', true);
            classification_result(n_session,n_dataset) = errTest;
            
            %save relevant variables for analysis 
            time_bin_size(n_dataset) = size(individual_runs{n_dataset}.fr_adapted,1);
                
            individual_runs{n_dataset}.task.hDebug= []; %remove saving debugger, created errors later when loading the file
            LabelNames_ind = {task.trialparams(:).action}';
            individual_runs{n_dataset}.Labels = LabelNames_ind;
            LabelNames = vertcat(LabelNames, LabelNames_ind(data_subset) );
            cueType_ind = {task.trialparams(:).cueType}';
            cueType = vertcat(cueType, cueType_ind(data_subset));
            session_date_val = {session_dates{n_session}};
            if flag8TrialBlocks || flag16TrialBlocks
                
                session_date_val{1,1} = [session_date_val{1,1} '_' num2str(dayIdx)];
            end
            
            session_date_ind =  repmat(session_date_val,numTrials , 1);
            session_date = vertcat(session_date, session_date_ind);
            time_phase_labels_ind = repmat({individual_runs{n_dataset}.phase_labels'},numTrials , 1);
            time_phase_labels = vertcat(time_phase_labels, time_phase_labels_ind);
            time_trial_ind = repmat({individual_runs{n_dataset}.relt},numTrials , 1);
            time_trial = vertcat(time_trial, time_trial_ind);
            
            %load audio data
            
            audio_file_name = fullfile(audio_data_pathway, [task.taskString '_errorTrials.mat']);
            
            if exist(audio_file_name)
                audio = load(fullfile(audio_data_pathway, [task.taskString '_errorTrials.mat']));
                Audio_raw = vertcat(Audio_raw, audio.audio_per_trial');
                Audio_envelope = vertcat(Audio_envelope, audio.envelope_per_trial');
                Audio_time = vertcat(Audio_time, audio.time_per_trial');
            else                              
                Audio_raw = vertcat(Audio_raw, cell(numTrials,1));
                Audio_envelope = vertcat(Audio_envelope, cell(numTrials,1));
                Audio_time = vertcat(Audio_time, cell(numTrials,1));
            end
    end
 
  
    %when combining datasets of different blocks, select features that are
    %present in each block. e.g. ratefilt might kick features our in block
    %1 but not block 2 -> allows combining datasets
    if (length(good_datasets) >1)
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
    if length(good_datasets) > 2  
        for kk = 3:(length(good_datasets))
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

    for n_dataset = 1:length(good_datasets) 
        
        featdef_ind = individual_runs{1,n_dataset}.featdef;
        feature_index_to_keep= ismember(featdef_ind(:,{ 'nsp', 'channel' ,'unit'}),combined_all); 
        %find index of features present in both datasets
        
        %adapt featureset and data according to features to keep
        featdef_ind_sub = featdef_ind(feature_index_to_keep,:);
        dataset_channel = featdef_ind_sub.dataset_channel;        
        data_ind = individual_runs{1, n_dataset}.fr_adapted;
        data_ind = data_ind(:,feature_index_to_keep,:); 
        
        pca_tmp = permute(data_ind,[2,1,3]);
    
         dPca_data_tmp = cat(3,dPca_data_tmp, pca_tmp);
        
        %separate channels according to brain area
        if strcmp(subject_id, 's2')
            SMG_idx = dataset_channel <= 96;
            PMV_idx = dataset_channel > 96 & dataset_channel <= 192;
            S1_idx = dataset_channel > 192;
            AIP_idx = dataset_channel < 0; %does not exist for s2
            M1_idx = dataset_channel < 0; %does not exist for s2
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
        Go_data = [array2table(TrialNumber') cell2table(LabelNames) cell2table(cueType)  cell2table(SMG_Go) cell2table(PMV_Go) cell2table(S1X_Go) cell2table(GoLabels) ...
        cell2table(session_date) cell2table(time_phase_labels) cell2table(time_trial) cell2table(Audio_raw) cell2table(Audio_envelope) cell2table(Audio_time)];

    elseif strcmp(subject_id, 's3')
        Go_data = [array2table(TrialNumber') cell2table(LabelNames) cell2table(cueType)  cell2table(SMG_Go) cell2table(PMV_Go) cell2table(S1X_Go) cell2table(AIP_Go) cell2table(M1_Go) cell2table(GoLabels) ...
        cell2table(session_date) cell2table(time_phase_labels) cell2table(time_trial) cell2table(Audio_raw) cell2table(Audio_envelope) cell2table(Audio_time)];
    end 
  
    Go_data = renamevars(Go_data,'Var1','TrialNumber');
      
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
 
   filename_save = [subject_id '_' session_dates{n_session} '_' spike_sorting_type '_' TaskCue '.mat'];
   
   if flag8TrialBlocks || flag16TrialBlocks
      filename_save = [subject_id '_' session_dates{n_session} '_' spike_sorting_type '_' TaskCue '_' num2str(dayIdx) '.mat'];
   end 
   filename_save = fullfile(save_data_pathway,filename_save);
   
    
   save(filename_save, 'Go_data', '-v7.3'); 
 

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
subject_id = 's3';  % s2 or p3 or n1
%spike_sorting_type = 'sorting_aligned_thr_-4.5';
%spike_sorting_type = 'sorting_aligned_noratefilt_4.5';

%spike_sorting_type = 'unsorted_aligned_noratefilt';
%spike_sorting_type = 'unsorted_aligned';

%save_data_pathway = ['D:\Users\Sarah\Documents\Saved_Data\InternalSpeechPaper\' subject_id '\Data\IndividualFiles\' spike_sorting_type];
save_data_pathway = ['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subject_id '\Data\IndividualFiles\' TaskCue '\' spike_sorting_type];
datafiles = dir([save_data_pathway '\*.mat']);

dPca_comb = cell(1,5);
flagdPCA = true;
for n_session = 1:length(datafiles)
    
    Data = load(fullfile(datafiles(n_session).folder, datafiles(n_session).name));
    if n_session == 1
        Go_data = Data.Go_data;       
    else
        if flagdPCA
            Go_data = vertcat(Go_data, Data.Go_data);
        else
            Go_data = vertcat(Go_data, Data.Go_data(:,1:15));
        end
    end 
    if flagdPCA
        dPca_comb(1:5) = arrayfun(@(x) cat(1,dPca_comb{x}, Data.Go_data.dPCA{x}), 1:5, 'UniformOutput', false);
    end 
end 

Go_data.TrialNumber = (1:size(Go_data,1))';

Go_data.GoLabels = (preproc.image2class_simple(Go_data.LabelNames))';
if flagdPCA
    Go_data.dPCA(1:5) = dPca_comb;
end

%keyboard

%%
flagIncludeAudio = false;
%filename 
filenameTmp = strsplit(datafiles(1).name, '.mat');
filenameTmp = strsplit(filenameTmp{1}, '_');
filenameTmp = filenameTmp(setdiff(1:length(filenameTmp),2));

if flagIncludeAudio
    filenameTmp{end+1} = 'Audio';
end 
filename_save1 = ['Table_' strjoin(filenameTmp, '_') '.mat'];

%remove audio data to save file space
if ~flagIncludeAudio 
    audioIdx = ismember(Go_data.Properties.VariableNames, {'Audio_raw','Audio_envelope','Audio_time'});
    Go_data = Go_data(:,~audioIdx);
end

save_data_pathway = ['C:\Users\Sarah\OneDrive - California Institute of Technology\Data\InternalSpeechPaper\' subject_id '\Data'];

filename_save = fullfile(save_data_pathway,filename_save1);

%if ~exist(filename_save)
    save(filename_save, 'Go_data', '-v7.3'); 
    disp([ 'saved ' filename_save1])
%end 
