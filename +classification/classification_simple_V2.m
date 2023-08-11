function [accTrain,accTest] = classification_simple_V2(Data,Labels, varargin)

%perform cross validation classification based on "Data"

[varargin,NbrRep] = Utilities.argkeyval('repetition',varargin, 1 );
[varargin,k] = Utilities.argkeyval('k',varargin, 8);
[varargin,classifier] = Utilities.argkeyval('classifierType',varargin,'diaglinear');
[varargin,flagErrorMatrix] = Utilities.argkeyval('ErrorMatrix',varargin,false);
[varargin,PCA_analysis] = Utilities.argkeyval('PCA',varargin,true);
[varargin,flagLeaveOneOut] = Utilities.argkeyval('LeaveOneOut',varargin,false);
[varargin,PCA_var_p] = Utilities.argkeyval('PCA_percentage',varargin,90);
[varargin,random_perm] = Utilities.argkeyval('randomPerm',varargin,false);
[varargin,Tuning_analysis] = Utilities.argkeyval('TunedUnits',varargin,false);

%[varargin,flagAverageBeforePCA] = Utilities.argkeyval('flagAverageBeforePCA',varargin, true);
Utilities.argempty(varargin);

NbrLabels = length(unique(Labels));

if Tuning_analysis && PCA_analysis
    error('Do you want to do both tuning and PCA?')
    %if on purpose remove 
end 

% Data in trials x features
% Labels in trials x 1 
    
%% decode cue - not working well at all, and very variable even if I average over 8 trials 

%figure(); gscatter(score(:,1), score(:,2), LabelsDecodeCue)
%PCA

if size(Data,1) ~= size(Labels,1)
    error('Wrong label size')
end

errTrain = zeros(NbrRep,k);
errTest = zeros(NbrRep,k);

for rep = 1:NbrRep
    
    if flagLeaveOneOut
        cv = cvpartition(size(Data,1),'LeaveOut');
    else
        cv = cvpartition(size(Data,1), 'KFold', k);
    end 
    
    %randomize trial labels
    if random_perm
        Labels = Labels(randperm(length(Labels)));
    end 

    for cvRun = 1:cv.NumTestSets %

        trIdx = find(cv.training(cvRun));
        teIdx = find(cv.test(cvRun));
        %TrainingSet
        DataTrain = Data(trIdx,:);
        
        %LabelsTrain(:,rep) = Labels(trIdx); 
        LabelsTrain = Labels(trIdx); 

        DataTest = Data(teIdx,:); 
        LabelsTestAll{:,cvRun,rep} = Labels(teIdx);
        LabelsTest = Labels(teIdx);


        if Tuning_analysis
            numChannels = size(DataTrain,2);
            p_val = nan(1,numChannels);
            for ChannelNbr = 1:numChannels

                DataPerTrial = DataTrain(:,ChannelNbr);
                 [p_val(ChannelNbr), ~, ~] = kruskalwallis(DataPerTrial,LabelsTrain,'off');

            end
%            warning('not controling for number of comparison yet')
%             [~,~,p_val] = Utilities.MultipleComparisonsCorrection(p_val,'method', 'fdr'); 
% 
%     
             TunedChannels = p_val < 0.05; 
            
            if nnz(TunedChannels) == 0
                TunedChannels = 1:numChannels;
                warning('no tuned channels available - taking all channels ')
            end 
            DataTrain = DataTrain(:,TunedChannels);
            DataTest = DataTest(TunedChannels); 
            
           
            
            
        end 
        
        
        if PCA_analysis
            [coeff, ~, ~,~, explained] = pca(DataTrain); 
            
            variance = cumsum(explained); 
            
            idx_90 = find(variance > PCA_var_p,1); 
            %idx_90
            PCA_DataTrain = DataTrain*coeff;
            PCA_DataTest = DataTest*coeff;
            DataTrain = PCA_DataTrain(:,1:idx_90);
            DataTest = PCA_DataTest(:,1:idx_90);
            %disp(['Performing PCA analysis: ' num2str(idx_90)])
            %figure(); gscatter(PCA_DataTrain(:,1), PCA_DataTrain(:,2), LabelsTrain)
            %figure(); gscatter(PCA_DataTrain(:,1), PCA_DataTrain(:,2), LabelsTrain)
        end 
        
        model = fitcdiscr(DataTrain, LabelsTrain, 'DiscrimType', classifier);
        %model = fitcnb(DataTrain, LabelsTrain);

       % model = fitcdiscr(DataTrain, LabelsTrain, 'DiscrimType', 'diaglinear');

        PredictedTrainAll{:,cvRun,rep} = predict(model, DataTrain);
        PredictedTestAll{:,cvRun,rep} = predict(model, DataTest);
        PredictedTrain = predict(model, DataTrain);
        PredictedTest = predict(model, DataTest);
        
        errTrain(rep,cvRun) = 1-(nnz(LabelsTrain == PredictedTrain)/numel(LabelsTrain));
        errTest(rep,cvRun) = 1-(nnz(LabelsTest == PredictedTest)/numel(LabelsTest));
        
    end

end 

accTrain = 1-mean(squeeze(errTrain));
accTest = 1-mean(squeeze(errTest));



if flagErrorMatrix
    
    LabelsTest = cell2mat(LabelsTestAll);
    LabelsTest = LabelsTest(:);
    
    PredictedTestAll = cell2mat(PredictedTestAll);
    PredictedTestAll = PredictedTestAll(:);
    %figure()
    [cm,gn] = confusionmat(LabelsTest, PredictedTestAll);
    Names = utile.label_number_to_grasp_name(gn);

    
    confusionchart(cm,Names)
    %disp(['Classification ITI correct: ' num2str(C(1)/sum(C(1,:))*100)])
    %disp(['Classification Speech correct: ' num2str(C(2,2)/sum(C(1,:))*100)])

    %ErrorMatrix = classification.misclassification(LabelsTest, PredictedTestAll);
    %utile.plot_error_matrix(ErrorMatrix, Names , '');
end


end

