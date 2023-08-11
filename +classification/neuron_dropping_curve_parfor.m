function [] = neuron_dropping_curve_parfor(phase_info,all_units_low_FR_removed,Figure_labels, varargin)

%I will be able to use this code to perform classification
%but I have to change a few things: 
%1) use not per trial dataset, use the one combining units as "features"
%2) average firing rate over phases -> 40 times x classification
%   - is this an issue...? I will be classifying 40 trials with more than
%   40 features which is kinda weird, how do I augment the number of trials
%   I have? 
%3) Maybe I should do PCA and use 


%Classification is performed by phase: we take one phase out of
%ITI, ImageCue, Delay and Action phase and we classify the rest of the data
%with it. There are two different ways of doing it: once by taking the
%whole data, and once per bin.

[varargin,testing_phase_name] = Utilities.argkeyval('testing_phase',varargin, 'Action'); %Inputs: ITI, ImageCue, Delay, Action
[varargin,num_of_iterations] = Utilities.argkeyval('iteration_number',varargin, 1000);
[varargin,k] = Utilities.argkeyval('number_of_folds',varargin, 4);
[varargin,classifier] = Utilities.argkeyval('classifier',varargin, 'diaglinear');


Utilities.argempty(varargin);

testing_phase = 2 %find(ismember(data.task.phaseNames, testing_phase_name));

%testing_phase = find(ismember(data.task.phaseNames, testing_phase_name));

%get labels and data from corresponding phase
%for repetition = 1:100
labels = phase_info.labels{testing_phase};
data_tot = phase_info.data{testing_phase};
data_tot = data_tot(:,logical(all_units_low_FR_removed)); %remove all units with too low firing rate

if size(data_tot,2) > 200 %fix max number of units at 150, not so useful otherwise
    max_feature_number = 200;
else
    max_feature_number = size(data_tot,2);
end 

result = cell(1,num_of_iterations);

% result2 = nan(1,num_of_iterations);
% parfor rep = 1:num_of_iterations
%     % processing....
%     result{rep} = true; % result of processing this iteration
%     result2(rep) = 1;
% end

num_phases = 4; %length(data.phase_time_idx)-1;
parfor rep = 1:num_of_iterations
%for rep = 1:num_of_iterations %for debugging
    err = nan(k,num_phases,max_feature_number);
    
    
    for num_features = 1:max_feature_number
        
        %completely randomize the available features each time
        features_rand = randperm(size(data_tot,2));
        %keep 1 feature in loop 1, 2 features in loop 2, etc. 
        features_idx = features_rand(1:num_features);
        %select the data with the selected features (= units)
        data_cv = data_tot(:,features_idx);
        neuralData_per_cell = arrayfun(@(x)data_cv(phase_info.trials{testing_phase}==x, :), 1:phase_info.trials{testing_phase}(end), 'UniformOutput', false);
        labels_cell = arrayfun(@(x)labels(phase_info.trials{testing_phase}==x), 1:phase_info.trials{testing_phase}(end), 'UniformOutput', false);
        
        cv = cvpartition(length(neuralData_per_cell), 'KFold', k);
        cv = repartition(cv); %pas forcement utile
        
        for run = 1:cv.NumTestSets %iterate over different data partitions
            %find training and testing data and labels indexes
            trIdx = find(cv.training(run));
            teIdx = find(cv.test(run));
            
            labels_train = labels_cell(trIdx);
            labels_train = cat(1,labels_train{:});
            
            labels_test = labels_cell(teIdx);
            labels_test = cat(1,labels_test{:});
            %test_labels = cat(1,labels_test{:});
            
            training = neuralData_per_cell(trIdx);
            training = cat(1, training{:}); %concatenate data for training
            
            data_testphase = neuralData_per_cell(teIdx);
            data_testphase =  cat(1, data_testphase{:});
            
            
            model2 = fitcdiscr(training, labels_train, 'DiscrimType', classifier);
            
            for current_phase = 1:num_phases
                
                if current_phase == testing_phase
                    predicted_labels = predict(model2, data_testphase);
                    
                    err(run,current_phase,num_features) = classerror(labels_test, predicted_labels);
                    %err2(run, current_phase, num_features) = classerror(labels_test, predicted_labels);
                    
                else
                    data_t = phase_info.data{current_phase};
                    data_t = data_t(:,logical(all_units_low_FR_removed)); % we have to put the correct units, and remove the ones that don-t have enough 
                                                                          % otherwise bad classification
                    data_t = data_t(:,features_idx);
                   
                    labels_t = phase_info.labels{current_phase};
                    
                    neuralData_per_cell_t = arrayfun(@(x)data_t(phase_info.trials{current_phase}==x, :), 1:phase_info.trials{current_phase}(end), 'UniformOutput', false);
                    labels_cell_t = arrayfun(@(x)labels_t(phase_info.trials{current_phase}==x), 1:phase_info.trials{current_phase}(end), 'UniformOutput', false);
                    
                    labels_test_t = labels_cell_t(teIdx);
                    labels_test_t = cat(1,labels_test_t{:});
                    %test_labels = cat(1,labels_test_t{:});
                    
                    data_testphase_t = neuralData_per_cell_t(teIdx);
                    data_testphase_t =  cat(1, data_testphase_t{:});
                    
                    predicted_labels = predict(model2, data_testphase_t);
                    err(run,current_phase,num_features) = classerror(labels_test_t, predicted_labels);
                    %  result{rep} = err2;
                end
                %err(run, current_phase, num_features) = classerror(test_labels, predicted_labels);
                
               % err(run, current_phase, num_features) = classerror(test_labels, predicted_labels);
                %   err(:,:,:,rep) = err2;
                
            end
            %sliding_window_classification(data, model2, 'bins', 10, 'step',5);
            
        end
    end
    
    result{rep} = err;
end


disp('here')



%end
%end
%One error per phase
save_name = [Figure_labels.subject 'NeuronDroppingCurve_3_grasps_MODEL_CUE_PHASE_Region_' Figure_labels.unit_region '_' date '.mat'];
%save('p3_NeuronDroppingCurve_AIPandPMV_02_11.mat', 'result', 'Figure_labels')
save(save_name,  'result', 'Figure_labels')
%%
PhaseNames = {' ITI  ', ' Image Cue ', ' Delay ' , ' Action '};
err = cat(4,result{:});
colors_plot = utile.get_color_rgb_codes(PhaseNames);
Fontsize = 15; 
%err = [cross validation x current_phase x number of features x repetition]
figure();
for phase = 1:size(err  ,2)
    
err_specific = squeeze(err(:,phase,:,:));

%first average over the 4 cross validations
%mean_cross_validations= squeeze(mean(err,1));
mean_specific = squeeze(mean(err_specific));
ci_specific = bootci(1e3, {@nanmean, mean_specific'});

mean_specific = mean(mean_specific,2); 

err_ci_specific(1,:) = ci_specific(2,:) - mean_specific'; 
err_ci_specific(2,:) = mean_specific' - ci_specific(1,:); 

%Check that ci is correct in errorbar
% figure()
% hold on 
% plot(abs(1-mean_specific)*100,'b')
% hold on 
% plot(abs(1-ci_specific(1,:))*100, 'r*')
% hold on 
% plot(abs(1-ci_specific(2,:))*100,'g*')
% 
%figure();
hold on 
%e = errorbar(1:length(mean_err),abs(1 -mean_err)*100, err_ci1*100, err_ci2*100);
e(phase) =errorbar(1:length(mean_specific),abs(1 -mean_specific)*100,err_ci_specific(1,:)*100, err_ci_specific(2,:)*100);
hold on 
%for i = 1:length(e)
e(phase).Color = colors_plot{phase};
%end 
%plot(abs(1-mean_err)*100,'.');
end 
title([Figure_labels.subject '  - Neuron Dropping Curve' ],'FontSize' , Fontsize);
    %- Training over Action phase'],'FontSize', Fontsize);
xlabel('Number of Neurons','FontSize', Fontsize);
ylabel('Mean Classification Accuracy [%]','FontSize', Fontsize);
legend(PhaseNames, 'FontSize', Fontsize);






% figure();
% 
% %first average over the 4 cross validations
% mean_cross_validations= squeeze(mean(err,1));
% %generate standard deviation over the 1000 iterations
% mean_err = squeeze(mean(mean_cross_validations,3))';
% err_ci = bootci(1e3, {@nanmean, mean_err'});
% std_dev = std(mean_cross_validations, [],3)';
% 
% 
% %figure();
% errorbar(abs(1 -mean_err)*100, std_dev*100);
% %plot(abs(1-mean_err)*100,'.');
% title([Figure_labels.unit_region ' Neuron Dropping Curve']);
% xlabel('Number of Neurons');
% ylabel('Mean Classification Accuracy [%]');
% legend(LegendNames);



% color = {'rp','bp','mp','bp'};
% chance_level = 1/(length(unique(labels)))*100;
% for kk = 4 %1:size(err,3)
%   % subplot(2,2,kk)
%     err_data = squeeze(err(:,:,kk,:));
%     err_data = permute(err_data, [1 3 2]);
%     err_data = reshape(err_data, [], size(err_data,1), 1);
%   %  mean_err = squeeze(mean(err(:,:,kk)))
%     mean_err = mean(err_data)
%     %blub = mean(squeeze(std(err)),2);
%     std_deviation = std(err_data)
%  %   blubi = plot(mean_err*100, color{kk});
%     hold on
%     ylim([0 80]);
%     xlim([0.8 4.2]);
%     errorbar(abs(1 -mean_err)*100, std_deviation*100, color{kk});
%     hold on
%     xL = get(gca,'XLim');
%    % l = line([1:size(err,2)],repmat(chance_level, [1, size(err,2)]),xL,'Color','k','LineStyle','--','Linewidth', 2);
%     l = line([0.8 4.2],[chance_level,chance_level],xL,'Color','r','LineStyle','--','Linewidth', 1);
%
%     legend(l, 'Chance level')
%     legend boxoff                   % Hides the legend's axes (legend border and background)
%
%     set(gca,'xtick',[1:size(err,2)],'xticklabel',PhaseNames)
%     xlabel('Phase');
%     ylabel('Classification Accuracy [%]');
%     truc = [data.region ': Classification using' PhaseNames{kk},'phase'];
%     title(truc);
% end
% %%
%
%
% %sliding window for each phase, obtaining several errors per phase
%
%
%
%


end

