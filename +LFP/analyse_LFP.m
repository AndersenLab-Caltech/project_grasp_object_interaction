clc
%clear all 
close all

session_dates = {'20210712','20210722','20210729','20210923','20210930','20211011','20211018','20211027','20211103','20211230'};
DataFolder = 'D:\Users\Sarah\Documents\Saved_Data\s2_LFP';
data_all = [];
labels_all = [];

for n_sess = 1:length(session_dates)
    disp(['Processing session ' num2str(n_sess)])
    data = load(fullfile(DataFolder, ['SpeechProcessingWrittenCue_' session_dates{n_sess}]));
    
    data_p = abs(data.data.p); 
    data_all = [data_all,data_p];
    labels_all = [labels_all;data.data.class];
end 

%%
%data_all: 1500 x 640 x 288 x 6: frequencies x trials x channels x phases
% I DID NOT KICK OUT WRONG TRIALS HERE

[s_,idx_] = sort(labels_all); 

data_SMG = data_all(:,:,1:96,:); 
data_SMG_sorted = data_all(:,idx_,1:96,:); 

figure(); 

plot(squeeze(mean(data_SMG(10:20,:,:,2))))

figure();
imagesc(squeeze(mean(data_SMG(300:1500,:,:,2))))
colorbar

figure();
imagesc(squeeze(data_SMG(10:20,:,1,6))')
colorbar


data_per_cond = arrayfun(@(x) squeeze(mean(data_all(:,labels_all == x, :,:),2)), unique(labels_all), 'UniformOutput',false);


figure();
imagesc(squeeze(data_per_cond{1}(300:1500,1:96,6)))
colorbar

%%

frequency_range = [4,8;8,12;12,30;30,70;70,150;150,300; 300,500;500,700; 700,1500];
frequency_range = [800,1000;1200,1500];

PhaseNames = {'ITI', 'Cue', 'Delay1', 'Imagined', 'Delay2', 'Action'};
Labels_session = reshape(repmat([1:10],8*8,1), [],1);
l_phases = length(PhaseNames); 
l_fre = length(frequency_range);
err_all = zeros(l_fre,l_phases);

data_all_average = {};
%SOMETHING IS WRONG WITH THE LAST SESSION? THE VALUES ARE WEIRD
data_SMG = data_SMG(:,1:576,:,:);
labels_all = labels_all(1:576); 
%%
for n_fr = 1:l_fre
    
   fr_tmp = frequency_range(n_fr,:); 
   data_tmp_all = [];
    for n_phase = 1:l_phases

        data_tmp = squeeze(mean(data_SMG(fr_tmp(1):fr_tmp(2),:,:, n_phase)));
        [accTrain,accTest] = classification.classification_simple(data_tmp,labels_all, 'LeaveOneOut', false,'PCA', true,'PCA_percentage', 90,'classifierType', 'diagQuadratic');
        err_all(n_fr, n_phase) = accTest*100;
        %disp([num2str((1-mean(errTest))*100)])
        data_tmp_all = cat(3,data_tmp_all,data_tmp);
    end 
    data_all_average{n_fr} = data_tmp_all;
end

figure(); 
bar(err_all', 'LineWidth', 2); 
legend(num2str(frequency_range))
ylabel('classification accuracy')
xticklabels(PhaseNames)
yline(1/8*100, 'r--', 'Chance','Linewidth', 2)

%save('D:\Users\Sarah\Code\Python\Paper2Speech\LFP_data', 'data_all_average', 'labels_all', '-v7.3');