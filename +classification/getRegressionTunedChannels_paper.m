function [tunedCombinedChannels, tunedChannelsPerPhase, tunedChannelsPerBin, sumPhase, sumBin,numTunedChannelsPerCategory, tunedChannelsPerPhasePerCategory,numTunedChannelsPerCategoryBin,p_per_phase, p_per_phaseOrig] = getRegressionTunedChannels_paper(Data,Labels, TimePhaseLabels,varargin)

% Takes as input a cell array and phase labels and returns a bool of tuned units

%default values
multipleCompare = true;
combineTunedChannels = true;
flagBinperBin = false;
flagShuffleTest = false; 

% Loading optional arguments
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'multcompare'
            multipleCompare = varargin{2};
        case 'combinetunedchannels'
            combineTunedChannels = varargin{2};
        case 'binperbintuning'
            flagBinperBin = varargin{2};                          
        case 'flagshuffletest'
            flagShuffleTest = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
      varargin(1:2) = [];
end

%the names of what is being analyzed
categoryNames = utile.label_number_to_grasp_name(unique(Labels))'; 

alpha = 0.05;
numBins = size(Data{1},1);

numConditions = numel(unique(Labels));

[groupCount,~] = groupcounts(Labels);

Condition_names = [{'ITI'}];
Condition = repmat(Condition_names, [max(groupCount),1]);

ConLabels =preproc.image2class_simple(Labels)';

Condition = [Condition; ConLabels];

%Presence vector X
for i = 1:numConditions
    X(:,i) = ismember(Condition, categoryNames{i});
end

numChannels = size(Data{1},2);
numPhases = numel(unique(TimePhaseLabels{1}));

p_per_phase = ones(numPhases,numConditions, numChannels);
p_per_phaseOrig = ones(numPhases, numConditions,numChannels);
p_per_bin = ones(numBins, numConditions, numChannels);
pFStat = ones(numPhases,numChannels);
r_per_phase = zeros(numPhases,numChannels); 
r_per_phase_adjusted = zeros(numPhases,numChannels); 

p_multcomp_per_phase = ones(numPhases,numConditions,numChannels);
p_multcomp_per_bin = ones(numBins, numConditions, numChannels);
try
    ITI_DataAll =  cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== 1,:),1),Data,TimePhaseLabels, 'UniformOutput', false));

catch
   keyboard

end 

for n_phase = 1:numPhases
    
    try
        DataPerPhase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== n_phase,:),1),Data,TimePhaseLabels, 'UniformOutput', false));
    catch
        keyboard
    end 


  for n_channel = 1:size(DataPerPhase,2)
      
      DataPerTrial = DataPerPhase(:,n_channel);

      if flagShuffleTest % randomize trials 
              DataPerTrial = DataPerTrial(randperm(length(DataPerTrial)));
      end
   
      ITI_Data = ones(max(groupCount),1)*mean(mean(ITI_DataAll(:,n_channel)));
      
      %Perform linear regression test for each channel for each phase 
          
      FR = [ITI_Data; DataPerTrial]; 
       
      mdl = fitlm(X,FR);
      pFStatTmp = anova(mdl,'summary');
      %save the pvalue of f-stastistic (is model better than constant coefficients?)
      pFStat(n_phase, n_channel) = pFStatTmp.pValue(2);
      p_per_phase(n_phase,:,n_channel) = mdl.Coefficients.pValue(2:end); 
      p_per_phaseOrig(n_phase,:,n_channel) =  p_per_phase(n_phase,:,n_channel);
      r_per_phase(n_phase, n_channel) = mdl.Rsquared.Ordinary;
      r_per_phase_adjusted(n_phase, n_channel) = mdl.Rsquared.Adjusted;


  end    
end 

if flagBinperBin

    for n_bin = 1:numBins
        disp(['Bin number ' num2str(n_bin)])
        DataPerBin = cell2mat(arrayfun(@(x,y) Data{x,1}(n_bin,:),1:size(Data,1), 'UniformOutput', false)');

        for n_channel = 1:numChannels

          DataPerBinTrial = DataPerBin(:,n_channel);
          DataPerBinTrialOrdered = DataPerBinTrial;
          ITI_Data = ones(max(groupCount),1)*mean(mean(ITI_DataAll(:,n_channel)));

          %perform linear regression for each channel for each phase for each of the 5 grasps 
          FR = [ITI_Data; DataPerBinTrialOrdered];  
          mdl = fitlm(X,FR);
          p_val_bin = mdl.Coefficients.pValue(2:end);
          
          if multipleCompare
              [~,~,p_val_bin] = utile.MultipleComparisonsCorrection(p_val_bin,'method', 'fdr');
          end 
      
          p_per_bin(n_bin,:,n_channel) = p_val_bin; 

        end
    end 
end 


if multipleCompare
    disp('Performing multiple comparison more strict')

    for n_channel = 1:numChannels
        [~,~,p_multcomp_per_phase(:,:,n_channel)] = utile.MultipleComparisonsCorrection(p_per_phase(:,:,n_channel),'method', 'fdr'); 
    end 
    %replace old p values by multcompare p values
    p_per_phase = p_multcomp_per_phase; 
end 


tunedChannelsPerPhasePerCategory = p_per_phase < alpha;
numTunedChannelsPerCategory = sum(tunedChannelsPerPhasePerCategory,3);
numTunedUnitsPhase = cell2mat(arrayfun(@(x) nnz(squeeze(sum(squeeze(tunedChannelsPerPhasePerCategory(x,:,:))))), 1:size(tunedChannelsPerPhasePerCategory,1), 'UniformOutput', false));

 

tunedChannelsPerPhase = logical(squeeze(sum(tunedChannelsPerPhasePerCategory,2)))';
sumPhase = sum(tunedChannelsPerPhase);
FStatSumPhase = sum(pFStat < 0.05,2); %not identical to sumPhase... maybe I should incorporate both to be more exact for selection of tuned channels?
warning('think about fstatistic to make better tuned units prediction?')
tunedChannelsPerBin = p_per_bin < 0.05;
numTunedChannelsPerCategoryBin = sum(tunedChannelsPerBin,3);

sumBin = sum(logical(squeeze(sum(tunedChannelsPerBin,2))),2); 

if combineTunedChannels
    tunedCombinedChannels = logical(sum(tunedChannelsPerPhase,2));    
end 




end

