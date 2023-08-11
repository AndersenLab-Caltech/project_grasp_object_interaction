function [TunedCombinedChannels, TunedChannelsPerPhase, TunedChannelsPerBin, sumPhase, sumBin,NbrTunedChannelsPerGrasp,NumberOfUnitsTunedToSeveralGrasps, TunedChannelsPerPhasePerGrasp,NbrTunedChannelsPerGraspBin,p_per_phase, p_per_phaseOrig] = getRegressionTunedChannels(Data,Labels, TimePhaseLabels,varargin)
% TunedChannelsPerBin, sumPhase, sumBin

% Takes as input a cell array and phase labels and returns a bool of tuned units

[varargin,multipleCompare] = Utilities.argkeyval('multcompare',varargin, true);
[varargin,multipleComparePhase] = Utilities.argkeyval('multcomparePhase',varargin, true);
[varargin,combineTunedChannels] = Utilities.argkeyval('combineTunedChannels',varargin, true);
[varargin,flagBinperBin] = Utilities.argkeyval('BinperBinTuning',varargin, true);
[varargin,removeITItuning] = Utilities.argkeyval('removeITItuning',varargin, false);
[varargin,flagAnalyzeGrasps] = Utilities.argkeyval('analyzeGrasps',varargin, false);
[varargin,flagAnalyzeDifferentCues] = Utilities.argkeyval('analyzeDifferentCues',varargin, false);
[varargin,flagShuffleTest] = Utilities.argkeyval('flagShuffleTest',varargin, false);

[varargin,unit_region] = Utilities.argkeyval('unit_region',varargin, '');
[varargin,cue_type] = Utilities.argkeyval('cue_type',varargin, '');


[varargin,SelectedGrasp] = Utilities.argkeyval('SelectedGrasp',varargin, 'all');


%the names of what is being analyzed
if ~ flagAnalyzeDifferentCues
    graspNames = utile.label_number_to_grasp_name(unique(Labels))';
    graspNames = [graspNames{:}];
else 
    CueNames = {'Image Cue',' Audio Cue',' Written Cue'};
    graspNames = CueNames(unique(Labels)');
end 
phaseNames = {'ITI', 'Cue', 'Delay','Action'};

Utilities.argempty(varargin);
alpha = 0.05;
numBins = size(Data{1},1);

%graspNames = {'Lateral','WritingTripod','MediumWrap','PalmarPinch','Sphere3Finger'};
numGrasps = length(graspNames);
fixed_trial_number = unique(histc(Labels, unique(Labels))); % the number of trials should be the same for all 

% %I don't want to add the ITI as baseline when I compare the different cues,
% %as the ITI is not necessarily different than the others 
if flagAnalyzeDifferentCues
    Condition_names = [graspNames];
    Condition = repmat(Condition_names, [fixed_trial_number,1]);
    Condition = reshape(Condition, [fixed_trial_number*(numGrasps),1]);

else
    Condition_names = ['ITI',graspNames];
    %Condition_names = [graspNames];
    Condition = repmat(Condition_names, [fixed_trial_number,1]);
    %Condition = reshape(Condition, [fixed_trial_number*(numGrasps),1]);
    %Condition = {'}
    Condition = reshape(Condition, [fixed_trial_number*(numGrasps+1),1]);

end 




%Presence vector X
for i = 1:numGrasps
    X(:,i) = ismember(Condition, graspNames{i});
end

%X = [zeros(1,5);X];

numChannels = size(Data{1},2);
UniqueTmp = unique(TimePhaseLabels{1});
numPhases = UniqueTmp(end);

p_per_phase = ones(numPhases,numGrasps, numChannels);
p_per_phaseOrig = ones(numPhases, numGrasps,numChannels);
p_per_bin = ones(numBins, numGrasps, numChannels);
pFStat = ones(numPhases,numChannels);
r_per_phase = zeros(numPhases,numChannels); 
r_per_phase_adjusted = zeros(numPhases,numChannels); 

p_multcomp_per_phase = ones(numPhases,numGrasps,numChannels);
p_multcomp_per_bin = ones(numBins, numGrasps, numChannels);
try
    ITI_DataAll =  cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== 1,:),1),Data,TimePhaseLabels, 'UniformOutput', false));

catch
    warning('Different length of phases - taking TimephaseLabels from 1st one')
    ITI_DataAll =  cell2mat(arrayfun(@(x) mean(x{1,1}(TimePhaseLabels{1}== 1,:),1),Data, 'UniformOutput', false));

end 
      warning('changed the ITI data to ALL THE SAME instead of random ITI for each')

for PhaseNbr = 1:numPhases
    try
        %DataPerPhase = arrayfun(@(x,y) x{1,1}(y{:}== PhaseNbr,:),Data,TimePhaseLabels, 'UniformOutput', false);

        DataPerPhase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== PhaseNbr,:),1),Data,TimePhaseLabels, 'UniformOutput', false));
    catch
        disp('Not the same number of units everywhere, how freaking annoying. Cannot even just take less because that will not work')
    end 
      
    %blub = cell2mat(arrayfun(@(x,y) x{1,1}(y{:}== 1,:),Data,TimePhaseLabels, 'UniformOutput', false));
    %blub2 = arrayfun(@(x,y) mean(x{1,1}(y{:}== 1,96)),Data,TimePhaseLabels, 'UniformOutput', false);


  for ChannelNbr = 1:size(DataPerPhase,2)
      
      DataPerTrial = DataPerPhase(:,ChannelNbr);
      %I DID NOT EXTRACT THE DATA CORRECTLY - dumbass. Data per trial is in
      %an random order, while the X that I predefined is in an EXACT order
      %- I just gave it the data in random order so it didn't find any
      %significant results -_-
      %DataPerTrialOrdered = arrayfun(@(x) DataPerTrial(Labels == x), unique(Labels), 'UniformOutput', false)
      DataPerTrialOrdered = cell2mat(arrayfun(@(x) DataPerTrial(Labels == x), unique(Labels), 'UniformOutput', false));

      
      if flagShuffleTest % randomize trials 
              DataPerTrialOrdered = DataPerTrialOrdered(randperm(length(DataPerTrialOrdered)));
      end
          
      %Compute average ITI over all the different grasps. Compare to
      %different phases per grasp
      %ITI_Data1 = cell2mat(arrayfun(@(x) ITI_DataAll(Labels == x,ChannelNbr),unique(Labels), 'UniformOutput', false)');
        
      %different ITI EVERYWHERE
      %ITI_Data = mean(cell2mat(arrayfun(@(x) ITI_DataAll(Labels == x,ChannelNbr),unique(Labels), 'UniformOutput', false)'),2);
      
      %same ITI EVERYWHERE
      ITI_Data = ones(fixed_trial_number,1)*mean(mean(ITI_DataAll(:,ChannelNbr)));
      
     
      %DataPerTrial
      %Perform linear regression test for each channel for each phase 
      
      if flagAnalyzeDifferentCues
          FR = [DataPerTrialOrdered]; %I don't want to add ITI as baseline when I compare differnet cues
      else
          
          FR = [ITI_Data; DataPerTrialOrdered]; % Right now: add the mean of the whole ITI data (2s) as baseline. 
      end 
      mdl = fitlm(X,FR);
      pFStatTmp = anova(mdl,'summary');
      %save the pvalue of f-stastistic (is model better than constant
      %coefficients)
      pFStat(PhaseNbr, ChannelNbr) = pFStatTmp.pValue(2);
      
      p_per_phase(PhaseNbr,:,ChannelNbr) = mdl.Coefficients.pValue(2:end); %Grasps in correct order
      p_per_phaseOrig(PhaseNbr,:,ChannelNbr) =  p_per_phase(PhaseNbr,:,ChannelNbr);
      r_per_phase(PhaseNbr, ChannelNbr) = mdl.Rsquared.Ordinary;
      r_per_phase_adjusted(PhaseNbr, ChannelNbr) = mdl.Rsquared.Adjusted;


  end
  
       
end 

if flagBinperBin

    for BinPhaseNbr = 1:numBins
        disp(['Bin number ' num2str(BinPhaseNbr)])
        DataPerBin = cell2mat(arrayfun(@(x,y) Data{x,1}(BinPhaseNbr,:),1:size(Data,1), 'UniformOutput', false)');

        for ChannelNbr = 1:numChannels

          DataPerBinTrial = DataPerBin(:,ChannelNbr);
          DataPerBinTrialOrdered = cell2mat(arrayfun(@(x) DataPerBinTrial(Labels == x), unique(Labels), 'UniformOutput', false));
          %ITI_Data = mean(cell2mat(arrayfun(@(x) ITI_DataAll(Labels == x,ChannelNbr),unique(Labels), 'UniformOutput', false)'),2);
          ITI_Data = ones(fixed_trial_number,1)*mean(mean(ITI_DataAll(:,ChannelNbr)));
          %perform linear regression for each channel for each phase for each
          %of the 5 grasps 
          FR = [ITI_Data; DataPerBinTrialOrdered]; % Right now: add the mean of the whole ITI data (2s) as baseline. 
          mdl = fitlm(X,FR);

          p_per_bin(BinPhaseNbr,:,ChannelNbr) = mdl.Coefficients.pValue(2:end); %Grasps in conrrect order


        end

    end 

end 


%Do multiple comparison testing for each Channel (we did 4 test per unit (one for each phase), so
%multcompare between 4
if multipleCompare
            disp('Performing multiple comparison')

    for ChannelNbr = 1:numChannels
        
        [~,~,p_multcomp_per_phase(:,:,ChannelNbr)] = Utilities.MultipleComparisonsCorrection(p_per_phase(:,:,ChannelNbr),'method', 'fdr'); 
        [~,~,p_multcomp_per_bin(:,:,ChannelNbr)] = Utilities.MultipleComparisonsCorrection(p_per_bin(:,:,ChannelNbr),'method', 'fdr'); 

    end 
    %replace old p values by multcompare p values
    p_per_phase = p_multcomp_per_phase; 
    p_per_bin = p_multcomp_per_bin;

end 

% %to only do multiple compare for the phase but not for each bin
% 
if multipleComparePhase
    disp('Performing multiple comparison more strict')

    for ChannelNbr = 1:numChannels
        [~,~,p_multcomp_per_phase(:,:,ChannelNbr)] = Utilities.MultipleComparisonsCorrection(p_per_phase(:,:,ChannelNbr),'method', 'fdr'); 
    end 
    %replace old p values by multcompare p values
    p_per_phase = p_multcomp_per_phase; 
end 

%to only do multiple compare for the phase but not for each bin, also
%separately for each phase

% if multipleComparePhase
%     disp('Performing multiple comparison less strict')
%   % warning('multiple comparison per phase, individually for each phase')
%     for ChannelNbr = 1:numChannels
%         for i = 1:numPhases
%             [~,~,p_multcomp_per_phase(i,:,ChannelNbr)] = Utilities.MultipleComparisonsCorrection(p_per_phase(i,:,ChannelNbr),'method', 'fdr'); 
% 
%         end 
%     end 
% %    replace old p values by multcompare p values
%     p_per_phase = p_multcomp_per_phase; 
% end 


TunedChannelsPerPhasePerGrasp = p_per_phase < 0.05;
NbrTunedChannelsPerGrasp = sum(TunedChannelsPerPhasePerGrasp,3);
NumberTunedUnitsPhase = cell2mat(arrayfun(@(x) nnz(squeeze(sum(squeeze(TunedChannelsPerPhasePerGrasp(x,:,:))))), 1:size(TunedChannelsPerPhasePerGrasp,1), 'UniformOutput', false));

%report every unit that is tuned (regardless to which grasp it is tuned to)

if strcmp(SelectedGrasp, 'all')
    GraspIdx = [1:size(TunedChannelsPerPhasePerGrasp,2)];
else
    GraspIdx = find(ismember(graspNames, SelectedGrasp));
end

TunedChannelsPerPhase = logical(squeeze(sum(TunedChannelsPerPhasePerGrasp(:,GraspIdx,:),2)))';
sumPhase = sum(TunedChannelsPerPhase);
FStatSumPhase = sum(pFStat < 0.05,2); %not identical to sumPhase... maybe I should incorporate both to be more exact for selection of tuned channels?
warning('think about fstatistic to make better tuned units prediction?')
TunedChannelsPerBin = p_per_bin < 0.05;
NbrTunedChannelsPerGraspBin = sum(TunedChannelsPerBin,3);

sumBin = sum(logical(squeeze(sum(TunedChannelsPerBin,2))),2); 


if flagAnalyzeGrasps
    Meta = struct; 
    Meta.number_phases = numPhases; 
    Meta.number_units = numChannels; 
    Meta.number_grasps = numGrasps;
    Meta.fixed_trial_number = fixed_trial_number; 
    Meta.Grasp_names = graspNames; 
    Meta.PhaseNames = phaseNames; 
    Figure_labels.unit_region = unit_region;
    Figure_labels.cue_type = cue_type; 
    [NumberOfUnitsTunedToSeveralGrasps] = analyze.analyze_grasp_tuning(TunedChannelsPerPhasePerGrasp, Meta,Figure_labels);
%Tuned to 0, 1, 2 , 3, 4, or 5 grasps

end 

if combineTunedChannels
    %if removeITItuning
    %    TunedCombinedChannels = logical(sum(p_per_phase(:,2:end) < 0.05,2));
    %else
        TunedCombinedChannels = logical(sum(TunedChannelsPerPhase,2));
    %end 
    
end 


% mean(r_per_phase(:, TunedCombinedChannels),2)
% mean(r_mdl(:, TunedCombinedChannels),2)
%NumberTunedUnitsPhase = cell2mat(arrayfun(@(x) nnz(squeeze(sum(squeeze(TunedChannelsPerPhase(x,:,:))))), 1:size(TunedChannelsPerPhase,1), 'UniformOutput', false));




% switchingUnits = 0;
% ActionUnits = 0;
% VisualUnits = 0;
% VisuoMotor = 0; 
% NotTuned = 0;
% 
% for nbrUnit = 1:size(TunedChannelsPerPhase,3)
%     disp(num2str(nbrUnit))
%    
%     %extract only the grasp with the lowest p value to better compare
%         %to anova tuned units 
%     
%     if nnz( p_per_phase(2,:,nbrUnit) < alpha) > 1
%         
%         CueTuning =  p_per_phase(2,:,nbrUnit) == min(p_per_phase(2,:,nbrUnit));
%     else
%         CueTuning = p_per_phase(2,:,nbrUnit) < alpha;
%         
%     end 
%     
%     if nnz( p_per_phase(4,:,nbrUnit) < alpha) > 1
%         %extract only the grasp with the lowest p value to better compare
%         %to anova tuned units 
%         ActionTuning =  p_per_phase(4,:,nbrUnit) == min(p_per_phase(4,:,nbrUnit));
%     else
%         ActionTuning = p_per_phase(4,:,nbrUnit) < alpha;
%         
%      end 
% 
%     %take into account units that are tuned to both 
%     CueTuning = TunedChannelsPerPhase(2,:,nbrUnit);
%     ActionTuning = TunedChannelsPerPhase(4,:,nbrUnit);
%     
%     if (isequal(CueTuning, ActionTuning) && nnz(CueTuning)~= 0 && nnz(ActionTuning)~=0)
%         VisuoMotor = VisuoMotor + 1; 
%     elseif (~isequal(CueTuning, ActionTuning) && nnz(CueTuning)~= 0 && nnz(ActionTuning)~=0)
%         CueTuning
%         ActionTuning
%         switchingUnits = switchingUnits +1;
%        
%     elseif  (~isequal(CueTuning, ActionTuning) && nnz(CueTuning)== 0 && nnz(ActionTuning)~=0)
%         ActionUnits = ActionUnits +1;
%     elseif  (~isequal(CueTuning, ActionTuning) && nnz(CueTuning)~= 0 && nnz(ActionTuning)==0)
%         VisualUnits = VisualUnits +1;
% 
%     elseif  (isequal(CueTuning, ActionTuning) && nnz(CueTuning)== 0 && nnz(ActionTuning)==0)
%         NotTuned = NotTuned +1; 
%     else
%         error('Should not occur')
%     end
%     
%     
% end 

%numbers_cue2 = [VisualUnits, ActionUnits, VisuoMotor, switchingUnits];
%Labels_pie_2 = {'Cue';'Action'; 'Visuomotor'; 'Switching'};
%utile.plot_pie_chart(numbers_cue2, Labels_pie_2, Figure_labels,sum(numbers_cue2)); 

%We combine the units that are tuned in either of the phases and use it for
%classification
if combineTunedChannels
    %if removeITItuning
    %    TunedCombinedChannels = logical(sum(p_per_phase(:,2:end) < 0.05,2));
    %else
        TunedCombinedChannels = logical(sum(TunedChannelsPerPhase,2));
    %end 
    
end 
% 
% BoolNames = {'False', 'True'};
% disp('Tuned channel parameters: '); 
% disp('-----------------------------------------------------------------')
% disp(['Multcompare: ' BoolNames{multipleCompare +1} ])
% disp(['combineTunedChannels: ' BoolNames{combineTunedChannels +1} ])
% disp(['Number of Tuned Combined Units: ' num2str(nnz(TunedCombinedChannels)) ])
% disp(['Number of ITI units: ' num2str(nnz(TunedChannelsPerPhase(:,1))) ])
% disp(['Number of Cue units: ' num2str(nnz(TunedChannelsPerPhase(:,2))) ])
% disp(['Number of Delay units: ' num2str(nnz(TunedChannelsPerPhase(:,3))) ])
% disp(['Number of Action units: ' num2str(nnz(TunedChannelsPerPhase(:,4))) ])

NumberOfUnitsTunedToSeveralGrasps = 0
end

