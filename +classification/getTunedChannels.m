function [TunedCombinedChannels, TunedChannelsPerPhase, TunedChannelsPerBin, sumPhase, sumBin] = getTunedChannels(Data,Labels, TimePhaseLabels,varargin)

% Takes as input a cell array and phase labels and returns a bool of tuned units


%default values
flagMultipleCompare = true;
flagBinperBin = false;
flagCombineTunedChannels = true;
flagRemoveITItuning = false;
flagShuffleTest = false;


% Loading optional arguments
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'multcompare'
            flagMultipleCompare = varargin{2};
        case 'binperbintuning'
            flagBinperBin = varargin{2};              
        case 'combinetunedchannels'
            flagCombineTunedChannels = varargin{2};             
        case 'removeitituning'
            flagRemoveITItuning = varargin{2};
        case 'flagshuffletest'
            flagShuffleTest = varargin{2};           
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
      varargin(1:2) = [];
end



numChannels = size(Data{1},2);
UniqueTmp = unique(TimePhaseLabels{1});
numPhases = length(UniqueTmp);
numBins = size(Data{1},1);
numUnits = size(Data{1},2);
p_per_phase = ones(numChannels, numPhases);
p_per_bin = ones(numChannels, numBins);

p_multcomp_per_phase = ones(numChannels, numPhases);
p_multcomp_per_bin = ones(numChannels, numBins);

for n_phase = 1:numPhases
    try
        DataPerPhase = cell2mat(arrayfun(@(x,y) mean(x{1,1}(y{:}== UniqueTmp(n_phase),:),1),Data,TimePhaseLabels, 'UniformOutput', false));
    catch
        keyboard
    end 
        
 if flagShuffleTest
     Labels = Labels(randperm(length(Labels)));
 end 
    
  for n_channel = 1:numChannels
      
      DataPerTrial = DataPerPhase(:,n_channel);
      %perform kruskal wallis test 
      [p_per_phase(n_channel,n_phase), ~, ~] = kruskalwallis(DataPerTrial,Labels,'off');

  end
  
       
end 

if flagBinperBin

    for n_bin = 1:numBins
        disp([ 'Bin nbr ' num2str(n_bin)]);
        DataPerBin = cell2mat(arrayfun(@(x,y) Data{x,1}(n_bin,:),1:size(Data,1), 'UniformOutput', false)');

        for n_channel = 1:numChannels

          DataPerBinTrial = DataPerBin(:,n_channel);
          %perform kruskal wallis test 
          [p_val_bin, ~, ~] = kruskalwallis(DataPerBinTrial,Labels,'off');
          
          p_per_bin(n_channel,n_bin) = p_val_bin;

        end

    end 

end 

if flagMultipleCompare
    for n_channel = 1:numChannels
        [~,~,p_multcomp_per_phase(n_channel,:)] = utile.MultipleComparisonsCorrection(p_per_phase(n_channel,:),'method', 'fdr'); 

    end 


%DELETE
  %  for n_unit = 1:numUnits
   %     [~,~,p_multcomp_per_bin(:,n_unit)] = utile.MultipleComparisonsCorrection(p_per_bin(:,n_unit)','method', 'fdr');
    %end 
    %replace old p values by multcompare p values
    p_per_phase = p_multcomp_per_phase; 
   % p_per_bin = p_multcomp_per_bin;


end 

TunedChannelsPerPhase = p_per_phase < 0.05;
sumPhase = sum(TunedChannelsPerPhase);

TunedChannelsPerBin = p_per_bin < 0.05;
sumBin = sum(TunedChannelsPerBin); 
%We combine the units that are tuned in either of the phases and use it for
%classification
if flagCombineTunedChannels
    if flagRemoveITItuning
        TunedCombinedChannels = logical(sum(p_per_phase(:,2:end) < 0.05,2));
    else
        TunedCombinedChannels = logical(sum(p_per_phase < 0.05,2));
    end 
    
end 

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


end

