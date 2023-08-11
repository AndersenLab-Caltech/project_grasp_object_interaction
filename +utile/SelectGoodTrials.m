function [TrialIdx] = SelectGoodTrials(Data,Labels,TimePhaseLabels,varargin)
% returns best trial idx of a cell of data 

%Select trials if they have a significant ttest between ITI and Action
%Phase. I want to order the trials from best to worst

%I select the 8 best trial indexes (the most channels have a significant
%ttest between ITI and action phase) and the 8 worst trial indexes (the
%lest channels that have a significant ttest between ITI and Action phase)

[varargin,flagBestTrial] = Utilities.argkeyval('BestTrial',varargin, true);
Utilities.argempty(varargin);

ITI_Data = cellfun(@(x,y) x(y == 1,:),Data, TimePhaseLabels, 'UniformOutput', false);
Action_Data = cellfun(@(x,y) x(y == 4,:),Data, TimePhaseLabels, 'UniformOutput', false);

fixedTrialNbr = unique(histc(Labels, unique(Labels))); 
H_Value = zeros(fixedTrialNbr,size(Data{1},2));
%H_Value = zeros(size(Data,1),size(Data{1},2));
UniqueLabels = unique(Labels);

for GraspNbr = 1:length(UniqueLabels)
    
    GraspType = UniqueLabels(GraspNbr);
    ITI_DataGrasp = ITI_Data(Labels == GraspType);
    Action_DataGrasp = Action_Data(Labels == GraspType);
    %keep the idx of the trials that we are extracting
    IdxGrasp = find(Labels == GraspType);
    
    for ChannelNbr = 1:size(Data{1},2)

        [H_Value(:,ChannelNbr) ]= cell2mat(cellfun(@(x,y) ttest2(x(:,ChannelNbr),y(:,ChannelNbr)),ITI_DataGrasp, Action_DataGrasp,'UniformOutput', false));

    end
    
    H_ValuePerTrial = nansum(H_Value,2); 
    [Val,Idx] = sort(H_ValuePerTrial);
    %Sort the idx of the trials from worst trial to best trial for that
    %specific grasp (I am doing it this way since I want to keep the same
    %number of best vs worst trials for each grasp to be able to compare
    %reliably)
    SortedIdx(:,GraspNbr) = IdxGrasp(Idx);
    
end 

%extract the worst trial index. the sort function sorts from smallest to
%biggest
TopSelection = 32; 

WorstTrialIdx = SortedIdx(1:TopSelection,:);
WorstTrialIdx = WorstTrialIdx(:);
WorstTrialsLabels = Labels(WorstTrialIdx);

%extract best trial index. 
BestTrialIdx = SortedIdx(end-(TopSelection-1):end,:);
BestTrialIdx = BestTrialIdx(:);

BestTrialLabels = Labels(BestTrialIdx);


if flagBestTrial %return best trial idx
    TrialIdx = 	BestTrialIdx; 
else %return worst trial idx
    TrialIdx = WorstTrialIdx; 
end 

end 