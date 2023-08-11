function [TableAll] = acousticOnsetAlignment2(TableAll, varargin)
%Align the neural data to acoustic onset
%I will call this function once per turn -> there should only be data in
%the table that is from that session day
%[varargin,threshold] = Utilities.argkeyval('AcousticThreshold',varargin, 0.8);

Utilities.argempty(varargin);

numTrials = size(TableAll.Trials{1,1}  ,2);
FsFramework = 30000; %frequency of framework and audio data
FsNeuralData = 1/0.05; %frequency of extracted neural data

unique_sessions = unique(TableAll.session_date);


time_phase_labels = TableAll.time_phase_labels{1};
%Find time start of cue phase (of 1 trial)
CuePhase = find(time_phase_labels == 2,1)*0.05;
%find time start of action phase 
ActionPhase = find(time_phase_labels == 6,1)*0.05;
%find indexes of action phase
ActionPhaseTimes = find(TableAll.time_phase_labels{1,1} == 6);
timeND =(1:length(TableAll.Trials{1,1}))/FsNeuralData;

%TimeOnsetActionPhase = timeND(ActionPhaseTimes(1));
for nbr_sessions = 1:length(unique_sessions)
    session_idx = ismember(TableAll.session_date,unique_sessions(nbr_sessions));
    
    Table = TableAll(session_idx,:);
    Labels = Table.Cue_labels{1};
    UniqueLabels = unique(Labels);
    
    AlignedActionPhase = cell(size(Table,1),8);
    AudioActionAligned = {};
    AudioAction = {};
    %figure();
    %loop through number of units classes
    figure();
    for classNbr =1:length(UniqueLabels)
        %extract all trials of one class
        IdxAll = find(Labels == classNbr);
        %select one class as the template 
        Template =  Table.Audio_envelope_subsampled{1,1}{IdxAll(1), 1};
        timeND = (1:length(Template))/FsNeuralData;

        %audio data of action phase
        AudioTemplateAction = Template(timeND > ActionPhase);
        
        AlignToStartTime = find((AudioTemplateAction' > 0.005),1)/FsNeuralData - 0.25;
        %detectSpeech(AudioTemplateAction', FsNeuralData,'Window',2)
        try
            AudioTemplateAction = delayseq(AudioTemplateAction',-AlignToStartTime, FsNeuralData);
        catch
            error('something is off')
            Template =  Table.Audio_envelope_subsampled{1,1}{IdxAll(3), 1};
            AudioTemplateAction = Template(timeND > ActionPhase);
            AlignToStartTime = find((AudioTemplateAction' > 0.005),1)/FsNeuralData - 0.25;
            AudioTemplateAction = delayseq(AudioTemplateAction',-AlignToStartTime, FsNeuralData);

        end 

%         subplot(8,1,classNbr)
%         %plot(AudioTemplateAction); hold on; plot((AudioTemplateAction' > 0.005)*max(AudioTemplateAction))
%         plot(Template); hold on; plot((AudioTemplateAction' > 0.005)*max(AudioTemplateAction))
%         for i = 1:6
%             %plot(Template)
%             xline(find(time_phase_labels == i,1)) %, 'k')
%         end 

   %end 
        %loop through each individual trial to align them to template
        for nbrRep = 1:length(IdxAll)
            %individual trial data
            AudioTmp = Table.Audio_envelope_subsampled{1,1}{IdxAll(nbrRep), 1};
            timeNDTmp = (1:length(AudioTmp))/FsNeuralData;
            AudioTmpAction = AudioTmp(timeNDTmp > ActionPhase);
                      
            %calculate the number of datapoints and the time (in s) at which the template and the current trial are the most similar 
            delayDataPoints(classNbr,nbrRep) = finddelay(AudioTemplateAction, AudioTmpAction);
            delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsNeuralData;
            
            AudioAction{classNbr,nbrRep} = AudioTmpAction'; %make sure all audio datasets are the same length

            %Calculate the delay for optimal alignment
            AudioActionAlignedTmp = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsNeuralData);
            NameIdx = classNbr + 4;

             for unitNbr = 1:size(Table,1)               
                NeuralAction{unitNbr, classNbr}(:,nbrRep) =  Table{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);
                %AlignedActionPhase{unitNbr,classNbr}(:,nbrRep) =delayseq(NeuralAction{unitNbr, classNbr}(:,nbrRep) , delayTime(classNbr,nbrRep), FsNeuralData);
                AlignedActionPhase{unitNbr,classNbr}(:,nbrRep) =delayseq(NeuralAction{unitNbr, classNbr}(:,nbrRep) , -delayTime(classNbr,nbrRep), FsNeuralData);

             end 
            
            AudioActionAligned{classNbr,nbrRep} = AudioActionAlignedTmp;
            
            %figure(); plot(AudioActionAlignedTmp)

        end 
        
    end 
    
    TableAll = [TableAll, AlignedActionPhase];
    
   %PLOT AUDIO ALIGNEMNT 

%     for i = 1:length(unique(Labels))
%             figure(); 
% 
%         %subplot(4,2,2*i-1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
%         %subplot(4,2,2*i); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')
%         subplot(1,2,1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
%         subplot(1,2,2); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')
%         
%         
%          
%     end 
    
   ClassNames = arrayfun(@(x) preproc.image2class_simple(x), unique(Labels), 'UniformOutput', false);
    Colors = utile.get_color_rgb_codes(ClassNames);
 
    fig = figure('units','normalized','outerposition',[0.35 0.05 0.4 0.8]);
    
    for i = 1:length(unique(Labels))

        subplot(8,2,2*i-1); plot(cell2mat(AudioAction(i,:))); title('Unaligned Audio data')
        subplot(8,2,2*i); plot(cell2mat(AudioActionAligned(i,:))); 
        title([ClassNames{i} ' aligned'])
        %subplot(1,2,1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
        %subplot(1,2,2); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')        
    end
    sgtitle(unique(TableAll.session_date))
    disp('Acoustic alignment finished')
    
    %speechIdx = detectSpeech(AudioActionAligned{1,1}, FsNeuralData)

    
    
    %Plot Units
    
%     for d = 1:nnz(Table.nsp ==1) % all SMG units, the other ones don't do much
%         figure();   
% 
%         FigTitles = {'Unaligned Neural Data', 'Aligned Neural Data'};
%         for j = 1:2
%                 subplot(2,1,j)
%             for i = 1:8      
%                 if j== 1
%                     data_tmp = NeuralAction{d,i}';
%                 else
%                     data_tmp = AlignedActionPhase{d,i}';
%                 end 
%                 ci = bootci(1000, {@mean,data_tmp});
%                 Mean_FR = squeeze(mean(data_tmp));
% 
%                 err_ci(1,:) = ci(2,:) - Mean_FR; 
%                 err_ci(2,:) = Mean_FR - ci(1,:);  
% 
%                 timer_period = 0.05; 
%                 time_idx = ([0:length(Mean_FR)-1]*timer_period); 
% 
%                 ER =utile.shadedErrorBar(time_idx,Mean_FR, err_ci,'lineprops',{'Color', Colors{i}});
%                 hold on         
%                 
%                 ylim([0 15])
%             end         
%             title(FigTitles{j})
%         end 
%     end 
                
     
    
    
    
end

 

