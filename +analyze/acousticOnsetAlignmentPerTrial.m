function [TableAll] = acousticOnsetAlignmentPerTrial(TableAll, varargin)
%Align the neural data to acoustic onset
%I will call this function once per turn -> there should only be data in
%the table that is from that session day
%[varargin,threshold] = Utilities.argkeyval('AcousticThreshold',varargin, 0.8);

%parameter that defines if we align or not align the data
[varargin,flagUnalign] = Utilities.argkeyval('unalign',varargin, 0 );

Utilities.argempty(varargin);

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
SMG_GoAlignedAll = [];
PMV_GoAlignedAll = [];
S1X_GoAlignedAll = [];
Audio_GoAlignedAll = [];
%TimeOnsetActionPhase = timeND(ActionPhaseTimes(1));
for nbr_sessions = 1:length(unique_sessions)
    session_idx = ismember(TableAll.session_date,unique_sessions(nbr_sessions));
    
    Table = TableAll(session_idx,:);
    Labels = Table.GoLabels;
    UniqueLabels = unique(Labels);
    SMG_GoAligned = cell(size(Table.SMG_Go));
    PMV_GoAligned = cell(size(Table.PMV_Go));
    S1X_GoAligned = cell(size(Table.S1X_Go));
    Audio_GoAlined = cell(size(Table.SMG_Go));

    AlignedActionPhase = cell(size(Table,1),8);
    AudioActionAligned = {};
    AudioAction = {};
    %figure();
    %loop through number of units classes
   % figure();
    for classNbr =1:length(UniqueLabels)
        %extract all trials of one class
        IdxAll = find(Labels == UniqueLabels(classNbr));
        %select one class as the template 
        Template =  Table.Audio_envelope_subsampled{IdxAll(1), 1};
        timeND = (1:length(Template))/FsNeuralData;

        %audio data of action phase
        AudioTemplateAction = Template(timeND > ActionPhase);
        
        %AlignToStartTime = find((normalize(AudioTemplateAction, 'range', [0 1])' > 0.8),1)/FsNeuralData - 0.25;
        AlignToStartTime = find((normalize(AudioTemplateAction, 'range', [0 1])' > 0.4),1)/FsNeuralData - 0.5;

        %AlignToStartTime = find(AudioTemplateAction,1 > 0.05)/FsNeuralData - 0.25;

        %detectSpeech(AudioTemplateAction', FsNeuralData,'Window',2)
        try
            AudioTemplateAction = delayseq(AudioTemplateAction',-AlignToStartTime, FsNeuralData);
        catch
            error('something is off')
            Template =  Table.Audio_envelope_subsampled{IdxAll(3), 1};
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
            AudioTmp = Table.Audio_envelope_subsampled{IdxAll(nbrRep), 1};
            timeNDTmp = (1:length(AudioTmp))/FsNeuralData;
            AudioTmpAction = AudioTmp(timeNDTmp > ActionPhase);
            AudioTmpAction = AudioTmpAction(1:30);
            %calculate the number of datapoints and the time (in s) at which the template and the current trial are the most similar 
            delayDataPoints(classNbr,nbrRep) = finddelay(AudioTemplateAction, AudioTmpAction);
            delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsNeuralData;
            
            AudioAction{classNbr,nbrRep} = AudioTmpAction'; %make sure all audio datasets are the same length
            
            %Calculate the delay for optimal alignment
            %AudioActionAlignedTmp = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsNeuralData);
            
            %we can unalign the data by a certain parameter. defined here
            if flagUnalign
                DelayPoint = delayDataPoints(classNbr,nbrRep);
            else
                DelayPoint = -delayDataPoints(classNbr,nbrRep);
            end
            AudioActionAlignedTmp = delayseq(AudioTmpAction',  DelayPoint);

             for unitNbr = 1:size(Table.SMG_Go{1,1},2)               
                NeuralAction =  Table.SMG_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
                %if Nan at the end it's okay. If Nan at the beginning ->
                %artefact. Pad with data from delay phase
                if delayDataPoints(classNbr,nbrRep) > 0
                    SMG_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
                else
                    %we pad the data with the last x values of the delay
                    %phase
                    NeuralDelay2 =  Table.SMG_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
                    Padding = NeuralDelay2(end-abs(delayDataPoints):end);
                    SMGTmp =utile.delayseq(NeuralAction, DelayPoint);
                    SMGTmp(1:length(Padding)) = Padding;
                    SMG_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =SMGTmp;
                    
                end 
             end 
             
             for unitNbr = 1:size(Table.PMV_Go{1,1},2)               
                NeuralAction =  Table.PMV_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
                
                if delayDataPoints(classNbr,nbrRep) > 0
                    PMV_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
                else
                    %we pad the data with the last x values of the delay
                    %phase
                    NeuralDelay2 =  Table.PMV_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
                    Padding = NeuralDelay2(end-abs(delayDataPoints):end);
                    PMVTmp =utile.delayseq(NeuralAction, DelayPoint);
                    PMVTmp(1:length(Padding)) = Padding;
                    PMV_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =PMVTmp;
                    
                end 
                
             end 
             
             for unitNbr = 1:size(Table.S1X_Go{1,1},2)               
                NeuralAction =  Table.S1X_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
                if delayDataPoints(classNbr,nbrRep) > 0
                    S1X_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
                else
                    %we pad the data with the last x values of the delay
                    %phase
                    NeuralDelay2 =  Table.S1X_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
                    Padding = NeuralDelay2(end-abs(delayDataPoints):end);
                    S1XTmp =utile.delayseq(NeuralAction, DelayPoint);
                    S1XTmp(1:length(Padding)) = Padding;
                    S1X_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =S1XTmp;
                    
                end
             end
            
            Audio_GoAlined{IdxAll(nbrRep),1} =  AudioActionAlignedTmp;
            AudioActionAligned{classNbr,nbrRep} = AudioActionAlignedTmp;
            
            %figure(); plot(AudioActionAlignedTmp)

        end 

    end 
    
    SMG_GoAlignedAll = [SMG_GoAlignedAll; SMG_GoAligned];
    PMV_GoAlignedAll = [PMV_GoAlignedAll; PMV_GoAligned];
    S1X_GoAlignedAll = [S1X_GoAlignedAll; S1X_GoAligned];
    Audio_GoAlignedAll = [Audio_GoAlignedAll; Audio_GoAlined];

    ClassNames = preproc.image2class_simple(unique(Labels));
    Colors = utile.get_color_rgb_codes(ClassNames);
 
    
    
%     %size test - turns out some are 30 and some 31 - shorten all
%     SizeAudio = cell2mat(cellfun(@(x) size(x,1), AudioAction, 'UniformOutput', false));
%     %if length(unique(SizeAudio)) ~= 1
%         %%% WORK FROM HERE%%%
%     minSize = min(min(SizeAudio));
%     warning(['shorten Audio data to have the same size :' num2str(minSize)])
%     try
% 
%         AudioAction = cellfun(@(x) x(1:30),AudioAction, 'UniformOutput', false);
%     catch
%         disp('blu')
%     end 
%     AudioActionAligned = cellfun(@(x) x(1:30), AudioActionAligned, 'UniformOutput', false);
        %% MAKE THEM ALL SHORTER 
        %e.g here issue is some are 30 and some are 31 
    %end 
    fig = figure('units','normalized','outerposition',[0.35 0.05 0.4 0.8]);

    for i = 1:length(unique(Labels))

        subplot(8,2,2*i-1); plot(cell2mat(AudioAction(i,:))); title('Unaligned Audio data')
        subplot(8,2,2*i); plot(cell2mat(AudioActionAligned(i,:))); 
        title([ClassNames{i} ' aligned'])
        %subplot(1,2,1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
        %subplot(1,2,2); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')        
    end
    sgtitle(unique_sessions(nbr_sessions))
    disp('Acoustic alignment finished')

    
    
    
end

SMG_GoAligned = SMG_GoAlignedAll;
PMV_GoAligned = PMV_GoAlignedAll;
S1X_GoAligned = S1X_GoAlignedAll;
Audio_SubAlinged = Audio_GoAlignedAll;

TableAll = [TableAll, cell2table(SMG_GoAligned),cell2table(PMV_GoAligned), cell2table(S1X_GoAligned), cell2table(Audio_SubAlinged)];
 
