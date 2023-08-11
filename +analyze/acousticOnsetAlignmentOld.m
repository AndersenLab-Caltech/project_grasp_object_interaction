function [Table] = acousticOnsetAlignment(Table, varargin)
%Align the neural data to acoustic onset
%I will call this function once per turn -> there should only be data in
%the table that is from that session day
[varargin,threshold] = Utilities.argkeyval('AcousticThreshold',varargin, 0.8);

Utilities.argempty(varargin);



numTrials = size(Table.Trials{1,1}  ,2);
FsFramework = 30000;
FsNeuralData = 1/0.05;

unique_sessions = unique(Table.session_date);

UnitCount_previousSession = 0; 

time_phase_labels = Table.time_phase_labels{1};
CuePhase = find(time_phase_labels == 2,1)*0.05;
ActionPhase = find(time_phase_labels == 6,1)*0.05;

for nbr_sessions = 1:length(unique_sessions)
    session_idx = ismember(Table.session_date,unique_sessions(nbr_sessions));
    
    Table_date = Table(session_idx,:);
    AcousticOnsetTimeIdx = zeros(numTrials,1);
    AcousticOnsetTime = zeros(numTrials,1);

    %extract the exact Auditory onset times
    for trialNbr = 1:numTrials

        audioEnvTmp = Table_date.Audio_envelope{1,1}{trialNbr, 1};
        %audioRawTmp = Table_date.Audio_raw{1,1}{trialNbr, 1};
        audioEnvTmp =  Table_date.Audio_raw{1,1}{trialNbr, 1};

        NeuralTmp = mean((Table_date.Trials{trialNbr,1}),2);
        timeFW = (1:length(audioEnvTmp))/FsFramework;

        timeND =(1:length(Table_date.Trials{1,1}))/FsNeuralData;
        %find idx of Action phase
        ActionPhaseTimes = find(Table_date.time_phase_labels{1,1} == 6);
        TimeOnsetActionPhase = timeND(ActionPhaseTimes(1));

%         figure();
%         plot(timeFW, audioEnvTmp);
%         hold on 
%     %     plot(timeND,NeuralTmp)
%     %     hold on 
%         yL = get(gca,'YLim');
%         line([TimeOnsetActionPhase, TimeOnsetActionPhase], yL,'Color','k','LineStyle','--','Linewidth', 2);
%         sound(audioRawTmp, FsFramework)

            %take the max after the action phase started
       
            
        BackgroundNoise = audioEnvTmp((timeFW > 2) & (timeFW < 4));
        minBackground = min(BackgroundNoise);
        maxBackground = max(BackgroundNoise);
        meanBackground = mean(BackgroundNoise);
        
       % figure();findchangepts(audioEnvTmp,'MaxNumChanges',10)
        
        figure(); plot(timeFW,(audioEnvTmp > maxBackground)*0.4);
        hold on;  plot(timeFW, audioEnvTmp); hold on;
        xline(ActionPhase, 'k'); hold on;
        xline(CuePhase, 'k'); hold on;

        
        %figure(); plot(timeFW, audioEnvTmp);    
            
%
    end 

    close all
    
%     audioTmp = Table_date.Audio_raw{1,1}{1, 1};
%     timeFW = (1:length(audioEnvTmp))/FsFramework;
%     ActionAudio = audioTmp(timeFW > ActionPhase);
%     timeActionFW = timeFW(timeFW > ActionPhase);
%     
    Labels = Table_date.Cue_labels{1};
    UniqueLabels = unique(Labels);
%     % write gui that lets you chose which trial to take for alignment to
    % get best results!
    
    %timeActionNeural = 
    AlignedActionPhase = cell(size(Table,1),8);
    for classNbr =1 %:length(UniqueLabels)
        IdxAll = find(Labels == classNbr);
        Template =  Table_date.Audio_raw{1,1}{IdxAll(1), 1};
        timeFW = (1:length(Template))/FsFramework;
        timeActionFW = timeFW(timeFW > ActionPhase);

        AudioTemplateAction = Template(timeFW > ActionPhase);

        for nbrRep = 1:length(IdxAll)
           % figure();

            AudioTmp = Table_date.Audio_raw{1,1}{IdxAll(nbrRep), 1};
            timeFWTmp = (1:length(AudioTmp))/FsFramework;
            AudioTmpAction = AudioTmp(timeFWTmp > ActionPhase);
            
            %subplot(1,2,1)
            %plot(timeActionFW,AudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal')
             
            delayDataPoints(classNbr,nbrRep) = finddelay(AudioTemplateAction, AudioTmpAction);
            delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsFramework;

%             if delayTime(classNbr,nbrRep) < 0
%                 s1 = alignsignals(AudioTmpAction,AudioTemplateAction); 
%                 subplot(1,2,2)
%                 plot(s1); hold on; plot(AudioTemplateAction + 0.75); title('Aligned signal');legend('AudioTmp', 'Template')
%             else
%                 s1 = alignsignals(AudioTemplateAction,AudioTmpAction);
%                 subplot(1,2,2)
%                 plot(AudioTmpAction); hold on; plot(s1 +0.75);title('Aligned signal'); legend('AudioTmp', 'Template')
%             end
            
            
%             figure();
%             subplot(1,2,1)
%             plot(timeActionFW,AudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal'); legend('tmp', 'template');
%             subplot(1,2,2)
%             AdaptedAudioTmpAction = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsFramework);
%             plot(timeActionFW,AdaptedAudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal')

            %correlate Neural data 
            for unitNbr = 1:size(Table_date,1)
                NameIdx = classNbr + 4;
                NeuralAction =  Table_date{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);
                
                AlignedActionPhase{unitNbr,classNbr}(:,nbrRep) =delayseq(NeuralAction, -delayTime(classNbr,nbrRep), 20);
                %NeuralDelay = delayseq(NeuralTemplateAction, -delayTime(classNbr,nbrRep), 20);

            end 

        end 
        
    end 
    
    %correlate neural data and look at calculated delay -> does not work
    %well lol 

    for unitNbr = 1:size(Table_date,1)

        for classNbr = 1:length(UniqueLabels)

            NameIdx = find(ismember(Table_date.Properties.VariableNames, preproc.image2class_simple(classNbr)));
            NeuralAction =  Table_date{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,1);

            for nbrRep = 1:8
                NeuralTmpAction = Table_date{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);                        
                delayTimeNeural(unitNbr,classNbr,nbrRep) =finddelay(NeuralAction, NeuralTmpAction)*0.05;

            end 

        end 

    end 
    
    MeanNeuralDelay  = squeeze(mean(delayTimeNeural,1));


%     audioEnvTmp =  Table_date.Audio_raw{1,1}{trialNbr, 1};
%     ActionTmp = audioEnvTmp(timeFW>ActionPhase);
%     
%     finddelay(ActionAudio, ActionTmp)
%     
%     s1 = alignsignals(ActionAudio,ActionTmp);
%     
%     
%     figure();
%     plot(ActionAudio); hold on; plot(ActionTmp);
%         
%     figure();
%     plot(s1); hold on; plot(ActionTmp);

    end
    
    
    
    
    
    %compute the IDx for neural data from acoustic onset
    [~, AcOnIdxNeuralData] = arrayfun(@(n) min(abs(timeND-n)), AcousticOnsetTime);

    ExtractionLength = 2*FsNeuralData; %extract 2 seconds of neural data after acoustic onset
    %Check: sometimes, if the subject said the world VERY late, there are
    %no 2 seconds left AFTER the end of the trial. In that case, we change
    %the Start Index of that trial to the end time - 2 seconds
    LateSpeakingTrials = find(AcOnIdxNeuralData + ExtractionLength > size(Table_date.Trials{1,1},1));   
    if ~isempty(LateSpeakingTrials)
        disp([' Trial(s) ' num2str(LateSpeakingTrials') ' started speaking late - adapt start time' ])
        AcOnIdxNeuralData(LateSpeakingTrials) = size(Table_date.Trials{1,1},1) - ExtractionLength;
    end 
    %extract 2 seconds from acoustinc onset for Action phase Data
    EndIdx = AcOnIdxNeuralData + ExtractionLength;

    newAction = {};
    Trials_AcousticOnset = {};
    %time_phase_labels_acousticOnset = {};

    ClassNames_c = arrayfun(@(x) preproc.image2class_simple(x), unique(Table.Cue_labels{1,1}),'UniformOutput',false); 

    for UnitNbr = 1:size(Table_date,1)
        %For each unit go through each trial and according to the right
        %Acoustic Onset extract new Action phase
        for trialNbr = 1:numTrials
            
            Trials_AcousticOnset{UnitNbr,1}(:,trialNbr) = Table_date.Trials{UnitNbr,1}(~(Table.time_phase_labels{1,1} == 4),trialNbr); %I'm taking on purpose the time phase labels of the first session so that it is the same for all of them, otherwise there are differences for different sessions... 
            newAction{UnitNbr,1}(:,trialNbr) = Table_date.Trials{UnitNbr,1}(AcOnIdxNeuralData(trialNbr):EndIdx(trialNbr),trialNbr);
            %Trials_new{UnitNbr,1}(:,trialNbr) = [];
        end 
        %Combine the new 'total phase' 
        Trials_AcousticOnset{UnitNbr,1} = [Trials_AcousticOnset{UnitNbr,1}; newAction{UnitNbr,1}];
        
        %Replace the time phase labels. In order to replace in the old
        %table, we need to know the total number of units of the previous
        %phase !! 
        
        %disp(['Processing unit ' num2str(UnitNbr + UnitCount_previousSession) ])
        Table(UnitNbr + UnitCount_previousSession, 'time_phase_labels') = {[Table.time_phase_labels{1,1}(~(Table.time_phase_labels{1,1} == 4)); 4*ones(size(newAction{UnitNbr,1},1),1)]};

        %Replace the grasp info
        Table(UnitNbr + UnitCount_previousSession,ClassNames_c') = arrayfun(@(n) Trials_AcousticOnset{UnitNbr,1}(:, Table.Cue_labels{1,1}  == n), unique(Table.Cue_labels{1,1}), 'UniformOutput', false)';
    end 

    UnitCount_previousSession = find(session_idx == 0,1) -1;
    %Classes = arrayfun(@(n) Trials_AcousticOnset(

    figure();
    for i = 1:trialNbr
        plot(timeFW, Table_date.Audio_envelope{1,1}{i, 1}, 'Color', cell2mat(utile.get_color_rgb_codes({preproc.image2class_simple(Table_date.Cue_labels{1,1}(i))})));
        hold on 
        yL = get(gca,'YLim');
        b = line([TimeOnsetActionPhase, TimeOnsetActionPhase], yL,'Color','k','LineStyle','--','Linewidth', 2);
    end 
    title('Acoustic envelope for each trial')
    legend(b, 'Action phase start')
    xlabel('Time')
    ylabel('Acoustic Sound');
end

 

