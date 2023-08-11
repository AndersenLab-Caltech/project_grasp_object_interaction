function [TableAll] = acousticOnsetAlignment(TableAll, varargin)
%Align the neural data to acoustic onset
%I will call this function once per turn -> there should only be data in
%the table that is from that session day
%[varargin,threshold] = Utilities.argkeyval('AcousticThreshold',varargin, 0.8);

Utilities.argempty(varargin);

numTrials = size(TableAll.Trials{1,1}  ,2);
FsFramework = 30000; %frequency of framework and audio data
FsNeuralData = 1/0.05; %frequency of extracted neural data

unique_sessions = unique(TableAll.session_date);

UnitCount_previousSession = 0; 

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
    %AcousticOnsetTimeIdx = zeros(numTrials,1);
    %AcousticOnsetTime = zeros(numTrials,1);

    Labels = Table.Cue_labels{1};
    UniqueLabels = unique(Labels);
    % write gui that lets you chose which trial to take for alignment to
    % get best results!
    
    %timeActionNeural = 
    AlignedActionPhase = cell(size(Table,1),8);
    AudioActionAligned = {};
    AudioAction = {};
    %loop through number of units classes
    for classNbr =1:length(UniqueLabels)
        %extract all trials of one class
        IdxAll = find(Labels == classNbr);
        %select one class as the template 
        %Template =  Table.Audio_raw{1,1}{IdxAll(1), 1};
        %Template =  Table.Audio_envelope{1,1}{IdxAll(1), 1};
        Template =  Table.Audio_envelope_subsampled{1,1}{IdxAll(1), 1};

        %timeFW = (1:length(Template))/FsFramework;
        timeFW = (1:length(Template))/FsNeuralData;

        %time of action phase
        timeActionFW = timeFW(timeFW > ActionPhase);
        %audio data of action phase
        AudioTemplateAction = Template(timeFW > ActionPhase);
        
        %loop through each individual trial to align them to template
        for nbrRep = 1:length(IdxAll)
           % individual trial data
            %AudioTmp = Table.Audio_raw{1,1}{IdxAll(nbrRep), 1};
            %AudioTmp = Table.Audio_envelope{1,1}{IdxAll(nbrRep), 1};
            AudioTmp = Table.Audio_envelope_subsampled{1,1}{IdxAll(nbrRep), 1};

            %timeFWTmp = (1:length(AudioTmp))/FsFramework;
            timeFWTmp = (1:length(AudioTmp))/FsNeuralData;

            AudioTmpAction = AudioTmp(timeFWTmp > ActionPhase);
            
            %subplot(1,2,1)
            %plot(timeActionFW,AudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal')
             
            %calculate the number of datapoints and the time (in s) at which the template and the current trial are the most similar 
            delayDataPoints(classNbr,nbrRep) = finddelay(AudioTemplateAction, AudioTmpAction);
            %delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsFramework;
            delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsNeuralData;

            
          %  AudioAction{classNbr,nbrRep} = AudioTmpAction(1:44281)'; %make sure all audio datasets are the same length
            AudioAction{classNbr,nbrRep} = AudioTmpAction'; %make sure all audio datasets are the same length

            %Calculate the delay for optimal alignment
            
            %AudioActionAlignedTmp = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsFramework);
            AudioActionAlignedTmp = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsNeuralData);

            %AudioActionAligned{classNbr,nbrRep} = AudioActionAlignedTmp(1:44281);      
            AudioActionAligned{classNbr,nbrRep} = AudioActionAlignedTmp;

%             if delayTime(classNbr,nbrRep) < 0
%                 s1 = alignsignals(AudioTmpAction,AudioTemplateAction); 
%                 subplot(1,2,2)
%                 plot(s1); hold on; plot(AudioTemplateAction + 0.75); title('Aligned signal');legend('AudioTmp', 'Template')
%             else
%                 s1 = alignsignals(AudioTemplateAction,AudioTmpAction);
%                 subplot(1,2,2)
%                 plot(AudioTmpAction); hold on; plot(s1 +0.75);title('Aligned signal'); legend('AudioTmp', 'Template')
%             end
            
            
            %figure(); detectSpeech(AudioTmpAction', FsFramework)
            
            %detect speech and returns start and end index of it; returns Lx2 array, where L is the number of differnet speech bubbles there are
            %speechIdx{nbrRep} = detectSpeech(AudioTmpAction', FsFramework);
            %speechIdx{nbrRep} = detectSpeech(AudioTmpAction', FsNeuralData);

            %align data according to start of speech index
            %AudioActionAligned2{classNbr,nbrRep} = [zeros( FsFramework/3,1);AudioTmpAction(speechIdx{nbrRep}(1):44281)'] ;
            

%             figure();
%             plot(timeActionFW(1:44281),AudioTmpAction(1:44281)); hold on; 
%             xline(timeActionFW(speechIdx{nbrRep}(1) - .05*FsFramework)); hold on;
%             %xline(timeActionFW(speechIdx{nbrRep}(2) - .05*FsFramework)); hold on;
%             xline(timeActionFW(speechIdx{nbrRep}(1) - .05*FsFramework)); hold on;
% 
%             try
%                 xline(timeActionFW(speechIdx{nbrRep}(1) + 0.5*FsFramework))
%             catch
%                 xline(timeActionFW(end))
% 
%             end
            
%              figure();
%              plot(timeActionFW(1:44281),AudioTmpAction(1:44281)); hold on; 
%              xline(timeActionFW(speechIdx{nbrRep}(1))); hold on;
%              legend('Audio', 'speech detection')
%             xline(timeActionFW(speechIdx{nbrRep}(2,2))); hold on;

         

%             figure();
%             subplot(1,2,1)
%             plot(timeActionFW,AudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal'); legend('tmp', 'template');
%             subplot(1,2,2)
%              AdaptedAudioTmpAction = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsFramework);
%             plot(timeActionFW,AdaptedAudioTmpAction); hold on; plot(timeActionFW,AudioTemplateAction +0.75); title('Original signal')
%            

            
            
            

            
            %correlate Neural data 
            
            

            NameIdx = classNbr + 4;
            for unitNbr = 1:size(Table,1)
                %NeuralAction(unitNbr, :) =  Table{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);
                NeuralAction{unitNbr, classNbr}(:,nbrRep) =  Table{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);

                AlignedActionPhase{unitNbr,classNbr}(:,nbrRep) =delayseq(NeuralAction{unitNbr, classNbr}(:,nbrRep) , -delayTime(classNbr,nbrRep), FsNeuralData);
                %NeuralDelay = delayseq(NeuralTemplateAction, -delayTime(classNbr,nbrRep), 20);

            end 

        end 
        
    end 
    
%     figure();
%     subplot(2,1,1)
%     plot(AlignedActionPhase{14, 1} )
%     subplot(2,1,2)
%     plot(AlignedActionPhase{14, 1} )

    TableAll = [TableAll, AlignedActionPhase];
    
    
    %PLOT ALIGNEMNT 
    
%     for i = 1:length(unique(Labels))
%         figure(); 
% 
%         %subplot(4,2,2*i-1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
%         %subplot(4,2,2*i); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')
%         subplot(2,2,1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
%         subplot(2,2,2); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')
%         
% %         subplot(2,2,3); plot(cell2mat(AudioAction(i,:))); title('Audio data')
% %         subplot(2,2,4); 
% %         for j = 1:8
% %             plot(AudioActionAligned2{i,j}); hold on; title('Audio data speech detection aligned')       
% %         end 
%          
%     end 
%     
    
    
%     for i = 1:8
%         figure()
%         for j = 1:8
%         plot(AudioActionAligned2{i,j}); hold on;        
%         end 
%                  
%     end 
    
%     %correlate neural data and look at calculated delay -> does not work
%     %well lol 
% 
%     for unitNbr = 1:size(Table,1)
% 
%         for classNbr = 1:length(UniqueLabels)
% 
%             NameIdx = find(ismember(Table.Properties.VariableNames, preproc.image2class_simple(classNbr)));
%             NeuralAction =  Table{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,1);
% 
%             for nbrRep = 1:8
%                 NeuralTmpAction = Table{unitNbr,NameIdx}{1,1}(ActionPhaseTimes,nbrRep);                        
%                 delayTimeNeural(unitNbr,classNbr,nbrRep) =finddelay(NeuralAction, NeuralTmpAction)*0.05;
% 
%             end 
% 
%         end 
% 
%     end 
%     
%    MeanNeuralDelay  = squeeze(mean(delayTimeNeural,1));


%     audioEnvTmp =  Table.Audio_raw{1,1}{trialNbr, 1};
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

 

