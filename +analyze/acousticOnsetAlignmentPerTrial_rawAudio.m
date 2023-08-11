function [TableAll] = acousticOnsetAlignmentPerTrial_rawAudio(TableAll, varargin)
%Align the neural data to acoustic onset
%I will call this function once per turn -> there should only be data in
%the table that is from that session day
%[varargin,threshold] = Utilities.argkeyval('AcousticThreshold',varargin, 0.8);

%parameter that defines if we align or not align the data
[varargin,flagUnalign] = Utilities.argkeyval('unalign',varargin, 0 );

Utilities.argempty(varargin);

FsFramework = 30000; %frequency of framework and audio data
FsNeuralData = 1/0.05; %frequency of extracted neural data

FsNeuralData = FsFramework;

unique_sessions = unique(TableAll.session_date);
flagAudioEnvelope = true; %if true, use envelope. if false, use raw audio

%Method:
% Aling data with subsampled audio envelope -> gives better results than
% using raw data
% Find delay sequences
% shift raw audio data based on these sequences


for nbr_sessions = 1:length(unique_sessions)
    session_idx = ismember(TableAll.session_date,unique_sessions(nbr_sessions));
    
    Table = TableAll(session_idx,:);

    timeNeuralData = (1:length(Table.time_phase_labels{1}))*0.05;
    ActionPhaseTimeStart = cellfun(@(x) timeNeuralData(find(x ==6,1)), Table.time_phase_labels, 'UniformOutput', false);

    timeAll = cellfun(@(x) (1:length(x))/FsNeuralData, Table.Audio_raw, 'UniformOutput', false);
    %timeActionIdx = cellfun(@(x) x > 6

    %find audio data that corresponds to action phase data. action phase start
    %is indicated by the variable ActionPhaseTimeStart. 
    %Time for audio data = timeAll. 
    
    %define which audio type to use for alignment
    if flagAudioEnvelope
        AudioDataType = Table.Audio_envelope;
    else
        AudioDataType = Table.Audio_raw;
    end
    
    AudioDataEnv = Table.Audio_envelope;
    AudioDataRaw = Table.Audio_raw;
    
    AudioActionEnv = cellfun(@(x,y,z) x(y > z), AudioDataEnv, timeAll,ActionPhaseTimeStart, 'UniformOutput', false);
    AudioActionRaw = cellfun(@(x,y,z) x(y > z), AudioDataRaw, timeAll,ActionPhaseTimeStart, 'UniformOutput', false);

    %find minimum size audio signal
    SizeAudioMin = min(cell2mat(cellfun(@(x) size(x,2), AudioActionEnv, 'UniformOutput', false)));

    % shorten audio signal to same length
    AudioActionEnv = cellfun(@(x) x(1:SizeAudioMin), AudioActionEnv, 'UniformOutput', false);
    AudioActionRaw = cellfun(@(x) x(1:SizeAudioMin), AudioActionRaw, 'UniformOutput', false);

    %downsample envelope audio signal to reduce complexity -> helps with
    %proper aligning
    ds_factor = 8; %downsample by 2 -> 15kHz 

    AudioActionDS = cellfun(@(x) downsample(x,ds_factor), AudioActionEnv, 'UniformOutput', false);

    %normalize signal to same maximum amplitude 
   % AudioActionDS_N = cellfun(@(x) x/max(abs(x)), AudioActionDS, 'UniformOutput', false);
    %zscore signal
    %AudioActionDS_N = cellfun(@(x) zscore(x), AudioActionDS, 'UniformOutput', false);
    AudioActionDS_N = AudioActionDS;
    
    %Downsampled Envelope data 
    AudioActionAlignedEnvDS = cell(size(AudioActionDS_N));
    %Envelope data 
    AudioActionAlignedEnv = cell(size(AudioActionEnv));
    %Raw data
    AudioActionAlignedRaw = cell(size(AudioActionRaw));

    %AudioRaw aligned to speech start
    AudioActionAlignedSpeechStart = cell(size(AudioActionRaw));

    
    %compute cross correlation between signals - perform this per class
    Labels = Table.GoLabels;
    ClassNamesAll = preproc.image2class_simple(Labels);
    UniqueLabels = unique(Labels);
    ClassNames = preproc.image2class_simple(unique(Labels));
    
    
    
    figure();

    n_samples = length(AudioActionRaw);
    
    %save('D:\Users\Sarah\HomonymPaper\AlignmentVariable\time_offset.mat', 'time_offset')

    %load annotated time offset for session 
    time_offset_data = load('D:\Users\Sarah\HomonymPaper\AlignmentVariable\time_offset_20230117_raw.mat');
    time_offset = time_offset_data.time_offset;
    
    
%     %define offset by hand: put timepoint at beginning of audio speech
%     %onset
%     for sample = 1:n_samples
%        time_offset = NaN(n_samples, 2);

%         audio_sample = AudioActionRaw{sample};
%         plot(1:length(audio_sample), audio_sample); % where time is in ms or seconds
%         
%         if exist('roi_point', 'var')
%             delete(roi_point);
%         end
%         
%         roi_point = drawpoint;
% 
%         % May have to round to nearest number depending on how you set this up
%         time_offset(sample,:) = roi_point.Position;
%     
%     end
    
    %keyboard
    
    AudioActionAlignedTest = cell(size(AudioActionDS_N));
    %align audio based on defined audio start in the ariable time_offset
    
    for sample = 1:n_samples
        audio_sample = AudioActionRaw{sample};
        AudioActionAlignedTest{sample} = delayseq(audio_sample',-time_offset(sample,1))';
       
    end
    
    
    

    
 
%     %plot to verify
%     figure(); 
%     a= 129;
%     b = a+15;
%     nbrToPlot = a:b; 
%     subplot(1,2,1)
%     for i = nbrToPlot
%         plot( AudioActionRaw{i} + i)
%         hold on;
%     end 
%     yticks(nbrToPlot)
%     yticklabels(ClassNamesAll(nbrToPlot))
%     title('unAligned')
%     
%      subplot(1,2,2)
%     for i = nbrToPlot
%         plot( AudioActionAlignedTest{i} +i)
%         hold on;
%     end 
%     title('Aligned')
    

    Audio_ActionRawAlinged = AudioActionAlignedTest;
    TableAll = [TableAll,  cell2table(Audio_ActionRawAlinged)];
%end 
    %%
    
    
    
%     for n_Audio = 1:10 %length(AudioActionRaw)
%         
%         figure();
%         %[a,b(n_Audio,:)] =detectSpeech(AudioActionRaw{n_Audio}', FsFramework)
%         detectSpeech(AudioActionRaw{n_Audio}', FsFramework,'Thresholds',[0.01, 0.1])
% 
%     end
    
   

    for classNbr = 1:length(UniqueLabels)
%        extract all trials of one class
        IdxAll = find(Labels == UniqueLabels(classNbr));
       % select one class as the template 
        Template =  AudioActionDS_N{IdxAll(1), 1};
        blub = {};

        blubAligned = {};
         for nbrRep = 1:length(IdxAll)
             AudioTmp = AudioActionDS_N{IdxAll(nbrRep), 1};

             
             
            % find delay using envelope of signal. Use delay to align audio
            % raw data. 
             delayIdx = finddelay(Template,  AudioTmp,FsFramework*0.1);
             
            % don't want to use cirshift as it just shifts the signal
            % (e.g. put datapoints from the end to the beginning, don't
           %  think it's ideal. Zero padding might be better? not sure.
             
            a = circshift(AudioTmp,-(delayIdx-length(AudioTmp)));
            b = delayseq(AudioTmp',-delayIdx);
       
             AudioActionAligned{IdxAll(nbrRep)} = circshift(AudioTmp,-(delayIdx-length(AudioTmp)));
             
            % Align downsamples audio envelope, envelope, and raw data 
             AudioActionAlignedEnvDS{IdxAll(nbrRep)} = delayseq(AudioTmp',-delayIdx)';
             AudioActionAlignedEnv{IdxAll(nbrRep)} = delayseq(AudioActionEnv{IdxAll(nbrRep), 1}',-delayIdx*ds_factor)';
             AudioActionAlignedRaw{IdxAll(nbrRep)} = delayseq(AudioActionRaw{IdxAll(nbrRep), 1}',-delayIdx*ds_factor)';

             
             blubAligned{nbrRep} = AudioActionAlignedRaw{IdxAll(nbrRep)};
             blub{nbrRep} = AudioActionRaw{IdxAll(nbrRep),1};
         end 
        
        figure();
        
        subplot(1,2,1); plot(cell2mat(blub')'); title('Unaligned Audio data')
        subplot(1,2,2); plot(cell2mat(blubAligned')');   title('Aligned Audio data')
        sgtitle(ClassNames{classNbr});
    

    end

end 
%AudioActionRaw
keyboard
 ouput_data = struct;
 output_data.audio_raw_aligned = cell2mat(AudioActionAlignedRaw);
 output_data.audio_env_aligned = cell2mat(AudioActionAlignedEnv);
 output_data.audio_raw = cell2mat(AudioActionRaw);
 output_data.audio_env = cell2mat(AudioActionEnv);
 output_data.audio_action_aligned_to_start = cell2mat(AudioActionAlignedTest);
 output_data.Fs = FsFramework;
 output_data.Labels = Labels;
 output_data.ClassNames =  ClassNames; 
 
%  save('D:\Users\Sarah\Code\Python\Paper3Homophones\Data.mat', 'output_data')
%  
%  
%  for i = 1:length(AudioActionAlignedTest)
%      tmp_sound = AudioActionAlignedTest{i};
%      
%      sound(tmp_sound, FsFramework)
%  
%  end 
 

%end 
% 
% for i = 1:24
%     figure()
%     plot(blub{i});
%     hold on;
%     plot(blubAligned{i} + 5 )
%     hold on;
%     plot(Template +10)
%     legend('NonAligned', 'Aligned', 'Template')
% end 
% 
% 
% %Find time start of cue phase (of 1 trial)
% CuePhase = find(time_phase_labels == 2,1)*0.05;
% %find time start of action phase 
% ActionPhase = find(time_phase_labels == 6,1)*0.05;
% %find indexes of action phase
% ActionPhaseTimes = find(TableAll.time_phase_labels{1,1} == 6);
% SMG_GoAlignedAll = [];
% PMV_GoAlignedAll = [];
% S1X_GoAlignedAll = [];
% Audio_GoAlignedAll = [];
% %TimeOnsetActionPhase = timeND(ActionPhaseTimes(1));
% for nbr_sessions = 1:length(unique_sessions)
%     session_idx = ismember(TableAll.session_date,unique_sessions(nbr_sessions));
%     
%     Table = TableAll(session_idx,:);
%     Labels = Table.GoLabels;
%     UniqueLabels = unique(Labels);
%     SMG_GoAligned = cell(size(Table.SMG_Go));
%     PMV_GoAligned = cell(size(Table.PMV_Go));
%     S1X_GoAligned = cell(size(Table.S1X_Go));
%     Audio_GoAlined = cell(size(Table.SMG_Go));
% 
%     AlignedActionPhase = cell(size(Table,1),8);
%     AudioActionAlignedEnvDS = {};
%     AudioActionEnv = {};
%     %figure();
%     %loop through number of units classes
%    % figure();
%     for classNbr =1:length(UniqueLabels)
%         %extract all trials of one class
%         IdxAll = find(Labels == UniqueLabels(classNbr));
%         %select one class as the template 
%         Template =  Table.Audio_raw{IdxAll(1), 1};
%         timeND = (1:length(Template))/FsNeuralData;
% 
%         %audio data of action phase
%         AudioTemplateAction = Template(timeND > ActionPhase);
%         
%         %AlignToStartTime = find((normalize(AudioTemplateAction, 'range', [0 1])' > 0.8),1)/FsNeuralData - 0.25;
%         AlignToStartTime = find((normalize(AudioTemplateAction, 'range', [0 1])' > 0.4),1)/FsNeuralData - 0.5;
% 
%         %AlignToStartTime = find(AudioTemplateAction,1 > 0.05)/FsNeuralData - 0.25;
% 
%         %detectSpeech(AudioTemplateAction', FsNeuralData,'Window',2)
%         try
%             AudioTemplateAction = delayseq(AudioTemplateAction',-AlignToStartTime, FsNeuralData);
%         catch
%             error('something is off')
%             Template =  Table.Audio_envelope_subsampled{IdxAll(3), 1};
%             AudioTemplateAction = Template(timeND > ActionPhase);
%             AlignToStartTime = find((AudioTemplateAction' > 0.005),1)/FsNeuralData - 0.25;
%             AudioTemplateAction = delayseq(AudioTemplateAction',-AlignToStartTime, FsNeuralData);
% 
%         end 
% 
% %         subplot(8,1,classNbr)
% %         %plot(AudioTemplateAction); hold on; plot((AudioTemplateAction' > 0.005)*max(AudioTemplateAction))
% %         plot(Template); hold on; plot((AudioTemplateAction' > 0.005)*max(AudioTemplateAction))
% %         for i = 1:6
% %             %plot(Template)
% %             xline(find(time_phase_labels == i,1)) %, 'k')
% %         end 
% 
%    %end 
%         %loop through each individual trial to align them to template
%         for nbrRep = 1:length(IdxAll)
%             %individual trial data
%             AudioTmp = Table.Audio_raw{IdxAll(nbrRep), 1};
%             timeNDTmp = (1:length(AudioTmp))/FsNeuralData;
%             AudioTmpAction = AudioTmp(timeNDTmp > ActionPhase);
%            % AudioTmpAction = AudioTmpAction(1:30);
%             %calculate the number of datapoints and the time (in s) at which the template and the current trial are the most similar 
%             delayDataPoints(classNbr,nbrRep) = finddelay(AudioTemplateAction, AudioTmpAction);
%             delayTime(classNbr,nbrRep) =delayDataPoints(classNbr,nbrRep)/FsNeuralData;
%             
%             AudioActionEnv{classNbr,nbrRep} = AudioTmpAction'; %make sure all audio datasets are the same length
%             
%             %Calculate the delay for optimal alignment
%             %AudioActionAlignedTmp = delayseq(AudioTmpAction', -delayTime(classNbr,nbrRep), FsNeuralData);
%             
%             %we can unalign the data by a certain parameter. defined here
%             if flagUnalign
%                 DelayPoint = delayDataPoints(classNbr,nbrRep);
%             else
%                 DelayPoint = -delayDataPoints(classNbr,nbrRep);
%             end
%             AudioActionAlignedTmp = delayseq(AudioTmpAction',  DelayPoint);
% 
% %              for unitNbr = 1:size(Table.SMG_Go{1,1},2)               
% %                 NeuralAction =  Table.SMG_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
% %                 %if Nan at the end it's okay. If Nan at the beginning ->
% %                 %artefact. Pad with data from delay phase
% %                 if delayDataPoints(classNbr,nbrRep) > 0
% %                     SMG_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
% %                 else
% %                     %we pad the data with the last x values of the delay
% %                     %phase
% %                     NeuralDelay2 =  Table.SMG_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
% %                     
% %                     Padding = NeuralDelay2(end-abs(delayDataPoints):end);
% %                     SMGTmp =utile.delayseq(NeuralAction, DelayPoint);
% %                     SMGTmp(1:length(Padding)) = Padding;
% %                     SMG_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =SMGTmp;
% %                     
% %                 end 
% %              end 
%              
% %              for unitNbr = 1:size(Table.PMV_Go{1,1},2)               
% %                 NeuralAction =  Table.PMV_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
% %                 
% %                 if delayDataPoints(classNbr,nbrRep) > 0
% %                     PMV_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
% %                 else
% %                     %we pad the data with the last x values of the delay
% %                     %phase
% %                     NeuralDelay2 =  Table.PMV_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
% %                     Padding = NeuralDelay2(end-abs(delayDataPoints):end);
% %                     PMVTmp =utile.delayseq(NeuralAction, DelayPoint);
% %                     PMVTmp(1:length(Padding)) = Padding;
% %                     PMV_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =PMVTmp;
% %                     
% %                 end 
% %                 
% %              end 
%              
% %              for unitNbr = 1:size(Table.S1X_Go{1,1},2)               
% %                 NeuralAction =  Table.S1X_Go{IdxAll(nbrRep), 1}(ActionPhaseTimes,unitNbr);
% %                 if delayDataPoints(classNbr,nbrRep) > 0
% %                     S1X_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =utile.delayseq(NeuralAction, DelayPoint);
% %                 else
% %                     %we pad the data with the last x values of the delay
% %                     %phase
% %                     NeuralDelay2 =  Table.S1X_Go{IdxAll(nbrRep), 1}(find(TableAll.time_phase_labels{1,1} == 5),unitNbr);
% %                     Padding = NeuralDelay2(end-abs(delayDataPoints):end);
% %                     S1XTmp =utile.delayseq(NeuralAction, DelayPoint);
% %                     S1XTmp(1:length(Padding)) = Padding;
% %                     S1X_GoAligned{IdxAll(nbrRep),1}(:,unitNbr) =S1XTmp;
% %                     
% %                 end
% %              end
%             
%             Audio_GoAlined{IdxAll(nbrRep),1} =  AudioActionAlignedTmp;
%             AudioActionAlignedEnvDS{classNbr,nbrRep} = AudioActionAlignedTmp;
%             
%             %figure(); plot(AudioActionAlignedTmp)
% 
%         end 
% 
%     end 
%     
%     SMG_GoAlignedAll = [SMG_GoAlignedAll; SMG_GoAligned];
%     PMV_GoAlignedAll = [PMV_GoAlignedAll; PMV_GoAligned];
%     S1X_GoAlignedAll = [S1X_GoAlignedAll; S1X_GoAligned];
%     Audio_GoAlignedAll = [Audio_GoAlignedAll; Audio_GoAlined];
% 
%     Colors = utile.get_color_rgb_codes(ClassNames);
%  
%     fig = figure('units','normalized','outerposition',[0.35 0.05 0.4 0.8]);
%     
%     
% %     %size test - turns out some are 30 and some 31 - shorten all
%      SizeAudio = cell2mat(cellfun(@(x) size(x,1), AudioActionEnv, 'UniformOutput', false));
% %     %if length(unique(SizeAudio)) ~= 1
% %         %%% WORK FROM HERE%%%
%     minSize = min(min(SizeAudio));
%     warning(['shorten Audio data to have the same size :' num2str(minSize)])
%     try
% 
%         AudioActionEnv = cellfun(@(x) x(1:minSize),AudioActionEnv, 'UniformOutput', false);
%     catch
%         disp('blu')
%     end 
%      AudioActionAlignedEnvDS = cellfun(@(x) x(1:minSize), AudioActionAlignedEnvDS, 'UniformOutput', false);
%         %% MAKE THEM ALL SHORTER 
%         %e.g here issue is some are 30 and some are 31 
%     %end 
% 
%     for i = 1:length(unique(Labels))
% 
%         subplot(8,2,2*i-1); plot(cell2mat(AudioActionEnv(i,:))); title('Unaligned Audio data')
%         subplot(8,2,2*i); plot(cell2mat(AudioActionAlignedEnvDS(i,:))); 
%         title([ClassNames{i} ' aligned'])
%         %subplot(1,2,1); plot(cell2mat(AudioAction(i,:))); title('Audio data')
%         %subplot(1,2,2); plot(cell2mat(AudioActionAligned(i,:))); title('Audio data aligned')        
%     end
%     sgtitle(unique_sessions(nbr_sessions))
%     disp('Acoustic alignment finished')
% 
%     
%     
%     
% end
% 
% SMG_GoAligned = SMG_GoAlignedAll;
% PMV_GoAligned = PMV_GoAlignedAll;
% S1X_GoAligned = S1X_GoAlignedAll;
% Audio_SubAlinged = Audio_GoAlignedAll;
% 
% TableAll = [TableAll, cell2table(SMG_GoAligned),cell2table(PMV_GoAligned), cell2table(S1X_GoAligned), cell2table(Audio_SubAlinged)];
%  
