%Checking and saving auditory files

clc
clear all 

%%
%%%%%To do for new session: %%%%
% - change session date
% - add session date to variable "all_sessions"
% - set flagControlDelayTime to true
% - run code
% - open variable DelayTable and verify that the delay is always the same
% for different trials
% - add the delay to the variable offset_timestamps 
% - set flagControlDelayTime to false -> run again

%COONTINUE HERE FOR 20211011!


%%
GraspCue = 'Speech';

save_data_pathway = 'D:\Users\Sarah\Documents\Saved_Data\s2_2021_speech_aligned_correctly';
save_data = true; 
%auditory files are saved in PMV array

subject_id = 's2';  % s2 or p3
subject = hst.Subject(subject_id);
%sessions_date = '20210712';
%sessions_date = '20210722';
%sessions_date = '20210729';
%sessions_date = '20210923';
%sessions_date = '20210930';
%sessions_date = '20211011';
%sessions_date = '20211018';
%sessions_date = '20211027';
sessions_date = '20211103';
%sessions_date = '20220323';

TaskName = 'Speech';
flagControlDelayTime = false; %PUT TO TRUE WHEN NEW SESSION AND ADD TO OFFSET_TIMESTAMPS
flagPlotFigures = false; 
filenames = dir(['D:\Users\Sarah\Documents\Saved_Data\s2_2021_speech_aligned_correctly\AudioFiles\FilteredAudio\' subject_id  '_' sessions_date '_' TaskName '*'] );
%DelayTable
offset_timestamps = [60527,59901;60061,59782;59869,59868; 59660,59617;59961,60258;60050,60019;60165,60197;59661,59658;59958,60232]; %given in timepoints; order:  20210712,20210722,20210729


if flagControlDelayTime
    TrialsToTest = 10;
else
    TrialsToTest = 1;
end 

maxTime =237781; %238610 %238611
% for 20210923 [59660;59617]  59660 59617
%[60527;59901]
%[60061;59782]
%[59869;59868]
%[59660;59617]

all_sessions = {'20210712','20210722','20210729','20210923','20210930','20211011','20211018','20211027','20211103'};
Offset_table = table( all_sessions', offset_timestamps);
session = hst.Session(sessions_date, subject);
taskfiles = session.getTaskFiles('Speech');

[numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech\ExcelFiles\' [subject_id '_good_trials_' GraspCue '.xlsx']]);     
session_date = str2double(sessions_date);
good_datasets = numb(session_date==numb(:,1), 2:end);
good_datasets(isnan(good_datasets)) = [];

if (length(good_datasets) <1)
    error ('Add to excel sheet, No good dataset present, skip or check for problem')
elseif good_datasets == 0
    error('No good dataset present for this day, has been checked, skip')  
elseif any(ismember(good_datasets, 1000))
    warning('There is a problem with one of the datasets, ask Spencer? Removing it for now');
    good_datasets(ismember(good_datasets,1000)) = []; 

elseif good_datasets == 100
    disp(['No auditory dataset for session ' num2str(sessions_date)])       
else
    disp('Everything seems ok');
end

removeErrorTrials = true; 

%% STEP  2 - Separate files per trial and save
 for task_n =1:length(filenames) 
  
        task = hst.Task(taskfiles{good_datasets(task_n)});
               
        if flagControlDelayTime
            offset_time = 0; %
        else
            offset_time = Offset_table{ismember(Offset_table{:,1}, sessions_date),2}(task_n);
        end

        if isempty(offset_time)
            error(['add offset time to session ' sessions_date ])
        end 

        start_times = task.trialTimes(:,1);
        average_trial_time = mean(task.trialTimes(:,2));
        end_times = task.trialTimes(:,1) + average_trial_time;
             
        if removeErrorTrials
                if task_n == 1 %audio datasets
                    data_name = 'Speech';
                elseif task_n == 2 %written datasets
                    data_name = 'Written';
                end 
                data_subset = setdiff(1:task.numTrials,  preproc.errorTrials(sessions_date,data_name));

            else 
                data_subset = 1:task.numTrials; 
        end 
       
        
        
        num_trials = length(data_subset); 
        audio_per_trial = cell(1,num_trials); 
        envelope_per_trial = cell(1,num_trials);
        envelope_per_trial_downsampled = cell(1,num_trials);
        %envelope_per_trial_cleaner
        time_per_trial = cell(1,num_trials); 
        
        %load the correct filtered audio file 
        [Recorded_sound,Fs] = audioread(fullfile(filenames(task_n).folder, filenames(task_n).name));
        Recorded_sound = Recorded_sound';
       
        
        time = (1:length(Recorded_sound))/Fs; 
         %Plot audio file 
        figure(); plot(time,Recorded_sound); 
        xlabel('Time [s]')
        ylabel('Audio signal'); 
        title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])
        %sound(Recorded_sound, Fs) 
        
        if flagControlDelayTime
            %load ns file to control if the delay time is the same
            file_dir = ['Z:\Data\s2\' sessions_date '\PMV'];
            audio_file = dir(fullfile(file_dir, [task.srcFile '*aligned.ns5']));
            ns = Blackrock.NSx(fullfile(audio_file.folder, audio_file.name));
            Recorded_sound = ns.read('Channel', 129); %microphone channel is 129
            Recorded_sound = Recorded_sound;%/max(abs(Recorded_sound));%normalize
            
        end 
        
        
        
         for trial_n = 1:TrialsToTest:num_trials
            
            trial = data_subset(trial_n);
            %Extracting correct time points by finding start and end idx
            TimeAudio = (1:length(Recorded_sound))/Fs;
            StartIdx = find(TimeAudio > start_times(trial) - offset_time/Fs,1);
            EndIdx = find(TimeAudio > end_times(trial) - offset_time/Fs,1);
            audio_tmp = Recorded_sound(StartIdx:EndIdx);
            %sound(audio_tmp, Fs);

            %For control purposes: read in raw audio data to calculate
            %delay time
            
            if flagControlDelayTime
                audio_tmp_raw = ns.read('time', [start_times(trial) end_times(trial)]);
                audio_tmp1 = audio_tmp_raw(1,:); %/ max(abs(audio_tmp_raw(1,:)));  % 1 corresponds to microphone array
                %verify that it is equal to the delay in the table 
                DelayTable(task_n,trial_n) = finddelay(audio_tmp,audio_tmp1);
                
                % [c,lags] = xcorr(audio_tmp,audio_tmp1);
                % stem(lags,c)
                
                if flagPlotFigures
                    figure(); 
                    plot(audio_tmp);
                    hold on;
                    plot(audio_tmp1);
                    legend()
                end 
            end 
            %Envelope of signal: 
            env = envelope(audio_tmp, Fs/60, 'peak'); 
            time_axis = (1:size(audio_tmp,2)) / Fs; 
     
            if flagPlotFigures
                figure();
                plot(time_axis, audio_tmp); 
                hold on; plot(time_axis, env, 'r');
                xlabel('Time [s]')
                ylabel('Audio signal'); 
                title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n) ' Trial ' num2str(trial)])
            end 
         
            DownsampleAudio = downsample(audio_tmp,Fs/20);
            DownsampleEnvelope = downsample(env,Fs/20);
           
            audio_per_trial{1,trial_n} = audio_tmp; 
            envelope_per_trial{1,trial_n} = env; 
            envelope_per_trial_downsampled{1,trial_n} = DownsampleEnvelope;
            time_per_trial{1,trial_n} = time_axis; 

            if flagPlotFigures
                figure();
                subplot(2,1,1)
                plot(time_axis,audio_tmp);
                hold on
                plot((0:length(DownsampleAudio)-1)*0.05,DownsampleAudio, 'LineWidth', 2);
                subplot(2,1,2)
                plot(time_axis,env);
                hold on
                plot((0:length(DownsampleEnvelope)-1)*0.05,DownsampleEnvelope, 'LineWidth', 2);
            end 
         end 
        
        if removeErrorTrials
            filename_save = [subject_id '_' sessions_date '_errorTrials_' GraspCue '_audio_file' num2str(task_n) '.mat'];
        else
            filename_save = [subject_id '_' sessions_date '_' GraspCue '_audio_file' num2str(task_n) '.mat'];
        end 
        
        MinLength = min(cell2mat(cellfun(@(x) size(x,2),envelope_per_trial,'UniformOutput', false)));
        envelope_per_trial = cellfun(@(x) x(:,1:MinLength), envelope_per_trial,'UniformOutput', false);
        audio_per_trial = cellfun(@(x) x(:,1:MinLength), audio_per_trial,'UniformOutput', false);
        LengthDS = unique(cell2mat(cellfun(@(x) size(x,2),envelope_per_trial_downsampled,'UniformOutput', false)));
        if ~flagControlDelayTime
            if length(LengthDS) ~= 1; error('Downsampled data does not have equal length');  end
        end

        if ~flagControlDelayTime %save data if we are not conrolling delay
            if save_data
                filename = [save_data_pathway '\AudioFiles\' filename_save ];
                
                disp(['save data under dataname: ' filename_save ' and subject: ' subject_id]); 

                save(filename, 'audio_per_trial', 'envelope_per_trial','envelope_per_trial_downsampled', 'time_per_trial'); 
                disp(['The data ' filename_save '  has been saved'])
            end 
        end 
    
        TaskName = {'AuditoryCue', 'WrittenCue'};
        Cue_labels = [task.trialparams.Image_code];
        Cue_labels = Cue_labels(data_subset); 
        [unique_cue_types] = unique(Cue_labels); 
        num_t = histc(Cue_labels, unique_cue_types);
        Label_names = arrayfun(@(x) preproc.image2class_simple(x), unique_cue_types,'UniformOutput', false);
        time =(1:length(envelope_per_trial{1,1}))/Fs;
        
        if ~flagControlDelayTime
            for blub = 1:length(unique_cue_types )

                %MinLength = min(cell2mat(cellfun(@(x) size(x,2),envelope_per_trial,'UniformOutput', false)));
                %envelope_per_trial = cellfun(@(x) x(:,1:MinLength), envelope_per_trial,'UniformOutput', false);

                Audio = cell2mat(audio_per_trial((Cue_labels == unique_cue_types(blub)))');
                figure(); 
                for i = 1:num_t(blub)
                plot(time,Audio(i,:) + 0.5*i)
                hold on
                end 
                plot(time,mean(Audio),'k', 'Linewidth',2)
                title([TaskName{task_n} ' ' Label_names{blub}])
                Audio_mean_per_cue{blub} = mean(Audio); 

            end 

%             time = (1:length(envelope_per_trial_downsampled{1,1}))/20;
%              for blub = 1:length(unique_cue_types )
% 
%                 %MinLength = min(cell2mat(cellfun(@(x) size(x,2),envelope_per_trial,'UniformOutput', false)));
%                 %envelope_per_trial = cellfun(@(x) x(:,1:MinLength), envelope_per_trial,'UniformOutput', false);
% 
%                 Audio = cell2mat(envelope_per_trial_downsampled((Cue_labels == unique_cue_types(blub)))');
%                 figure(); 
%                 for i = 1:num_t(blub)
%                 plot(time,Audio(i,:) + 0.5*i)
%                 hold on
%                 end 
%                 plot(time,mean(Audio), 'Linewidth',2)
%                 title([TaskName{task_n} ' ' Label_names{blub}])
%                 Audio_mean_per_cue{blub} = mean(Audio); 
% 
%              end 
        end 
 end  
 
 
%sound(audio_per_trial{15},Fs)



 
