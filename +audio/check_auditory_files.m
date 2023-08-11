%Checking and saving auditory files

clc
clear all 
close all 

%save data 


GraspCue = 'Speech';

save_data_pathway = 'D:\Users\Sarah\Documents\Saved_Data\s2_2020_aligned_correctly';
save_data = true; 
%auditory files are saved in PMV array

subject_id = 's2';  % s2 or p3
subject = hst.Subject(subject_id);
sessions_date = '20210712';
%sessions_date = '20210722';


session = hst.Session(sessions_date, subject);
taskfiles = session.getTaskFiles('Speech');

nsp_audio = 'PMV';
%nsp_audio = 'S1X_S1'; %should not need that -> won't work because means
%the wrong config file was loaded. 
%nsp_audio = 'SMG_AIP';

data_files = env.get('data');

file_dir = fullfile(data_files{1},subject_id,sessions_date, nsp_audio );
filename = [subject_id '_good_trials_' GraspCue '.xlsx'];

[numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_grasps\ExcelFiles\' filename]);     
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


 for task_n = 1:length(good_datasets) 
  
        task = hst.Task(taskfiles{good_datasets(task_n)});
        %Find the correct audio file: it's in PMV and it's a ns5 file
        audio_file = dir(fullfile(file_dir, [task.srcFile '*aligned.ns5']));
        %Load the ns file
        ns = Blackrock.NSx(fullfile(file_dir , audio_file.name));
        %Find the channel number which contains the recording of microphone
        ChannelLabels = {ns.ChannelInfo.Label};
        %find nonempty channels to be able to use ismember
        non_empty_channels = ~cellfun(@isempty,ChannelLabels);
        %find microphone idx
        microphone_idx = find(ismember(ChannelLabels(non_empty_channels),'microphone'));
        %find associated channel
        microphone_channel = find(non_empty_channels,microphone_idx);
        %load sound file
        
        %I can use ns read to extract the time that I want (provided by the
        %task file!) e.g. ns.read('Channel', microphone_channel, 'time',
        %[start_time, length])
        
        start_times = task.trialTimes(:,1);
        average_trial_time = mean(task.trialTimes(:,2));
        end_times = task.trialTimes(:,1) + average_trial_time;
        
        
        
        Fs = ns.Fs; 
        num_trials = length(start_times); 
        audio_per_trial = cell(1,num_trials); 
        envelope_per_trial = cell(1,num_trials);
        %envelope_per_trial_cleaner
        time_per_trial = cell(1,num_trials); 
        %extract audio files per trial 
       
        Recorded_sound = ns.read('Channel',microphone_channel);
        %normalize data
        Recorded_sound = Recorded_sound/ max(abs(Recorded_sound)); 
        time = (1:length(Recorded_sound))/Fs; 
         %Plot audio file 
        figure(); plot(time,Recorded_sound); 
        xlabel('Time [s]')
        ylabel('Audio signal'); 
        title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])
        %sound(Recorded_sound, Fs) 
        
         filename = 'C:\Users\Sarah\Documents\Audio\ex2.wav';
        audiowrite(filename,Recorded_sound,Fs);

            
        
        for trial = 1:num_trials
            
            %Error - not working
            %audio_tmp = ns.read('Channel',microphone_channel, 'time', [start_times(trial) end_times(trial)]);
            
            %Extracting correct time points by finding start and end idx
            TimeAudio = (1:length(Recorded_sound))/Fs;
            StartIdx = find(TimeAudio > start_times(trial),1);
            EndIdx = find(TimeAudio > end_times(trial),1);
            audio_tmp = Recorded_sound(StartIdx:EndIdx);
            
            %AudioBeginning = Recorded_sound(1:EndIdx);
            %filename = 'C:\Users\Sarah\Documents\Audio\AudioBeginning.wav';
            %audiowrite(filename,AudioBeginning,Fs);

            %audio_tmp = ns.read('Channel',microphone_channel);

            %normalize the data to make it sound better!
            audio_tmp  = audio_tmp / max(abs(audio_tmp));
            
         
            
            %Envelope of signal: 
            env = envelope(audio_tmp, Fs/6, 'peak'); 
            time_axis = (1:size(audio_tmp,2)) / Fs; 
           
            
            %Plot audio file 
            
%             LMax = logical(zeros(size(env))); 
%             LMin = logical(zeros(size(env)));
%             interest_time = 8.9*Fs:10.8*Fs; 
%             LMax(interest_time) = islocalmax(env(interest_time));
%             LMin(interest_time) = islocalmin(env(interest_time));
            %audio_final = zeros(
            %put local min and max to 0 if not between 8 and 10.5 seconds            
            figure();
            plot(time_axis, audio_tmp); 
            hold on; plot(time_axis, env, 'r');
%             hold on; plot(time_axis(LMax), env(LMax), 'k*');
%             hold on;  plot(time_axis(LMin), env(LMin), 'k*');
            xlabel('Time [s]')
            ylabel('Audio signal'); 
            title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n) ' Trial ' num2str(trial)])
            %hold on 
            %[max_val, max_idx] = max(env);
            %Play audio file 
            
%             figure();
%             subplot(1,2,1)
%             plot(time_axis, audio_tmp)
%             subplot(1,2,2)
%             %plot(time_axis, wdenoise(audio_tmp,4))
%             plot(time_axis, wdenoise(audio_tmp,10))



            
            sound(audio_tmp, Fs); 
            audio_per_trial{1,trial} = audio_tmp; 
            envelope_per_trial{1,trial} = env; 
            time_per_trial{1,trial} = time_axis; 

        end 
        
        filename_save = [subject_id '_' sessions_date '_' GraspCue '_audio_file.mat'];
      

    if save_data
        filename = [save_data_pathway '\AudioFiles\' filename_save ];
        disp(['save data under dataname: ' filename_save ' and subject: ' subject_id]); 
       
        save(filename, 'audio_per_trial', 'envelope_per_trial', 'time_per_trial'); 
        disp(['The data ' filename_save '  has been saved'])
    end 
    
 end  
 
 
 
 
Cue_labels = [task.trialparams.Image_code];
unique_cue_types = unique(Cue_labels); 
Label_names = arrayfun(@(x) preproc.image2class_simple(x), unique_cue_types,'UniformOutput', false);
time =(1:length(envelope_per_trial{1,1}))/Fs;
for blub = 1:length(unique_cue_types )
    
    Audio = cell2mat(envelope_per_trial((Cue_labels == unique_cue_types(blub)))');
    figure(); 
    for i = 1:8
        plot(time,Audio(i,:))
        hold on
    end 
    plot(time,mean(Audio), 'Linewidth',2)
    title(Label_names{blub})
    Audio_mean_per_cue{blub} = mean(Audio); 
end 


% timeDiff =  -23.5690; 
%  disp('----------------------------------------------------')
%         
%         disp(['Time difference: ' num2str(timeDiff) ' s ']);
%         %disp(['Length of Recording: ' num2str(length_recorded_sound) ' s ']);
%         
%         
%         disp('---------Trial start times in video:   -------------------------------------------')
%         
%         video_trial_start_time = task.trialTimes(:,1) + timeDiff;
%         trial_start_time = arrayfun(@(x) datestr(seconds(x),'MM:SS'), video_trial_start_time,'UniformOutput', false);
% 
        



%         %Plot audio file 
%         figure(); plot(Recorded_sound); 
%         xlabel('Time [s]')
%         ylabel('Audio signal'); 
%         title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])


%  for i = unique(labels)
%      
%      blub2 = envelope_per_trial(find(labels == i));
%      blub = cell2mat(envelope_per_trial(find(labels == i))');
%      blub_normalize = normr(blub); 
%      %figure();
%      %for j = 1:size(blub,1)
%      %    plot(time_per_trial{1}, blub(j,:));
%      %    hold on
%      %end 
%      
%      figure();
%      for j = 1:size(blub,1)
%          plot(time_per_trial{1}, blub_normalize(j,:));
%          hold on
%      end 
%      plot(time_per_trial{1}, mean(blub_normalize), 'k', 'LineWidth', 2);
%      ylabel('Normalized Audio signal'); 
%      xlabel('Time [s]'); 
%      title([preproc.image2class_simple(i)]);
%  end 
%  
%  
%  
% Recorded_sound = ns.read('Channel',microphone_channel);
% Recorded_sound = Recorded_sound/ max(abs(Recorded_sound)); 
% time_axis = (1:size(Recorded_sound,2)) / ns.Fs; 
% 
% %Plot audio file 
% figure(); plot(time_axis(:,250*Fs:350*Fs), Recorded_sound(:,250*Fs:350*Fs)); 
% xlabel('Time [s]')
% ylabel('Audio signal'); 
% title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])
% sound(Recorded_sound(:,270*Fs:450*Fs), Fs)


%         Recorded_sound = ns.read('Channel',microphone_channel);
%         time_axis = (1:size(Recorded_sound,2)) / ns.Fs; 
% 
%         

%   
%         
%         %Plot audio file 
%         start_time = task.trialTimes(1,1); 
%         end_time = task.trialTimes(end,1) + task.trialTimes(end,2);
%         Recorded_sound_subset = Recorded_sound(1,start_time*Fs:end_time*Fs);
%         time_subset = time_axis(1,start_time*Fs:end_time*Fs);
%         figure(); plot(time_subset, Recorded_sound_subset); 
%         xlabel('Time [s]')
%         ylabel('Audio signal'); 
%         title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])
%         
%        
%         figure(); plot(time_axis, env, 'r'); hold on; plot(time_axis, Recorded_sound);
%         
%         %Play audio file 
%         sound(Recorded_sound_subset, Fs); 
%         %spectrogram(Recorded_sound)    
% 
% 


 
%         %Play audio file 
%         %sound(Recorded_sound_subset, ns.Fs); 
%         spectrogram(Recorded_sound)    
%         
%         T =1/Fs;
%         L = length(Recorded_sound_subset);  
%         Y = fft(Recorded_sound_subset);
%         P2 = abs(Y/L); 
%         P1 = P2(1:L/2+1); 
%         P1(2:end-1) = 2*P1(2:end -1);
%         f = Fs*(0:(L/2))/L;
%         figure();
%         plot(f,P1) 
%         title('Single-Sided Amplitude Spectrum of X(t)')
%         xlabel('f (Hz)')
%         ylabel('|P1(f)|')
% 
%         %bandpass filtering sound 
%         Fn = Fs/2;                                              % Nyquist Frequency (Hz)
%         Wp = 1000/Fn;                                           % Passband Frequency (Normalised)
%         Ws = 1010/Fn;                                           % Stopband Frequency (Normalised)
%         Rp =   1;                                               % Passband Ripple (dB)
%         Rs = 150;                                               % Stopband Ripple (dB)
%         [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
%         [z,p,k] = cheby2(n,Rs,Ws,'low');                        % Filter Design
%         [soslp,glp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability
%         figure(3)
%         freqz(soslp, 2^16, Fs)                                  % Filter Bode Plot
%         filtered_sound = filtfilt(soslp, glp, Recorded_sound_subset);
%         sound(filtered_sound, Fs)
%         
%         
%         bfil=fft(Recorded_sound_subset); %fft of input signal
%         fs = Fs;
%         wn=[4000 8000]/(fs/2);   %bandpass
%         [b,a]=butter(6,wn);
%         fvtool(b,a);
%         f=filter(b,a,z);
%         afil=fft(f);
%         subplot(2,1,1);plot(real(bfil));
%         title('frequency respones of input signal');
%         xlabel('frequency');ylabel('magnitude');
%         subplot(2,1,2);plot(real(afil));
%         title('frequency respones of filtered signal');
%         xlabel('frequency');ylabel('magnitude');
%         
%         
%         Fn  = Fs/2;                                 % Nyquist Frequency
%         Fco =   2000;                                 % Passband (Cutoff) Frequency
%         Fsb =   3000;                                 % Stopband Frequency
%         Rp  =    100;                                 % Passband Ripple (dB)
%         Rs  =   1000;                                 % Stopband Ripple (dB)
%         [n,Wn]  = buttord(Fco/Fn, Fsb/Fn, Rp, Rs);  % Filter Order & Wco
%         [b,a]   = butter(n,Wn, 'high');             % Lowpass Is Default Design
%         [sos,g] = tf2sos(b,a);                      % Second-Order-Section For STability
%         figure(1)
%         freqz(sos, 2048, Fs)
%         filtered_sound = filtfilt(sos, g, Recorded_sound_subset);
%         
%         [A,B,C,D] = butter(10,[3000 Fn-1]/Fn);
%         %d = designfilt('bandpassiir','FilterOrder',20, ...
%         %  'HalfPowerFrequency1',2000,'HalfPowerFrequency2',3000, ...
%         %'SampleRate',Fs);
%         [sos,g] = ss2sos(A,B,C,D);
%        filtered_sound2 = filtfilt(sos, g, Recorded_sound_subset);
% 
%         
%         %fvt = fvtool(sos,d,'Fs',1500);
%         %legend(fvt,'butter','designfilt')
%         
%         
%         time_subset = time_ax  is(1,70*Fs:106*Fs);
%         figure(); plot(time_subset, filtered_sound2); 
%         xlabel('Time [s]')
%         ylabel('Audio signal'); 
%         title(['Audio signal - ' session_date ' - Block ' num2str(task_n)])
%         sound(filtered_sound2)
%         %sound(filtered_sound)
%         %Plot audio file 
%         time_subset = time_axis(1,70*Fs:106*Fs);
%         figure(); plot(time_subset, Recorded_sound_subset); 
%         xlabel('Time [s]')
%         ylabel('Audio signal'); 
%         title(['Audio signal - ' session_date ' - Block ' num2str(task_n)])
        
 

