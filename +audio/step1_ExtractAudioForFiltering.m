%Checking and saving auditory files

clc
clear all 
close all 

%save data 



save_data_pathway = 'D:\Users\Sarah\Documents\Saved_Data\s2_2021_speech_aligned_correctly\AudioFiles\Tmp';

save_data = true; 
%auditory files are saved in PMV array

subject_id = 's2';  % s2 or p3
subject = hst.Subject(subject_id);
sessions_date = '20210712';
%sessions_date = '20210729';
%sessions_date = '20210923';
%sessions_date = '20210930';
%sessions_date = '20211018';
%sessions_date = '20211027';
%sessions_date = '20211103';
%sessions_date = '20220323'; % missing microphone audio :( 
%sessions_date = '20211230';
%sessions_date = '20230113';
%sessions_date = '20230106';
%sessions_date = '20230117';
%sessions_date ='20211011';
sessions_date ='20210722';

GraspCue = 'Speech';
%GraspCue = 'SpeechTrainingInternal';

%GraspCue = 'Homonym_task';

%sessions_date = '20220523';
session = hst.Session(sessions_date, subject);
taskfiles = session.getTaskFiles('Speech');

nsp_audio = 'PMV';

if strcmp(subject_id, 's2')

    if str2double(sessions_date) > 20220300
        flagNewRig = true; 
        nsp_audio = 'APX';
        env.set('arrays', {'APX', 'S1X_S1'});
    else
        flagNewRig = false; 
    end 
    
elseif  strcmp(subject_id, 's3')
    flagNewRig = true; 
    nsp_audio = 'MPX';
    env.set('arrays', {'MPX', 'S1X_S1'});
end 
    
%nsp_audio = 'S1X_S1'; %should not need that -> won't work because means
%the wrong config file was loaded. 
%nsp_audio = 'SMG_AIP';

data_files = env.get('data');

file_dir = fullfile(data_files{1},subject_id,sessions_date, nsp_audio );
filename = [subject_id '_good_trials_' GraspCue '.xlsx'];

if strcmp(subject_id, 's2')
    [numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech\ExcelFiles\' filename]);     

else
    [numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech - s3\ExcelFiles\' filename]);     

end 
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


%% STEP  1 - Save Files needed for extraction 
 for task_n = 1:length(good_datasets) 
  
        task = hst.Task(taskfiles{good_datasets(task_n)});
        %Find the correct audio file: it's in PMV and it's a ns5 file
        if flagNewRig
            audio_file = dir(fullfile(file_dir, [task.srcFile '*.ns5']));
        else
            audio_file = dir(fullfile(file_dir, [task.srcFile '*aligned.ns5']));
        end 
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
        
        %Recorded_sound = ns.read('Channel',10);

        %Recorded_sound2 = ns.read(
        %normalize data
        Recorded_sound = Recorded_sound/ max(abs(Recorded_sound)); 
        time = (1:length(Recorded_sound))/Fs; 
         %Plot audio file 
        figure(); plot(time,Recorded_sound); 
        xlabel('Time [s]')
        ylabel('Audio signal'); 
        title(['Audio signal - ' sessions_date ' - Block ' num2str(task_n)])
        %sound(Recorded_sound, Fs) 
        
        filename_save = [subject_id '_' sessions_date '_' GraspCue '_audio_file' num2str(task_n)];       
        filename = fullfile(save_data_pathway,[filename_save '.wav']); %'C:\Users\Sarah\Documents\Audio\ex2.wav';
      
        %change file names - if ever needed again
        %filename_save = [subject_id '_' sessions_date '_errorTrials_' GraspCue '_audio_file' num2str(task_n)];
        %filename = fullfile(save_data_pathway,[filename_save '.wav']); %'C:\Users\Sarah\Documents\Audio\ex2.wav';
        %filenameOld = fullfile('D:\Users\Sarah\Documents\Saved_Data\InternalSpeechPaper\s2_InternalSpeech_Paper\AudioFiles',[filename_save '.mat']);
        %filename_save_new = [task.taskString '_errorTrials'];
        %filenameNew= fullfile('D:\Users\Sarah\Documents\Saved_Data\InternalSpeechPaper\s2_InternalSpeech_Paper\AudioFiles',[filename_save_new '.mat']);       
        %if ~exist(filenameNew)
        %    movefile(filenameOld, filenameNew,'f');
        %end 
        
         audiowrite(filename,Recorded_sound,Fs);
         disp(['Audio data saved in ' save_data_pathway])
%   
        for trial = 1
            
            TimeAudio = (1:length(Recorded_sound))/Fs;
            StartIdx = find(TimeAudio > start_times(trial),1);
            EndIdx = find(TimeAudio > end_times(trial),1);           

            AudioBeginning = Recorded_sound(1:EndIdx);
            filename = fullfile(save_data_pathway, [filename_save '_beginning.wav']);
            audiowrite(filename,AudioBeginning,Fs);


        end 
        
 end  
 
 
 
