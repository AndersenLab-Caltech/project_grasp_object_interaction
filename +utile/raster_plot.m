%Raster plot for GraspShape experiment - adapted from Spencer

clc
clear all 
close all 

%save data 

GraspCue = 'Speech'

sessions_date = '20210923';

save_data_pathway = 'D:\Users\Sarah\Documents\Saved_Data\s2_2021_speech_aligned_correctly';
save_data = true; 
flagSortedData = true; 
flagPlotFigures = true; 
flagSpikeAnalysis = false; 
flagThresholdData = false; %takes sorted data (right now all ,should maybe change to a specific value)
%auditory files are saved in PMV array

subject_id = 's2';  % s2 or p3
subject = hst.Subject(subject_id);
%sessions_date = '20190911';
%sessions_date = '20190610';
%sessions_date = '20190417';
%sessions_date = '20191028';
%sessions_date = '20191016';

%brain_region = 'S1X_S1'; %SMG_AIP %S1X_S1 %PMV
%brain_region = 'PMV'; %SMG_AIP %S1X_S1 %PMV
%brain_region = 'SMG_AIP'; %SMG_AIP %S1X_S1 %PMV
brain_regions = {'SMG_AIP', 'PMV', 'S1X_S1'};
%brain_regions = {'SMG_AIP'};

save_fig_folder = ['D:\Users\Sarah\Pictures\InternalSpeech\Raster_Plot\' sessions_date]; 

 
session = hst.Session(sessions_date, subject);

taskfiles = session.getTaskFiles('Speech');

dataZ = env.get('data');

filename = [subject_id '_good_trials_' GraspCue '.xlsx'];

[numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech\ExcelFiles\' filename]);     
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

flagPlotAuditory = true; %adds the auditory signal on top (only if present)
DatesRemoveTrials = []; %add if there are dates to remove
 for task_n = 1:length(good_datasets) 
     %Channels_per_trial_timestemps = cell([num_electrodes, num_trials, length(brain_regions)]);
     Channels_per_trial_timestemps = {}; 
     electrode_idx = [];
     Waveforms = {};
     
     if flagPlotAuditory
           % AuditoryFile = dir(fullfile(save_data_pathway,'AudioFiles', [subject_id '_' num2str(session_date) '_' GraspCue '_audio_file.mat']));
            AuditoryFile = dir(fullfile(save_data_pathway,'AudioFiles', [subject_id '_' num2str(session_date) '_' GraspCue '_audio_file' num2str(task_n) '.mat']));

            
            if ~isempty(AuditoryFile)
                AudioData_tmp = load(fullfile(AuditoryFile.folder, AuditoryFile.name));
                AudioData = AudioData_tmp.envelope_per_trial; %could change to other version of Audio file

            else
                warning('No Auditory File present');  
                AudioData = [];
            end 
      end 

     for brain_r = 1:length(brain_regions)
         

        brain_region = brain_regions{brain_r};
        task = hst.Task(taskfiles{good_datasets(task_n)});

        %Compute the approximate start and end time of a trial (based on
        %the average trial time)
        start_times = task.trialTimes(:,1);
        average_trial_time = mean(task.trialTimes(:,2));
        end_times = task.trialTimes(:,1) + average_trial_time;
        
        %find the correct files 
        file_dir = fullfile(dataZ{1},subject_id,sessions_date, brain_region );
        if flagSortedData
            
            dataFiles = dir(fullfile(file_dir, 'SortedObjects',[task.srcFile '*aligned*.nev']));
            save_fig_folder_new = fullfile(save_fig_folder,'SortedSpikes' );
            file_dir = fullfile(file_dir, 'SortedObjects');
            disp('Using Sorted Data');
        elseif flagThresholdData
            dataFiles = dir(fullfile(file_dir, [task.srcFile '*aligned*3.5.nev']));
            save_fig_folder_new = fullfile(save_fig_folder,'Threshold_-3.5' );
            disp('Using thresholded Data at -3.5'); 
        else 
            dataFiles = dir(fullfile(file_dir, [task.srcFile '*aligned.nev']));
            save_fig_folder_new = save_fig_folder; 
            disp('Using unsorted unthresholded data'); 
        end 
        
       
          %non aligned data 
         dataFiles = dir(fullfile(file_dir, [task.srcFile '*.nev']));
        
        %Load the nev file
        for nbr_f = 1%:length(dataFiles)
            
            file_name = dataFiles(nbr_f).name; 
            save_name = fullfile(save_fig_folder_new, extractBefore(file_name,'.'));
            nev = Blackrock.NEV(fullfile(file_dir , file_name));
            Fs = 30000;
            data = nev.read('WAVEFORMS', 'time', [start_times(1), end_times(1)]);
            
            %ns = task.getNeuralDataObject('SMG', 'ns6');
            %data2 = ns.read
            disp(['Number of spikes that have been recorded during the trial: ' num2str(size(data.Timestamps,1))])

            unique_channels = unique(data.Channels);         
            electrode_idx = [electrode_idx; brain_r*(ones(length(unique_channels),1))]; 
            
            if ismember(sessions_date, DatesRemoveTrials)
               trialsToKeep = utile.getTrialsToKeep(sessions_date, GraspCue); 
            else
                trialsToKeep = 1:length(start_times);
            end
            
            group_data = [task.trialparams([trialsToKeep]).Image_code];
            %labels_name  = arrayfun(@(x) preproc.image2class_simple(x), unique(group_data), 'UniformOutput', false); 
            labels_name  = preproc.image2class_simple(unique(group_data));

            labels_per_trial = arrayfun(@(x) preproc.image2class_simple(x), group_data, 'UniformOutput', false); 
            data_per_group = arrayfun(@(x) find(group_data == x), unique(group_data), 'UniformOutput', false)';

            average_phase_start_time = mean(task.phaseTimes - (task.phaseTimes(:,1)));
            phase_start_time_timestamps = average_phase_start_time*Fs; 


            %phase_time_idx = arrayfun(@(x) find(time_phase_labels == x,1), unique(time_phase_labels));
            %Phase_names = {'ITI', 'Image Cue', 'Delay' 'Action'};
            Waveforms_tmp = {};
            Channels_per_trial_timestemps_tmp = {}; 
            
            
            
            for trialNbr = 1:length(trialsToKeep)
                data_tmp = nev.read('WAVEFORMS', 'time', [start_times(trialNbr), end_times(trialNbr)]);
                time_tmp = data_tmp.Timestamps;
                
                if flagSortedData
                    isNotNoise = ~(data_tmp.Units == 255);
                else
                    isNotNoise = logical(ones(length(data_tmp.Units),1));
                end 
                time_tmp = time_tmp(isNotNoise);
                time_adapted = time_tmp - (time_tmp(1)-1); %subtract the starting timetrial to get them all under same frame 
                
                disp(['Check: last time point is ' num2str(time_adapted(end)/ Fs)])

                Channels_timestamps_tmp = arrayfun(@(x) time_adapted(data_tmp.Channels(isNotNoise) == x), unique_channels,'UniformOutput', false);
                Channels_per_trial_timestemps_tmp(:,trialNbr) = Channels_timestamps_tmp; 
                             
                Waveforms_tmp_tmp = arrayfun(@(x) data_tmp.Waveforms(:,data_tmp.Channels(isNotNoise) == x), unique_channels,'UniformOutput', false);
                Waveforms_tmp(:,trialNbr) = Waveforms_tmp_tmp; 
                %plot.raster(Channels_timestamps', 'FS', Fs, 'PSTH',false)
            end
            Channels_per_trial_timestemps = [Channels_per_trial_timestemps; Channels_per_trial_timestemps_tmp];
            Waveforms = [Waveforms; Waveforms_tmp];
            clear Channels_per_trial_timestemps_tmp
            clear Waveforms_tmp

        end 
     end
        
    data_per_brain_region  = arrayfun(@(x) find(electrode_idx == x), unique(electrode_idx), 'UniformOutput', false)';
     brain_region_2 = {'SMG', 'PMV', 'S1X'};

     %brain_region_2 = {'SMG'};

     if flagPlotFigures
           %Plot Trials per Electrode
           for b = 1:size(Channels_per_trial_timestemps,1)
               disp([' Figure ' num2str(b) '/' num2str(size(Channels_per_trial_timestemps,1)) ])
               save_figure_name = [save_name '_Electrode_' num2str(b) '.jpg' ];
               
                PlottingData = Channels_per_trial_timestemps(b,:);
                if ~isempty(PlottingData{1,1})
                TSTR = [' Region: ' brain_regions{electrode_idx(b)} ' Electrode ' num2str(b)]; 
                utile.raster_adapted(Channels_per_trial_timestemps(b,:), 'FS', Fs, 'PSTH',true, 'GROUPS', data_per_group, 'SORT','NONE','TITLE',...
                   TSTR,'grouplabels',labels_name', 'phase_time_changes',phase_start_time_timestamps, 'saveFigureName', save_figure_name ); %, 'GROUPLABELSS',labels_name')
                end 
            end 

           %Plot Electrodes per Trial
           for trialNbr = 1:1:size(Channels_per_trial_timestemps,2)
               TSTR = [' Trial ' num2str(trialNbr) ' - ' labels_per_trial{trialNbr}]; 
               save_figure_name = [save_name '_Trial_' num2str(trialNbr) '.jpg'];
                brain_region_2 = {'SMG'};
                data_per_brain_region = data_per_brain_region(1);
               if isempty(AudioData)
                    utile.raster_adapted(Channels_per_trial_timestemps(:,trialNbr), 'FS', Fs, 'PSTH',true, 'GROUPS', data_per_brain_region, 'SORT','NONE', 'TITLE', TSTR,...
                       'phase_time_changes',phase_start_time_timestamps, 'grouplabels', brain_region_2', 'saveFigureName', save_figure_name, 'AudioData',[] );
               else
                  utile.raster_adapted(Channels_per_trial_timestemps(:,trialNbr), 'FS', Fs, 'PSTH',true, 'GROUPS', data_per_brain_region, 'SORT','NONE', 'TITLE', TSTR,...
                       'phase_time_changes',phase_start_time_timestamps, 'grouplabels', brain_region_2', 'saveFigureName', save_figure_name, 'AudioData',AudioData{trialNbr} );  
               end 
              
           end 
     end 
  
      %%     
           %units_interesting = 73;
           
           %trials_colors = [1,2,3,7,9, 10,12,13, 17, 19,20, 27,28, 37, 38,39 ];
           
       if flagSpikeAnalysis
           trials_colors = 1:40;
           trial_end_time = time_adapted(end)/ Fs; 
           
           for units_interesting = 97 %1:size(Channels_per_trial_timestemps,1)
               for trial_n = 1:length(trials_colors)

                   trial_nbr = trials_colors(trial_n); %trials_colors(1); 

                   timestemps_c = Channels_per_trial_timestemps(units_interesting  ,trial_nbr);
                   waveform_c = cell2mat(Waveforms(units_interesting, trial_nbr)); 
                   all_points = timestemps_c{1,1};

                   uv = unique(all_points);
                   [a, n]  = histcounts(all_points,uv);

                   time = n(1:(end-1))/ Fs; 
                   hold on
                   sorted_all = all_points;

                   %start_time_normal = 7*Fs; 
                   start_time_normal = 1; 

                   end_time_normal = 10*Fs; 
                   %start_time_int = 10*Fs; 


                   %normal_waveform = waveform_c(:,sorted_all < start_idx_Fr); 
                   %interesting_waveform = waveform_c(:,sorted_all > start_idx_Fr); 
                   waveform_normal_idx = logical((sorted_all > start_time_normal).*(sorted_all < end_time_normal));
                   waveform_interesting_idx = ~waveform_normal_idx; 
                   normal_waveform = waveform_c(:,waveform_normal_idx ); 
                   interesting_waveform = waveform_c(:,waveform_interesting_idx); 

                   average_normal_waveform = mean(normal_waveform,2); 
                   zero_corr = zeros(1, size(interesting_waveform,2));
                   for num_int = 1:size(interesting_waveform,2)
                        [C1,lag1] = xcorr(average_normal_waveform,interesting_waveform(:,num_int));        
                        zero_corr(num_int) = C1(lag1 == 0);
                   end 
                    high_corr_idx = (zero_corr > mean(zero_corr));

                   waveFig = figure('units','normalized','outerposition',[0 0 1 1]);
                   subplot(2,2,1); 

                   plot(mean(normal_waveform,2), 'g', 'LineWidth', 2)
                   hold on 
                   plot(mean(interesting_waveform(:,high_corr_idx),2), 'r', 'LineWidth', 2)
                   hold on 
                   plot(mean(interesting_waveform(:,~high_corr_idx),2), 'm', 'LineWidth', 2)

                  legend(['"Normal" Mean Waveform (' num2str(size(normal_waveform,2)) ' Spikes)'], ...
                  ['"High Cross - Correlation Interesting" Mean Waveform (' num2str(size(interesting_waveform(:,high_corr_idx),2)) ' Spikes)'],...
                  ['"Lower Correlation Interesting" Mean Waveform (' num2str(size(interesting_waveform(:,~high_corr_idx),2)) ' Spikes)'],...
                    'Location', 'southoutside');
                   title(['Unit ' num2str(units_interesting), ' Trial ' num2str(trial_nbr)])


                   subplot(2,2,2)
                   plot(time', a, '*');

                   ylim([0 inf])
                   xlim([0 trial_end_time])
                   xlabel('Time')
                   ylabel('Number of Spikes at timepoint')
                   title('Firing Rates over trial')

                   l = line([start_time_normal/Fs,start_time_normal/Fs],[0, max(a)], 'Color', 'black','LineStyle','--');
                   hold on 
                   l2 = line([end_time_normal/Fs,end_time_normal/Fs],[0, max(a)], 'Color', 'black','LineStyle','--');
                   legend(l, 'Cut of time for "normal" vs "interesting" wavelength')

                   subplot(2,2,3); 
                   title('"Normal" Waveform')
                   for i = 1:size(normal_waveform,2)
                       hold on
                       plot(normal_waveform(:,i))
                   end

                   hold on 
                   plot(mean(normal_waveform,2), 'k', 'LineWidth', 2)

                   subplot(2,2,4); 
                   title('"Interesting" Waveform')

                   for i = 1:size(interesting_waveform,2)
                       hold on
                      plot(interesting_waveform(:,i))
                   end 
                   hold on 
                   plot(mean(interesting_waveform,2), 'k', 'LineWidth', 2)


                 %If figure name is provided, save figure
                 blub = ['Unit ' num2str(units_interesting), ' Trial ' num2str(trial_nbr)];

                 save_figure_name_wave = [save_name '_Waveform_Unit_' num2str(units_interesting), '_Trial_' num2str(trial_nbr) '.jpg' ];
                   if ischar(save_figure_name_wave)
                        saveas(waveFig, save_figure_name_wave); 
                        disp(['Figure  '  blub ' saved']);    
                  end 

               end
           end 
           
       end 
 end
        
  
