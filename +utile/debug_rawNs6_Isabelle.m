%clear all
%close all
% task resource

%%
GraspCue = 'Speech';

session_dates = {'20210923', '20210930', '20211011'};

session_date = '20210923';
session_date = '20210930';
session_date = '20211011';

session_dates = {'20210923', '20210930'};


%%


save_data_pathway = 'D:\Users\Sarah\Documents\Saved_Data\s2_2021_speech_aligned_correctly';

subject_id = 's2';  % s2 or p3
subject = hst.Subject(subject_id);

brain_regions = {'SMG_AIP'};

for n_sess = 1:length(session_dates)
    
    session_date = session_dates{n_sess};
    ChannelsToPlot = get_channels_to_plot(session_date);

    save_fig_folder = ['D:\Users\Sarah\Pictures\RasterPlots\Richard\' session_date '\RawData2\']; 


    if ~exist(save_fig_folder, 'file')
        mkdir(save_fig_folder);
    end 

    session = hst.Session(session_date, subject);

    taskfiles = session.getTaskFiles('Speech');

    dataZ = env.get('data');

    filename = [subject_id '_good_trials_' GraspCue '.xlsx'];

    [numb,txt,raw] = xlsread(['C:\Users\Sarah\Dropbox\Code\project_speech\ExcelFiles\' filename]);     
    session_date = str2double(session_date);

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
        disp(['No auditory dataset for session ' num2str(session_date)])       
    else
        disp('Everything seems ok');
    end

    DatesRemoveTrials = {'20191016'};

    cueTypes = {'Auditory', 'Written'};
    flagPlotAuditory = true; %adds the auditory signal on top (only if present)

    %%
    debug = Debug.Debugger('Speech');
    subject = hst.Subject('S2');
    map = Blackrock.ArrayMap(subject.hArray{1}.MapFile);

    taskfiles = session.getTaskFiles('Speech');

    for n_task = 1:length(good_datasets)

        taskfile_nev = taskfiles{good_datasets(n_task)};

        task_nev = hst.Task(taskfile_nev,debug);
        task = hst.Task(taskfile_nev,debug);

        ns = task_nev.getNeuralDataObject('SMG','ns6'); 
        nev2 = task_nev.getNeuralDataObject('SMG','nev'); 

        ns = ns{1};
        % get relevant trial timing 



        % create spike filter
        %this is the blackrock filter
        filtobj{1} = Blackrock.getSpikeContentFilter('butter_o4_bp_250-5000');
        [b,a] = iirnotch(60/(ns.Fs/2),500/(ns.Fs/2)); % 60 notch filter with q factor 10

        nTrials = length([task.trialdata.et_trialStart]);


        Fs = 30000;

        average_phase_start_time = mean(task.phaseTimes - (task.phaseTimes(:,1)));

        %average_phase_start_time = [ 0, 2.0874,3.6502,4.2151, 5.7896, 6.3534];


       PhaseNames = {'ITI' ,'Cue', 'D1', 'Internal', 'D2', 'Speech'};
        for n_trial = 1:nTrials
            t_trialStart = task.data.neuralTime(task.trialdata(n_trial).et_trialStart);
            t_trialEnd = task.data.neuralTime(task.trialdata(n_trial).et_trialCompleted);

            % get raw ns6 data with 10s buffer on either end
            dt_raw = ns.read('Time',[t_trialStart (t_trialEnd)])';
            %get spike waveforms 
            nev2_waveform = nev2{1}.read('WAVEFORMS', 'time', [(t_trialStart +average_phase_start_time(4)),(t_trialStart +average_phase_start_time(5))])';




            dt_filt = filter(filtobj{1},dt_raw);
            time = (1:size(dt_filt,1))/Fs;

            for nChannel = 1:length(ChannelsToPlot)
                channel = ChannelsToPlot(nChannel);
                fig = figure('units','normalized','outerposition',[0 0 0.7 0.4]);

                subplot(1,2,1)

                plot(time,dt_filt(:,channel));

                for n_phase = 1:length(average_phase_start_time)
                    title([num2str(session_date) ' - Channel ' num2str(channel) ' Trial # ' num2str(n_trial) ' ' cueTypes{n_task}])
                    xline(average_phase_start_time(n_phase), 'k--',PhaseNames{n_phase},'linewidth', 1, 'fontsize',12)
                    hold on
                    xline((average_phase_start_time(4)), 'r','linewidth', 1);
                    hold on 
                    xline((average_phase_start_time(5)), 'r','linewidth', 1);

                end 
                ylabel('Filtered signal')

                subplot(1,2,2)
                 plot(nev2_waveform.Waveforms(:,nev2_waveform.Channels == channel))
                ylabel('Waveforms')
                title(['InternalSpeech waveforms - ' task.task.TrialParams(n_trial).action])

        %         blub = dt_filt(:,channel);
        %         blub((dt_filt(:,channel) < 18)) = 0;
        %         figure(); plot(time,blub);
        %         


                save_figure_name = [save_fig_folder 'FilteredData_Channel_' num2str(channel) '_Trial_' num2str(n_trial) '_' cueTypes{n_task} '.jpg'];
                saveas(fig, save_figure_name); 
            end 
             close all

            %xline(
        end 
        close all

    end 
end 

%%
function channels_to_plot = get_channels_to_plot(session_date)

    if strcmp(session_date, '20210923')
        channels_to_plot = [15,81,83,85,88,90,93,96];
    elseif strcmp(session_date, '20210930')
         channels_to_plot = [35,63,77,79,88,90,93];
    elseif strcmp(session_date, '20211011')
         channels_to_plot =[71,74,75,83,85,88];
    end 
end 



% for i = 1:6
%     subplot(2,3,i)
%     
%     plot(time,dt_filt(:,i));
% 
%     
%     
% end 


















% % dt_filt2 = filter(filtobj{2}{1}, filtobj{2}{2}, dt_filt);
%   dt_filt2 = filter(filtobj{3}{1}, filtobj{3}{2}, dt_filt); % 180
% %  dt_filt2 = filter(filtobj{4}{1}, filtobj{4}{2}, dt_filt2); % 300
%   
%  % dt_filt2 = filter(filtobj{5}{1}, filtobj{5}{2}, dt_filt2); % 420
% % dt_filtfilt = filtfilt(filtobj,dt_raw);
% 
% t_samples_rel = 0:((2*ns.Fs)-1);
% 
% ch = 8
% figure
% subplot(1,2,1)
% plot(t_samples_rel,dt_filt(:,ch))
% hold on
% plot(t_samples_rel,dt_filt2(:,ch))
% subplot(1,2,2)
% dt = dt_filt(~isnan(dt_filt(:,ch)),ch);
% dt2 = dt_filt2(~isnan(dt_filt2(:,ch)),ch);
% [p1,f1] = pspectrum(dt,ns.Fs,'power','FrequencyLimits',[0 1000]);
% plot(f1,p1)
% hold on
% [p2,f2] = pspectrum(dt2,ns.Fs,'power','FrequencyLimits',[0 1000]);
% plot(f2,p2)
% xlabel('freq')
% ylabel('power')
% 
% %% collate trials
% numTr = size(task.trialdata,2);
% t_samples_rel = 0:((2*ns.Fs)-1);
% for t = 1:10%numTr
%     t_trialStart = task.data.neuralTime(task.trialdata(t).et_trialStart);
%     t_trialEnd = task.data.neuralTime(task.trialdata(t).et_trialCompleted);
%     
%     % get raw ns6 data with 10s buffer on either end
%     dt_raw = ns.read('Time',[t_trialStart (t_trialStart+2)])';
%     dt_filt = filter(filtobj{1},dt_raw);
% %     dt_filt2 = filter(filtobj{3}{1}, filtobj{3}{2}, dt_filt);
% %     dt_filt2 = filter(filtobj{4}{1}, filtobj{4}{2}, dt_filt2);
%  %   dt_filt2 = filter(filtobj{4}{1}, filtobj{4}{2}, dt_filt2);
%     
%     for ch = 1:96
%         dt = dt_filt(~isnan(dt_filt(:,ch)),ch);
%         [p1,f1] = pspectrum(dt,ns.Fs,'power','FrequencyLimits',[0 1000]);
%         chP(ch,t) = max(p1);
%         chF(ch,t) = f1(p1==max(p1));
%     end
% end
% medP = median(chP'); medF = median(chF');
% figure('Position',[400 500 1250 450])
% bar(medP)
% ylim([0 200])
% ylabel('median max(power)')
% yyaxis right
% plot(medF,'ok', 'MarkerFaceColor','k')
% ylim([0 500])
% ylabel('median frequency with highest power')
% 
% 
% 
% % for each channel
% % go over first 10 trials
% % get the power for the first 2 seconds of the trial
% % if at any point the power goes over 20, round the frequency to the
% % nearest multiple of 5 and log
% 
% % end result: for each session, for each channel, median max power for each
% % trial, and median max freq
% 
% %the % trials where power
% % went too high, and the median power where this occurred
% 
% % 180 Hz artifact
% % could try to eliminate any channels with artifact above 50 power?
% % 
% %caxis([0 1])