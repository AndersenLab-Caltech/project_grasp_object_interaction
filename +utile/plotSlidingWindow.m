
clc
clear all
%close all

%Plot figures that evaluate which average bin number is best

%AverageBinNumbers = [2 3 4 5 6 7 8 9 15 20];
%AverageBinNumbers = [2  6 15 20];
AverageBinNumbers = 15
unit_region = 'PMV';
chanceLevel = 1/5*100;
flagOldData = false; 
DataFolder = ['C:\Users\Sarah\Documents\Saved_Data\Data_aligned_correctly\CueComparison\DecodingMethods\' unit_region '\PCA_90_percent_bin_per_bin_']; 
yLimVal = 100; 



for blub = 1:length(AverageBinNumbers)
    AverageBinNumber = AverageBinNumbers(blub); % +1  
    
    if ~flagOldData
        load([DataFolder num2str(AverageBinNumber) '_average.mat'])
       [numb,NameTmp,raw] = xlsread('C:\Users\Sarah\Dropbox\Code\project_analysis_grasps\ExcelFilesNew\s2_good_trials_ImageCue_AuditoryCue_WrittenCue_2s'); 


    else
        load([DataFolder num2str(AverageBinNumber) '_average_OldData.mat'])

    end

    size(DataAll(1).TrainError)
    AverageBinNumber = AverageBinNumber +1; % +1  

    %% Figure over time for individual cues, mean over all sessions + individual
    CommonSession = DataAll(1).Sessions(ismember(DataAll(1).Sessions, DataAll(2).Sessions));
    if ~flagOldData %old data does not have written cue 
        CommonSession = CommonSession(ismember(CommonSession, DataAll(3).Sessions));
    end 
    numPhases = size(DataAll(1).TrainError,3);

    DataPerSession = []; 
    clear DataPerSession
    clear DataMean

    for i = 1:size(DataAll,2)
        DataName= DataAll(i).Name;
        CommonIdx = ismember(DataAll(i).Sessions,CommonSession);
        TestingData = DataAll(i).TestError(:,:,:,CommonIdx);
        TimePhaseLabels = DataAll(i).TimePhase(1);
        numSessions = size(TestingData,4);
        figure();
        BinPhaseChanges = find(diff(TimePhaseLabels{1,1})) - AverageBinNumber ;
        DataPlot = zeros(numPhases, numSessions);

        %figure();
        for SessionNbr = 1:numSessions   
            %subplot(2,5,SessionNbr)
            time = (1:length(DataPlot))*0.05; %rate = 0.05
            time = time + AverageBinNumber*0.05;
            DataPlot(:,SessionNbr) = (1-squeeze(mean(mean(TestingData(:,:,:,SessionNbr)))))*100;   
            plot(time,DataPlot(:,SessionNbr))
            hold on
            DataPerSession(i,:,SessionNbr) = DataPlot(:,SessionNbr);
        end 
        hold on
        DataMean(:,i) = mean(DataPlot,2);
        a = plot(time,DataMean(:,i), 'k', 'LineWidth', 2);
        ylim([0, yLimVal])

        for i = 1:length(BinPhaseChanges) %(data.phase_time(2) - data.phase_time(1))+1
            yL = get(gca,'YLim');
            l(i) = line([time(BinPhaseChanges(i)) time(BinPhaseChanges(i))],yL,'Color','k','LineStyle','--','Linewidth', 2);
        end
        
        c = line([0 time(end) ],[chanceLevel,chanceLevel],'Color','r','LineStyle','--','Linewidth', 0.75);

            
        title([DataName ' - Averaged over ' num2str(AverageBinNumber*0.05) 's'])
        legend([a l(1) c],[{DataAll.Name}, 'Phase Change','Chance Level'] )
        xlabel('Time [s]')
        ylabel('Classification accuracy [%]')
        xlim([0, time(end)])
        clear TestingError

    end 

    %% Figure mean over session for each individual cue 

    figure(); 

    for i = 1:size(DataMean,2)
        a(i) =plot(time,DataMean(:,i), 'LineWidth', 2);
        hold on

    end
    ylim([0 yLimVal]);

    for i = 1:length(BinPhaseChanges) %(data.phase_time(2) - data.phase_time(1))+1
            yL = get(gca,'YLim');
            l(i) = line([time(BinPhaseChanges(i)) time(BinPhaseChanges(i))],yL,'Color','k','LineStyle','--','Linewidth', 2);
    end
    c = line([0 time(end) ],[chanceLevel,chanceLevel],'Color','r','LineStyle','--','Linewidth', 0.75);

    title(['Compare onset of decoding for all three cues - Averaged over ' num2str(AverageBinNumber*0.05) 's'])
    legend([a l(1) c],[{DataAll.Name}, 'Phase Change', 'Chance Level'] )
    xlabel('Time [s]')
    ylabel('Classification accuracy [%]')
    xlim([0, time(end)])
    %% for each session, each different cue type + number of recording 

    figure()

    for i = 1:size(DataPerSession,3)
        subplot(2,4,i)
        ylim([0, 100])

        for j = 1:size(DataPerSession,1)
            hold on
             a(j) = plot(time,DataPerSession(j,:,i), 'LineWidth' ,1.5);
        end 

        for j = 1:length(BinPhaseChanges) %(data.phase_time(2) - data.phase_time(1))+1
            yL = get(gca,'YLim');
            l(j) = line([time(BinPhaseChanges(j)) time(BinPhaseChanges(j))],yL,'Color','k','LineStyle','--','Linewidth', 2);
        end
        
        c = line([0 time(end) ],[chanceLevel,chanceLevel],'Color','r','LineStyle','--','Linewidth', 0.75);

        title(CommonSession(i))
        
        if ~flagOldData
           labels2 = arrayfun(@(x) [DataAll(x).Name ' (ReNbr: ' num2str(numb(2*i -1,x+1)) ')'],1:size(DataAll,2) , 'UniformOutput', false);
        else
           labels2 = arrayfun(@(x) [DataAll(x).Name ],1:size(DataAll,2) , 'UniformOutput', false);
        end 
        
        legend(a,labels2 )
        xlim([0, time(end)])

    end 
    sgtitle(['Compare onset of decoding for all three cues per individual Session - Averaged over ' num2str(AverageBinNumber*0.05) 's'])
    xlabel('Time [s]')
    ylabel('Classification accuracy [%]')
    
    %% compare decoding by recording number
    if ~flagOldData

        figure();
        ColorsTest = {'g','b','r', 'r'};

        for blub = 1:size(DataPerSession,1)
            subplot(3,1,blub)
            RecordNbr = numb(:,1+blub);
            RecordNbr = RecordNbr(~isnan(RecordNbr));
            DataTmp = squeeze(DataPerSession(blub,:,:));
            DataPerRecNbr = arrayfun(@(x) DataTmp(:,RecordNbr == x), unique(RecordNbr), 'UniformOutput', false);

           for i = 1:size(DataTmp,2)
               hold on 
               plot(time', DataTmp(:,i),[ ColorsTest{RecordNbr(i)} ':'],'LineWidth', 1 );

           end 

           for j = 1:length(BinPhaseChanges) %(data.phase_time(2) - data.phase_time(1))+1
                yL = get(gca,'YLim');
                l(j) = line([time(BinPhaseChanges(j)) time(BinPhaseChanges(j))],yL,'Color','k','LineStyle','--','Linewidth', 2);
           end


            c1 = plot(time', mean(DataPerRecNbr{1},2),'g', 'LineWidth', 2 )
            hold on
            c2 = plot(time', mean(DataPerRecNbr{2},2),'b', 'LineWidth', 2)
            hold on
            c3 = plot(time', mean(DataPerRecNbr{3},2),'r', 'LineWidth', 2)
            title(DataAll(blub).Name)
            xlim([0, time(end)])

        end 
        sgtitle(['Average decoding accuracy based on which recording it was on the specific day'])
        xlabel('Time [s]')
        ylabel('Classification accuracy [%]')

        legend([c1 c2 c3 ] , {'First recording', 'Second recording', 'Third recording' })
    
    end 
end 