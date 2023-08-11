function [Output] = ZScoreCell(Input,Labels,varargin)
%function that Returns the Z Score of a matrix in cell form. The input argument can
%chose which type of Z scoring is chosed

%Z-scoring: compute the mean of the entire data of one channel/unit over one block and substract it from 
%each activity, and then divide by the standard deviation of the entire block.
    
%Chose normalization method. median or mean
[varargin,ZscoreMethod] = Utilities.argkeyval('ZScoreMethod',varargin, 'mean'); 
[varargin,flagPlotFigures] = Utilities.argkeyval('PlotFigures',varargin, 0); 

Utilities.argempty(varargin);
U_Labels= unique(Labels);
colors = utile.get_color_rgb_codes(utile.label_number_to_grasp_name(U_Labels));





if isa(Input, 'cell')
    
    %Test if empty - if yes, do not apply
    if isempty(Input) || isempty(find(cellfun(@isempty, Input) ==0))
        disp('Empty Input')
        Output = Input; 
    return 
    end

    %MCellNew = cell(size(Input));
    MCellNew = Input;

        for channelNbr = 1:size(Input{1},2)
            M = cell2mat(cellfun(@(x) x(:,channelNbr), Input, 'UniformOutput', false)');
            if strcmp( ZscoreMethod, 'median')
                Z = ( M - median(M(:)) ) / mad(M(:));

            elseif strcmp(ZscoreMethod, 'mean')
                Z = ( M - mean(M(:)) ) / std(M(:));
            else
                error('Not recognized method')
            end
                
            for trialNbr = 1:size(MCellNew,1)
                try
                    MCellNew{trialNbr}(:,channelNbr) = Z(:,trialNbr);
                catch
                    disp('hello')
                end 
            end 

            for i = 1:length(U_Labels)
                MAll(:,i,channelNbr) = mean(M(:,Labels == U_Labels(i)),2);
                ZAll(:,i,channelNbr) = mean(Z(:,Labels == U_Labels(i)),2);

            end
        end 
 
    %plot the average over all channel to compare 
    figure();
    subplot(2,1,1)

    time = 0.05*(1:length(M));
    for i = 1:length(U_Labels)
        hold on
        plot(time,mean(MAll(:,i,:),3), 'Color', colors{i});
    end
    title('Original')
    
    subplot(2,1,2)
    for i = 1:length(U_Labels)
        hold on
        plot(time,nanmean(ZAll(:,i,:),3), 'Color', colors{i});
    end
    title('Z-scored')
    sgtitle('Firing rate averaged over entire session')
    
    Output = MCellNew;
elseif isa(Input, 'double')
    %Z-scoring: compute the mean of the entire data of one channel/unit over one block and substract it from 
    %each activity, and then divide by the standard deviation of the entire block.
    M = Input;
    
    if strcmp( ZscoreMethod, 'median')
        Z = ( M - median(M(:)) ) / mad(M(:));

    elseif strcmp(ZscoreMethod, 'mean')
        Z = ( M - mean(M(:)) ) / std(M(:));
    else
        error('Not recognized method')
    end
    Output = Z;
    %plot the average over the channel to compare
    
    if flagPlotFigures
        figure();
        subplot(2,1,1)

        time = 0.05*(1:length(M));
        for i = 1:length(U_Labels)
            hold on
            plot(time,mean(M(:,Labels == U_Labels(i)),2), 'Color', colors{i});
        end

        subplot(2,1,2)
        for i = 1:length(U_Labels)
            hold on
            plot(time,mean(Z(:,Labels == U_Labels(i)),2), 'Color', colors{i});
        end   
    end 
end

end

