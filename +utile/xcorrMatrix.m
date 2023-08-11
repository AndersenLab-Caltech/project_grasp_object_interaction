function [outputArg1] = xcorrMatrix(AudioSignal)
%returns correlation matrix of audio signals
%** not needed, becaue the function corrcoef will immediately give me a
%very similar value and can take into account matrices, so it's easier
outputArg1 = [];
outputArg2 = [];

for i = 1:length(AudioSignal)
    
    for j = 1:length(AudioSignal)
            outputArg1(i,j) = xcorr(AudioSignal{i}, AudioSignal{j},0,'normalized');
            %outputArg2(i,j) = corrcoef(AudioSignal{i}, AudioSignal{j});
    end
    
end
mean(mean(outputArg1))

blub = cell2mat(AudioSignal);
blub2 = corrcoef(blub);

outputArg1 - blub2;

end

