function [yCI95] = calculate_CI(data)
%calculate 95% confidence interval of provided data
    N = size(data, 1);
    sem = std(data) / sqrt(N);  % standard error of the mean
    %sem = cellfun((@(x) std(x)), dataTmp, 'UniformOutput', false); % this just calcs the SEM within each cell, trying to go between cells?
    CI95 = tinv([0.025 0.975], N-1);  % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, sem, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
end