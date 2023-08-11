function [ err ] = classerror( y,yhat )

% DO NOT USE IT. THE CLASSIFICATION ERROR CAN VARY WHEN THE CLASSES ARE NOT
% EQUALLY DISTRIBUTED, AS IT IS NOT A WEIGHTED ERROR. IDK what is best tbh
% but do not use


%CLASSERROR computes the classification error for each class and averages
%it
%   
%   Input:
%       y:      the labels
%       yhat:   the output of the classifier
%
%   output:
%       err:    the class-averaged classification error
classes = unique(y);
err_ = zeros(1,length(classes));
rep = histc(y, unique(y))
for c=1:length(classes)
    err_(c) = sum((y~=yhat) & (y == classes(c)))./sum(y==classes(c));% + sum((y~=yhat) & (y == 1))./sum(y==1));
end


%err_
err = mean(err_);
end


