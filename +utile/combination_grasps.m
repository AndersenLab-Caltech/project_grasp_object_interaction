function [uni_row_comb_tester,Figure_labels] = combination_grasps( Figure_labels, varargin )

[varargin,n] = Utilities.argkeyval('n',varargin, 26); %Inputs: ITI, ImageCue, Delay, Action
[varargin,k] = Utilities.argkeyval('k',varargin, 5); %Inputs: ITI, ImageCue, Delay, Action

Utilities.argempty(varargin);



Figure_labels.Figure_Names = cell(n,1);
uni_row_comb_tester = zeros(n,6);
count = 0; 

for i = 2:size(uni_row_comb_tester,2) 
    C = nchoosek(1:k,i);
    for j = 1:size(C,1)
        count = count +1;
        Figure_labels.Figure_Names{count} = [num2str(i) 'Grasps_Fig' num2str(j) '.mat'];
        uni_row_comb_tester(count, C(j,:)) = 1;
    end  
end

end

