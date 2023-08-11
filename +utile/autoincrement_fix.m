% remove values of script 
session_date = '20230124';
array_names = {'MPx', 'S1X_S1'};

for n_array = 1:length(array_names)
    
    data_file = ['Z:\Data\S3\' session_date '\' array_names{n_array}];

    taskfiles = dir(data_file);

    for id = 1:length(taskfiles)

        if taskfiles(id).bytes > 0
            % Get the file name 
            [~, f,ext] = fileparts(fullfile(data_file,taskfiles(id).name));
            %split file based on character '-'
            C = strsplit(f,'-');
            %check length of last element 
            number_elements =  C{end};
            %if it is 6 -> autoincrement was on
            if length(number_elements) == 6
                %rename file by removing last 3 elements 
                rename = strcat(f(1:(end-3)),ext);
                movefile(fullfile(data_file,taskfiles(id).name), fullfile(data_file,rename)); 
            end 
        end 
    end

end 
