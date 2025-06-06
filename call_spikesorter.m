%Grasp project Sarah David Mac
 
clc 
clear all
close all

%subject_id = 's4'; %GB
%subject_id = 's3';  %AN
subject_id = 's2'; %FG

subject = hst.Subject(subject_id);

%threshold session using blackrock filter. You can define the threshold
%here. 
% Blackrock.thresholdSession('20230721', 's3', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
% Blackrock.thresholdSession('20230803', 's3', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
% Blackrock.thresholdSession('20230724', 's3', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')

%Blackrock.thresholdSession('20250517', subject_id, 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
% Blackrock.thresholdSession('20250519', subject_id, 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
Blackrock.thresholdSession('20250520', subject_id, 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')

% Blackrock.thresholdSession('20230725', 's2', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
% Blackrock.thresholdSession('20230803', 's2', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
% Blackrock.thresholdSession('20230810', 's2', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')

%%
  keyboard
%%

tasktype = {};
end_comment = {};
session_dates = {'20250519'};
taskfileTest = 2;

idxToRemove = [];
for n_session_dates= 1:length(session_dates)
    session_date = session_dates{n_session_dates};
    session = hst.Session(session_date, subject);
    taskfiles = session.getTaskFiles('GraspObject');
    
    for i = 1:length(taskfiles)

        try
            task = hst.Task(taskfiles{i});
            tasktype{n_session_dates,i} = task.task.params.user.cue_type;
            %end_comment{blub,i} = task.userEndComment;  
            end_comment{n_session_dates,i} = task.parameterName;
            
            if size(task.trialparams,2) ~= task.numTrials
                disp('probably aborted trial')
                idxToRemove = [idxToRemove, i];
            end 
        catch ME
            disp('???')
            idxToRemove = [idxToRemove, i];
        end 
    end
    end_comment;
end  

good_blocks = setdiff(1:length(taskfiles), idxToRemove);
good_blocks
taskfiles(good_blocks)
%%
 keyboard
 
h = SpikeSorter.GUI(fileparts(taskfiles{1}))

task = hst.Task(taskfiles{1});



n_session_dates = [session_dates' tasktype end_comment]

blub2 = 'Z:\Data\s2\20190423\PMV\20190423-144935-145018-GraspShape.mat'
%Check the endcomments to find out which sessions to combine while spike
%sorting
%end
h = SpikeSorter.GUI(fileparts(blub2))

%h = SpikeSorter.GUI(fileparts(taskfiles{1}), 'subject', 's3')
h = SpikeSorter.GUI(fileparts(taskfiles{1}));

%for new rig: seems like 

session = hst.Session('20171003', 's2');
taskfiles = session.getTaskFiles('GraspShape');
task = hst.Task(taskfiles{1});

%first and second are 'good' 
 