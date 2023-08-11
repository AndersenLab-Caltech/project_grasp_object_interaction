%Grasp project Sarah David Mac

clc 
clear all
close all

subject_id = 's3';  %AN
%subject_id = 's2'; %FG

subject = hst.Subject(subject_id);

%threshold session using blackrock filter. You can define the threshold
%here. 
Blackrock.thresholdSession('20230721', 's3', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')
Blackrock.thresholdSession('20230803', 's3', 'THRESHOLD', -4.5, 'noise_model', 'Blackrock')

%%

tasktype = {};
end_comment = {};
session_dates = {'20230810'};
taskfileTest = 2;

idxToRemove = [];
for blub= 1:length(session_dates)
    session_date = session_dates{blub};
    session = hst.Session(session_date, subject);
    taskfiles = session.getTaskFiles('Speech');
    %taskfiles = session.getTaskFiles('SpeechTrainingInternal');

    %taskfiles = session.getTaskFiles('RecordOnly');
    %h = SpikeSorter.GUI(fileparts(taskfiles{1})) %Select the nev files to import them. If several sessions the same day, select them all at the same time


    for i = 1:length(taskfiles)

        try
            task = hst.Task(taskfiles{i});
            tasktype{blub,i} = task.task.params.user.cue_type;
            %end_comment{blub,i} = task.userEndComment;  
            end_comment{blub,i} = task.parameterName;
            
            if size(task.trialparams,2) ~= task.numTrials
                disp('probably aborded trial')
                idxToRemove = [idxToRemove, i];
            end 
        catch
            disp('???')
            %probably need to remove this one too 
            idxToRemove = [idxToRemove, i];
           % keyboard
           % error('Probabily need to change the file names and remove the last three numbers')
        end 
    end
    end_comment;
end  
%end 

good_blocks = setdiff(1:length(taskfiles), idxToRemove);
good_blocks
taskfiles(good_blocks)

 keyboard
h = SpikeSorter.GUI(fileparts(taskfiles{1}))

task = hst.Task(taskfiles{1});



blub = [session_dates' tasktype end_comment]

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
 