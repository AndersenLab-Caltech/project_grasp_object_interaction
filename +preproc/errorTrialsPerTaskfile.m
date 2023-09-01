function wrongTrials = errorTrialsPerTaskfile(task,subject)
%Returns the indexes of wrong trials for each session day

taskfile = task.taskString;
 

if strcmp(subject, 's2')
    if strcmp(taskfile, '20230720-164128-164400-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230720-165606-165608-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230720-170055-170229-GraspObject')
        wrongTrials = [18,19,20,21,27,28];
    elseif strcmp(taskfile, '20230720-170740-170814-GraspObject')
        wrongTrials = [13,14,20,21];
    elseif strcmp(taskfile, '20230720-170740-171307-GraspObject')
        wrongTrials = [7,8,9,16,17,22];
    elseif strcmp(taskfile, '20230720-171816-171832-GraspObject')
        wrongTrials = [11,12,13,14,21,22,23,29];

    elseif strcmp(taskfile, '20230725-154852-155152-GraspObject')
        wrongTrials = [13,14,30,34,35,36];
    elseif strcmp(taskfile, '20230725-155855-155923-GraspObject')
        wrongTrials = [6,7,21,25,31];
    elseif strcmp(taskfile, '20230725-160457-160518-GraspObject')
        wrongTrials = [23,24,36];

    elseif strcmp(taskfile, '20230803-165329-165754-GraspObject')
        wrongTrials = [7,8,11,16,17,18,19,20];
    elseif strcmp(taskfile, '20230803-165329-170230-GraspObject')
        wrongTrials = [8,11,12];
    elseif strcmp(taskfile, '20230803-170617-170806-GraspObject')
        wrongTrials = [4,5,11,12,18];
    elseif strcmp(taskfile, '20230803-170617-171133-GraspObject')
        wrongTrials = [4,7,8,9,16,17];
    elseif strcmp(taskfile, '20230803-171500-171634-GraspObject')
        wrongTrials = [3,13,17,18];
    elseif strcmp(taskfile, '20230803-171500-172007-GraspObject')
        wrongTrials = [5,12,17];

    elseif strcmp(taskfile, '20230810-163923-163940-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230810-163923-164343-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230810-164655-164723-GraspObject')
        wrongTrials = [14];
    elseif strcmp(taskfile, '20230810-164655-165027-GraspObject')
        wrongTrials = [9,10];
    elseif strcmp(taskfile, '20230810-165337-165355-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230810-165337-165705-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20230817-154012-154020-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230817-154012-154432-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230817-154806-154814-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230817-154806-155138-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230817-155456-155505-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230817-155456-155910-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20230824-155938-155946-GraspObject')
        wrongTrials = [7,14];
    elseif strcmp(taskfile, '20230824-155938-160341-GraspObject')
        wrongTrials = [18];
    elseif strcmp(taskfile, '20230824-160653-160701-GraspObject')
        wrongTrials = [7,10];
    elseif strcmp(taskfile, '20230824-160653-161048-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230824-161418-161426-GraspObject')
        wrongTrials = [19];
    elseif strcmp(taskfile, '20230824-161418-161834-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20230831-164149-164156-GraspObject')
        wrongTrials = [18];
    elseif strcmp(taskfile, '20230831-164149-164554-GraspObject')
        wrongTrials = [8,14,15,17,20];
    elseif strcmp(taskfile, '20230831-164956-165130-GraspObject')
        wrongTrials = [10,17,18];
    elseif strcmp(taskfile, '20230831-164956-165520-GraspObject')
        wrongTrials = [4,8,13];
    elseif strcmp(taskfile, '20230831-165918-170037-GraspObject')
        wrongTrials = [17,18];
    elseif strcmp(taskfile, '20230831-165918-170420-GraspObject')
        wrongTrials = [];

 %  elseif strcmp(taskfile, '')
 %       wrongTrials = [];

    else
        
        disp(task.task.params.user.cue_type)
        disp(taskfile)
        keyboard
        error([' Unknown session day - add to errorTrials'])    
      
    end

elseif strcmp(subject, 's3')
    
     if strcmp(taskfile,   '20230721-095625-100057-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20230721-095625-100644-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20230721-101242-101247-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20230721-101242-101923-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20230721-102408-102434-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20230721-102408-103034-GraspObject')
        wrongTrials = []; 

    elseif strcmp(taskfile,'20230724-111804-111856-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230724-111804-112515-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230724-113033-113037-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230724-113033-113625-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230724-114426-114430-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230724-114426-115214-GraspObject')
        wrongTrials = []; 
   
    elseif strcmp(taskfile,'20230803-095402-095720-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230803-095402-100940-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230803-101445-101511-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230803-101445-102359-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230803-102859-103015-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230803-102859-103540-GraspObject')
        wrongTrials = []; 

        % elseif strcmp(taskfile,'')
    %    wrongTrials = [];   
        
     else
        keyboard
        error([taskfile ' Unknown session day - add to errorTrials'])    
      
       
    end

    
    
    
    
end






