function wrongTrials = errorTrialsPerTaskfile(task,subject)
%Returns the indexes of wrong trials for each session day

taskfile = task.taskString;
 

if strcmp(subject, 's2')
    if strcmp(taskfile, '')
        wrongTrials = [];
   
  %  elseif strcmp(taskfile, '')
 %       wrongTrials = [];
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






