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

    elseif strcmp(taskfile,   '20230721-101242-101923-GraspObject')
        wrongTrials = [];   
        
     
    elseif strcmp(taskfile,'20230721-102408-102434-GraspObject')
        wrongTrials = [];   


    elseif strcmp(taskfile,'20230721-102408-103034-GraspObject')
        wrongTrials = [];      


        % elseif strcmp(taskfile,'')
    %    wrongTrials = [];   
        
     else
        keyboard
        error([taskfile ' Unknown session day - add to errorTrials'])    
      
       
    end

    
    
    
    
end






