function wrongTrials = errorTrialsPerTaskfile(task,subject)
%Returns the indexes of wrong trials for each session day
%cueforSpeech = 'Written';
%cueforSpeech = 'Speech';
%cueforSpeech = 'All';



taskfile = task.taskString;
 

if strcmp(subject, 's2')
    if strcmp(taskfile, '20210712-135346-135354-Speech')
        wrongTrials = [18,62];
    elseif strcmp(taskfile, '20210712-135346-141239-Speech')
        wrongTrials = [];
    elseif strcmp(taskfile, '20210722-145830-145837-Speech')
        wrongTrials = [];        
    elseif strcmp(taskfile, '20210722-144833-144835-Speech')
        wrongTrials = [11];
        
    elseif strcmp(taskfile, '20210729-120048-120627-Speech')
        wrongTrials = [57];
        
  elseif strcmp(taskfile, '20210729-114822-114949-Speech')
        wrongTrials = [15,33];
   elseif strcmp(taskfile, '20210923-204418-204521-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile,  '20210923-210053-210206-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20210930-201806-201816-Speech')
        wrongTrials = [16,33,45,58];
   elseif strcmp(taskfile, '20210930-200549-200623-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20211011-200846-200951-Speech')
        wrongTrials = [52];
   elseif strcmp(taskfile, '20211011-195707-195721-Speech')
        wrongTrials = [51,60];
   elseif strcmp(taskfile, '20211018-194758-194912-Speech')
        wrongTrials = [48];
   elseif strcmp(taskfile, '20211018-200142-200212-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20211027-200643-200830-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20211027-202500-203440-Speech')
        wrongTrials = [23,30,31,37,42,45,52,55];
   elseif strcmp(taskfile, '20211103-193909-194312-Speech')
        wrongTrials = [60];
   elseif strcmp(taskfile, '20211103-192911-192924-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20220323-203042-203053-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20220323-201236-202005-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20211230-194323-194559-Speech')
       wrongTrials = [];
   elseif strcmp(taskfile, '20220523-200117-200505-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20220523-201328-201410-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20220523-202637-202645-Speech')
        wrongTrials = [];
   elseif strcmp(taskfile, '20220523-203755-203811-Speech')
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
    
     if strcmp(taskfile,'20230106-093139-093221-SpeechTrainingInternal')
        wrongTrials = [];   
        
     elseif strcmp(taskfile,'20230106-094536-094652-SpeechTrainingInternal')
        wrongTrials = [];   
        
     elseif strcmp(taskfile,'20230106-095152-095224-SpeechTrainingInternal')
        wrongTrials = [];   
        
    elseif strcmp(taskfile, '20230706-095711-095719-Speech')
        wrongTrials = [49];
        
    elseif strcmp(taskfile, '20230706-101104-101115-Speech')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230706-102545-102549-Speech')
        wrongTrials = [40];
    elseif strcmp(taskfile,'20230706-103527-103555-Speech')
        wrongTrials = [58];   
        
    elseif strcmp(taskfile,'20230712-095444-095450-Speech')
        wrongTrials = [55];   
            
    elseif strcmp(taskfile,'20230712-100419-100427-Speech')
        wrongTrials = [8];   
                
    elseif strcmp(taskfile,'20230712-112848-112851-Speech')
        wrongTrials = [50];   
                
    elseif strcmp(taskfile,'20230713-094451-094549-Speech')
        wrongTrials = [4];   
                
    elseif strcmp(taskfile,'20230713-100153-100201-Speech')
        wrongTrials = [];   
                
    elseif strcmp(taskfile,'20230713-113252-113309-Speech')
        wrongTrials = [];   
        
     elseif strcmp(taskfile,'20230713-114225-114245-Speech')
         wrongTrials = [];   
         
 elseif strcmp(taskfile,'20230810-102731-102739-SpeechTrainingInternal')
         wrongTrials = [];
               
  elseif strcmp(taskfile,'20230810-103657-104513-SpeechTrainingInternal')
        wrongTrials = [];
   % elseif strcmp(taskfile,'')
    %    wrongTrials = [];   
        
    else
        error([taskfile ' Unknown session day - add to errorTrials'])    
      
       
    end

    
    
    
    
end

% if strcmp(session_date, '')
%     if strcmp(speechCue,'Speech')
%         wrongTrials = [];
%     elseif strcmp(speechCue,'Written')
%         wrongTrials = [];       
%     elseif strcmp(speechCue,'All')
%         wrongTrials = [];
%     end     
% end 




