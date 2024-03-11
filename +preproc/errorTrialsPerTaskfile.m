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

    elseif strcmp(taskfile, '20230907-164131-164600-GraspObject')
        wrongTrials = [18];
    elseif strcmp(taskfile, '20230907-164131-164942-GraspObject')
        wrongTrials = [8,14,15,17,20];
    elseif strcmp(taskfile, '20230907-165333-165440-GraspObject')
        wrongTrials = [10,17,18];
    elseif strcmp(taskfile, '20230907-165333-165822-GraspObject')
        wrongTrials = [4,8,13];
    elseif strcmp(taskfile, '20230907-170205-170303-GraspObject')
        wrongTrials = [17,18];
    elseif strcmp(taskfile, '20230907-170205-170643-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20230921-163332-163522-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230921-163332-163909-GraspObject')
        wrongTrials = [4,5];
    elseif strcmp(taskfile, '20230921-164256-164356-GraspObject')
        wrongTrials = [10,13,14];
    elseif strcmp(taskfile, '20230921-164256-164739-GraspObject')
        wrongTrials = [14,18];
    elseif strcmp(taskfile, '20230921-165134-165232-GraspObject')
        wrongTrials = [9,20];
    elseif strcmp(taskfile, '20230921-165134-165621-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20230922-153857-154413-GraspObject')
        wrongTrials = [14,18,19,20];
    elseif strcmp(taskfile, '20230922-153857-154807-GraspObject')
        wrongTrials = [17,18,19,20];
    elseif strcmp(taskfile, '20230922-155144-155235-GraspObject')
        wrongTrials = [3];
    elseif strcmp(taskfile, '20230922-155144-155614-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230922-155952-160047-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20230922-155952-160420-GraspObject')
        wrongTrials = [13];

    elseif strcmp(taskfile, '20231101-154538-154724-GraspObject')
        wrongTrials = [13,16,17,18,19];
    elseif strcmp(taskfile, '20231101-155157-155307-GraspObject')
        wrongTrials = [8,9,18];
    elseif strcmp(taskfile, '20231101-155657-155720-GraspObject')
        wrongTrials = [19,20];
    elseif strcmp(taskfile, '20231101-160116-160225-GraspObject')
        wrongTrials = [9,18,20];
    elseif strcmp(taskfile, '20231101-160633-160640-GraspObject')
        wrongTrials = [18,19];
    elseif strcmp(taskfile, '20231101-161052-161104-GraspObject')
        wrongTrials = [5,16,17];

    elseif strcmp(taskfile, '20231103-153057-153202-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20231103-153601-153646-GraspObject')
        wrongTrials = [3,18];
    elseif strcmp(taskfile, '20231103-154819-154826-GraspObject')
        wrongTrials = [9,20];
    elseif strcmp(taskfile, '20231103-160743-160800-GraspObject')
        wrongTrials = [7,18];
    elseif strcmp(taskfile, '20231103-161216-161236-GraspObject')
        wrongTrials = [12,13,14];
    elseif strcmp(taskfile, '20231103-161713-161733-GraspObject')
        wrongTrials = [7];

    elseif strcmp(taskfile, '20231201-143639-143900-GraspObject')
        wrongTrials = [17,18,19,20];
    elseif strcmp(taskfile, '20231201-144327-144421-GraspObject')
        wrongTrials = [6,13,14,19];
    elseif strcmp(taskfile, '20231201-144810-144912-GraspObject')
        wrongTrials = [3,7,12,13,17,18,19];
    elseif strcmp(taskfile, '20231201-145309-145324-GraspObject')
        wrongTrials = [13,19,20];
    elseif strcmp(taskfile, '20231201-145714-145733-GraspObject')
        wrongTrials = [7,13];
    elseif strcmp(taskfile, '20231201-150129-150147-GraspObject')
        wrongTrials = [2,4,5,8,9,10,15,16];

    % SHUFFLED IMAGES SESSIONS
    elseif strcmp(taskfile, '20240130-160210-160331-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20240130-160826-160903-GraspObject')
        wrongTrials = [18];
    elseif strcmp(taskfile, '20240130-161349-161358-GraspObject')
        wrongTrials = [7,8,14,18,19,20,23];
    elseif strcmp(taskfile, '20240130-161846-161902-GraspObject')
        wrongTrials = [15,23];
    elseif strcmp(taskfile, '20240130-162325-162349-GraspObject')
        wrongTrials = [12,18,19,22,24];
    elseif strcmp(taskfile, '20240130-162816-162829-GraspObject')
        wrongTrials = [20];

    elseif strcmp(taskfile, '20240216-154209-154222-GraspObject')
        wrongTrials = [11,15,18];
    elseif strcmp(taskfile, '20240216-154659-154726-GraspObject')
        wrongTrials = [7,18,24];
    elseif strcmp(taskfile, '20240216-155234-155302-GraspObject')
        wrongTrials = [12];
    elseif strcmp(taskfile, '20240216-155713-155731-GraspObject')
        wrongTrials = [9,10,16];
    elseif strcmp(taskfile, '20240216-160152-160210-GraspObject')
        wrongTrials = [8,10,18];
    elseif strcmp(taskfile, '20240216-160643-160659-GraspObject')
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

    elseif strcmp(taskfile,'20230830-101127-101544-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230830-101127-102232-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230830-102835-103038-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230830-102835-103700-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230830-104300-104459-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230830-104300-105138-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20230921-111034-111552-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230921-111034-112448-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230921-113100-113731-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230921-113100-114358-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230921-115000-115126-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230921-115000-115741-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20230929-114653-114831-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230929-115538-115603-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230929-120223-120502-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230929-121108-121111-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230929-121731-121812-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20230929-122435-122441-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20231005-111022-111131-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231005-111809-112005-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231005-112620-112638-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231005-113414-113418-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231005-114035-114040-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231005-114924-114931-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20231030-110135-111056-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231030-111725-111811-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231030-112509-112542-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231030-113219-113239-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231030-113942-114002-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231030-114723-114736-GraspObject')
        wrongTrials = [];

    % SHUFFLED IMAGES SESSIONS
    elseif strcmp(taskfile,'20231207-112847-112859-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231207-112847-113437-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231207-113902-113941-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231207-114506-114514-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20231212-105541-105711-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231212-105541-110209-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231212-110837-110847-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231212-111920-111928-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20231212-112429-112441-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20231212-112429-113231-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240112-100249-100306-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240112-100249-100838-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240112-101247-101256-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240112-101247-101759-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240112-102209-102216-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240112-102209-102704-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240119-095536-095924-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240119-100435-100447-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240119-100919-100929-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240119-101450-101518-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240119-101939-101947-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240119-102416-102431-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240208-113355-113411-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240208-113830-113856-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240208-114421-114511-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240208-114959-115017-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240208-115504-120100-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240208-120553-120611-GraspObject')
        wrongTrials = [14,15];

    elseif strcmp(taskfile,'20240214-100633-100800-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240214-101247-101257-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240214-101723-101733-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240214-102155-102351-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240214-102836-102845-GraspObject')
        wrongTrials = [8]; 
    elseif strcmp(taskfile,'20240214-103305-103318-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240229-101108-101128-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240229-101616-101639-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240229-102049-102103-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240229-102518-102539-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240229-102946-103003-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240229-103417-103438-GraspObject')
        wrongTrials = [14,15];
        % elseif strcmp(taskfile,'')
    %    wrongTrials = [];   
        
     else
        keyboard
        error([taskfile ' Unknown session day - add to errorTrials'])    
      
       
    end

    
    
    
    
end






