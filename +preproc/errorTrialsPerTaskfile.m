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

    % SHUFFLED IMAGES SESSIONS %

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

    elseif strcmp(taskfile, '20240625-154025-154614-GraspObject')
        wrongTrials = [21];
    elseif strcmp(taskfile, '20240625-155110-155212-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20240625-155624-155730-GraspObject')
        wrongTrials = [22];
    elseif strcmp(taskfile, '20240625-160154-160220-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20240625-160634-160656-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20240625-161112-161122-GraspObject')
        wrongTrials = [19];

    elseif strcmp(taskfile, '20240703-153433-153443-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20240703-153912-154011-GraspObject')
        wrongTrials = [4,5,6,12,16,19];
    elseif strcmp(taskfile, '20240703-154449-154627-GraspObject')
        wrongTrials = [10,11,12,18,22,23,24];
    elseif strcmp(taskfile, '20240703-155312-155331-GraspObject')
        wrongTrials = [6,13,20,23];
    elseif strcmp(taskfile, '20240703-155815-155824-GraspObject')
        wrongTrials = [7,11,14,18,19,22,23,24];
    elseif strcmp(taskfile, '20240703-160309-160317-GraspObject')
        wrongTrials = [6,7,10,11,15,18,19,22,23,24];

    elseif strcmp(taskfile, '20240710-153111-153122-GraspObject')
        wrongTrials = [12,22,23];
    elseif strcmp(taskfile, '20240710-153542-153634-GraspObject')
        wrongTrials = [8,20];
    elseif strcmp(taskfile, '20240710-154047-154056-GraspObject')
        wrongTrials = [22];
    elseif strcmp(taskfile, '20240710-154519-154612-GraspObject')
        wrongTrials = [6,7];
    elseif strcmp(taskfile, '20240710-155039-155053-GraspObject')
        wrongTrials = [15,24];
    elseif strcmp(taskfile, '20240710-155511-155521-GraspObject')
        wrongTrials = [2,6,10,17,22];

    elseif strcmp(taskfile, '20240716-153015-153134-GraspObject')
        wrongTrials = [19,23];
    elseif strcmp(taskfile, '20240716-153602-153709-GraspObject')
        wrongTrials = [17,18,23];
    elseif strcmp(taskfile, '20240716-154138-154147-GraspObject')
        wrongTrials = [11,17,21,23,24];
    elseif strcmp(taskfile, '20240716-154613-154714-GraspObject')
        wrongTrials = [16,22];
    elseif strcmp(taskfile, '20240716-155131-155139-GraspObject')
        wrongTrials = [5,9,14];
    elseif strcmp(taskfile, '20240716-155558-155608-GraspObject')
        wrongTrials = [10,24];

    elseif strcmp(taskfile, '20241023-155048-155133-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241023-155636-155716-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241023-160132-160142-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241023-160613-160640-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241023-161055-161110-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241023-161530-161539-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20241028-152905-153037-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241028-153455-153533-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241028-153952-154001-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241028-154437-154449-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241028-154918-154928-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241028-155346-155402-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20241111-151559-151608-GraspObject')
        wrongTrials = [1,2,3];
    elseif strcmp(taskfile, '20241111-152048-152128-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241111-152549-152559-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241111-153034-153046-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241111-153511-153522-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20241111-153946-153959-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20250120-161115-161137-GraspObject')
        wrongTrials = [23];
    elseif strcmp(taskfile, '20250120-161616-161630-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250120-162044-162057-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250120-162513-162522-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250120-162940-162949-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250120-163402-163411-GraspObject')
        wrongTrials = [];

    % VARIED SIZES SESSIONS %

    elseif strcmp(taskfile, '20250305-155049-155103-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250305-160306-160316-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250305-160833-160851-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250305-161353-161529-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250305-162100-162109-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250305-162601-162609-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20250312-153338-153355-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250312-154008-154211-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250312-154731-154835-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250312-155330-155344-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250312-155824-155845-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250312-160338-160349-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20250328-154150-154157-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250328-154654-154807-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250328-155333-155419-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250328-155909-155930-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250328-160455-160505-GraspObject')
        wrongTrials = [17];
    elseif strcmp(taskfile, '20250328-161032-161048-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile, '20250402-153558-153613-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250402-154115-154153-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250402-154641-154710-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250402-155158-155215-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250402-155657-155716-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250402-160157-160210-GraspObject')
        wrongTrials = [15,19,22,27];

    elseif strcmp(taskfile, '20250501-112005-112012-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250501-112506-112538-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250501-113027-113038-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250501-113535-113601-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250501-114044-114055-GraspObject')
        wrongTrials = [28];
    elseif strcmp(taskfile, '20250501-114545-114606-GraspObject')
        wrongTrials = [5,12,19,25];

    elseif strcmp(taskfile, '20250517-163354-163411-GraspObject')
        wrongTrials = [7,11,23];
    elseif strcmp(taskfile, '20250517-163925-163953-GraspObject')
        wrongTrials = [10,13,16,20];
    elseif strcmp(taskfile, '20250517-164448-164500-GraspObject')
        wrongTrials = [1,3,9,10,11,12,22,26];
    elseif strcmp(taskfile, '20250517-165007-165022-GraspObject')
        wrongTrials = [ 4,5,11,14,18,19,25,26];
    elseif strcmp(taskfile, '20250517-165547-165602-GraspObject')
        wrongTrials = [6,11,14,16,18,22,26];
    elseif strcmp(taskfile, '20250517-170056-170109-GraspObject')
        wrongTrials = [2,7,9,10,11,18,22,23,24];

    elseif strcmp(taskfile, '20250519-105808-105820-GraspObject')
        wrongTrials = [9,10,11,12,22,23,28];
    elseif strcmp(taskfile, '20250519-110425-110517-GraspObject')
        wrongTrials = [15,23,24,25,26];
    elseif strcmp(taskfile, '20250519-111003-111028-GraspObject')
        wrongTrials = [3,4,5,9,12,13,14,19,22,23,24,25];
    elseif strcmp(taskfile, '20250519-111559-111615-GraspObject')
        wrongTrials = [9,10,15,21,22,23,24];
    elseif strcmp(taskfile, '20250519-112127-112139-GraspObject')
        wrongTrials = [1,2,3,8,13,15,18,19,20,26,27,28];
    elseif strcmp(taskfile, '20250519-112643-112657-GraspObject')
        wrongTrials = [1,2,3,9,10,12,17,21,25,26];

    elseif strcmp(taskfile, '20250520-105649-105843-GraspObject')
        wrongTrials = [23];
    elseif strcmp(taskfile, '20250520-110403-110434-GraspObject')
        wrongTrials = [18,23,27];
    elseif strcmp(taskfile, '20250520-110921-111000-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile, '20250520-111443-111502-GraspObject')
        wrongTrials = [6,12,15,16,20,22,27,28];
    elseif strcmp(taskfile, '20250520-112002-112010-GraspObject')
        wrongTrials = [16,23];
    elseif strcmp(taskfile, '20250520-112456-112509-GraspObject')
        wrongTrials = [9,14,21];



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

    % SHUFFLED IMAGES SESSIONS %

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

    elseif strcmp(taskfile,'20240829-155905-155913-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240829-160425-160441-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240829-160900-160911-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240829-161422-161451-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240829-161913-161925-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240829-162353-162404-GraspObject')
        wrongTrials = [];

    % VARIED_SIZE SESSIONS %

    elseif strcmp(taskfile,'20240521-110839-110853-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240521-110839-111606-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240521-112503-112651-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240521-113734-114016-GraspObject')
        wrongTrials = [18,26,27,28,32];
    elseif strcmp(taskfile,'20240521-114722-114746-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240521-115446-115456-GraspObject')
        wrongTrials = [14,15];

    elseif strcmp(taskfile,'20240523-094246-094429-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240523-095408-095418-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240523-100124-100136-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240523-100734-100748-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240523-101809-101819-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240523-102422-102436-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240621-094059-094342-GraspObject')
        wrongTrials = [6,31]; 
    elseif strcmp(taskfile,'20240621-095207-095218-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240621-100010-100029-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240621-100705-100718-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240621-101327-101338-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240621-101929-101940-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240718-140204-140212-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-141437-141447-GraspObject')
        wrongTrials = [9,10,11,16,17,18]; 
    elseif strcmp(taskfile,'20240718-142439-142447-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-143149-143201-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240718-144542-144553-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-145448-145459-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240719-094048-094352-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-095355-095419-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-100008-100029-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-100702-100725-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240719-101325-101336-GraspObject')
        wrongTrials = [8]; 
    elseif strcmp(taskfile,'20240719-101940-101955-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240719-143531-143639-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-144329-144345-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-145009-145032-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-145704-145721-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240719-150313-150327-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240719-150927-150941-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240808-152825-152911-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240808-153657-153720-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240808-154557-154616-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240808-155602-155611-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240808-160314-160324-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240808-160928-160945-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240809-105639-110749-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240809-111340-111401-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240809-111951-112013-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240809-112602-112621-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240809-113218-113234-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240809-113858-113911-GraspObject')
        wrongTrials = [];

    % 50/50 GO-NOGO SESSIONS %

    elseif strcmp(taskfile,'20241016-113418-113429-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241016-114054-114103-GraspObject')
        wrongTrials = [23,24]; 
    elseif strcmp(taskfile,'20241016-114626-114646-GraspObject')
        wrongTrials = [11,17,21,24,28]; 
    elseif strcmp(taskfile,'20241016-115221-115302-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20241016-115815-115826-GraspObject')
        wrongTrials = [7,14]; 
    elseif strcmp(taskfile,'20241016-120359-120413-GraspObject')
        wrongTrials = [3];

    elseif strcmp(taskfile,'20241024-101048-101059-GraspObject')
        wrongTrials = [29]; 
    elseif strcmp(taskfile,'20241024-101651-101730-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241024-102314-102334-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241024-102934-103007-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20241024-103555-103617-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241024-104200-104215-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20241120-114331-114339-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241120-114915-114953-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241120-115509-115537-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241120-120101-120111-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20241120-120627-120643-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241120-121155-121206-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20241121-105635-105651-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241121-110206-110231-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241121-110803-110815-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241121-111333-111353-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20241121-111913-111928-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20241121-112448-112459-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250123-095705-095744-GraspObject')
        wrongTrials = [15]; 
    elseif strcmp(taskfile,'20250123-100334-100403-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250123-100919-100933-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250123-101616-101630-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250123-102309-102319-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250123-102836-102923-GraspObject')
        wrongTrials = [];

    % COMBINATIONS SESSIONS %

    elseif strcmp(taskfile,'20250211-125620-125632-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250211-130241-130256-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250211-130853-130910-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250211-132727-132737-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250211-132727-133707-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250211-132727-134502-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250212-113556-113605-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250212-114206-114241-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250212-114831-114839-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250212-115755-115804-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250212-115755-120459-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250212-115755-121127-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250408-105620-105629-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-110421-110431-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-111021-111109-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-111847-111911-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250408-112726-112736-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-113431-113444-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250408-144506-144518-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-145216-145232-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-145823-145832-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250408-150422-150501-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250408-151139-151150-GraspObject')
        wrongTrials = [10]; 
    elseif strcmp(taskfile,'20250408-151812-151826-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250409-115658-115707-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-120724-120749-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-121336-121353-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-122025-122129-GraspObject')
        wrongTrials = [16];
    elseif strcmp(taskfile,'20250409-122825-122835-GraspObject')
        wrongTrials = [35]; 
    elseif strcmp(taskfile,'20250409-123515-123522-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20250409-150257-150435-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-151031-151105-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-151736-151802-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-152358-152427-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250409-153053-153103-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250409-153727-153746-GraspObject')
        wrongTrials = [8,39];

    elseif strcmp(taskfile,'20250424-094649-094659-GraspObject')
        wrongTrials = [22]; 
    elseif strcmp(taskfile,'20250424-095333-095354-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250424-100005-100025-GraspObject')
        wrongTrials = [6]; 
    elseif strcmp(taskfile,'20250424-100745-100757-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20250424-101453-101548-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20250424-102211-102226-GraspObject')
        wrongTrials = [];

        % elseif strcmp(taskfile,'')
    %    wrongTrials = [];   
        
    else
        
        disp(task.task.params.user.cue_type)
        disp(taskfile)
        keyboard
        error([' Unknown session day - add to errorTrials'])    
      
     end

     
elseif strcmp(subject, 's4')
    
     if strcmp(taskfile,   '20240307-133830-133839-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240307-133830-134703-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240307-135309-135446-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240307-135309-140133-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240307-140857-140906-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240307-140857-141641-GraspObject')
        wrongTrials = []; 

    elseif strcmp(taskfile,'20240415-162730-162737-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240415-162730-163715-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240415-164610-164631-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240415-164610-165405-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240415-170112-170122-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240415-170112-170756-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240418-122926-123043-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240418-124055-124106-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240418-124737-124820-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240418-130618-130628-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240418-131311-131330-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240418-132028-132039-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240429-144113-144417-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240429-145020-145117-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240429-145759-145841-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240429-150505-150520-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240429-151129-151339-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240429-151935-151948-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240509-145935-150256-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240509-151005-151059-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240509-152331-152436-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240509-153249-153259-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240509-153939-154059-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240509-154650-154704-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240523-151547-151646-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240523-152629-152731-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240523-153502-153513-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240523-154652-154708-GraspObject')
        wrongTrials = [27];  
    elseif strcmp(taskfile,'20240523-155544-155557-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240523-160147-160200-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240610-143532-143846-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240610-144505-144542-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240610-145213-145226-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240610-145819-145924-GraspObject')
        wrongTrials = [27];  
    elseif strcmp(taskfile,'20240610-150627-150641-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240610-151237-151247-GraspObject')
        wrongTrials = [];

    % REAL IMAGES    

    elseif strcmp(taskfile,'20240613-143221-143547-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240613-143221-145119-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240613-150115-150245-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240613-150115-150959-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240613-151818-151827-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240613-151818-153040-GraspObject')
        wrongTrials = [9];

    elseif strcmp(taskfile,'20240617-160042-160057-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240617-161105-161116-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240617-162104-162117-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240617-163622-163637-GraspObject')
        wrongTrials = [];  
    elseif strcmp(taskfile,'20240617-164325-164355-GraspObject')
        wrongTrials = [];   
    elseif strcmp(taskfile,'20240617-165646-165659-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240620-155701-155736-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240620-160401-160449-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240620-161055-161126-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240620-161801-161824-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240620-162434-162449-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240620-163058-163112-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240627-142509-142633-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240627-143258-143334-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240627-144024-144055-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240627-144713-145248-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240627-145849-145902-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240627-150511-150525-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240715-143149-143354-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240715-144056-144138-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240715-144745-144804-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240715-145404-145427-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240715-150039-150053-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240715-150700-150714-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240718-150705-150739-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-151417-151438-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-152047-152115-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-152737-152831-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240718-153444-153501-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240718-154116-154129-GraspObject')
        wrongTrials = [];

    elseif strcmp(taskfile,'20240822-155251-155611-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240822-160345-160401-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240822-161051-161106-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240822-161908-161946-GraspObject')
        wrongTrials = [];
    elseif strcmp(taskfile,'20240822-162623-162637-GraspObject')
        wrongTrials = []; 
    elseif strcmp(taskfile,'20240822-163238-163253-GraspObject')
        wrongTrials = [];
        % elseif strcmp(taskfile,'')
    %    wrongTrials = []; 

     else
        keyboard
        error([taskfile ' Unknown session day - add to errorTrials'])    
      
       
    end

    
    
    
    
end






