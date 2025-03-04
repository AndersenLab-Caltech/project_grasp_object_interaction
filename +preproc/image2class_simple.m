function [ class_all ] = image2class_simple( imagename_all )

%Returns the name of a class if provided a number, or the number of the
%class if provided a name. 

if ischar(imagename_all)
    image_n = 1;
else
    image_n = length(imagename_all); 
end

class_all = {};

for n_words = 1:image_n

    if ~isnumeric(imagename_all) 

        if image_n ~= 1
            imagename = imagename_all{n_words};
        else
            imagename = imagename_all;
        end 

        if strcmp('MediumWrap_Hand', imagename)
             class =1;

        elseif strcmp('PalmarPinch_Hand', imagename)
             class =2;

        elseif strcmp('Lateral_Hand', imagename)
             class =3;

        elseif strcmp('Sphere3Finger_Hand', imagename)
              class =4;

        elseif strcmp('MediumWrap_Hand_Object', imagename)
               class =5;

        elseif strcmp('PalmarPinch_Hand_Object', imagename)
               class =6;

        elseif strcmp('Lateral_Hand_Object', imagename)
               class =7;

        elseif strcmp('Sphere3Finger_Hand_Object', imagename)
                class =8;

        elseif strcmp('MediumWrap_Object', imagename)
                class =9;

        elseif strcmp('PalmarPinch_Object', imagename)
                class =10;

        elseif strcmp('Lateral_Object', imagename)
                class =11;

        elseif strcmp('Sphere3Finger_Object', imagename)
                class =12;

        elseif strcmp('MediumWrap_Hand_Shuffled', imagename)
                class =13;

        elseif strcmp('PalmarPinch_Hand_Shuffled', imagename)
                class =14;

        elseif strcmp('Lateral_Hand_Shuffled', imagename)
                class =15;

        elseif strcmp('Sphere3Finger_Hand_Shuffled', imagename)
                class =16;

        elseif strcmp('MediumWrap_Hand_Object_Shuffled', imagename)
               class =17;

        elseif strcmp('PalmarPinch_Hand_Object_Shuffled', imagename)
               class =18;

        elseif strcmp('Lateral_Hand_Object_Shuffled', imagename)
               class =19;

        elseif strcmp('Sphere3Finger_Hand_Object_Shuffled', imagename)
                class =20;

        elseif strcmp('MediumWrap_Object_Shuffled', imagename)
                class =21;

        elseif strcmp('PalmarPinch_Object_Shuffled', imagename)
                class =22;

        elseif strcmp('Lateral_Object_Shuffled', imagename)
                class =23;

        elseif strcmp('Sphere3Finger_Object_Shuffled', imagename)
                class =24;

        elseif strcmp('Sphere3Finger_Object_Large', imagename)
                class =25;
        
        elseif strcmp('Sphere3Finger_Object_Medium', imagename)
                class =26;
        
        elseif strcmp('Sphere3Finger_Object_Small', imagename)
                class =27;
        
        elseif strcmp('Sphere3Finger_Hand_Object_Large', imagename)
                class =28;
        
        elseif strcmp('Sphere3Finger_Hand_Object_Medium', imagename)
                class =29;
        
        elseif strcmp('Sphere3Finger_Hand_Object_Small', imagename)
                class =30;

        elseif strcmp('Sphere3Finger_Hand_Large', imagename)
                class =31;
        
        elseif strcmp('Sphere3Finger_Hand_Medium', imagename)
                class =32;
        
        elseif strcmp('Sphere3Finger_Hand_Small', imagename)
                class =33;

        elseif strcmp('PalmarPinch_Object_Large', imagename)
                class =34;
        
        elseif strcmp('PalmarPinch_Object_Medium', imagename)
                class =35;
        
        elseif strcmp('PalmarPinch_Object_Small', imagename)
                class =36;
        
        elseif strcmp('PalmarPinch_Hand_Object_Large', imagename)
                class =37;
        
        elseif strcmp('PalmarPinch_Hand_Object_Medium', imagename)
                class =38;
        
        elseif strcmp('PalmarPinch_Hand_Object_Small', imagename)
                class =39;

        elseif strcmp('PalmarPinch_Hand_Large', imagename)
                class =40;
        
        elseif strcmp('PalmarPinch_Hand_Medium', imagename)
                class =41;
        
        elseif strcmp('PalmarPinch_Hand_Small', imagename)
                class =42;

        elseif strcmp('MediumWrap_Object_Large', imagename)
                class =43;
        
        elseif strcmp('MediumWrap_Object_Medium', imagename)
                class =44;
        
        elseif strcmp('MediumWrap_Object_Small', imagename)
                class =45;
        
        elseif strcmp('MediumWrap_Hand_Object_Large', imagename)
                class =46;
        
        elseif strcmp('MediumWrap_Hand_Object_Medium', imagename)
                class =47;
        
        elseif strcmp('MediumWrap_Hand_Object_Small', imagename)
                class =48;

        elseif strcmp('MediumWrap_Hand_Large', imagename)
                class =49;
        
        elseif strcmp('MediumWrap_Hand_Medium', imagename)
                class =50;
        
        elseif strcmp('MediumWrap_Hand_Small', imagename)
                class =51;

        elseif strcmp('Lateral_Object_Large', imagename)
                class =52;
        
        elseif strcmp('Lateral_Object_Medium', imagename)
                class =53;
        
        elseif strcmp('Lateral_Object_Small', imagename)
                class =54;
        
        elseif strcmp('Lateral_Hand_Object_Large', imagename)
                class =55;
        
        elseif strcmp('Lateral_Hand_Object_Medium', imagename)
                class =56;
        
        elseif strcmp('Lateral_Hand_Object_Small', imagename)
                class =57;

        elseif strcmp('Lateral_Hand_Large', imagename)
                class =58;
        
        elseif strcmp('Lateral_Hand_Medium', imagename)
                class =59;
        
        elseif strcmp('Lateral_Hand_Small', imagename)
                class =60;

        elseif strcmp('PalmarPinch_Hnad_Object', imagename) % misspelled 'Hand' and unsure how to correct it at this stage
                class = 61;

        elseif strcmp('Sphere3Finger_Combination_deck', imagename)
                class =62;

        elseif strcmp('Sphere3Finger_Combination_block', imagename)
                class =63;

        elseif strcmp('Sphere3Finger_Combination_rod', imagename)
                class =64;

        elseif strcmp('Sphere3Finger_Combination_ball', imagename)
                class =65;

        elseif strcmp('Lateral_Combination_deck', imagename)
                class =66;

        elseif strcmp('Lateral_Combination_block', imagename)
                class =67;

        elseif strcmp('Lateral_Combination_rod', imagename)
                class =68;

        elseif strcmp('Lateral_Combination_ball', imagename)
                class =69;

        elseif strcmp('MediumWrap_Combination_deck', imagename)
                class =70;

        elseif strcmp('MediumWrap_Combination_block', imagename)
                class =71;

        elseif strcmp('MediumWrap_Combination_rod', imagename)
                class =72;
               
        elseif strcmp('MediumWrap_Combination_ball', imagename)
                class =73;

        elseif strcmp('PalmarPinch_Combination_deck', imagename)
                class =74;

        elseif strcmp('PalmarPinch_Combination_block', imagename)
                class =75;

        elseif strcmp('PalmarPinch_Combination_rod', imagename)
                class =76;

        elseif strcmp('PalmarPinch_Combination_ball', imagename)
                class =77;

        else

            error([ imagename ' - Unknown label, add it to list']);
        end
        class_all{n_words} = class; 
        
        
    elseif(isnumeric(imagename_all))
        imagename = imagename_all(n_words);

        if imagename ==1 
             class = 'MediumWrap_Hand';

        elseif imagename == 2
             class ='PalmarPinch_Hand';

        elseif  imagename == 3
             class = 'Lateral_Hand';

        elseif imagename ==4
              class ='Sphere3Finger_Hand'; 

        elseif imagename== 5
               class = 'MediumWrap_Hand_Object';

        elseif imagename ==6
               class ='PalmarPinch_Hand_Object';

        elseif imagename == 7
               class ='Lateral_Hand_Object';

        elseif imagename == 8
                class ='Sphere3Finger_Hand_Object'; 

        elseif imagename == 9
                class ='MediumWrap_Object'; 
              
        elseif imagename == 10
                class ='PalmarPinch_Object'; 

        elseif imagename == 11
                class ='Lateral_Object'; 

        elseif imagename == 12
                class ='Sphere3Finger_Object';

        elseif imagename ==13 
             class = 'MediumWrap_Hand_Shuffled';

        elseif imagename == 14
             class ='PalmarPinch_Hand_Shuffled';

        elseif  imagename == 15
             class = 'Lateral_Hand_Shuffled';

        elseif imagename ==16
              class ='Sphere3Finger_Hand_Shuffled'; 

        elseif imagename== 17
               class = 'MediumWrap_Hand_Object_Shuffled';

        elseif imagename ==18
               class ='PalmarPinch_Hand_Object_Shuffled';

        elseif imagename == 19
               class ='Lateral_Hand_Object_Shuffled';

        elseif imagename == 20
                class ='Sphere3Finger_Hand_Object_Shuffled'; 

        elseif imagename == 21
                class ='MediumWrap_Object_Shuffled'; 
              
        elseif imagename == 22
                class ='PalmarPinch_Object_Shuffled'; 

        elseif imagename == 23
                class ='Lateral_Object_Shuffled'; 

        elseif imagename == 24
                class ='Sphere3Finger_Object_Shuffled';

        elseif imagename == 25
                class ='Sphere3Finger_Object_Large';

        elseif imagename == 26
                class ='Sphere3Finger_Object_Medium';

        elseif imagename == 27
                class ='Sphere3Finger_Object_Small';

        elseif imagename == 28
                class ='Sphere3Finger_Hand_Object_Large';

        elseif imagename == 29
                class ='Sphere3Finger_Hand_Object_Medium';

        elseif imagename == 30
                class ='Sphere3Finger_Hand_Object_Small';

        elseif imagename == 31
                class ='Sphere3Finger_Hand_Large';

        elseif imagename == 32
                class ='Sphere3Finger_Hand_Medium';

        elseif imagename == 33
                class ='Sphere3Finger_Hand_Small';

        elseif imagename == 34
                class ='PalmarPinch_Object_Large';

        elseif imagename == 35
                class ='PalmarPinch_Object_Medium';

        elseif imagename == 36
                class ='PalmarPinch_Object_Small';

        elseif imagename == 37
                class ='PalmarPinch_Hand_Object_Large';

        elseif imagename == 38
                class ='PalmarPinch_Hand_Object_Medium';

        elseif imagename == 39
                class ='PalmarPinch_Hand_Object_Small';

        elseif imagename == 40
                class ='PalmarPinch_Hand_Large';

        elseif imagename == 41
                class ='PalmarPinch_Hand_Medium';

        elseif imagename == 42
                class ='PalmarPinch_Hand_Small';

        elseif imagename == 43
                class ='MediumWrap_Object_Large';

        elseif imagename == 44
                class ='MediumWrap_Object_Medium';

        elseif imagename == 45
                class ='MediumWrap_Object_Small';

        elseif imagename == 46
                class ='MediumWrap_Hand_Object_Large';

        elseif imagename == 47
                class ='MediumWrap_Hand_Object_Medium';

        elseif imagename == 48
                class ='MediumWrap_Hand_Object_Small';

        elseif imagename == 49
                class ='MediumWrap_Hand_Large';

        elseif imagename == 50
                class ='MediumWrap_Hand_Medium';

        elseif imagename == 51
                class ='MediumWrap_Hand_Small';

        elseif imagename == 52
                class ='Lateral_Object_Large';

        elseif imagename == 53
                class ='Lateral_Object_Medium';

        elseif imagename == 54
                class ='Lateral_Object_Small';

        elseif imagename == 55
                class ='Lateral_Hand_Object_Large';

        elseif imagename == 56
                class ='Lateral_Hand_Object_Medium';

        elseif imagename == 57
                class ='Lateral_Hand_Object_Small';

        elseif imagename == 58
                class ='Lateral_Hand_Large';

        elseif imagename == 59
                class ='Lateral_Hand_Medium';

        elseif imagename == 60
                class ='Lateral_Hand_Small';

        elseif imagename == 61
                class ='PalmarPinch_Hnad_Object'; % misspelled 'Hand' and unsure how to correct it at this stage

        elseif imagename == 62
                class ='Sphere3Finger_Combination_deck';

        elseif imagename == 63
                class ='Sphere3Finger_Combination_block';

        elseif imagename == 64
                class ='Sphere3Finger_Combination_rod';

        elseif imagename == 65
                class ='Sphere3Finger_Combination_ball';

        elseif imagename == 66
                class ='Lateral_Combination_deck';

        elseif imagename == 67
                class ='Lateral_Combination_block';

        elseif imagename == 68
                class ='Lateral_Combination_rod';

        elseif imagename == 69
                class ='Lateral_Combination_ball';

        elseif imagename == 70
                class ='MediumWrap_Combination_deck';

        elseif imagename == 71
                class ='MediumWrap_Combination_block';

        elseif imagename == 72
                class ='MediumWrap_Combination_rod';

        elseif imagename == 73
                class ='MediumWrap_Combination_ball';

        elseif imagename == 74
                class ='PalmarPinch_Combination_deck';

        elseif imagename == 75
                class ='PalmarPinch_Combination_block';

        elseif imagename == 76
                class ='PalmarPinch_Combination_rod';

        elseif imagename == 77
                class ='PalmarPinch_Combination_ball';

        else
            error([ imagename 'Unknown grasp, add it to list pls']);
        end
        class_all{n_words} = class; 

    end 

end 


if isnumeric(class_all{1})
    class_all = cell2mat(class_all);
end




end

