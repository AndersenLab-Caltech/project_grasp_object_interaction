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

        else
            error([ imagename 'Unknown grasp, add it to list']);
        end
        class_all{n_words} = class; 

    end 

end 


if isnumeric(class_all{1})
    class_all = cell2mat(class_all);
end




end

