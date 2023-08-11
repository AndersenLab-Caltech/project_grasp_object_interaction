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

