function [colors] = get_color_rgb_codes(imagename)
colors = cell(1, length(imagename));

%Input the labels of what you want to plot, the correct color will be
%output. E.g. labels of different graps, or different phases. 
for i =1:length(imagename)
    
     if strcmp('Lateral', imagename{i}) ==1 
        colors{i} = utile.rgb('Red');
    elseif strcmp('WritingTripod',  imagename{i}) ==1 
        colors{i} = utile.rgb('LawnGreen');
    elseif strcmp('MediumWrap',  imagename{i}) ==1 
        colors{i} = utile.rgb('Cyan');
    elseif strcmp('PalmarPinch',  imagename{i}) ==1 
        colors{i} = utile.rgb('Gold');
    elseif strcmp('Sphere3Finger',  imagename{i}) ==1 
        colors{i} = utile.rgb('Magenta');
    elseif strcmp('redballoon',  imagename{i}) ==1 
        colors{i} = utile.rgb('Silver');    
        
    elseif strcmp('Lateral No Go',  imagename{i}) ==1 
        colors{i} = utile.rgb('Tomato');
    elseif strcmp('WritingTripod No Go',  imagename{i}) ==1 
        colors{i} = utile.rgb('LightGreen');
    elseif strcmp('MediumWrap No Go',  imagename{i}) ==1 
        colors{i} = utile.rgb('LightCyan');
    elseif strcmp('PalmarPinch No Go',  imagename{i}) ==1 
        colors{i} = utile.rgb('Yellow');
    elseif strcmp('Sphere3Finger No Go',  imagename{i}) ==1 
        colors{i} = utile.rgb('Violet');
 
     elseif strcmp('Lateral Image', imagename{i}) ==1 
        colors{i} = utile.rgb('Green');
    elseif strcmp('WritingTripod Image',  imagename{i}) ==1 
        colors{i} = utile.rgb('LawnGreen');
    elseif strcmp('MediumWrap Image',  imagename{i}) ==1 
        colors{i} = utile.rgb('MediumAquamarine');
    elseif strcmp('PalmarPinch Image',  imagename{i}) ==1 
        colors{i} = utile.rgb('SeaGreen');
    elseif strcmp('Sphere3Finger Image',  imagename{i}) ==1 
        colors{i} = utile.rgb('ForestGreen'); 
        
    elseif strcmp('Lateral Auditory',  imagename{i}) ==1 
        colors{i} = utile.rgb('Blue');
    elseif strcmp('WritingTripod Auditory',  imagename{i}) ==1 
        colors{i} = utile.rgb('DarkCyan');
    elseif strcmp('MediumWrap Auditory',  imagename{i}) ==1 
        colors{i} = utile.rgb('LightCyan');
    elseif strcmp('PalmarPinch Auditory',  imagename{i}) ==1 
        colors{i} = utile.rgb('RoyalBlue');
    elseif strcmp('Sphere3Finger Auditory',  imagename{i}) ==1 
        colors{i} = utile.rgb('Navy');  
        
        
     elseif strcmp('Lateral Written',  imagename{i}) ==1 
        colors{i} = utile.rgb('Yellow');
    elseif strcmp('WritingTripod Written',  imagename{i}) ==1 
        colors{i} = utile.rgb('LemonChiffon');
    elseif strcmp('MediumWrap Written',  imagename{i}) ==1 
        colors{i} = utile.rgb('Moccasin');
    elseif strcmp('PalmarPinch Written',  imagename{i}) ==1 
        colors{i} = utile.rgb('Khaki');
    elseif strcmp('Sphere3Finger Written',  imagename{i}) ==1 
        colors{i} = utile.rgb('Gold');  
        
        
        
    elseif  strcmp('Action',  imagename{i}) ==1 
        colors{i} = utile.rgb('SteelBlue');       
    elseif strcmp('Cue',  imagename{i}) ==1 
        colors{i} = utile.rgb('Orange');  
    elseif strcmp('Image Cue',  imagename{i}) ==1 
        colors{i} = utile.rgb('Orange');  
    elseif strcmp('CuePhase',  imagename{i}) ==1 
        colors{i} = utile.rgb('Orange');  
    elseif strcmp('Both',  imagename{i}) ==1 
        colors{i} = utile.rgb('LightSlateGray');       
    elseif strcmp('PMV',  imagename{i}) ==1 
        colors{i} = utile.rgb('Red');
    elseif strcmp('SMG',  imagename{i}) ==1 
        colors{i} = utile.rgb('Blue');        
    elseif strcmp('SMG&PMV',  imagename{i}) ==1 
        colors{i} = utile.rgb('Purple');      
    elseif strcmp('V',  imagename{i}) ==1 
        colors{i} = utile.rgb('LightSlateGray');       
    elseif strcmp('S',  imagename{i}) 
        colors{i} = utile.rgb('LightSteelBlue');
           

    elseif  strcmp('ITI',  imagename{i}) ==1 
        colors{i} = utile.rgb('Silver');    
        
    elseif  strcmp('Delay',  imagename{i}) ==1 
        colors{i} = utile.rgb('Cyan'); 
        
     elseif  strcmp('Tuned',  imagename{i}) ==1 
        colors{i} = utile.rgb('Green'); 
        
     elseif  strcmp('Non-tuned',  imagename{i}) ==1 
        colors{i} = utile.rgb('DarkGrey'); 
        
     elseif  strcmp('Removed',  imagename{i}) ==1 
        colors{i} = utile.rgb('DarkRed'); 
        
     elseif  strcmp('Group Mean',  imagename{i}) ==1 
        colors{i} = utile.rgb('Black'); 
        
     elseif  strcmp('MotorImageryLeft',  imagename{i}) ==1 
        colors{i} = utile.rgb('Fire'); 
        
    elseif  strcmp('MotorImageryRight',  imagename{i}) ==1 
        colors{i} = utile.rgb('Green'); 
        
     elseif  strcmp('MotorImagery',  imagename{i}) ==1 
        colors{i} = utile.rgb('Green'); 
    
%      elseif  strcmp('Image',  imagename{i}) ==1 
%         colors{i} = utile.rgb('FireBrick'); 
%         
 
        
    elseif  strcmp('Speaking',  imagename{i}) 
        colors{i} = utile.rgb('Blue'); 
        
     elseif  strcmp('Audio',  imagename{i}) 
        colors{i} = utile.rgb('SeaGreen'); 
        
     elseif  strcmp('Go Trial',  imagename{i}) 
        colors{i} = utile.rgb('Green'); 
        
     elseif  strcmp('No Go Trial',  imagename{i}) 
        colors{i} = utile.rgb('Red'); 
        
     elseif  strcmp('Blue',  imagename{i}) 
        colors{i} = utile.rgb('Blue'); 
        
     elseif  strcmp('Green',  imagename{i}) 
        colors{i} = utile.rgb('Green'); 
        
     elseif  strcmp('Yellow',  imagename{i}) 
        colors{i} = utile.rgb('Fire'); 
        
    
     elseif  strcmp('Brown',  imagename{i}) 
        colors{i} = utile.rgb('Brown'); 
        
    
     elseif  strcmp('Gray',  imagename{i}) 
        colors{i} = utile.rgb('SlateGray'); 
        
    
     elseif  strcmp('Grasps',  imagename{i}) 
        colors{i} = utile.rgb('Blue'); 
        
     elseif  strcmp('Colors',  imagename{i}) 
        colors{i} = utile.rgb('Red'); 
        
     elseif  strcmp('Colors',  imagename{i}) 
        colors{i} = utile.rgb('Red'); 
    
     elseif strcmp('ImageCue', imagename{i})
         colors{i} = utile.rgb('Green');
         
     elseif strcmp('AuditoryCue', imagename{i})
         colors{i} = utile.rgb('Blue');
         
      elseif strcmp('WrittenCue', imagename{i})
         colors{i} = utile.rgb('Fire');
         
      elseif strcmp('Image', imagename{i})
         colors{i} = utile.rgb('Green');
         
     elseif strcmp('Auditory', imagename{i})
         colors{i} = utile.rgb('Blue');
         
      elseif strcmp('Written', imagename{i})
         colors{i} = utile.rgb('Fire');
         
      elseif strcmp('BothCues', imagename{i})
         colors{i} = utile.rgb('Cyan');
         
     elseif strcmp('Training', imagename{i})
     colors{i} = utile.rgb('Green');
     
     elseif strcmp('Testing', imagename{i})
     colors{i} = utile.rgb('Blue');
     
    elseif strcmp('Visuomotor', imagename{i}) 
        colors{i} = utile.rgb('LightSlateGray');       
        
     elseif strcmp('Action&Cue', imagename{i}) 
        colors{i} = utile.rgb('LightSlateGray');       
        
    elseif strcmp('Switching', imagename{i})
        colors{i} = utile.rgb('LightSteelBlue');
        
       
     elseif strcmp('Plan', imagename{i})
        colors{i} = utile.rgb('LightSteelBlue');
        
    else  
        error(['Condition ' imagename{i} ' not present, add color'])
        
    end 
end 

end

