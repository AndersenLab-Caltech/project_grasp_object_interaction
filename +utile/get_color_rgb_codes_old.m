function [colors] = get_color_rgb_codes(imagename)
colors = cell(1, length(imagename));

%Input the labels of what you want to plot, the correct color will be
%output. E.g. labels of different graps, or different phases. 
for i =1:length(imagename)
    
     if mean(ismember('Lateral', imagename{i})) ==1 
        colors{i} = utile.rgb('Red');
    elseif mean(ismember('WritingTripod',  imagename{i})) ==1 
        colors{i} = utile.rgb('LawnGreen');
    elseif mean(ismember('MediumWrap',  imagename{i})) ==1 
        colors{i} = utile.rgb('Cyan');
    elseif mean(ismember('PalmarPinch',  imagename{i})) ==1 
        colors{i} = utile.rgb('Gold');
    elseif mean(ismember('Sphere3Finger',  imagename{i})) ==1 
        colors{i} = utile.rgb('Magenta');
    elseif  mean(ismember('redballoon',  imagename{i})) ==1 
        colors{i} = utile.rgb('Silver');    
        
        
    elseif  mean(ismember('Action',  imagename{i})) ==1 
        colors{i} = utile.rgb('SteelBlue');       
    elseif mean(ismember('Cue',  imagename{i})) ==1 
        colors{i} = utile.rgb('Orange');   
     elseif mean(ismember('Cue',  imagename{i})) ==1 
        colors{i} = utile.rgb('Orange');
    elseif mean(ismember('Both',  imagename{i})) ==1 
        colors{i} = utile.rgb('LightSlateGray');       
     elseif mean(ismember('PMV',  imagename{i})) ==1 
        colors{i} = utile.rgb('Red');
    elseif mean(ismember('SMG',  imagename{i})) ==1 
        colors{i} = utile.rgb('Blue');        
    elseif mean(ismember('SMG&PMV',  imagename{i})) ==1 
        colors{i} = utile.rgb('Purple');      
    elseif mean(ismember('V',  imagename{i})) ==1 
        colors{i} = utile.rgb('LightSlateGray');       
    elseif contains('S',  imagename{i}) 
        colors{i} = utile.rgb('LightSteelBlue');
           

    elseif  mean(ismember('ITI',  imagename{i})) ==1 
        colors{i} = utile.rgb('Silver');    
        
    elseif  mean(ismember('Delay',  imagename{i})) ==1 
        colors{i} = utile.rgb('Cyan'); 
        
     elseif  mean(ismember('Tuned',  imagename{i})) ==1 
        colors{i} = utile.rgb('Green'); 
        
     elseif  mean(ismember('Non-tuned',  imagename{i})) ==1 
        colors{i} = utile.rgb('DarkGrey'); 
        
     elseif  mean(ismember('Removed',  imagename{i})) ==1 
        colors{i} = utile.rgb('DarkRed'); 
        
     elseif  mean(ismember('Group Mean',  imagename{i})) ==1 
        colors{i} = utile.rgb('Black'); 
        
     elseif  mean(ismember('MotorImageryLeft',  imagename{i})) ==1 
        colors{i} = utile.rgb('Fire'); 
        
    elseif  mean(ismember('MotorImageryRight',  imagename{i})) ==1 
        colors{i} = utile.rgb('Green'); 
        
     elseif  mean(ismember('MotorImagery',  imagename{i})) ==1 
        colors{i} = utile.rgb('Green'); 
    
     elseif  mean(ismember('Image',  imagename{i})) ==1 
        colors{i} = utile.rgb('FireBrick'); 
        
 
        
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
         
      elseif strcmp('BothCues', imagename{i})
         colors{i} = utile.rgb('Cyan');
         
    else  
        error('Grasp not present, add color')
        
    end 
end 

end

