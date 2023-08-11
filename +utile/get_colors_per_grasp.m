function [colors] = get_colors_per_grasp(grasp_names)
colors = cell(1, length(grasp_names));

for i =1:length(grasp_names)
    
    if strcmp( 'Lateral' ,grasp_names{i})
        colors{i} = 'r';
    elseif strcmp( 'WritingTripod' ,grasp_names{i})
        colors{i} = 'g';
        %colors{i} = utile.rgb('ForestGreen');
    elseif strcmp( 'MediumWrap' ,grasp_names{i})
        colors{i} = 'b';
    elseif strcmp( 'PalmarPinch' ,grasp_names{i})
        colors{i} = 'c';
    elseif strcmp( 'Sphere3Finger' ,grasp_names{i})
        colors{i} = 'm';
    elseif strcmp('redballoon', grasp_names{i})
        colors{i} = 'k';
        
    elseif strcmp('Group Mean', grasp_names{i})
        colors{i} = 'k';
     elseif strcmp('Image', grasp_names{i})
        colors{i} = 'k';
     elseif strcmp('Speaking', grasp_names{i})
        colors{i} = 'b';
     elseif strcmp('Audio', grasp_names{i})
        colors{i} = 'm';
    elseif strcmp('Go Trial', grasp_names{i})
        colors{i} = 'g';
        
    elseif strcmp('No Go Trial', grasp_names{i})
        colors{i} = 'r';
        
    elseif strcmp('MotorImagery', grasp_names{i})
        colors{i} = 'g';
        
  %  elseif strcmp('Speaking', grasp_names{i})
  %      colors{i} = 'b';
     elseif strcmp('Yellow', grasp_names{i})
        colors{i} = utile.rgb('Fire'); %Fire
        
     elseif strcmp('Brown', grasp_names{i})
        colors{i} = utile.rgb('Brown');
        
     elseif strcmp('Gray', grasp_names{i})
        colors{i} = utile.rgb('SlateGray');
        
    elseif strcmp('Blue', grasp_names{i})
        colors{i} = utile.rgb('Blue');
        
    elseif strcmp('Green', grasp_names{i})
        colors{i} = utile.rgb('Green');
        
     elseif strcmp('Grasps', grasp_names{i})
        colors{i} = 'g';
        
     elseif strcmp('Colors', grasp_names{i})
        colors{i} = 'm';
     
    elseif strcmp('SMG', grasp_names{i})
        colors{i} = utile.rgb('Blue');
        
    elseif strcmp('PMV', grasp_names{i})
        colors{i} = utile.rgb('Green');
    
    elseif strcmp('MotorImageryLeft', grasp_names{i})
        colors{i} = utile.rgb('Fire');
        
    elseif strcmp('MotorImageryRight', grasp_names{i})
        colors{i} = utile.rgb('Green');
    
    elseif strcmp('S1X', grasp_names{i})
        colors{i} = utile.rgb('Fire');
    else
        
     
        
        error(['Grasp ' grasp_names ' not present, add color'])
        
    end 
end 

end

