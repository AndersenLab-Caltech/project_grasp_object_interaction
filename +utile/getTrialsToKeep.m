function [subset] = getTrialsToKeep(date, Cue)
%Returns the trials that do not have artefacts for specific session days.
%Depending on the session, the trials with artefacts are different. 
%subset will contain the indexes of the trials that do not have artefacts
%in them


if strcmp(date, '20191016')
    
    if strcmp(Cue, 'Speaking_Colors_ActionPhase_8s')
        subset = [2,4,5,6,8,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26,29,30,31,34,35,36,38,40];
        
    elseif strcmp(Cue,'MotorImagery_vs_GraspSpeaking_vs_ColorSpeaking_8s_ActionPhase')
        subset = [2,4,5,6,10,12,13,14,16,17,18,19,20,21,22,23,24,25,26,29,30,31,34,35,38,40];
    end
    
    
end   


end

