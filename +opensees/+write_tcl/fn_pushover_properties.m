function [ control_node, control_dof, max_displacement, step_size ] = fn_pushover_properties( first_story_node, analysis, story )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

control_node = first_story_node(end);
if strcmp(analysis.pushover_direction,'x')
    control_dof = 1;
    max_displacement = analysis.pushover_drift_x*sum(story.story_ht);
elseif strcmp(analysis.pushover_direction,'z')
    control_dof = 3;
    max_displacement = analysis.pushover_drift_z*sum(story.story_ht);
end
step_size = max_displacement / analysis.pushover_num_steps;
        
end

