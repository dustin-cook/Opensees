function [ control_node, control_dof, max_displacement, step_size ] = fn_pushover_properties( first_story_node, analysis, story )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

control_node = first_story_node(end);
if strcmp(analysis.pushover_direction,'x')
    control_dof = 1;
elseif strcmp(analysis.pushover_direction,'z')
    control_dof = 3;
end
max_displacement = analysis.pushover_drift*sum(story.story_ht);
step_size = max_displacement / analysis.pushover_num_steps;
        
end

