function [ ] = fn_assemble_story(input_dir, grid_id, start, story_group_id, direction)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(grid_id)
    story_group.id(i) = i;
    story_group.story_group_id(i) = story_group_id;
    story_group.grid_line_id(i) = grid_id(i);
    story_group.orientation(i) = direction(i);
    if story_group.orientation(i) == 1
        story_group.x_start(i) = 0;
        story_group.z_start(i) = start(i);
    elseif story_group.orientation(i) == 3
        story_group.x_start(i) = start(i);
        story_group.z_start(i) = 0;
    else
        error('Grid Orientation not Recognized')
    end

end
T = table(story_group.id', story_group.story_group_id', story_group.grid_line_id', story_group.orientation', story_group.x_start', story_group.z_start');
T.Properties.VariableNames = {'id' 'story_group_id' 'grid_line_id' 'orientation' 'x_start' 'z_start'};
writetable(T,[input_dir filesep 'story_group.csv'],'WriteVariableNames',true);
end

