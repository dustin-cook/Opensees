function [ ] = fn_assemble_story(input_dir, story_group_id)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
id = 0;

for j = 1:length(story_group_id)
    for i = 1:length(story_group_id{j}.grid_id)
        id = id + 1;
        story_group.id(id) = i;
        story_group.story_group_id(id) = j;
        story_group.grid_line_id(id) = story_group_id{j}.grid_id(i);
        story_group.orientation(id) = story_group_id{j}.direction(i);
        if story_group.orientation(id) == 1
            story_group.x_start(id) = 0;
            story_group.z_start(id) = story_group_id{j}.start(i);
        elseif story_group.orientation(i) == 3
            story_group.x_start(id) = story_group_id{j}.start(i);
            story_group.z_start(id) = 0;
        else
            error('Grid Orientation not Recognized')
        end
    end
end

T = table(story_group.id', story_group.story_group_id', story_group.grid_line_id', story_group.orientation', story_group.x_start', story_group.z_start');
T.Properties.VariableNames = {'id' 'story_group_id' 'grid_line_id' 'orientation' 'x_start' 'z_start'};
writetable(T,[input_dir filesep 'story_group.csv'],'WriteVariableNames',true);

end

