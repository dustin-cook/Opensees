function [ ] = fn_stack_stories(input_dir, story_ht, story_dead_load, story_live_load, story_groups, model_id)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(story_groups)
    story.id(i) = i;
    story.model_id(i) = model_id;
    story.story_group_id(i) = story_groups(i);
    story.story_ht(i) = story_ht(i);
    story.story_dead_load(i) = story_dead_load(i);
    story.story_live_load(i) = story_live_load(i);
    if i == 1
        story.y_offset(i) = 0;
    else
        story.y_offset(i) = sum(story_ht(1:(i-1)));
    end
    
end

T = table(story.id', story.model_id', story.story_group_id', story.y_offset', story.story_ht', story.story_dead_load', story.story_live_load');
T.Properties.VariableNames = {'id','model_id','story_group_id','y_offset','story_ht', 'story_dead_load', 'story_live_load'};
writetable(T,[input_dir filesep 'story.csv'],'WriteVariableNames',true);
end

