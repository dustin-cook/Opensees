function [ max_drift_profile ] = fn_drift_profile( disp, story, node )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

max_drift_profile = zeros(length(story.id),1);
for i = 1:length(story.id)
    if i == 1
        nodal_drifts_x = disp(node.id(node.story == i),:)/story.story_ht(i);
    else
        nodal_drifts_x = (disp(node.id(node.story == i & node.mass > 0),:) - disp(node.id(node.story == i-1 & node.mass > 0),:))/story.story_ht(i);
    end
    max_drifts_all_slab_x = max(abs(nodal_drifts_x),[],2);
    max_drift_profile(i) = max(max_drifts_all_slab_x);
end

end

