function [ max_drift_profile ] = fn_drift_profile( disp, story )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

max_drift_profile = zeros(length(story.id),1);
for i = 1:length(story.id)
    if i == 1
        nodal_drifts_x = disp(story.nodes_on_slab{i},:)/story.story_ht(i);
    else
        nodal_drifts_x = (disp(story.nodes_on_slab{i},:) - disp(story.nodes_on_slab{i-1},:))/story.story_ht(i);
    end
    max_drifts_all_slab_x = max(abs(nodal_drifts_x),[],2);
    max_drift_profile(i) = max(max_drifts_all_slab_x);
end

end

