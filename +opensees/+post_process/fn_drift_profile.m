function [ max_drift_profile ] = fn_drift_profile( disp, story, node )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

max_drift_profile = zeros(length(story.id),1);
for i = 1:length(story.id)
    nodal_drifts_x = [];
    max_nodal_drifts = [];
    if i == 1
        nodal_drifts_x = disp(node.story == i,:)/story.story_ht(i);
        max_nodal_drifts = max(abs(nodal_drifts_x),[],2);
        max_drift_profile(i) = max(max_nodal_drifts);
    else
        nodes_below = node(node.story == i-1,:);
        count = 0;
        for j = 1:length(nodes_below.id)
            node_this_story = (node.story == i & node.x == nodes_below.x(j) & node.z == nodes_below.z(j));
            if sum(node_this_story) > 0
                count = count +1;
                nodal_drifts_x(count,:) = (disp(node_this_story,:) - disp(node.id == nodes_below.id(j),:))/story.story_ht(i);
            end
        end
        max_nodal_drifts = max(abs(nodal_drifts_x),[],2);
        max_drift_profile(i) = max(max_nodal_drifts);
    end
end

end

