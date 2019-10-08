function [ max_twist_profile ] = fn_twist_profile( disp_TH, story, node, direction )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

record_nodes = node(node.record_disp == 1,:);
max_twist_profile = zeros(length(story.id),1);
for s = 1:height(story)
    story_nodes = record_nodes(record_nodes.story == story.id(s),:);
    node_center = story_nodes(story_nodes.center == 1,:);
    eastern_nodes = story_nodes(story_nodes.x == max(story_nodes.x),:);
    [~, idx] = min(abs(eastern_nodes.z - node_center.z));
    node_east = eastern_nodes(idx,:);
    nodal_twist_TH = (disp_TH.(['node_' num2str(node_east.id) '_TH']).(['disp_' direction '_TH']) - disp_TH.(['node_' num2str(node_center.id) '_TH']).(['disp_' direction '_TH']));
    max_twist_profile(s) = max(abs(nodal_twist_TH));
end

end

