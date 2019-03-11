function [ max_drift_profile ] = fn_drift_profile( disp_TH, story, node, direction )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

max_drift_profile = zeros(length(story.id),1);
for i = 1:length(story.id)
    nodes_this_story = node(node.story == i,:);
    if i == 1
        for n = 1:height(nodes_this_story)
            if nodes_this_story.record_disp(n)
                node_id = nodes_this_story.id(n);
                nodal_drifts_TH = disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_' direction '_TH'])/story.story_ht(i);
                max_nodal_drifts(n) = max(abs(nodal_drifts_TH));
            else
                max_nodal_drifts(n) = NaN;
            end
        end
        max_drift_profile(i) = max(max_nodal_drifts);
    else
        for n = 1:height(nodes_this_story)
            node_below = node(node.story == i-1 & node.x == nodes_this_story.x(n) & node.z == nodes_this_story.z(n),:);
            if ~isempty(node_below) && nodes_this_story.record_disp(n)
                node_id = nodes_this_story.id(n);
                nodal_drifts_TH = (disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_' direction '_TH']) - disp_TH.(['node_' num2str(node_below.id) '_TH']).(['disp_' direction '_TH']))/story.story_ht(i);
                max_nodal_drifts(n) = max(abs(nodal_drifts_TH));
            else
                max_nodal_drifts(n) = NaN;
            end
            
        end
        max_drift_profile(i) = max(max_nodal_drifts);
    end
end

end

