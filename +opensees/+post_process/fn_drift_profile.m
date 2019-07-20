function [ max_drift_profile ] = fn_drift_profile( disp_TH, story, node, direction )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

max_drift_profile = zeros(length(story.id),1);
for s = 1:height(story)
    this_story = story.id(s);
    if this_story > 0
        nodes_this_story = node(node.story == this_story & node.on_slab == 1,:);
        if this_story == 1 && s == 1
            for n = 1:height(nodes_this_story)
                if nodes_this_story.record_disp(n)
                    node_id = nodes_this_story.id(n);
                    nodal_drifts_TH = disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_' direction '_TH'])/story.story_ht(s);
                    max_nodal_drifts(n) = max(abs(nodal_drifts_TH));
                else
                    max_nodal_drifts(n) = NaN;
                end
            end
            max_drift_profile(s) = max(max_nodal_drifts);
        else
            for n = 1:height(nodes_this_story)
                node_below = node(node.story == this_story-1 & node.x == nodes_this_story.x(n) & node.z == nodes_this_story.z(n) & node.on_slab == 1,:);
                if ~isempty(node_below) && nodes_this_story.record_disp(n)
                    node_id = nodes_this_story.id(n);
                    nodal_drifts_TH = (disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_' direction '_TH']) - disp_TH.(['node_' num2str(node_below.id) '_TH']).(['disp_' direction '_TH']))/story.story_ht(s);
                    max_nodal_drifts(n) = max(abs(nodal_drifts_TH));
                else
                    max_nodal_drifts(n) = NaN;
                end

            end
            max_drift_profile(s) = max(max_nodal_drifts);
        end
    end
end

end

