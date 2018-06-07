function [ edp_profile ] = fn_calc_max_repsonse_profile( edp, story, ave, unit )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

max_edp_all_nodes = max(abs(edp),[],2)/unit;
edp_profile = [max_edp_all_nodes(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    if ave == 1 % take the mean for all nodes at level
        edp_profile(i+1) = mean(max_edp_all_nodes(story.nodes_on_slab{i}));
    else % Take the Max of all nodes at level
        edp_profile(i+1) = max(max_edp_all_nodes(story.nodes_on_slab{i}));
    end
end


end

