function [ edp_profile ] = fn_calc_max_repsonse_profile( max_edp, story, node, ave )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

edp_profile = zeros(length(story.id),1);
for s = 1:height(story)
    if ave == 1 % take the mean for all nodes at level
        edp_profile(s) = mean(max_edp(node.story == story.id(s) & node.on_slab == 1 & ~isnan(max_edp)));
    else % Take the Max of all nodes at level
        edp_profile(s) = max(abs(max_edp(node.story == story.id(s) & node.on_slab == 1 & ~isnan(max_edp))));
    end
end


end

