function [ edp_profile ] = fn_calc_max_repsonse_profile( max_edp, story, ave )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

edp_profile = zeros(length(story.id),1);
for i = 1:length(story.id)
    if ave == 1 % take the mean for all nodes at level
        edp_profile(i) = mean(max_edp(story.nodes_on_slab{i}));
    else % Take the Max of all nodes at level
        edp_profile(i) = max(max_edp(story.nodes_on_slab{i}));
    end
end


end

