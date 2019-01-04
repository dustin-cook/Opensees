function [ model, element ] = fn_basic_analysis_properties( model, story, element )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Building Properties
model.num_stories = length(story.id);

% Default DCR
model.DCR_raw_max = 4; % Set to middle of interp if not already defined (ie create a k factor of 0.85)
element.DCR_raw_max_V = model.DCR_raw_max*ones(height(element),1);

end

