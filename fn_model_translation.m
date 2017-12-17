function [ I_bm, A_bm, E_bm, I_col, A_col, E_col, damp_ratio, story_weight_k, foundation_fix, story_ht_in, bay_width_in, num_stories, num_bays, floor_ht, story_mass, building_ht, building_mass ] = fn_model_translation( model )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Object discritization
num_stories = model.num_stories;
num_bays = model.num_bays;
bay_width_in = model.bay_width_in;
damp_ratio = model.damp_ratio;


story_ht_in = str2num(strrep(strrep(model.story_ht_in{1},'[',''),']',''));
foundation_fix = str2num(strrep(strrep(model.foundation_fix{1},'[',''),']',''));
story_weight_k = str2num(strrep(strrep(model.story_weight_k{1},'[',''),']',''));
E_col = str2num(strrep(strrep(model.E_col{1},'[',''),']',''));
A_col = str2num(strrep(strrep(model.A_col{1},'[',''),']',''));
I_col = str2num(strrep(strrep(model.I_col{1},'[',''),']',''));
E_bm = str2num(strrep(strrep(model.E_bm{1},'[',''),']',''));
A_bm = str2num(strrep(strrep(model.A_bm{1},'[',''),']',''));
I_bm = str2num(strrep(strrep(model.I_bm{1},'[',''),']',''));

% Object Methods

floor_ht = zeros(length(story_ht_in+1),1);
for i = 1:length(story_ht_in)
    floor_ht(i+1) = sum(story_ht_in(1:i));
end
story_mass = story_weight_k/386;
building_ht = sum(story_ht_in);
building_mass = sum(story_mass);
end

