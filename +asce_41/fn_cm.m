function [ c_m ] = fn_cm( num_stories, T, hazus_class )
%UNTITLED Summary of this function goes here
%   Calculates Cm factor based on Table 7-4 from ASCE 41-13
%   Created by: Dustin Cook 5/4/18

if num_stories >= 3 && T < 1
    build_class_table = readtable(['inputs' filesep 'building_class.csv'],'ReadVariableNames',true);
    build_class = build_class_table(strcmp(build_class_table.hazus_class,hazus_class),:);
    c_m = build_class.asce_41_13_cm(1);
else
    c_m = 1;
end

end

