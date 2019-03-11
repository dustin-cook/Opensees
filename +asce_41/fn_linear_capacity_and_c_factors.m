function [ model, element, joint ] = fn_linear_capacity_and_c_factors( model, story, ele_prop_table, element, analysis, joint_table, read_dir, write_dir )
% Description: Main script that post process an ASCE 41 analysis

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:

%% Import Packages
import asce_41.fn_calc_dcr
import asce_41.fn_c_factors
import asce_41.main_element_capacity

%% Calculation of C factors
% Calculate Element Capacities
[ element, ~ ] = main_element_capacity( story, ele_prop_table, element, analysis, joint_table, read_dir, write_dir );

% Default C1 and C2 Factors to 1.0
element.c1 = ones(height(element),1);
element.c2 = ones(height(element),1);

% Calculate the DCR
[ element, model.DCR_raw_max ] = fn_calc_dcr( element, 'NA', read_dir );

% Caclulate C1 and C2 Factors (Uses equation C7-3 for U-strength, and therefore doesn't need to iterate. However, I think it only applies to LSP)
for i = 1:height(element)
    ele_dir = element.direction{i};
    [ element.c1(i), element.c2(i) ] = fn_c_factors( model.site_class{1}, model.num_stories, model.(['T1_' ele_dir]), model.(['hazus_class_' ele_dir]), model.DCR_raw_max );
end

%% Amplify Final Element Forces and Displacements
% Calculate Element Capacities again
[ element, joint ] = main_element_capacity( story, ele_prop_table, element, analysis, joint_table, read_dir, write_dir );

% Forces
element.Pmax = element.Pmax.*element.c1.*element.c2;
element.Vmax = element.Vmax.*element.c1.*element.c2;
element.Mmax = element.Mmax.*element.c1.*element.c2;

% X Displacements
filter = strcmp(element.direction,'x');
c_factor_modifier = max(element.c1(filter))*max(element.c2(filter));
story.max_disp_x = story.max_disp_x*c_factor_modifier;
story.max_drift_x = story.max_drift_x*c_factor_modifier;

% Z Displacements
if isfield(story,'max_disp_z')
    filter = strcmp(element.direction,'z');
    c_factor_modifier = max(element.c1(filter))*max(element.c2(filter));
    story.max_disp_z = story.max_disp_z*c_factor_modifier;
    story.max_drift_z = story.max_drift_z*c_factor_modifier;
end

end

