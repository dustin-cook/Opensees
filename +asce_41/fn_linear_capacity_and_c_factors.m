function [ model, element, element_TH, element_PM ] = fn_linear_capacity_and_c_factors( model, story, ele_prop_table, element, element_TH, analysis )
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

c1c2 = 0;
percent_error = inf;
idx = 1;
while percent_error > 0.01 % Iterate on c factors until change is less than 1%
    c1c2_old = c1c2;

    % Calculate Element Capacities
    [ element, element_TH, element_PM ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis );
    vn_columns = element.Vn(strcmp(element.type,'column'));
    vn_progression(idx) = vn_columns(1);

    % Caclulate C1 and C2 Factors
    for i = 1:height(element)
        ele_dir = element.direction{i};
        [ element.c1(i), element.c2(i) ] = fn_c_factors( model.site_class{1}, model.num_stories, model.(['T1_' ele_dir]), model.(['hazus_class_' ele_dir]), model.DCR_raw_max );
    end

    % Calculate the DCR
    [ element, model.DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'NA' );

    % Calculate new iteration error
    c1c2 = element.c1 + element.c2;
    percent_error = max(abs((c1c2 - c1c2_old) ./ c1c2_old));
    idx = idx + 1;
    if idx > 100
        error('C factors will not converge')
    end
end
%% Amplify Final Element Forces and Displacements
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

