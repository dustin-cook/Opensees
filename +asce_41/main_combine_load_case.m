function [ ] = main_combine_load_case( analysis, ele_prop_table )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

import asce_41.*
write_dir = [analysis.out_dir filesep 'asce_41_data'];
load([write_dir filesep 'node_analysis.mat'])
load([write_dir filesep 'model_analysis.mat'])

% Load Analysis Data
for i = 1:3
    read_dir = [analysis.out_dir filesep 'load_case_' num2str(i)];
    load([read_dir filesep 'story_analysis.mat']);
    load([read_dir filesep 'element_analysis.mat']);
    story_cases.(['load_case_' num2str(i)]) = story;
    element_cases.(['load_case_' num2str(i)]) = element;
end

%% Calculate Envelopes
% EDP Profiles
story.max_accel_x = max([story_cases.load_case_1.max_accel_x,story_cases.load_case_2.max_accel_x],[],2);
story.max_disp_x = max([story_cases.load_case_1.max_disp_x,story_cases.load_case_2.max_disp_x],[],2);
story.max_disp_center_x = max([story_cases.load_case_1.max_disp_center_x,story_cases.load_case_2.max_disp_center_x],[],2);
story.max_drift_x = max([story_cases.load_case_1.max_drift_x,story_cases.load_case_2.max_drift_x],[],2);
if isfield(story,'max_drift_z')
    story.max_accel_z = max([story_cases.load_case_1.max_accel_z,story_cases.load_case_2.max_accel_z],[],2);
    story.max_disp_z = max([story_cases.load_case_1.max_disp_z,story_cases.load_case_2.max_disp_z],[],2);
    story.max_disp_center_z = max([story_cases.load_case_1.max_disp_center_z,story_cases.load_case_2.max_disp_center_z],[],2);
    story.max_drift_z = max([story_cases.load_case_1.max_drift_z,story_cases.load_case_2.max_drift_z],[],2);
end

%% Torsional Amplification
story.A_x = ones(height(story),1);
story.A_z = ones(height(story),1);

% Calculate the actual torsional amplification multiplier
TAR_x = story_cases.load_case_1.torsional_factor_x;
if sum(strcmp('torsional_factor_z',story_cases.load_case_1.Properties.VariableNames)) > 0
    TAR_z = story_cases.load_case_1.torsional_factor_z;
end

% X direction
if sum(TAR_x > 1.2) > 0
    % For linear analysis calcuate the accidental torional amplification factor
    % (this only applies to forces and displacements caused by accidental
    % torsion)
    for i = 1:length(TAR_x)
        story.A_x(i) = min([max([TAR_x(i)/1.2;1])^2,3]);
    end
end
% Amplify Displacements
story.max_disp_x_mod = max(story.max_disp_x, story.max_disp_x + story_cases.load_case_1.max_disp_x*(max(story.A_x)) - story_cases.load_case_3.max_disp_x*max(story.A_x)); 
story.max_disp_center_x_mod = max(story.max_disp_center_x, story.max_disp_center_x + story_cases.load_case_1.max_disp_center_x*(max(story.A_x)) - story_cases.load_case_3.max_disp_center_x*max(story.A_x)); 
story.max_drift_x_mod = max(story.max_drift_x, story.max_drift_x + story_cases.load_case_1.max_drift_x*(max(story.A_x)) - story_cases.load_case_3.max_drift_x*max(story.A_x)); 

% Z direction
if exist('TAR_z','var')
    if sum(TAR_z > 1.2) > 0
        % For linear analysis calcuate the accidental torional amplification factor
        % (this only applies to forces and displacements caused by accidental
        % torsion)
        for i = 1:length(TAR_z)
            story.A_z(i) = min([max([TAR_z(i)/1.2;1])^2,3]);
        end
    end
    % Amplify Displacements
    story.max_disp_z_mod = max(story.max_disp_z, story.max_disp_z + story_cases.load_case_1.max_disp_z*(max(story.A_z)) - story_cases.load_case_3.max_disp_z*max(story.A_z)); 
    story.max_disp_center_z_mod = max(story.max_disp_center_z, story.max_disp_center_z + story_cases.load_case_1.max_disp_center_z*(max(story.A_z)) - story_cases.load_case_3.max_disp_center_z*max(story.A_z)); 
    story.max_drift_z_mod = max(story.max_drift_z, story.max_drift_z + story_cases.load_case_1.max_drift_z*(max(story.A_z)) - story_cases.load_case_3.max_drift_z*max(story.A_z)); 
end

% SRSS for element demands
story.A_max = max(story.A_x,story.A_z);

%% Element Forces
element.Mmax_1_mod = element.Mmax_1 + element_cases.load_case_1.Mmax_1*(max(story.A_max)) - element_cases.load_case_3.Mmax_1*max(story.A_max); 
element.Mmax_2_mod = element.Mmax_2 + element_cases.load_case_1.Mmax_2*(max(story.A_max)) - element_cases.load_case_3.Mmax_2*max(story.A_max); 
element.Vmax_1_mod = element.Vmax_1 + element_cases.load_case_1.Vmax_1*(max(story.A_max)) - element_cases.load_case_3.Vmax_1*max(story.A_max); 
element.Vmax_2_mod = element.Vmax_2 + element_cases.load_case_1.Vmax_2*(max(story.A_max)) - element_cases.load_case_3.Vmax_2*max(story.A_max); 
element.Pmax_mod = element.Pmax + element_cases.load_case_1.Pmax*(max(story.A_max)) - element_cases.load_case_3.Pmax*max(story.A_max); 

%% Acceptance
% Determine amplificaton factors
[ model, element, story ] = fn_linear_capacity_and_c_factors( model, story, element );

% Determine acceptance factors
[ element ] = main_m_factors( ele_prop_table, element );
[ hinge ] = fn_linear_hinge_accept(element, ele_prop_table, node);

%% Save Data
save([write_dir filesep 'model_analysis.mat'],'model')
save([write_dir filesep 'story_analysis.mat'],'story')
save([write_dir filesep 'element_analysis.mat'],'element')
save([write_dir filesep 'hinge_analysis.mat'],'hinge')

end

