function [ hinge ] = fn_joint_hinge( joint )
% Description: Filters the hinge properties to find nonlinear parameters
% for a single joint

% Created By: Dustin Cook
% Date Created: 1/9/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import asce_41.fn_filter_asce41_table

% Load Beam Hinge Table 10-7 from ASCE 41-17
hinge_table = readtable(['+asce_41' filesep 'joint_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id

%% Begin Method
% Determine Condition
if strcmp(joint.class,'a') || strcmp(joint.class,'b')
    condition = 1; % Interior Joint
else
    condition = 2; % All other Joints
end
hinge_filt = hinge_table(hinge_table.condition == condition,:);

% Axial Load Ratio
p_ratio = joint.Pmax/(joint.a*joint.fc_e); % Currently assuming area and concrete strength comes from the colomn above
[ hinge_filt ] = fn_filter_asce41_table( hinge_filt, p_ratio, 'p_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

% Transverse Reinforcement
hinge_filt = hinge_filt(strcmp(hinge_filt.trans_rien,joint.trans_rien),:);

% Shear Ratio
v_ratio = joint.Vmax/joint.Vj;
[ hinge_filt ] = fn_filter_asce41_table( hinge_filt, v_ratio, 'v_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

%% Double Check only 1 row of the hinge table remains
if length(hinge_filt.a_hinge) ~= 1
    error('Hinge table filtering failed to find unique result')
end

%% Save filtered hinge to final table structure
hinge = hinge_filt;
end

