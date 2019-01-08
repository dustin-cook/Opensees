function [ element, element_TH, element_PM, joint ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis, joint_table )
% Description: Main script that calculates the strength of each element in 
% the model according to ASCE 41 

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:


%% Import Packages
import asce_41.*

%% Begin Method
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_prop_table.id == ele_id,:);
    ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
    
    %% Calculate Element Capacties
    [ ele, element_TH.(['ele_' num2str(element.id(i))]), element_PM.(['ele_' num2str(element.id(i))]) ] = fn_element_capacity( story, ele, ele_prop, ele_TH, analysis.nonlinear );
    disp([num2str(i), ' out of ', num2str(length(element.id)) ' elements complete' ])
    
    %% Caculate required development length and make sure there is enough
    [ ele.pass_aci_dev_length ] = fn_development_check( ele, ele_prop );
    ele_to_save(i,:) = ele;
end

element = ele_to_save;

%% Calculate Beam Column Strength 
for i =1:length(joint_table.id)
    joint = joint_table(i,:);
    beam_left = element(element.id == joint.beam_left,:); % Maximum of beam pos and neg nominal bending strength
    beam_right = element(element.id == joint.beam_right,:); % Maximum of beam pos and neg nominal bending strength 
    column_low = element(element.id == joint.column_low,:); % Minimum of column pos and negative nominal moment strength
    column_high = element(element.id == joint.column_high,:); % Minimum of column pos and negative nominal moment strength
    beam_strength_1 = sum([beam_left.Mn_pos,beam_right.Mn_neg]); % case 1: 1 pos and 1 neg bending
    beam_strength_2 = sum([beam_left.Mn_neg,beam_right.Mn_pos]); % case 2: 1 pos and 1 neg bending the other way
    joint.beam_strength = min([beam_strength_1,beam_strength_2]); % Min of two cases
    col_strength_1 = sum([column_low.Mn_pos,column_high.Mn_pos]); % case 1: positive bending for both columns
    col_strength_2 = sum([column_low.Mn_neg,column_high.Mn_neg]); % case 2: negative bending for both columns
    joint.column_strength = min([col_strength_1,col_strength_2]); % Min of two cases
    joint.col_bm_ratio = joint.column_strength/joint.beam_strength;
    
    joint_to_save(i,:) = joint;
end

joint = joint_to_save;

end

