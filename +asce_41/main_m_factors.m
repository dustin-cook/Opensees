function [ element ] = main_m_factors( ele_prop_table, element )
% Description: Main script for calculating ASCE 41 m-factors

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Import Packages
import asce_41.fn_m_factors

%% Go through each element and calculate the M factors
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    [ element_temp(i,:) ] = fn_m_factors( ele, ele_props );
end
element = element_temp;
end

