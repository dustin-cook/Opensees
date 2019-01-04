function [ element, element_TH, element_PM ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis )
% Description: Main script that calculates the strength of each element in 
% the model according to ASCE 41 

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:


%% Import Packages
import asce_41.*

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

end

