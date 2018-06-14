clear
close
clc

%% Import Packages
import asce_41.*

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = 'test';

%% Read in element and hinge data tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
m_table = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);

%% Go through each element and calculate the hinge properties
% CURRENTLY JUST FOR COLUMS, UPDATE TO BE FOR ALL MEMBERS
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    [ element(i,:) ] = fn_m_factors( m_table, ele, ele_props );
end

%% Save capacities to element database
writetable(element,[output_dir filesep 'element.csv'])