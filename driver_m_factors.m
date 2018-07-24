clear
close
clc

%% Import Packages
import asce_41.*

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '11DL11LL';

%% Read in element and m factor tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'element_analysis.mat'])
m_table.col = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.beam = readtable(['+asce_41' filesep 'linear_beam_m.csv'],'ReadVariableNames',true);

%% Go through each element and calculate the M factors
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    [ element_temp(i,:) ] = fn_m_factors( m_table, ele, ele_props );
end
element = element_temp;

%% Save capacities to element database
save([output_dir filesep 'element_analysis.mat'],'element')