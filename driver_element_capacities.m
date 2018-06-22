%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*

%% Read in element table
analysis.model_id = 3;
analysis.name = 'NL_10DL10LL';
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable([output_dir filesep 'element_linear.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    
    [ element_temp(i,:) ] = fn_element_capacity( ele, ele_prop );
    disp([num2str(i), ' out of ', num2str(length(element.id)) ' elements complete' ])
end

element = element_temp;
%% Save capacities to element database
writetable(element,[output_dir filesep 'element_linear.csv'])