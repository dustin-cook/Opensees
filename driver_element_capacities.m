%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*


%% Read in element table
analysis.model_id = 3;
analysis.name = '11DL11LL';
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'element_TH.mat'])
load([output_dir filesep 'element_analysis.mat'])

for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
    
    %% Calculate Element Capacties
    [ ele, element_TH.(['ele_' num2str(element.id(i))]), element_PM.(['ele_' num2str(element.id(i))]) ] = fn_element_capacity( ele, ele_prop, ele_TH );
    disp([num2str(i), ' out of ', num2str(length(element.id)) ' elements complete' ])
    
    %% Caculate required development length and make sure there is enough
    [ ele.pass_aci_dev_length ] = fn_development_check( ele, ele_prop );
    ele_to_save(i,:) = ele;
end

element = ele_to_save;
%% Save capacities to element database
save([output_dir filesep 'element_analysis.mat'],'element')
save([output_dir filesep 'element_TH.mat'],'element_TH')
save([output_dir filesep 'element_PM.mat'],'element_PM')