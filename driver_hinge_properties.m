clear
close
clc

%% Import Packages
import asce_41.*
import plotting_tools.*

%% Define Analysis and Model parameters
analysis.model_id = 11;
analysis.gm_id = 6;
analysis.name = 'test';

%% Read in element and hinge data tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'element_analysis.mat'])

%% Go through each element and calculate the hinge properties
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    if strcmp(ele.type,'beam')
        [ hinge ] = fn_beam_hinge( ele, ele_props );
    elseif strcmp(ele.type,'column')
        [ hinge ] = fn_col_hinge( ele, ele_props );
    elseif strcmp(ele.type,'wall')
        [ hinge ] = fn_wall_hinge( ele, ele_props );
    end
    
    % Plot Hinges
    plot_name = ['element_' num2str(ele.id)];
    fn_plot_backbone( ele, ele_props, hinge, output_dir, plot_name, 1)

    % save as element hinge table
    if strcmp(ele.type,'wall') && strcmp(ele.critical_mode,'shear')
        element.c_hinge(i,1) = hinge.c_hinge;
        element.d_hinge(i,1) = hinge.d_hinge;
        element.e_hinge(i,1) = hinge.e_hinge;
        element.f_hinge(i,1) = hinge.f_hinge;
        element.g_hinge(i,1) = hinge.g_hinge;
    else
        element.a_hinge(i,1) = hinge.a_hinge;
        element.b_hinge(i,1) = hinge.b_hinge;
        element.c_hinge(i,1) = hinge.c_hinge;
    end
    element.io(i) = hinge.io;
    element.ls(i) = hinge.ls;
    element.cp(i) = hinge.cp;
end

%% Save capacities to element database
save([output_dir filesep 'element_analysis.mat'],'element')
writetable(element,[output_dir filesep 'element_linear.csv'])