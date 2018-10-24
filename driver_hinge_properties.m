clear
close
clc

%% Import Packages
import asce_41.*
import plotting_tools.*

%% Define Analysis and Model parameters
analysis.model_id = 8;
analysis.gm_id = 6;
analysis.name = 'nonlinear';

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
    theta_yeild = ele.Mn_aci_pos*(ele.length)/(4*ele_props.e*ele_props.iz);
    Q_y = ele.Mn_aci_pos;
    Q_ult = ele.Mp_pos;
    post_yeild_slope = min([((Q_ult-Q_y)/Q_y)/hinge.a_hinge,0.1*(1/theta_yeild)]);
    force_vector = [0,1,post_yeild_slope*hinge.a_hinge+1,hinge.c_hinge,hinge.c_hinge];
    disp_vector = [0, theta_yeild, theta_yeild+hinge.a_hinge, theta_yeild+hinge.a_hinge+(hinge.b_hinge-hinge.a_hinge)/2, theta_yeild+hinge.b_hinge]; % ASSUMING y = d/2 NEED TO UPDATE
    plot(disp_vector,force_vector)
    ylabel('Q/Qy')
    xlabel('Total Rotation of Hinge')
    fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , ['element_' num2str(ele.id)], 2 )

    % save as element hinge table
    element.a_hinge(i) = hinge.a_hinge;
    element.b_hinge(i) = hinge.b_hinge;
    element.c_hinge(i) = hinge.c_hinge;
    element.io(i) = hinge.io;
    element.ls(i) = hinge.ls;
    element.cp(i) = hinge.cp;
end

%% Save capacities to element database
save([output_dir filesep 'element_analysis.mat'],'element')
writetable(element,[output_dir filesep 'element_linear.csv'])