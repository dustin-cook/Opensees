clear
close
clc

%% Import Packages
import asce_41.*
import plotting_tools.*

%% Define Analysis and Model parameters
analysis.model_id = 12;
analysis.gm_id = 8;
analysis.name = 'output_fix_polly';

%% Read in element and hinge data tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'element_analysis.mat'])

% element = element(:,1:41);

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
    
%     % Plot Hinges
%     if strcmp(ele.type,'beam') || strcmp(ele.type,'column')
%         theta_yeild = ele.Mn_aci_pos*(ele.length)/(4*ele_props.e*ele_props.iz);
%         Q_y = ele.Mn_aci_pos;
%         Q_ult = ele.Mp_pos;
%         post_yeild_slope = min([((Q_ult-Q_y)/Q_y)/hinge.a_hinge,0.1*(1/theta_yeild)]);
%         force_vector = [0,1,post_yeild_slope*hinge.a_hinge+1,hinge.c_hinge,hinge.c_hinge];
%         disp_vector = [0, theta_yeild, theta_yeild+hinge.a_hinge, theta_yeild+hinge.a_hinge+(hinge.b_hinge-hinge.a_hinge)/2, theta_yeild+hinge.b_hinge];
%         plot(disp_vector,force_vector)
%         xlabel('Rotation')
%     elseif strcmp(ele.type,'wall')
%         elastic_shear_stiffness = ele_props.g*ele_props.a/ele.length; 
%         f1 = hinge.f_hinge*ele.Vn_aci;
%         u1 = f1/elastic_shear_stiffness;
%         f2 = ele.Vn_aci;
%         u2 = (hinge.g_hinge/100)*ele.length;
%         f3 = ele.Vn_aci*1.001;
%         u3 = (hinge.d_hinge/100)*ele.length;
%         f4 = ele.Vn_aci*ele.c_hinge;
%         u4 = (hinge.e_hinge/100)*ele.length;
%         force_vector = [0,f1,f2,f3,f4]/f1;
%         disp_vector = [0,u1,u2,u3,u4];
%         plot(disp_vector,force_vector)
%         xlabel('Drift')
%     end
%     ylabel('Q/Qy')
%     fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , ['element_' num2str(ele.id)], 2 )
    % save as element hinge table
    if strcmp(ele.type,'wall')
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