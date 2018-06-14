%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*

%% Read in element table
analysis.model_id = 3;
analysis.name = 'test';
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    
%     if ele_id == 1 || ele_id == 2
%         element.Vn(i) = 999999999;
%         element.V0(i) = 999999999;
%         element.Vn_aci(i) = 999999999;
%         element.Mn_aci(i) = 999999999;
%         element.Mp(i) = 999999999;
%         element.Pn(i) = 999999999;
%     else
        [ element(i,:) ] = fn_element_capacity( ele, ele_prop );
        
%         % Shear capacity per ASCE 41
%         [ element.Vn(i), element.V0(i) ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele.Mmax, ele.Vmax, ele.Pmax );
% 
%         % Shear Capacity per ACI
%         [ ~, element.Vn_aci(i) ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele.Pmax, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a );
% 
%         % Moment Capcity per ACI
%         [ ~, element.Mn_aci(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );
% 
%         % Probable Moment Capcity
%         [ ~, element.Mp(i) ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );
% 
%         % Axial Capacity per ACI
%         [ ~, element.Pn(i) ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );
%     end
    
end

%% Save capacities to element database
writetable(element,[output_dir filesep 'element.csv'])