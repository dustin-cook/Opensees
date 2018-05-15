%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*

%% Read in element table
ele_file = ['element.csv'];
element_table = readtable(ele_file,'ReadVariableNames',true);

for i = 1:length(element_table.id)
    ele = element_table(i,:);
    % Shear capacity per ASCE 41
    [ element_table.Vn(i), element_table.V0(i) ] = fn_shear_capacity( ele.Av, ele.fy_e, ele.d, ele.S, ele.lambda, ele.fc_e, ele.Ag, ele.M, ele.V, ele.P );
    
    % Shear Capacity per ACI
    [ element_table.Vu_aci(i), element_table.Vn_aci(i) ] = fn_aci_shear_capacity( ele.fc_e, ele.b, ele.d, ele.P, ele.Av, ele.fy_e, ele.S, ele.lambda, ele.Ag );
    
    % Moment Capcity per ACI
    [ element_table.Mu_aci(i), element_table.Mn_aci(i) ] = fn_aci_moment_capacity( ele.fc_e, ele.b, ele.d, ele.As, ele.As_d, ele.fy_e, ele.clear_cover, ele.Es );
    
    % Probable Moment Capcity
    [ ~, element_table.Mp(i) ] = fn_aci_moment_capacity( ele.fc_e*1.15, ele.b, ele.d, ele.As, ele.As_d, ele.fy_e*1.15, ele.clear_cover, ele.Es );
    
    % Axial Capacity per ACI
    [ element_table.Pu(i), element_table.Pn(i) ] = fn_aci_axial_capacity( ele.fc_e, ele.Ag, ele.As, ele.fy_e );
    
end

%% Save capacities to element database
writetable(element_table,ele_file)