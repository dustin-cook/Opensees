function [ ele ] = fn_element_capacity( ele, ele_prop )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

% Rigid Elements
if ele.ele_id == 1 || ele.ele_id == 2
    ele.Vn = inf;
    ele.V0 = inf;
    ele.Vn_aci = inf;
    ele.Mn_aci = inf;
    ele.Mp = inf;
    ele.Pn = inf;
    
% Flexible Elements
else
    % Shear capacity per ASCE 41
    [ ele.Vn, ele.V0 ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele.Mmax, ele.Vmax, ele.Pmax );

    % Shear Capacity per ACI
    [ ~, ele.Vn_aci ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele.Pmax, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a );

    % Moment Capcity per ACI
    [ ~, ele.Mn_aci ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );

    % Probable Moment Capcity
    [ ~, ele.Mp ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );

    % Axial Capacity per ACI
    [ ~, ele.Pn ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );
end
end

