function [ ele ] = fn_element_capacity( ele, ele_prop )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

% Assign type to element table
ele.type = ele_prop.type;

% Rigid Elements
if ele.ele_id == 1 || ele.ele_id == 2
    ele.Vn = inf;
    ele.V0 = inf;
    ele.Vn_aci = inf;
    ele.Mn_aci = inf;
    ele.Mp = inf;
    ele.Pn = inf;
    ele.DCR_total_raw = 4.0;
    
% Flexible Elements
else
    if ~isfield(ele, 'DCR_total_raw')
        ele.DCR_total_raw = 4.0; % Set to middle of interp if not already defined (ie create a k factor of 0.85)
    end
    % Shear capacity per ASCE 41 (use lower bound strength since assuming shear is force controlled)
    [ ele.Vn, ele.V0 ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_n, ele_prop.d, ele_prop.S, ele_prop.lambda, ele_prop.fc_n, ele_prop.a, ele.Mmax, ele.Vmax, ele.Pmax, ele.DCR_total_raw );

    % Shear Capacity per ACI (use lower bound strength since assuming shear is force controlled)
    [ ~, ele.Vn_aci ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele.Pmax, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a );

    % Moment Capcity per ACI
    [ ~, ele.Mn_aci ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );

    % Probable Moment Capcity
    [ ~, ele.Mp ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.clear_cover, ele_prop.Es, ele.Pmax );

    % Axial Capacity per ACI (use lower bound strength since assuming axial is force controlled)
    [ ~, ele.Pn ] = fn_aci_axial_capacity( ele_prop.fc_n, ele_prop.a, ele_prop.As, ele_prop.fy_n );
end
end

