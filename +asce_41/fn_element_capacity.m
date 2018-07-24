function [ ele, ele_TH, ele_PM ] = fn_element_capacity( ele, ele_prop, ele_TH )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*
import aci_318.*

if ~ismember('DCR_total_raw', ele.Properties.VariableNames)
    ele.DCR_total_raw = 4.0; % Set to middle of interp if not already defined (ie create a k factor of 0.85)
end

% Axial Compression Capacity per ACI (use lower bound strength since assuming axial is force controlled)
[ ~, ele.Pn_aci_c, ~, ~ ] = fn_aci_axial_capacity( ele_prop.fc_n, ele_prop.a, ele_prop.As, ele_prop.fy_n );

% Axial Compression Capacity per ACI (use lower bound strength since assuming axial is force controlled)
[ ~, ~, ~, ele.Pn_aci_t ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );

% Shear capacity per ASCE 41
[ ele.Vn, ele.V0 ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele.Mmax, ele.Vmax, ele.Pmax, ele.DCR_total_raw );

% Shear Capacity per ACI
[ ~, ele.Vn_aci, ele.Vs_aci ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele.Pmax, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a );

% Modify Axial Force as Force Controlled Action (EQ 7-35)
[ Pmax_force_cotrolled ] = fn_force_controlled_action( ele.Pmax, ele.P_grav, 'cp', 'high', 1, 1 );

if strcmp(ele.type,'beam')
    % Moment Capcity per ACI
    [ ~, ele.Mn_aci_c ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, 0 );
    ele.Mn_aci_t = ele.Mn_aci_c;
    % Probable Moment Capcity
    [ ~, ele.Mp_c ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, 0 );
    ele.Mp_t = ele.Mp_c;
else
    % Moment Capcity per ACI
    [ ~, ele.Mn_aci_c ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, Pmax_force_cotrolled );

    % Moment Capcity per ACI
    [ ~, ele.Mn_aci_t ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, ele.Pmin );

    % Probable Moment Capcity
    [ ~, ele.Mp_c ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, ele.Pmax );

    % Probable Moment Capcity
    [ ~, ele.Mp_t ] = fn_aci_moment_capacity( ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, ele.Pmin );

    % Calculate PM Diagram Vectors
    % Percent Pmax
    vector_P = linspace(-0.9*ele.Pn_aci_t,ele.Pn_aci_c,25);
    for i = 1:length(vector_P)
        [ ~, vector_M(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, vector_P(i) );
    end
end

% Balanced Moment Capcity and Reinforcement Ratio
[ ~, ele.row_bal ] = fn_balanced_moment( ele_prop.fc_e, ele_prop.fy_e, ele_prop.w, ele_prop.d, ele_prop.As_d );

% Determin Flexure v Shear Critical
[ ele ] = fn_element_critical_mode( ele );

%% Calculate Capacity Time Histories
for i = 1:length(ele_TH.P_TH_1) %% ASSUMING P is uniform throughout member
    % Axial Capacity
    if ele_TH.P_TH_1(i) >= 0 
        ele_TH.Pn(i) = ele.Pn_aci_c; % Compressive Capacity
        [ ele_TH.P_force_controlled(i) ] = fn_force_controlled_action( ele_TH.P_TH_1(i), ele.P_grav, 'cp', 'high', 1, 1 );
    else
        ele_TH.Pn(i) = ele.Pn_aci_t; % Tensile Capacity
        ele_TH.P_force_controlled(i) = ele_TH.P_TH_1(i); % Tensile Capacity
    end
    % Shear Capacity
    if strcmp(ele.type,'column')
        [ ele_TH.Vn(i), ~ ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele.Mmax, ele.Vmax, ele_TH.P_force_controlled(i), ele.DCR_total_raw );
    else
        [ ~, ele_TH.Vn(i), ~ ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_TH.P_force_controlled(i), ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a );
    end  
end
    
if strcmp(ele.type,'beam')
    % Moment Capacity 
    ele_TH.Mn = ones(1,length(ele_TH.P_force_controlled))*ele.Mn_aci_c;
    ele_PM = [];
else
    % Moment Capacity 
    ele_TH.Mn = interp1(vector_P,vector_M,ele_TH.P_force_controlled);
    % Save PM Structure
    ele_PM.vector_P = [-ele.Pn_aci_t, vector_P, ele.Pn_aci_c];
    ele_PM.vector_M = [0, vector_M, 0];
end

% End Function
end

