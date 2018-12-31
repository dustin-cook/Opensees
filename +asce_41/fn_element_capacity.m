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

% Axial Tension Capacity per ACI (use lower bound strength since assuming axial is force controlled)
[ ~, ~, ~, ele.Pn_aci_t ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );

% Shear capacity per ASCE 41
[ ele.Vn, ele.V0 ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.As_d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele.Mmax, ele.Vmax, ele.Pmax, ele.DCR_total_raw );

% Check if walls
if strcmp(ele.type,'wall')
    % Is it force Controlled?
    rho_n = ele_prop.Av/(ele_prop.S*ele_prop.w);
    if rho_n < 0.0015 % ASCE 41-17 10.7.2.3
        error('Wall is force controlled, modify element')
    end
    
    % Are the axial loads too high for lateral resistance
    if ele.Pmax > 0.35*ele.Pn_aci_c 
        error('Wall has too much axial load to take lateral force, modify model')
    end
else
    rho_n = nan;
end

% Shear Capacity per ACI
[ ele.Vu_aci, ele.Vn_aci, ele.Vs_aci ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele.Pmax, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a, ele_prop.hw, ele.type, rho_n, ele_prop.As_d );

% Modify Axial Force as Force Controlled Action (EQ 7-35)
[ Pmax_force_cotrolled ] = fn_force_controlled_action( ele.Pmax, ele.P_grav, 'cp', 'high', 1, 1 );

if contains(ele_prop.description,'rigid')
    ele.Mn_aci_pos = inf;
    ele.Mn_aci_neg = inf;
    ele.Mp_pos = inf;
    ele.Mp_neg = inf;
else
    if strcmp(ele.type,'beam') 
        % Moment Capcity per ACI
        [ ~, ele.Mn_aci_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        [ ~, ele.Mn_aci_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        % Probable Moment Capcity
        [ ~, ele.Mp_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        [ ~, ele.Mp_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
    else
        % Moment Capcity per ACI
        [ ~, ele.Mn_aci_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, ele.P_grav, 0, 0 );
        [ ~, ele.Mn_aci_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, ele.P_grav, 0, 0 );

        % Probable Moment Capcity
        [ ~, ele.Mp_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, ele.P_grav, 0, 0 );
        [ ~, ele.Mp_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, ele.P_grav, 0, 0 );

        % Calculate PM Diagram Vectors
        % Percent Pmax
        vector_P = linspace(-0.9*ele.Pn_aci_t,ele.Pn_aci_c,25);
        for i = 1:length(vector_P) % Currently Assumes Column has uniform strength in each directions (ie symmetric layout)
            [ ~, vector_M(i) ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, vector_P(i), 0, 0 );
        end
    end
end

% Balanced Moment Capcity and Reinforcement Ratio
[ ele.row_bal ] = fn_balanced_moment( ele_prop.fc_e, ele_prop.fy_e );

% Determine Flexure v Shear Critical
[ ele ] = fn_element_critical_mode( ele, ele_prop );

% Check if walls are force controlled
if strcmp(ele.type,'wall')
    % For walls controlled by shear, if 
    if strcmp(ele.critical_mode,'shear') && ele.Pmax > 0.15*ele_prop.a*ele_prop.fc_e % ASCE 41-17 table 10-20 note b
        error('Wall is force controlled, too much axial load')
    end
end

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
        [ ele_TH.Vn(i), ~ ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.As_d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, max(abs([ele_TH.M_TH_1(i),ele_TH.M_TH_2(i)])), abs(ele_TH.V_TH_1(i)), ele_TH.P_force_controlled(i), ele.DCR_total_raw );
    else
        [ ~, ele_TH.Vn(i), ~ ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_TH.P_force_controlled(i), ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a, ele_prop.hw, ele.type, rho_n, ele_prop.As_d );
    end  
end

if strcmp(ele.type,'beam')
    % Moment Capacity 
    ele_TH.Mn_pos = ones(1,length(ele_TH.P_force_controlled))*ele.Mn_aci_pos;
    ele_TH.Mn_neg = ones(1,length(ele_TH.P_force_controlled))*ele.Mn_aci_neg;
    ele_PM = [];
else
    % Moment Capacity 
    ele_TH.Mn_pos = interp1(vector_P,vector_M,ele_TH.P_force_controlled);
    ele_TH.Mn_neg = ele_TH.Mn_pos; %assumes columns are the same in both directions
    % Save PM Structure
    ele_PM.vector_P = [-ele.Pn_aci_t, vector_P, ele.Pn_aci_c];
    ele_PM.vector_M = [0, vector_M, 0];
end

% End Function
end

