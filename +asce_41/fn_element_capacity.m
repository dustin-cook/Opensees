function [ ele, ele_TH, ele_PM ] = fn_element_capacity( story, ele, ele_prop, ele_TH, nonlinear )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*
import aci_318.*

%% Calc Axial Capacity
% Axial Compression Capacity per ACI (use lower bound strength since assuming axial is force controlled)
[ ~, ele.Pn_c, ~, ~ ] = fn_aci_axial_capacity( ele_prop.fc_n, ele_prop.a, ele_prop.As, ele_prop.fy_n );

% Axial Tension Capacity per ACI (use expected strength since tension is deformation controlled)
[ ~, ~, ~, ele.Pn_t ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );

% Axial Capacity Time History
for i = 1:length(ele_TH.P_TH_1)
    % Axial Capacity
    if ele_TH.P_TH_1(i) >= 0 
        ele_TH.Pn(i) = ele.Pn_c; % Compressive Capacity is the same for each timestep
        [ ele_TH.P_TH_linear(i) ] = fn_force_controlled_action( ele_TH.P_TH_1(i), ele.P_grav, 'cp', 'high', 1, 1 ); % Use force controlled axial loads for linear procedures
    else
        ele_TH.Pn(i) = ele.Pn_t; % Tensile Capacity is the same for each timestep
        ele_TH.P_TH_linear(i) = ele_TH.P_TH_1(i); % Tensile Capacity for linear procedures is not force controlled
    end
end

%% Calc Moment Capacity
if contains(ele_prop.description,'rigid')
    ele.Mn_pos = inf;
    ele.Mn_neg = inf;
    ele.Mp_pos = inf;
    ele.Mp_neg = inf;
    ele_PM = [];
    ele_TH.Mn_pos = inf;
    ele_TH.Mn_neg = inf;
    ele_TH.Mp_pos = inf;
    ele_TH.Mp_neg = inf;
    ele_TH.Mn_pos_linear = inf;
    ele_TH.Mn_neg_linear = inf;
else
    if strcmp(ele.type,'beam') % beams
        % Moment Capcity per ACI (assume no axial loads for beams)
        [ ~, ele.Mn_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        [ ~, ele.Mn_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        [ ~, ele.Mp_pos ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        [ ~, ele.Mp_neg ] = fn_aci_moment_capacity( 'neg', ele_prop.fc_e*1.15, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e*1.15, ele_prop.Es, 0, ele_prop.slab_depth, ele_prop.b_eff );
        % Moment Capacity Time History
        ele_TH.Mn_pos = ones(1,length(ele_TH.P_TH_1))*ele.Mn_pos;
        ele_TH.Mn_neg = ones(1,length(ele_TH.P_TH_1))*ele.Mn_neg;
        ele_TH.Mp_pos = ones(1,length(ele_TH.P_TH_1))*ele.Mp_pos;
        ele_TH.Mp_neg = ones(1,length(ele_TH.P_TH_1))*ele.Mp_neg;
        ele_TH.Mn_pos_linear = ele_TH.Mn_pos;
        ele_TH.Mn_neg_linear = ele_TH.Mn_neg;
        ele_PM = [];
    else % columns and walls
        % PM Interactions
        vector_P = linspace(-0.9*ele.Pn_t,ele.Pn_c,25); % axial force range
        for i = 1:length(vector_P) % Currently Assumes Column has uniform strength in each directions (ie symmetric layout)
            [ ~, vector_M(i) ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, vector_P(i), 0, 0 );
            [ ~, vector_Mp(i) ] = fn_aci_moment_capacity( 'pos', ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, 1.15*ele_prop.fy_e, ele_prop.Es, vector_P(i), 0, 0 );
        end
        % Save PM Structure
        ele_PM.vector_P = [-ele.Pn_t, vector_P, ele.Pn_c];
        ele_PM.vector_M = [0, vector_M, 0];
        % Moment Capacity Time History
        ele_TH.Mn_pos = interp1(vector_P,vector_M,ele_TH.P_TH_1);
        ele_TH.Mn_neg = ele_TH.Mn_pos; % assumes columns are the same in both directions
        ele_TH.Mp_pos = interp1(vector_P,vector_Mp,ele_TH.P_TH_1);
        ele_TH.Mp_neg = ele_TH.Mn_pos; % assumes columns are the same in both directions
        ele_TH.Mn_pos_linear = interp1(vector_P,vector_M,ele_TH.P_TH_linear);
        ele_TH.Mn_neg_linear = ele_TH.Mn_pos_linear; % assumes columns are the same in both directions
        % Moment Capcity
        ele.Mn_pos = min(ele_TH.Mn_pos); % Best estimate moment is the minumum from the time history (most conservative)
        ele.Mn_neg = min(ele_TH.Mn_pos);
        ele.Mp_pos = min(ele_TH.Mn_pos);
        ele.Mp_neg = min(ele_TH.Mn_pos);
    end
end

%% Shear Capacity
% Vye and Diplacement Ductility
[ ele ] = fn_disp_ductility( ele, ele_prop, story );

% Shear capacity is not a function of time
if strcmp(ele.type,'column')
    % Determine Ductility Factors
    if nonlinear % for nonlinear use the displacement ductility
        ductility_factor = ele.disp_duct; 
    else % for linear use max DCR
        ductility_factor = ele.DCR_raw_max_V; 
    end
    
    % The yield displacement is the lateral displacement of the column, determined using the effective rigidities 
    % from Table 10-5, at a shear demand resulting in flexural yielding of the plastic hinges, VyE.
    [ ele.Vn, ele.V0 ] = fn_shear_capacity( ele_prop.Av, ele_prop.fy_e, ele_prop.As_d, ele_prop.S, ele_prop.lambda, ele_prop.fc_e, ele_prop.a, ele_TH.M_TH_1, ele_TH.V_TH_1, ele.P_grav, ductility_factor );
    ele.Vs = NaN;
else
    [ ~, ele.Vn, ele.Vs ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.Av, ele_prop.fy_e, ele_prop.S, ele_prop.lambda, ele_prop.a, ele_prop.hw, ele.type, ele_prop.As_d );
    ele.V0 = NaN;
end 

% Shear capacity time history is uniform
ele_TH.Vn = ones(1,length(ele_TH.P_TH_1))*ele.Vn;


%% Perform Checks
% Balanced Moment Capcity and Reinforcement Ratio
[ ele.row_bal ] = fn_balanced_moment( ele_prop.fc_e, ele_prop.fy_e );

% Determine Flexure v Shear Critical
[ ele ] = fn_element_critical_mode( ele, ele_prop );

% Wall Checks
if strcmp(ele.type,'wall')
    % For walls controlled by shear, if axial loads are too high
    if strcmp(ele.critical_mode,'shear') && ele.Pmax > 0.15*ele_prop.a*ele_prop.fc_e % ASCE 41-17 table 10-20 note b
        error('Wall is force controlled, too much axial load')
    end
    % Are the axial loads too high for lateral resistance
    if ele.Pmax > 0.35*ele.Pn_c 
        error('Wall has too much axial load to take lateral force, modify model')
    end
end

% End Function
end

