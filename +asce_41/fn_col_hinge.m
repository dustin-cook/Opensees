function [ hinge ] = fn_col_hinge( ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Calculate terms
p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);

v_ratio = max([ele.Vn_aci/ele.V0 , 0.2]); % ASSUMING Vye = V from ACI, UPDATE

row_t = ele_props.Av/(ele_props.w*ele_props.S); % Should this be D or S?
As_tot = sum(str2num(strrep(strrep(ele_props.As{1},']',''),'[','')));
row_l = As_tot/(ele_props.w*ele_props.d);
if row_t < 0.0005
    error('Equations in table not valid: not enough transverse reinforcement')
end

% Calculate condition
condition = 1; % ASSUMPTION FOR NOW, UPDATE

% Caclulate Hinge Terms based on Table 10-8 of ASCE 41-17
if condition == 1
    row_t = min([row_t , 0.0175]);
    hinge.a_hinge = max([0.042 - 0.043*p_ratio + 0.63*row_t - 0.23*v_ratio , 0]);
    if p_ratio <= 0.5
        hinge.b_hinge = max([0.5/(5 + (p_ratio/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , hinge.a_hinge]);
    else
        b_at_5 = max([0.5/(5 + (0.5/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , hinge.a_hinge]);
        hinge.b_hinge = max([interp1([0.5,0.7],[b_at_5,0],p_ratio) , hinge.a_hinge]);
    end
    hinge.c_hinge = max([0.24 - 0.4*p_ratio, 0]);
    hinge.io = min([0.15*hinge.a_hinge , 0.005]);
    
    % don't let p ratio go below 0.1 for LS and CP criteria
    p_ratio_ls_cp = max([p_ratio,0.1]);
    a_ls_cp = max([0.042 - 0.043*p_ratio_ls_cp + 0.63*row_t - 0.23*v_ratio , 0]);
    if p_ratio_ls_cp <= 0.5
        b_ls_cp = max([0.5/(5 + (p_ratio_ls_cp/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , a_ls_cp]);
    else
        b_at_5 = max([0.5/(5 + (0.5/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , a_ls_cp]);
        b_ls_cp = max([interp1([0.5,0.7],[b_at_5,0],p_ratio_ls_cp) , a_ls_cp]);
    end
    hinge.ls = 0.5*b_ls_cp;
    hinge.cp = 0.7*b_ls_cp;
elseif condition == 2
    row_t = min([row_t , 0.0075]);
    fy_e_trans = ele_props.fy_e; % Assumes transverse rien is the same strength as long
    hinge.a_hinge = min([(row_t*fy_e_trans)/(8*row_l*ele_props.fy_e) , 0.025]);
    hinge.b_hinge = min([max([0.012 - 0.085*p_ratio + 12*row_t , hinge.a_hinge]) , 0.06]);
    hinge.c_hinge = min([0.15 + 36*row_t , 0.4]);
end

hinge = struct2table(hinge);

% Double Check only 1 row of the hinge table remains
if length(hinge.a_hinge) ~= 1
    error('Hinge table filtering failed to find unique result')
end


end

