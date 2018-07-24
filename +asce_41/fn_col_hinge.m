function [ hinge ] = fn_col_hinge( ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Assmptions
% 1. Splice region contains at least two ties
% 2. 


%% Calculate terms
p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);

v_ratio = max([ele.Vn_aci/ele.V0 , 0.2]); % ASSUMING Vye = V from ACI, UPDATE

row_t = ele_props.Av/(ele_props.w*ele_props.S); % Should this be D or S?
As_tot = sum(str2num(strrep(strrep(ele_props.As{1},']',''),'[','')));
row_l = As_tot/(ele_props.w*ele_props.d);
if row_t < 0.0005
    error('Equations in table not valid: not enough transverse reinforcement')
end

%% Calculate condition
if ele.pass_aci_dev_length == 1
    condition = 1; % Not Controlled by development length
else
    condition = [1,2]; % Controlled by development length (but run through both to find min
end

for i = 1:length(condition)
    %% Caclulate Hinge Terms based on Table 10-8 of ASCE 41-17
    if condition(i) == 1
        row_t = min([row_t , 0.0175]);
        hinge_filt.a_hinge = max([0.042 - 0.043*p_ratio + 0.63*row_t - 0.23*v_ratio , 0]);
        if p_ratio <= 0.5
            hinge_filt.b_hinge = max([0.5/(5 + (p_ratio/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , hinge_filt.a_hinge]);
        else
            b_at_5 = max([0.5/(5 + (0.5/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , hinge_filt.a_hinge]);
            hinge_filt.b_hinge = max([interp1([0.5,0.7],[b_at_5,0],p_ratio) , hinge_filt.a_hinge]);
        end
        hinge_filt.c_hinge = max([0.24 - 0.4*p_ratio, 0]);
        hinge_filt.io = min([0.15*hinge_filt.a_hinge , 0.005]);

        % don't let p ratio go below 0.1 for LS and CP criteria
        p_ratio_ls_cp = max([p_ratio,0.1]);
        a_ls_cp = max([0.042 - 0.043*p_ratio_ls_cp + 0.63*row_t - 0.23*v_ratio , 0]);
        if p_ratio_ls_cp <= 0.5
            b_ls_cp = max([0.5/(5 + (p_ratio_ls_cp/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , a_ls_cp]);
        else
            b_at_5 = max([0.5/(5 + (0.5/0.8)*(1/row_t)*(ele_props.fc_e/ele_props.fy_e)) - 0.01 , a_ls_cp]);
            b_ls_cp = max([interp1([0.5,0.7],[b_at_5,0],p_ratio_ls_cp) , a_ls_cp]);
        end
        hinge_filt.ls = 0.5*b_ls_cp;
        hinge_filt.cp = 0.7*b_ls_cp;
    elseif condition(i) == 2
        row_t = min([row_t , 0.0075]);
        fy_e_trans = ele_props.fy_e; % Assumes transverse rien is the same strength as long
        hinge_filt.a_hinge = min([(row_t*fy_e_trans)/(8*row_l*ele_props.fy_e) , 0.025]);
        hinge_filt.b_hinge = min([max([0.012 - 0.085*p_ratio + 12*row_t , hinge_filt.a_hinge]) , 0.06]);
        hinge_filt.c_hinge = min([0.15 + 36*row_t , 0.4]);
    end

    %% Double Check only 1 row of the hinge table remains
    if length(hinge_filt.a_hinge) ~= 1
        error('Hinge table filtering failed to find unique result')
    end
    
    % Save to temp struct
    hinge_temp.a_hinge(i) = hinge_filt.a_hinge;
    hinge_temp.b_hinge(i) = hinge_filt.b_hinge;
    hinge_temp.c_hinge(i) = hinge_filt.c_hinge;
    hinge_temp.io(i) = hinge_filt.io;
    hinge_temp.ls(i) = hinge_filt.ls;
    hinge_temp.cp(i) = hinge_filt.cp;
end

% Find Min values of cases run
hinge.a_hinge = min(hinge_temp.a_hinge);
hinge.b_hinge = min(hinge_temp.b_hinge);
hinge.c_hinge = min(hinge_temp.c_hinge);
hinge.io = min(hinge_temp.io);
hinge.ls = min(hinge_temp.ls);
hinge.cp = min(hinge_temp.cp);
    
% Save to table
hinge = struct2table(hinge);
end

