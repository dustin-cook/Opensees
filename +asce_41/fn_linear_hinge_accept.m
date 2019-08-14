function [ hinge_table ] = fn_linear_hinge_accept(element, ele_prop_table, node)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.fn_calc_dcr

%% Define DCRs for each acceptance level
[ele_cp, ~] = fn_calc_dcr( element, 'cp' );
[ele_ls, ~] = fn_calc_dcr( element, 'ls' );
[ele_io, ~] = fn_calc_dcr( element, 'io' );

%% Go through each element twice and define the hinges
count = 0;
for e = 1:height(element)
    if ~element.rigid(e) && ~element.elastic(e)
        ele_prop = ele_prop_table(ele_prop_table.id == element.ele_id(e),:);
        for i = 1:2
            % Identification
            count = count + 1;
            hinge.id(count,1) = count;
            hinge.element_id(count,1) = element.id(e);
            hinge.ele_side(count,1) = i;
            hinge.ele_direction{count,1} = element.direction{e};
            hinge.story(count,1) = element.story(e);
            hinge.node(count,1) = element.(['node_' num2str(i)])(e);
            node_x = node.x(node.id == element.(['node_' num2str(i)])(e));
            if i == 1 && node_x == 1571
                hinge.damage_recorded(count,1) = ele_prop.(['damage_' num2str(i) '_east']);
            else
                hinge.damage_recorded(count,1) = ele_prop.(['damage_' num2str(i)]);
            end
            hinge.direction{count,1} = 'primary';

            % Demands
            hinge.V_ratio(count,1) = element.(['DCR_raw_max_V_' num2str(i)])(e);
            hinge.M_ratio(count,1) = element.(['DCR_raw_max_M_' num2str(i)])(e);
            hinge.P_ratio_asce41(count,1) = element.Pmax(e) / element.Pn_c(e); % Assumes compression
            hinge.P_ratio_nominal(count,1) = element.Pmax(e) / (ele_prop.a*ele_prop.fc_n);
            hinge.P_ratio_expected(count,1) = element.Pmax(e) / (ele_prop.a*ele_prop.fc_e);
            
            % Modified Demands
            hinge.V_ratio_mod(count,1) = element.(['DCR_mod_max_V_' num2str(i)])(e);
            hinge.M_ratio_mod(count,1) = element.(['DCR_mod_max_M_' num2str(i)])(e);
            hinge.P_ratio_asce41_mod(count,1) = element.Pmax_mod(e) / element.Pn_c(e); % Assumes compression
            hinge.P_ratio_nominal_mod(count,1) = element.Pmax_mod(e) / (ele_prop.a*ele_prop.fc_n);
            hinge.P_ratio_expected_mod(count,1) = element.Pmax_mod(e) / (ele_prop.a*ele_prop.fc_e);

            % Factors
            hinge.m_cp(count,1) = element.(['m_cp_' num2str(i)])(e);
            hinge.m_ls(count,1) = element.(['m_ls_' num2str(i)])(e);
            hinge.m_io(count,1) = element.(['m_io_' num2str(i)])(e);
            
            % Accepance DCR
            hinge.V_ratio_accept_cp(count,1) = ele_cp.(['DCR_max_V_' num2str(i)])(e);
            hinge.M_ratio_accept_cp(count,1) = ele_cp.(['DCR_max_M_' num2str(i)])(e);
            hinge.P_ratio_accept_cp(count,1) = ele_cp.DCR_max_P(e); % Assumes compression
            hinge.V_ratio_accept_ls(count,1) = ele_ls.(['DCR_max_V_' num2str(i)])(e);
            hinge.M_ratio_accept_ls(count,1) = ele_ls.(['DCR_max_M_' num2str(i)])(e);
            hinge.P_ratio_accept_ls(count,1) = ele_ls.DCR_max_P(e); % Assumes compression
            hinge.V_ratio_accept_io(count,1) = ele_io.(['DCR_max_V_' num2str(i)])(e);
            hinge.M_ratio_accept_io(count,1) = ele_io.(['DCR_max_M_' num2str(i)])(e);
            hinge.P_ratio_accept_io(count,1) = ele_io.DCR_max_P(e); % Assumes compression

            % Acceptance
            if  ele_io.DCR_max_all(e) <= 1
                hinge.accept(count,1) = 1; % Passes IO
            elseif ele_ls.DCR_max_all(e) <= 1
                hinge.accept(count,1) = 2; % Passes LS
            elseif ele_cp.DCR_max_all(e) <= 1
                hinge.accept(count,1) = 3; % Passes CP
            else
                hinge.accept(count,1) = 4; % Fails all performance levels
            end
        end
    end
end

% Hinge to table
hinge_table = struct2table(hinge);

end

