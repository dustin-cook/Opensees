function [ element, DCR_raw_max ] = fn_calc_dcr( element, perform_level )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Assumuptions
% 1. M and V are deformation Controlled
% 2. Axial compression is force controlled while axial tension is
% deformation controlled

import asce_41.fn_force_controlled_action

for i = 1:height(element)
    ele = element(i,:);
    
    if ~ele.elastic && ~ele.rigid
        %% Calculate raw DCRs (straight from analysis)
        element.DCR_raw_max_P(i) = ele.Pmax / ele.Pn_c; % Assumes compression contolled
        element.DCR_raw_max_V_1(i) = ele.Vmax_1 / ele.Vn_1;
        element.DCR_raw_max_V_2(i) = ele.Vmax_2 / ele.Vn_2;
        element.DCR_raw_max_M_1(i) = ele.Mmax_1 ./ min(ele.Mn_pos_1,ele.Mn_neg_1);
        element.DCR_raw_max_M_2(i) = ele.Mmax_2 ./ min(ele.Mn_pos_2,ele.Mn_neg_2);
        element.DCR_raw_max_all(i) = max([element.DCR_raw_max_M_1(i),element.DCR_raw_max_M_2(i),element.DCR_raw_max_V_1(i),element.DCR_raw_max_V_2(i),element.DCR_raw_max_P(i)]);

        %% Calculated modified DCRs (c1c2 and torsion modifications)
        element.DCR_mod_max_P(i) = ele.Pmax_mod / ele.Pn_c; % Assumes compression contolled
        element.DCR_mod_max_V_1(i) = ele.Vmax_1_mod / ele.Vn_1;
        element.DCR_mod_max_V_2(i) = ele.Vmax_2_mod / ele.Vn_2;
        element.DCR_mod_max_M_1(i) = ele.Mmax_1_mod ./ min(ele.Mn_pos_1,ele.Mn_neg_1);
        element.DCR_mod_max_M_2(i) = ele.Mmax_2_mod ./ min(ele.Mn_pos_2,ele.Mn_neg_2);
        element.DCR_mod_max_all(i) = max([element.DCR_raw_max_M_1(i),element.DCR_raw_max_M_2(i),element.DCR_raw_max_V_1(i),element.DCR_raw_max_V_2(i),element.DCR_raw_max_P(i)]);

        if ~strcmp(perform_level,'NA')
            %% Modify DCR by m factor for deformation controlled actions
            element.DCR_max_M_1(i) = element.DCR_mod_max_M_1(i) / ele.(['m_' perform_level '_1']);
            element.DCR_max_M_2(i) = element.DCR_mod_max_M_2(i) / ele.(['m_' perform_level '_2']);
            element.DCR_max_V_1(i) = element.DCR_mod_max_V_1(i) / ele.(['m_' perform_level '_1']);
            element.DCR_max_V_2(i) = element.DCR_mod_max_V_2(i) / ele.(['m_' perform_level '_2']);
            [ element.Pmax_fc(i) ] = fn_force_controlled_action( ele.Pmax_mod, ele.P_grav, perform_level, 'high', ele.c1, ele.c2 );
            element.DCR_max_P(i) = element.Pmax_fc(i) / ele.Pn_c; % assumes compression
            element.DCR_max_all(i) = max([element.DCR_max_M_1(i),element.DCR_max_M_2(i),element.DCR_max_V_1(i),element.DCR_max_V_2(i),element.DCR_max_P(i)]);
        end
    end
end

% Total DCR
DCR_raw_max = max(element.DCR_raw_max_all);



end

