function [ element, DCR_max_raw ] = fn_calc_dcr( element, element_TH, perform_level )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Assumuptions
% 1. M and V are deformation Controlled
% 2. Axial compression is force controlled while axial tension is
% deformation controlled

% Initialize Variable
DCR_M_raw_TH_1 = ones(length(element.id),length(element_TH.ele_1.M_TH_1));
DCR_M_raw_TH_2 = ones(length(element.id),length(element_TH.ele_1.M_TH_2));

for i = 1:length(element.id)
    ele_TH = element_TH.(['ele_' num2str(element.id(i))]);

    %% Calculate raw DCRs
    DCR_P_raw_TH(i,:) = ele_TH.P_TH_1 ./ ele_TH.Pn;
    DCR_V_raw_TH(i,:) = ele_TH.V_TH_1 ./ ele_TH.Vn(i);
    filter = (ele_TH.M_TH_1 >= 0);
    DCR_M_raw_TH_1(i,filter) = abs(ele_TH.M_TH_1(filter)) ./ ele_TH.Mn_pos_linear(filter);
    DCR_M_raw_TH_1(i,~filter) = abs(ele_TH.M_TH_1(~filter)) ./ ele_TH.Mn_neg_linear(~filter);
    filter = (ele_TH.M_TH_2 >= 0);
    DCR_M_raw_TH_2(i,filter) = abs(ele_TH.M_TH_2(filter)) ./ ele_TH.Mn_pos_linear(filter);
    DCR_M_raw_TH_2(i,~filter) = abs(ele_TH.M_TH_2(~filter)) ./ ele_TH.Mn_neg_linear(~filter);

    %% Modify DCR by m factor for deformation controlled actions
    DCR_M_TH_1(i,:) = DCR_M_raw_TH_1(i,:) / element.(['m_' perform_level])(i);
    DCR_M_TH_2(i,:) = DCR_M_raw_TH_2(i,:) / element.(['m_' perform_level])(i);
    DCR_V_TH(i,:) = DCR_V_raw_TH(i,:) / element.(['m_' perform_level])(i);

    %% Calculate Quf for force controlled actions
    DCR_P_TH(i,:) = ele_TH.P_force_controlled ./ ele_TH.Pn;

    %% Calc Max DCR for Each Element
    element.DCR_raw_max_M(i) = max([DCR_M_raw_TH_1(i,:),DCR_M_raw_TH_2(i,:)]);
    element.DCR_raw_max_V(i) = max(DCR_V_raw_TH(i,:));
    element.DCR_raw_max_P(i) = max(DCR_P_raw_TH(i,:));
    element.DCR_raw_max_all(i) = max([element.DCR_raw_max_M(i),element.DCR_raw_max_V(i),element.DCR_raw_max_P(i)]);
    element.DCR_max_M(i) = max([DCR_M_TH_1(i,:),DCR_M_TH_2(i,:)]);
    element.DCR_max_V(i) = max(DCR_V_TH(i,:));
    element.DCR_max_P(i) = max(DCR_P_TH(i,:));
    element.DCR_max_all(i) = max([element.DCR_max_M(i),element.DCR_max_V(i),element.DCR_max_P(i)]);
end

% Total DCR
DCR_max_raw = max(element.DCR_raw_max_all);



end

