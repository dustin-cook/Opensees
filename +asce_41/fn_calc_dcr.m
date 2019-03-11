function [ element, DCR_raw_max ] = fn_calc_dcr( element, perform_level, read_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Assumuptions
% 1. M and V are deformation Controlled
% 2. Axial compression is force controlled while axial tension is
% deformation controlled

% Initialize Variable
DCR_P_raw_TH = ones(length(element.id),length(element_TH.ele_1.P_TH_1));
DCR_V_raw_TH = ones(length(element.id),length(element_TH.ele_1.V_TH_1));
DCR_M_raw_TH_1 = ones(length(element.id),length(element_TH.ele_1.M_TH_1));
DCR_M_raw_TH_2 = ones(length(element.id),length(element_TH.ele_1.M_TH_2));

for i = 1:height(element)
    load([read_dir filesep 'element_TH_' num2str(element.id(i)) '.mat'])

    %% Calculate raw DCRs (multipled by C factors)
    DCR_P_raw_TH(i,:) = element.c1(i) * element.c2(i) * ele_TH.P_TH_1 ./ ele_TH.Pn;
    DCR_V_raw_TH(i,:) = element.c1(i) * element.c2(i) * ele_TH.V_TH_1 ./ ele_TH.Vn(i);
    filter = (ele_TH.M_TH_1 >= 0);
    DCR_M_raw_TH_1(i,filter) = element.c1(i) * element.c2(i) * abs(ele_TH.M_TH_1(filter)) ./ ele_TH.Mn_pos_linear(filter);
    DCR_M_raw_TH_1(i,~filter) = element.c1(i) * element.c2(i) * abs(ele_TH.M_TH_1(~filter)) ./ ele_TH.Mn_neg_linear(~filter);
    filter = (ele_TH.M_TH_2 >= 0);
    DCR_M_raw_TH_2(i,filter) = element.c1(i) * element.c2(i) * abs(ele_TH.M_TH_2(filter)) ./ ele_TH.Mn_pos_linear(filter);
    DCR_M_raw_TH_2(i,~filter) = element.c1(i) * element.c2(i) * abs(ele_TH.M_TH_2(~filter)) ./ ele_TH.Mn_neg_linear(~filter);

    %% Calc Max DCR for Each Element
    element.DCR_raw_max_M(i) = max([DCR_M_raw_TH_1(i,:),DCR_M_raw_TH_2(i,:)]);
    element.DCR_raw_max_V(i) = max(DCR_V_raw_TH(i,:));
    element.DCR_raw_max_P(i) = max(DCR_P_raw_TH(i,:));
    element.DCR_raw_max_all(i) = max([element.DCR_raw_max_M(i),element.DCR_raw_max_V(i),element.DCR_raw_max_P(i)]);
    
    if ~strcmp(perform_level,'NA')
        %% Modify DCR by m factor for deformation controlled actions
        DCR_M_TH_1(i,:) = DCR_M_raw_TH_1(i,:) / element.(['m_' perform_level])(i);
        DCR_M_TH_2(i,:) = DCR_M_raw_TH_2(i,:) / element.(['m_' perform_level])(i);
        DCR_V_TH(i,:) = DCR_V_raw_TH(i,:) / element.(['m_' perform_level])(i);

        %% Calculate Quf for force controlled actions
        DCR_P_TH(i,:) = ele_TH.P_TH_linear ./ ele_TH.Pn;

        %% Calc Max DCR for Each Element
        element.DCR_max_M(i) = max([DCR_M_TH_1(i,:),DCR_M_TH_2(i,:)]);
        element.DCR_max_V(i) = max(DCR_V_TH(i,:));
        element.DCR_max_P(i) = max(DCR_P_TH(i,:));
        element.DCR_max_all(i) = max([element.DCR_max_M(i),element.DCR_max_V(i),element.DCR_max_P(i)]);
    end
end

% Total DCR
DCR_raw_max = max(element.DCR_raw_max_all);



end

