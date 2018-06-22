function [ element, DCR_max_raw ] = fn_calc_dcr( element, perform_level, c1, c2, seismicity )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

DCR_P_raw = element.Pmax ./ element.Pn;
DCR_M_raw = element.Mmax ./ element.Mn_aci;
for i = 1:length(element.id)
    if strcmp(element.type{i},'column')
        DCR_V_raw(i,1) = element.Vmax(i) / element.Vn(i);
    elseif strcmp(element.type{i},'beam')
        DCR_V_raw(i,1) = element.Vmax(i) / element.Vn_aci(i);
    else
        DCR_V_raw(i,1) = 1;
    end
end
DCR_max_raw = max([DCR_P_raw; DCR_V_raw; DCR_M_raw]);

% Modify DCR by m factor for deformation controlled (ASSUME ALL MOMENTS ARE DEFORMATION CONTROLLED)
element.DCR_M = DCR_M_raw ./ element.(['m_' perform_level]);

% Calculate Quf for force controlled actions (ASSUME ALL AXIAL AND SHEAR ARE FORCE CONTROLLED)
if strcmp(perform_level,'cp')
    x = 1;
else
    x = 1.3;
end
if strcmp(seismicity,'high')
    j = 2;
elseif strcmp(seismicity,'high')
    j = 1.5;
else
    j = 1;
end
P_quf_factors = x/(c1*c2*j); % ONLY SUPPOSED TO APPLY TO EQ DEMANDS, NOT GRAVITY (UPDATE)
v_quf_factors = x/(c1*c2*j);

element.DCR_P = P_quf_factors*DCR_P_raw;
element.DCR_V = v_quf_factors*DCR_V_raw;


end

