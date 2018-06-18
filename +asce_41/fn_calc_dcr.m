function [ element, DCR_max_raw ] = fn_calc_dcr( element, perform_level, c1, c2, seismicity )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

DCR_P_raw = element.Pmax ./ element.Pn;
DCR_V_raw = element.Vmax ./ element.Vn_aci;
DCR_M_raw = element.Mmax ./ element.Mn_aci;
DCR_max_raw = max([DCR_P_raw; DCR_V_raw; DCR_M_raw]);

% Modify DCR by m factor for deformation controlled (ASSUME ALL MOMENTS ARE DEFORMATION CONTROLLED)
element.DCR_M = element.Mmax ./ (element.Mn_aci .* element.(['m_' perform_level]));

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
P_quf = element.Pmax*x/(c1*c2*j); % ONLY SUPPOSED TO APPLY TO EQ DEMANDS, NOT GRAVITY (UPDATE)
v_quf = element.Vmax*x/(c1*c2*j);

element.DCR_P = P_quf ./ element.Pn;
element.DCR_V = v_quf ./ element.Vn_aci;


end

