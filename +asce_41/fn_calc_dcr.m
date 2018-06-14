function [ element, DCR_max_raw ] = fn_calc_dcr( element, perform_level )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

DCR_P_raw = element.Pmax ./ element.Pn;
DCR_V_raw = element.Vmax ./ element.Vn_aci;
DCR_M_raw = element.Mmax ./ element.Mn_aci;
DCR_max_raw = max([DCR_P_raw; DCR_V_raw; DCR_M_raw]);

% Modify DCR by m factor (CURRENTLY ASSUME ALL ARE DEFORMATION CONTROLLED)
element.DCR_P = element.Pmax ./ (element.Pn .* element.(['m_' perform_level]));
element.DCR_V = element.Vmax ./ (element.Vn_aci .* element.(['m_' perform_level]));
element.DCR_M = element.Mmax ./ (element.Mn_aci .* element.(['m_' perform_level]));

end

