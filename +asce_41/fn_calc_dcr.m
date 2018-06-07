function [ element, DCR_max ] = fn_calc_dcr( element )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

element.DCR_P = element.Pmax ./ element.Pn;
element.DCR_V = element.Vmax ./ element.Vn_aci;
element.DCR_M = element.Mmax ./ element.Mn_aci;
DCR_max = max([element.DCR_P; element.DCR_V; element.DCR_M]);

end

