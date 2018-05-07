function [ Pu, Pn ] = fn_aci_axial_capacity( fc, Ag, As, fy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Pn = 0.8*(0.85*fc*(Ag-As)+fy*As);
phi = 0.65;
Pu = phi*Pn;

end

