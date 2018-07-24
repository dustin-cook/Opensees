function [ Pu_c, Pn_c, Pu_t, Pn_t ] = fn_aci_axial_capacity( fc, Ag, As, fy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Scope
% 1. Non prestressed members only

%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));

%% Solve for Axial Compression Capacity
Pn_c = 0.8*(0.85*fc*(Ag-sum(As))+fy*sum(As));
phi = 0.65;
Pu_c = phi*Pn_c;

%% Solve for Axial Tension Capacity
Pn_t = fy*sum(As);
phi = 0.9;
Pu_t = phi*Pn_t;

end

