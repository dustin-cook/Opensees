function [ Pu, Pn ] = fn_aci_axial_capacity( fc, Ag, As, fy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));

%% Solve for Axial Capacity
Pn = 0.8*(0.85*fc*(Ag-sum(As))+fy*sum(As));
phi = 0.65;
Pu = phi*Pn;

end

