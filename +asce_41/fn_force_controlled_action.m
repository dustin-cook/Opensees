function [ mod_force ] = fn_force_controlled_action( raw_force, grav_force, perform_level, seismicity, c1, c2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Calculate Quf for force controlled actions (ASSUME ALL AXIAL IS FORCE CONTROLLED)
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

eq_force = raw_force - grav_force;
mod_force = grav_force + x*eq_force/(c1*c2*j);
end

