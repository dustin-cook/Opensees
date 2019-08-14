function [ mod_force ] = fn_force_controlled_action( raw_force, grav_force, perform_level, seismicity, c1, c2 )
% Decription: Function to modify demands for force controlled actions 
% according to ASCE 41-17

% Created By: Dustin Cook
% Date Created: 1/2/2019

%% Begin Method
% Calculate Quf for force controlled actions (ASSUME ALL AXIAL IS FORCE CONTROLLED)
if strcmp(perform_level,'cp')
    x = 1;
else
    x = 1.3;
end
if strcmp(seismicity,'high')
    j = 2;
elseif strcmp(seismicity,'moderate')
    j = 1.5;
else
    j = 1;
end

eq_force = raw_force - grav_force;
mod_force = grav_force + x*eq_force/(c1*c2*j); % ASCE 41-17 eq 7-35
end

