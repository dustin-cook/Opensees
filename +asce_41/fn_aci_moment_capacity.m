function [ Mu, Mn ] = fn_aci_moment_capacity( fc, b, d, P, As, As_d, fy, clear_cover )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Calculate intermediate/defualt Values
d_c = d - clear_cover;

%% Begin Method
% Find Location of Neutral Axis
c = 0;

% Calc compression zone
a = (As*fy + P)/(0.85*fc*b); % Assuming all steel is tension just tension steel and that steel is yielding

% Moment Capacity
Mn = a*0.85*fc*b*(d_c-a/2);
phi = 0.9; % Assuming tension steel strain is greater than 0.005
Mu = phi*Mn;


end

