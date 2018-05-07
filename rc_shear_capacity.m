%% Calculate Concrete Element Shear Capacity Based on ACI 318-11
clear all
close all
clc

%% Inputs
% fc = 6000; % psi
% b = 7.5; % in
% d = 900; % in
% p = 242144; % pounds
% 
% ast = (5/16)^2*pi;
% fy = 40000;
% s = 12;

fc = 6000; % psi
b = 12; % in
d = 300; % in
p = 310000; % pounds

ast = 0.6169*2;
fy = 40000;
s = 16;

%% Calculate intermediate/defualt Values
ag = b*d; % Assuming Rectangular Section
lambda = 1.0; % Default to Normal Weight Concrete

%% Begin Method
% Concrete Capacity
vc = 2*(1 + p/(2000*ag))*lambda*b*d*sqrt(fc); 

% Reinforcement Capacity
if ast > 0
    vs = ast*fy*d/s;
else
    vs = 0;
end

% Total Capacity
phi = 0.75;
vu = phi*(vc + vs);
disp(['Vu = ' num2str(round(vu/1000)) ' kips'])