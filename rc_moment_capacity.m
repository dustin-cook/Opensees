%% Calculate Concrete Element Flexural Capacity Based on ACI 318-11
clear all
close all
clc

%% Inputs
fc = 6000; % psi
b = 7.5; % in
d = 900; % in
p = 242144; % pounds

as = 0.6169*[1 1 1 1];
as_d = [d d-16 d-16*2 d-16*3];
fy = 40000;

% fc = 6000; % psi
% b = 12; % in
% d = 300; % in
% p = 310000; % pounds
% 
% as = 0.6169*4;
% fy = 40000;

%% Calculate intermediate/defualt Values

%% Begin Method
% Find Location of Neutral Axis
c = 0
a = (as*fy + p)/(0.85*fc*b); % Assuming just tension steel and that steel is yielding

% Moment Capacity
phi = 0.9; % Assuming tension steel strain is greater than 0.005
mu = phi*a*0.85*fc*b*(d-a/2);
disp(['Mu = ' num2str(round(mu/1000)) ' kip-in'])