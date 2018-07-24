function [ Vu, Vn, Vs ] = fn_aci_shear_capacity( fc, b, d, P, Av, fy, S, lambda, Ag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Inputs

%% Begin Method
% Concrete Capacity
Vc = 2*(1 + P/(2000*Ag))*lambda*b*d*sqrt(fc); 

% Reinforcement Capacity
if Av > 0
    Vs = Av*fy*d/S;
else
    Vs = 0;
end

% Total Capacity
phi = 0.75;
Vn = Vc + Vs;
Vu = phi*Vn;

end

