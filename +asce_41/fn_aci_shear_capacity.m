function [ Vu, Vn ] = fn_aci_shear_capacity( fc, b, d, P, Av, fy, S, lambda, Ag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Inputs

%% Begin Method
% Concrete Capacity
vc = 2*(1 + P/(2000*Ag))*lambda*b*d*sqrt(fc); 

% Reinforcement Capacity
if Av > 0
    vs = Av*fy*d/S;
else
    vs = 0;
end

% Total Capacity
phi = 0.75;
Vn = vc + vs;
Vu = phi*Vn;

end

