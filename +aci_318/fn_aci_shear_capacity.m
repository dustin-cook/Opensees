function [ Vu, Vn, Vs ] = fn_aci_shear_capacity( fc, b, d, P, Av, fy, S, lambda, Ag, hw, type, rho_t, As_d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Shear Walls
if strcmp(type,'wall')
    Acv = b*d;
    aspect_ratio = hw/d;
    alpha_c = min(max((-2)*(aspect_ratio-2) + 2, 2), 3);
    Vn = Acv*(alpha_c*lambda*sqrt(fc) + rho_t*fy); % ACI 318-14 eq 18.10.4.1
    Vs = nan;

%% Beams and Columns
else
    % Reformat steel depth
    As_d = str2double(strsplit(strrep(strrep(As_d{1},']',''),'[',''),','));
    d_eff = max(As_d);

    % Concrete Capacity
    Vc = 2*(1 + P/(2000*Ag))*lambda*b*d_eff*sqrt(fc); 

    % Reinforcement Capacity
    if Av > 0
        Vs = Av*fy*d_eff/S;
    else
        Vs = 0;
    end

    % Total Capacity
    Vn = Vc + Vs;

end
%% Ultimate Capacity
phi = 0.75;
Vu = phi*Vn;
end

