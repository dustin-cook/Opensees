function [ Vu, Vn, Vs ] = fn_aci_shear_capacity( fc, b, d, Av, fy, S, lambda, Ag, hw, type, As_d, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Assumptions:

%% Begin Method
if strcmp(type,'wall') % Shear Walls
    rho_t = Av/(S*b);
    if rho_t < 0.0015 % ASCE 41-17 10.7.2.3
        error('Wall is force controlled, modify element')
    end
    
    Acv = b*d;
    aspect_ratio = hw/d;
    alpha_c = min(max((-2)*(aspect_ratio-2) + 2, 2), 3);
    Vn = Acv*(alpha_c*lambda*sqrt(fc) + rho_t*fy); % ACI 318-14 eq 18.10.4.1
    Vs = nan;

else % Beams and Columns
    % Reformat steel depth
    As_d = str2double(strsplit(strrep(strrep(As_d{1},']',''),'[',''),','));
    d_eff = max(As_d);

    % Concrete Capacity
    if strcmp(type,'beam') % Beams
        Vc = 2*lambda*b*d_eff*sqrt(fc); %eq 22.5.5.1 in ACI 318-14
    elseif strcmp(type,'column') % Columns
        Vc = 2*(1 + P/(2000*Ag))*lambda*b*d_eff*sqrt(fc); %eq 22.5.6.1 in ACI 318-14, Nu is the Axial force coresponding to the V (so time based), taken as max axial for now.
    end

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

