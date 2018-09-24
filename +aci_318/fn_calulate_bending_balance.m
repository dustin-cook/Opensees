function [ balance_eq, fs ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% Assumes Beta1 = 0.85

%% Calculate Stress in stell
e_s = abs(0.003*(As_d-c)/c);
fs = min(e_s,fy/Es)*Es;

%% Calculate Balance Differential (perfect balance = 0)
if slab_depth > 0 % T-beam Section
    if beta_1*c > slab_depth
        balance_eq = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - (0.85*fc*(b_eff*slab_depth) + 0.85*fc*(b*(beta_1*c-slab_depth))) + P;
    else
        balance_eq = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - 0.85*fc*(b_eff*beta_1*c) + P;
    end
else % Rectangular Section
    balance_eq = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - 0.85*fc*(b*0.85*c) + P;
end

end

