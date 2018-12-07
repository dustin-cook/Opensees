function [ balance_eq, fs, e_s ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% Uses whitney stress block

%% Calculate Stress in steel
e_s = 0.003*(c-As_d)/c; % Positive = compression strain, negative = tension strain
fy_eff = fy - 0.85*fc; % Effective compresison steel stress which accounts for displaced concrete.
fs = max(min(e_s*Es,fy_eff),-fy); % Check both positive and negative yield stress

%% Whitney stress block depth
a = beta_1*c;

%% Calculate Balance Differential (perfect balance = 0)
if slab_depth > 0 % T-beam Section
    if a > slab_depth % Compresion zone is in the web
        balance_eq = sum(As.*fs) + 0.85*fc*b_eff*slab_depth + 0.85*fc*b*(a-slab_depth) - P;
    else % Compresion zone is all in the flange
        balance_eq = sum(As.*fs) + 0.85*fc*b_eff*a - P;
    end
else % Rectangular Section
    balance_eq = sum(As.*fs) + 0.85*fc*b*a - P;
end

end

