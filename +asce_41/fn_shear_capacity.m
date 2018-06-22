function [ Vn, V0 ] = fn_shear_capacity( Av, fy, d, s, lambda, fc, Ag, M, V, Nu )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d_eff = 0.8*d;
k = 1;
mv_ratio = min([max([M/V*d_eff,2]),4]);

if s/d_eff <= 0.75
    alpha = 1;
elseif s/d_eff >= 1
    alpha = 0;
else
    alpha = interp1([0.75,1],[1,0],s/d_eff);
end

V0 = alpha*Av*fy*d_eff/s + lambda*(6*sqrt(fc)/mv_ratio)*sqrt(1+Nu/(6*sqrt(fc)*Ag))*0.8*Ag;

Vn = k*V0;
end

