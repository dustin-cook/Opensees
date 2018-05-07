function [ Vn, V0 ] = fn_shear_capacity( Av, fy, d, s, lambda, fc, Ag, M, V, Nu )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d_eff = 0.8*d;
k = 1;
mv_ratio = min([max([M/V*d_eff,2]),4]);

V0 = Av*fy*d_eff/s + lambda*(6*sqrt(fc)/mv_ratio)*sqrt(1+Nu/(6*sqrt(fc)*Ag))*0.8*Ag;

Vn = k*V0;
end

