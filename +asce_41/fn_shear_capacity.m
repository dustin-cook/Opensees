function [ Vn, V0 ] = fn_shear_capacity( Av, fy, As_d, s, lambda, fc, Ag, M, V, Nu, DCR_max )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

As_d = str2double(strsplit(strrep(strrep(As_d{1},']',''),'[',''),','));
d_eff = max(As_d);
mv_ratio = min([max([M/(V*d_eff),2]),4]);

% Set alpha factor
if s/d_eff <= 0.75
    alpha = 1;
elseif s/d_eff >= 1
    alpha = 0;
else
    alpha = interp1([0.75,1],[1,0],s/d_eff);
end

% Set K factor
if DCR_max <= 2
    k = 1;
elseif DCR_max >= 6
    k = 0.7;
else
    k = interp1([2,6],[1,0.7],DCR_max);
end

% Calculate shear based on EQ 10-3 from ASCE 41-17
V0 = alpha*Av*fy*d_eff/s + lambda*(6*sqrt(fc)/mv_ratio)*sqrt(1+Nu/(6*sqrt(fc)*Ag))*0.8*Ag;
Vn = k*V0;

end

