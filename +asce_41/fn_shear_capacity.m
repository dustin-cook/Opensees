function [ Vn, V0 ] = fn_shear_capacity( Av, fy, d_eff, s, lambda, fc, Ag, M_TH, V_TH, Nug, ductility_factor )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mv_ratio = min([max([max(abs(M_TH./(V_TH*d_eff))),2]),4]);

% If Nu is in tension, set to zero
Nug = max([Nug,0]); % based on gravity load of asce 41-17 eq 7-3

% Set alpha factor
if s/d_eff <= 0.75
    alpha = 1;
elseif s/d_eff >= 1
    alpha = 0;
else
    alpha = interp1([0.75,1],[1,0],s/d_eff);
end

% Set K factor
if ductility_factor <= 2
    k = 1;
elseif ductility_factor >= 6
    k = 0.7;
else
    k = interp1([2,6],[1,0.7],ductility_factor); % The Ductility Factor This is max DCR for linear and displacement ductility for nonlinear.
end
% k = 0.9; % hardcode if you want

% Calculate shear based on EQ 10-3 from ASCE 41-17
V0 = alpha*Av*fy*d_eff/s + lambda*(6*sqrt(fc)/mv_ratio)*sqrt(1+Nug/(6*sqrt(fc)*Ag))*0.8*Ag;
Vn = k*V0;

end

