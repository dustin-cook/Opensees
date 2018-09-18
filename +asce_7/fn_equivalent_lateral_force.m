function [ V, Cs, Fx, Cvx ] = fn_equivalent_lateral_force( w, h )
% Equvalent Lateral Force Method as outlined in ASCE 7-16

%% Assumptions
% 1) T-long = 8 sec
% 2) Use Ta for period
% 3) Assume concrete moment frame
% 4) Assume modern site parameters for Imperial Valley
% 5) Importance factor = 1

%% Site Parameters
S1 = 0.667;
Sds = 1.21;
Sd1 = 0.667;

%% Concrete Moment Frame Properties (remove later to be more general)
R = 5; % assume IMF
Ie = 1.0;
Ct = 0.028; % table 12.8-1
hn = max(h);
x = 0.8;

%% Period Properties
T = Ct*hn^x;
T_L = 8;

%% Caclulate seismic response coefficient: Cs
if T <= T_L
    cs_max = Sd1/(T*(R/Ie)); % eq 12.8-3
else
    cs_max = Sd1*T_L/(T^2*(R/Ie)); % eq 12.8-4
end

if S1 >= 0.6
    cs_min = 0.5*S1/(R/Ie); % eq 12.8-6
else
    cs_min = max([0.044*Sds*Ie, 0.01]); % eq 12.8-5
end

Cs = min([max([Sds/(R/Ie),cs_min]),cs_max]); % eq 12.8-2

%% Calculate Base Shear: V
W = sum(w); % Total effective seismic weight
V = Cs*W; % Total Base shear eq 12.8-1

%% Calculate vertical distribution factor: Cvx
k = min([max([interp1([0.5,2.5],[1,2],T,'linear','extrap'),1]),2]);
Cvx = (w .* h.^k)/sum(w .* h.^k); % eq 12.8-12

%% Calculate Vertical Forces
Fx = Cvx*V; % eq 12.8-11

end

