function [ V, Cs, Fx, Cvx ] = fn_equivalent_lateral_force( w, h, T_elastic, S1, Sds, Sd1, Ie, run_drifts )
% Equvalent Lateral Force Method as outlined in ASCE 7-16

%% Assumptions
% 1) T-long = 8 sec
% 2) Use Ta for period
% 3) Assume concrete moment frame
% 4) Assume modern site parameters for Imperial Valley
% 5) Importance factor = 1

%% Site Parameters
% ICSB
% S1 = 0.667;
% Sds = 1.21;
% Sd1 = 0.667;

%% Concrete Moment Frame Properties (remove later to be more general)
% ICSB
% R = 5; % assume IMF
% Ie = 1.0;
% Ct = 0.016; % table 12.8-2
% x = 0.9; % table 12.8-2
% hn = max(h);

% Fixed Assumptions
T_L = 8;

% RC Frame Archetypes
R = 8; % assume IMF
% Ie = 1.0;
Ct = 0.016; % table 12.8-2
x = 0.9; % table 12.8-2
hn = max(h)/12;

%% Period Properties
if run_drifts
    T = T_elastic;
else
    % Coefficient for upper limit on calculated period
    % table 12.8-1
    Cu = interp1([0.1, 0.15, 0.2, 0.3],...
                 [1.7,1.6,1.5,1.4],...
                 min(max(Sds,0.1),0.3));
    Ta = Ct*hn^x;
    CuTa = Cu*Ta;
    T = min(T_elastic,CuTa);
end

%% Caclulate seismic response coefficient: Cs
Cs_short = Sds/(R/Ie); % eq 12.8-2

if T <= T_L
    Cs_T = Sd1/(T*(R/Ie)); % eq 12.8-3
else
    Cs_T = Sd1*T_L/(T^2*(R/Ie)); % eq 12.8-4
end

Cs = min(Cs_short,Cs_T); % Cs need not exceed Cs_T

% eq 12.8-5
if run_drifts
    cs_min_1 = 0; % not requred for drift assessment
else
    cs_min_1 = max([0.044*Sds*Ie, 0.01]); 
end  

% eq 12.8-6
if S1 >= 0.6
    cs_min_2 = 0.5*S1/(R/Ie);
else
    cs_min_2 = 0;
end

Cs = max([Cs, cs_min_1, cs_min_2]); % eq 12.8-2

%% Calculate Base Shear: V
W = sum(w); % Total effective seismic weight
V = Cs*W; % Total Base shear eq 12.8-1

%% Calculate vertical distribution factor: Cvx
k = min([max([interp1([0.5,2.5],[1,2],T,'linear','extrap'),1]),2]);
Cvx = (w .* h.^k)/sum(w .* h.^k); % eq 12.8-12

%% Calculate Vertical Forces
Fx = Cvx*V; % eq 12.8-11

end

