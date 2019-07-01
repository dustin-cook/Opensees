function [th_NC] = fn_eurocode_rotation_acceptance(L,h,b,s,Av,cov,fc,fy,P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Notes
% All acceptance criteria are in total element chord rotation

%% EUROCODE
% Convert Units
covm = cov/39.37;   % inches to m
hm = h/39.37;       % inches to m
bm = b/39.37;       % inches to m
sm = s/39.37;       % inches to m
Lm = L/39.37;       % inches to m
fcm = fc/145.04;    % PSI tp MPa
fym = fy/145.04;    % PSI tp MPa
Pm = P*4.4482/1000; % lbs to kN

% Fixed Parameters
gam_el = 1.5;
rho_d = 0;          % diagonal reinforcement 
sum_bi = 0.006;
wc = 1;
wt = 1;

% Dependant Parameters
ho = hm - 2*covm;
bo = bm - 2*covm;
rho_t = min(Av/(b*s),0.0075);
p_ratio_EU = (Pm/1000/(bm*hm*fcm));

% EQ A.2
alpha = (1-sm/(2*bo))*(1-sm/(2*ho))*(1-(sum_bi^2)/(6*ho*bo));

% EQN A.1
theta_NC = (1/gam_el)*0.016*(0.3^p_ratio_EU)*...
            ((max(0.01,wc)*fcm/max(0.01,wt))^0.225)*...
            ((Lm/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
            (1.25^(100*rho_d));




% EQN A.3
th_NC = (1/gam_el)*0.0145*(0.25^p_ratio_EU)*...
        ((max(0.01,wc)/max(0.01,wt))^0.3)*...
        (fcm^0.2)*((Lm/2/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
        (1.275^(100*rho_d));
    
% Polly Updates
    % sum Bi
    % Lv /== length
end

