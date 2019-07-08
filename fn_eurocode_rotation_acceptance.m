function [theta_NC, theta_SD, theta_DL] = fn_eurocode_rotation_acceptance(h,b,d,d_prm,s,As,Av,db,cov,fc,fy,P,M,V)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Notes
% All acceptance criteria are in total element chord rotation

%% EUROCODE
% Convert Units
covm = cov/39.37;     % inches to m
hm = h/39.37;         % inches to m
bm = b/39.37;         % inches to m
dm = d/39.37;         % inches to m
dm_prm = d_prm/39.37; % inches to m
sm = s/39.37;         % inches to m
dbm = db/39.37;         % inches to m
fcm = fc/145.04;      % PSI tp MPa
fym = fy/145.04;      % PSI tp MPa
Pm = P*4.4482/1e6;   % lbs to MN

% Fixed Parameters
gam_el = 1.5;
rho_d = 0;                            % diagonal reinforcement 
bi = [10,10,10,10,10,10,10,10]/39.37; % Hardcoded from ICSB column D-S11

% Dependant Parameters
wc = 0.5*As/(h*d);  % compression rein ratio (assume half of total reinformcent)
wt = 0.5*As/(h*d);  % tension rein ratio (assume half of total reinformcent)
Lv = (M/V)/39.37;    % Moment to shear demand ratio (convertedc to meters)
ho = hm - 2*covm;
bo = bm - 2*covm;
rho_t = Av/(b*s);
p_ratio_EU = (Pm/(bm*hm*fcm));

% EQ A.2
alpha = (1-sm/(2*bo))*(1-sm/(2*ho))*(1-(sum(bi.^2))/(6*ho*bo));

% EQN A.1
theta_NC = (1/gam_el)*0.016*(0.3^p_ratio_EU)*...
            ((max(0.01,wc)*fcm/max(0.01,wt))^0.225)*...
            ((Lv/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
            (1.25^(100*rho_d));

% A.3.2.3
theta_SD = 0.75*theta_NC;

% EQN A.10b
zm = dm-dm_prm;
av = 1; % shear cracking expected before flexure
phi_y = 0.003 / (hm/2);   % Yeild Curvature (assume is same as concrete crushing over half of the member depth)
theta_DL = phi_y*(Lv + av*zm)/3 + 0.0013*(1+1.5*(hm/Lv)) + 0.13*phi_y*dbm*fym/sqrt(fcm); % can just use the calculated yeild rotation instead of this is wanted

% % EQN A.3
% th_NC = (1/gam_el)*0.0145*(0.25^p_ratio_EU)*...
%         ((max(0.01,wc)/max(0.01,wt))^0.3)*...
%         (fcm^0.2)*((Lm/2/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
%         (1.275^(100*rho_d));
    
% Polly Updates
    % sum Bi
    % Lv /== length
    % wc and wt /== 1
end

