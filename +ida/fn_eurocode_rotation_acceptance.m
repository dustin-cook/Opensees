function [theta_NC, theta_SD, theta_DL] = fn_eurocode_rotation_acceptance(h,b,d,d_prm,cov,s,As_comp,As_ten,Av,db,bi,fc,fy,P,M,V,av,phi_y)
% Function to calculate the Eurocode rotational limit states for a
% reinforced concrete column. Concrete column is input using standard units
% of inches and lbs. Limit state rotation outputs defined in terms of total 
% rotaion. Methods are based on Annex A of Eurocode 8.

%   Created By: Polly Murray
%   Modified By: Dustin Cook
%   Date Last Modified: September 10th, 2019

% INPTUS
% h = Depth of the column gross cross section (inches)
% b = Width of the column gross cross section (inches)
% d = Depth to the center of the bottom longitudinal bars (inches)
% d_prm = Depth to the center of the top longitudinal bars (inches)
% s = Spacing of the transverse reinforcement (inches)
% As_comp = Total Area of compressive longitudinal reinforcement (inches^2)
% As_ten = Total Area of tensile longitudinal reinforcement (inches^2)
% Av = Cross sectional area of transverse reinforcement (inches^2)
% db = Average diameter of tensile longitudinal bars (inches)
% cov = Depth of Concrete cover to the centerline of the hoop(inches)
% bi = vector of centerline spacing of longitudinal bars laterally restrained
%      by a stirrup corner or a cross-tie along the perimeter of the 
%      cross-section.(inches)
% fc = Concrete compressive strength (psi)
% fy = Steel strength (psi)
% P = Maximum axial demand from analysis (lbs)
% M = Maximum moment demand from analysis (lbs-in)
% V = Maximum shear demand from analysis (lbs)
% av = 1 if shear cracking expected before flexure, av = 0 othewise 
% phi_y = Yeild Curvature of the end section (in/in^2)

% OUTPUS
% theta_NC = Near Collapse Rotational Limit State (rads of total rotation)
% theta_SD = Significant Damage Rotational Limit State (rads of total rotation)
% theta_DL = Damage Limitation Rotational Limit State (rads of total rotation)

% ASSUMPTIONS/NOTES
% 1) All acceptance criteria are in total element chord rotation
% 2) Only applies to primary seismic elements
% 3) Assumes there is no diagonal reinforcement

%% Convert Units from standard to metric
covm = cov/39.37;      % inches to m
hm = h/39.37;          % inches to m
bm = b/39.37;          % inches to m
dm = d/39.37;          % inches to m
dm_prm = d_prm/39.37;  % inches to m
sm = s/39.37;          % inches to m
dbm = db/39.37;        % inches to m
fcm = fc/145.04;       % PSI tp MPa
fym = fy/145.04;       % PSI tp MPa
bim = bi/39.37;        % inches to m
phi_y_m = phi_y*39.37; % in/in^2 tot m/m^2

%% Defined Fixed Parameters
gam_el = 1.5;                         % Defined as 1.5 for primary seismic elements according to A.3.2.2
rho_d = 0;                            % ratio of diagonal reinforcement 

%% Calculate Dependant Parameters
wc = As_comp/(h*d);        % compression rein ratio
wt = As_ten/(h*d);         % tension rein ratio
Lv = (M/V)/39.37;          % Moment to shear demand ratio (converted to meters)
ho = hm - 2*covm;          % depth of the concrete core, to the centerline of the hoop
bo = bm - 2*covm;          % width of the concrete core, to the centerline of the hoop
rho_t = Av/(b*s);          % transverse steel rein ratio
p_ratio = (P/(b*h*fc)); % axial load ratio

% EQ A.2
alpha = (1-sm/(2*bo))*(1-sm/(2*ho))*(1-(sum(bim.^2))/(6*ho*bo));

%% Calculate Near Collapse Limit State (EQN A.1)
theta_NC = (1/gam_el)*0.016*(0.3^p_ratio)*...
            ((max(0.01,wc)*fcm/max(0.01,wt))^0.225)*...
            ((Lv/hm)^0.35)*(25^(alpha*rho_t*(fym/fcm)))*...
            (1.25^(100*rho_d));

%% Calculate Significant Damage Limit State (A.3.2.3)
theta_SD = 0.75*theta_NC;

%% Calculate Damage Limitation Limit State (EQN A.10b)
zm = dm-dm_prm;
theta_DL = phi_y_m*(Lv + av*zm)/3 + 0.0013*(1+1.5*(hm/Lv)) + 0.13*phi_y_m*dbm*fym/sqrt(fcm); % can just use the calculated yeild rotation instead of this is wanted

end
