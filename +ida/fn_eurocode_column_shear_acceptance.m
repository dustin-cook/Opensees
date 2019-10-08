function[V_NC] = fn_eurocode_column_shear_acceptance(h,b,d,d_prm,s,As,Av,db,fc,fy,P,M,V,deform_pl,phi_y,d_NA)    
% Function to calculate the Eurocode shear strength limit states for a
% reinforced concrete column. Concrete column is input using standard units
% of inches and lbs. Limit state output is defined as the shear demand
% limit. Methods are based on Annex A of Eurocode 8.

%   Created By: Dustin Cook
%   Date Last Modified: September 10th, 2019

% INPTUS
% h = Depth of the column gross cross section (inches)
% b = Width of the column gross cross section (inches)
% d = Depth to the center of the bottom longitudinal bars (inches)
% d_prm = Depth to the center of the top longitudinal bars (inches)
% s = Spacing of the transverse reinforcement (inches)
% As = Total Area of longitudinal reinforcement (inches^2)
% Av = Cross sectional area of transverse reinforcement (inches^2)
% db = Average diameter of tensile longitudinal bars (inches)
% cov = Depth of Concrete cover to the centerline of the hoop(inches)
% fc = Concrete compressive strength (psi)
% fy = Steel strength (psi)
% P = Maximum axial demand from analysis (lbs)
% M = Maximum moment demand from analysis (lbs-in)
% V = Maximum shear demand from analysis (lbs)
% deform_pl = plastic deformation in terms of rotation
% phi_y = Yeild Curvature of the end section (in/in^2)
% d_NA = Depth of Nuetral Axis (inches)

% OUTPUS
% V_NC = Near Collapse Shear Resistance Limit State (lbs(

% ASSUMPTIONS/NOTES
% 1) Only applies to primary seismic elements

%% Convert Standard to metric
hm = h/39.37;           % m
bm = b/39.37;           % m
dm = d/39.37;           % m
dm_prm = d_prm/39.37;   % m
fcm = fc/145.04;        % MPa
fym = fy/145.04;        % MPa
N = (P*4.4482)/1000000; % MN
dbm = db/(39.37^2);     % m^2
phi_y_m = phi_y*39.37;  % in/in^2 tot m/m^2
d_NA_m = d_NA/39.37;    % m

%% Define Fixed Parameters
gam_el = 1.15; % 1.15 for primary seismic elements

%% Calculate Dependant Parameters
Lv = (M/V)/39.37;    % Moment to shear demand ratio (convertedc to meters)
rho_tot = As/(b*d);  % Total longitudinal reinforcement
rho_w = Av/(b*s);    % Transverse reinforcement
z = dm-dm_prm;       % length of internal lever arm (m)
Vw = rho_w*b*z*fym;  % EQ A.13 - contribution of transverse reinforcement to shear resistance
av = 1;              % Assume shear cracking to preceede flexure cracking
TH_y = phi_y_m * (Lv+av*z)/3 + 0.0013*(1+1.5*hm/Lv) + 0.13*phi_y_m*dbm*fym/sqrt(fcm); % EQ A.10b
u_delta_pl = deform_pl / TH_y; % plastic ductility demand

%% Calculate the shear strength acceptance criteria
V_NC_m = (1/gam_el)*(...% EQ A.12
         ((hm-d_NA_m)/(2*Lv))*min([N,0.55*dm*bm*fcm])+...
         (1-0.05*min([5,u_delta_pl]))*(...
         0.16*max([0.5,100*rho_tot])*...
         (1-0.16*min([5,Lv/h]))*sqrt(fc)*dm*bm + Vw...
         )...
         );

if Lv/hm <= 2
    delta = atan(hm/(2*Lv));
    Vr_max = ((4/7)*(1-0.02*min([5,u_delta_pl]))/gam_el) * ... % EQ A.16
                (1+1.35*(N/(bm*dm*fcm))) * ...
                (1+0.45*100*rho_tot) * ...
                sqrt(min([40,fcm]))*bm*z*sin(2*delta);
    V_NC_m = min([V_NC_m,Vr_max]);
end

%% Convert back to Standard Units
V_NC = (V_NC_m*1000000)/4.4482;

end