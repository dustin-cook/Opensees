function[V_NC] = fn_eurocode_column_shear_acceptance(h,b,d,d_prm,s,As,Av,db,fc,fy,P,M,V,deform_pl)    
% EUROCODE Limit state rotation determiniation

% Only for primary seismic columns controlled by shear

% Convert Standard to metric
hm = h/39.37;         % m
bm = b/39.37;         % m
dm = d/39.37;         % m
dm_prm = d_prm/39.37; % m
fcm = fc/145.04;      % MPa
fym = fy/145.04;      % MPa
N = (P*4.4482)/1000000; % MN
dbm = db/(39.37^2);   % m^2

% Fixed quantities
gam_el = 1.15; % 1.15 for primary seismic elements
x = hm/2;      % Compression Zone Depth (assume equal to half of section for now)

% Calculated Quantities
Lv = (M/V)/39.37;    % Moment to shear demand ratio (convertedc to meters)
rho_tot = As/(b*d);  % Total longitudinal reinforcement
rho_w = Av/(b*s);    % Transverse reinforcement
z = dm-dm_prm;       % length of internal lever arm (m)
Vw = rho_w*b*z*fym;  % EQ A.13 - contribution of transverse reinforcement to shear resistance
av = 1;              % Assume shear cracking to preceede flexure cracking
phi_y = 0.003 / x;   % Yeild Curvature (assume is same as concrete crushing over compresion zone depth) 
TH_y = phi_y * (Lv+av*z)/3 + 0.0013*(1+1.5*hm/Lv) + 0.13*phi_y*dbm*fym/sqrt(fcm); % EQ A.10b
u_delta_pl = deform_pl / TH_y;

% Calculate the shear strength acceptance criteria
V_NC_m = (1/gam_el)*(...% EQ A.12
         ((hm-x)/(2*Lv))*min([N,0.55*dm*bm*fcm])+...
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

% Convert back to Standard Units
V_NC = (V_NC_m*1000000)/4.4482;

end