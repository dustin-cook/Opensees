function [ l_d, l_dt, l_dt_raw, l_dht ] = fn_aci_development_length( fy, fc, A_tr, s, d_b, c_b, lambda, psi_t, psi_e, psi_s, psi_c, psi_r, psi_rc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Only considers deformed bars (ie no wires and no plain bars)
% 2. Does not consider headed bar development length
% 3. Does not consider reduction in dev length for As provided > As req
% 4. Assume only one bar at a time (no bundles, n = 1)

%% Intermediate Calculations
n = 1; % 1 Bar being developed at a time
psi_t_e = min([psi_t*psi_e, 1.7]); % Table 25.4.2.4
K_tr = 40*A_tr/(s*n); % EQ 25.4.2.3b
fc = min([fc,10000]); % 25.4.1.4

%% Calulate Tension Develpment Length of Deformed Bars
confinement = min([(c_b + K_tr)/d_b, K_tr, 2.5]);
l_dt_raw = (3*fy*psi_t_e*psi_s*d_b)/(40*lambda*sqrt(fc)*confinement); % EQ 25.4.2.3a
l_dt = max([l_dt_raw,12]); % 25.4.2.1

%% Calulate Tension Develpment Length of Standard Hooks
l_dht = (fy*psi_e*psi_c*psi_r*d_b)/(50*lambda*sqrt(fc)); % 25.4.3.1
l_dht = max([l_dht,8*d_b,6]); % 25.4.3.1

%% Calulate Compression Develpment Length of Bars
l_dc = (fy*psi_rc*d_b)/(50*lambda*sqrt(fc)); % 25.4.9.2
l_dc = max([l_dc,0.0003*fy*psi_rc*d_b,8]); % 25.4.9.2

%% Total Development Length
l_d = max([l_dt, l_dc]);

end

