% Calculate Development Length according to chapter 25.4 ACI 318-14

%% Define Inputs
psi_t = 1.3; % for top bars in members deeper than 12 in, Table 25.4.2.4
psi_e = 1.0; % for non coated bars, Table 25.4.2.4
psi_s = 1.0; % for no. 7 bars or greater, Table 25.4.2.4
lambda = 1.0; % for normal weigth concrete, Table 25.4.2.4
fy = 50000; % in psi
fc = 75000; % in psi
A_tr = 0.73; % total tranverse reinforcement area
s = 2; % transvers reinforcement spacing
n = 3; % number of bars being developed or spliced
c_b = 2; % full distance to edge of section or half distance to center of next bar
d_b = 0.5; % Diameter of bar

%% Intermediate Calculations
psi_t_e = min([psi_t*psi_e, 1.7]); % Table 25.4.2.4
K_tr = 40*A_tr/(s*n); % EQ 25.4.2.3b
fc = min([fc,10000]); % 25.4.1.4

%% Calulate Develpment Length
confinement = min([(c_b + K_tr)/d_b, K_tr, 2.5]);
l_d = (3*fy*psi_t_e*psi_s*d_b)/(40*lambda*sqrt(fc)*confinement); % EQ 25.4.2.3a
l_d = max([l_d,12]); % 25.4.2.1