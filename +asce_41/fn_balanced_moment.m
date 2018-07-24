function [ M_bal, row_bal ] = fn_balanced_moment( fc, fy, b, d, As_d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Does not consider axial loads (ie only for beams)
% 2. Assumes only 1 row of tension steel
% 3. Concrete yeild at 0.003 and steel yeild at 0.00207 strain

%% Maniluptae String Vector Properties to Arrays
As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));

%% Calculate depth of nuetral axis
d_prime = max(As_d);
c = d_prime*(0.003)/(0.003 + 0.00207);
a = 0.85*c;

%% Calculate total tenstion steel for balanced moment
As_bal = (0.85)*a*b*(fc/fy);
row_bal = As_bal/(b*d);

%% Calculate Balanced Moment capacity
M_bal = As_bal*fy*(d_prime - a/2);

end

