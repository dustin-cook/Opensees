function [ Mu, Mn ] = fn_aci_moment_capacity( fc, b, d, As, As_d, fy, clear_cover, Es )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));
As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));

%% Calculate intermediate/defualt Values
d_c = d - clear_cover;

%% Begin Method
% Find Location of Neutral Axis
balance_eq = 0.85*fc*b*0.85*d/2;
y = 0;
step = 0.01;
tolerance = (0.85*fc*b*0.85*d/2)/1000;
count = 1;
while abs(balance_eq(count)) > tolerance
    count = count + 1;
    y(count) = y(count-1) + step;
    if y >= d_c/2
        error('Nuetral Axis of Concrete Section Not Found')
    end
    c = d_c/2-y(count);
    e_s = abs(0.003*(As_d-c)/c);
    fs = min(e_s,fy/Es)*Es;
    balance_eq(count) = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - 0.85*fc*b*0.85*c;
end

% plot(y,balance_eq)
% grid on
% xlabel('Distance from Center')
% ylabel('Balance Equation')

% Moment Capacity
Mn = sum(As.*fs.*abs(As_d-c)) + 0.85*fc*b*0.85*c*(c*1.15/2);
phi = 0.9; % Assuming tension steel strain is greater than 0.005
Mu = phi*Mn;


end

