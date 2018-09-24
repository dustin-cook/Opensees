function [ Mu, Mn ] = fn_aci_moment_capacity( orientation, fc, b, d, As, As_d, fy, Es, P, slab_depth, b_eff )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Assuming tension steel strain is greater than 0.005 (tension controlled and phi = 0.9)

%% Import Packages
import aci_318.*

%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));
As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));
if strcmp(orientation,'neg')
    As = fliplr(As);
    As_d = d - fliplr(As_d);
end
[ beta_1 ] = fn_beta_1( fc );

%% Begin Method
% Find Location of Neutral Axis
y_prev = -d/2;
step = 0.1;
tolerance = 50; % lbs %(0.85*fc*b*0.85*d/2+abs(P))/5000;
balance_found = 0;
count = 0;
while balance_found == 0
    count = count + 1;
    y(count) = y_prev + step;
    y_prev = y(count);
    if y(count) >= d/2
        error('Nuetral Axis of Concrete Section Not Found')
    end
    c = d/2-y(count);
    [ balance_eq(count), fs ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
    if abs(balance_eq(count)) < tolerance
        balance_found = 1;
    elseif count > 1 && (sign(balance_eq(count)) ~= sign(balance_eq(count-1)))
        step = 0.00002;
        switch_index = count;
        while balance_found == 0
            count = count + 1;
            y(count) = y_prev - step;
            y_prev = y(count);
            if count > 5000 + switch_index
                error('Nuetral Axis of Concrete Section Not Found')
            end
            c = d/2-y(count);
            [ balance_eq(count), fs ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
            if abs(balance_eq(count)) < tolerance
                balance_found = 1;
            end
        end
    end
end

% plot(y,balance_eq)
% grid on
% xlabel('Distance from Center')
% ylabel('Balance Equation')

% Moment Capacity
if slab_depth > 0 % T-beam Section
    if beta_1*c > slab_depth
        Mn = sum(As.*fs.*abs(As_d-c)) + 0.85*fc*(b_eff*slab_depth)*(c-slab_depth/2) + 0.85*fc*(b*(beta_1*c-slab_depth))*(c*(1-beta_1/2)-slab_depth/2) + P*y(end);
    else
        Mn = sum(As.*fs.*abs(As_d-c)) + 0.85*fc*(b_eff*beta_1*c)*c*(1-beta_1/2) + P*y(end);
    end
else % Rectangular Section
    Mn = sum(As.*fs.*abs(As_d-c)) + 0.85*fc*(b*beta_1*c)*c*(1-beta_1/2) + P*y(end);
end
phi = 0.9; % Assuming tension steel strain is greater than 0.005
Mu = phi*Mn;


end

