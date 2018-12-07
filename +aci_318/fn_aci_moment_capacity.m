function [ Mu, Mn ] = fn_aci_moment_capacity( orientation, fc, b, h, As, As_d, fy, Es, P, slab_depth, b_eff )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Whitney stress block
% 2. for Phi assumes columns are tied and not spiral

%% Import Packages
import aci_318.*

%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));
As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));
if strcmp(orientation,'neg')
    As = fliplr(As);
    As_d = h - fliplr(As_d);
end
[ beta_1 ] = fn_beta_1( fc );

%% Begin Method
% Find Location of Neutral Axis
y_prev = -h/2;
step = 0.1;
tolerance = 50; % lbs %(0.85*fc*b*0.85*d/2+abs(P))/5000;
balance_found = 0;
count = 0;
while balance_found == 0
    count = count + 1;
    y(count) = y_prev + step;
    y_prev = y(count);
    if y(count) >= h/2
        error('Nuetral Axis of Concrete Section Not Found')
    end
    c = h/2-y(count);
    [ balance_eq(count), fs, e_s ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
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
            c = h/2-y(count);
            [ balance_eq(count), fs, e_s ] = fn_calulate_bending_balance( c, P, As, As_d, b, b_eff, slab_depth, fy, Es, fc, beta_1 );
            if abs(balance_eq(count)) < tolerance
                balance_found = 1;
            end
        end
    end
end

% Nominal Moment Capacity
a = beta_1*c;
if slab_depth > 0 % T-beam Section
    if a > slab_depth
        Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b_eff*slab_depth)*(h-slab_depth)/2 + (0.85*fc*b*(a-slab_depth))*((h - 3*a - slab_depth)/2);
    else
        Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b_eff*a)*(h-a)/2;
    end
else % Rectangular Section
    Mn = sum(As.*fs.*(h/2 - As_d)) + (0.85*fc*b*a)*(h-a)/2;
end

%% Phi and Design Capacity
if min(e_s) < -0.005 % check largest tension strain (tension is negative)
    phi = 0.9;
else
    phi = max(min((abs(min(e_s))-0.002)*(0.9-0.65)/(0.005-0.002)+0.65,0.9),0.65); % Interpolate for phi
end
Mu = phi*Mn;

%% Validation Plots 
if 1 == 0
    plot(y,balance_eq)
    grid on
    xlabel('Distance from Center')
    ylabel('Balance Equation')
end


end

