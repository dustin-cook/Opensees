function [ Mu, Mn ] = fn_aci_moment_capacity( fc, b, d, As, As_d, fy, clear_cover, Es, P )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Inital Setup
As = str2double(strsplit(strrep(strrep(As{1},'[',''),']',''),','));
As_d = str2double(strsplit(strrep(strrep(As_d{1},'[',''),']',''),','));

%% Calculate intermediate/defualt Values
d_c = d - clear_cover;

%% Begin Method
% Find Location of Neutral Axis
y_prev = -d/2;
step = 0.5;
tolerance = (0.85*fc*b*0.85*d/2+abs(P))/1000;
balance_found = 0;
count = 0;
while balance_found == 0
    count = count + 1;
    y(count) = y_prev + step;
    y_prev = y(count);
    if y(count) >= d_c/2
        error('Nuetral Axis of Concrete Section Not Found')
    end
    c = d_c/2-y(count);
    e_s = abs(0.003*(As_d-c)/c);
    fs = min(e_s,fy/Es)*Es;
    balance_eq(count) = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - 0.85*fc*b*0.85*c + P;
    if abs(balance_eq(count)) < tolerance
        balance_found = 1;
    elseif count > 1 && (sign(balance_eq(count)) ~= sign(balance_eq(count-1)))
        step = 0.01;
        switch_index = count;
        while balance_found == 0
            count = count + 1;
            y(count) = y_prev - step;
            y_prev = y(count);
            if count > 50 + switch_index
                error('Nuetral Axis of Concrete Section Not Found')
            end
            c = d_c/2-y(count);
            e_s = abs(0.003*(As_d-c)/c);
            fs = min(e_s,fy/Es)*Es;
            balance_eq(count) = sum(As.*fs.*((As_d-c)./abs(As_d-c))) - 0.85*fc*b*0.85*c + P;
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
Mn = sum(As.*fs.*abs(As_d-c)) + 0.85*fc*b*0.85*c*(c*1.15/2) + P*y(end);
phi = 0.9; % Assuming tension steel strain is greater than 0.005
Mu = phi*Mn;


end

