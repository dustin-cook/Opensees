clear all 
close all
clc

% Load Accel TH
NS_ground_floor_east = dlmread('chan_11_accel.tcl','\n');
NS_ground_floor_west = dlmread('chan_10_accel.tcl','\n');

% Time vector
x = 0.01:0.01:30;

% Plot Accel TH
figure
% subplot(2,1,1)
hold on
plot(x,NS_ground_floor_east(1:length(x)),'LineWidth',1.4,'DisplayName','Ground Floor East')
plot(x,NS_ground_floor_west(1:length(x)),'LineWidth',1.4,'DisplayName','Ground Floor West')
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('EW Frame Direction')
legend('location','northeast')
grid on
box on

% subplot(2,1,2)
% hold on
% plot(x,NS_building(1:length(x)),'LineWidth',1.4,'DisplayName','Ground Floor')
% % plot(x,EW_ground(1:length(x)),'LineWidth',1.4,'DisplayName','Free Field')
% xlabel('Time (s)')
% ylabel('Acceleration (g)')
% title('NS Wall Direction')
% grid on
% box on
