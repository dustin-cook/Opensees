clear all
close all
clc

%% Assumptions and Description

%% Import Packages
import aci_318.*

%% User Inputs
orientation = 'pos';
fc = 4000;
b = 14;
h = 24;
As = {'[3, 3]'};
As_d = {'[2.5,21.5]'};
fy = 60000;
Es = 29000000;
P = linspace(-350000,1150000,100);
% P = [1100000,623700,504400,0,-350000];
slab_depth = 0;
b_eff = 0;

%% Call Moment Capacity Function
for i = 1:length(P)
    [ Mu(i), Mn(i) ] = fn_aci_moment_capacity( orientation, fc, b, h, As, As_d, fy, Es, P(i), slab_depth, b_eff );
end

%% Plots results
hold on
M_book = [0,521.,559.7,297,0];
P_book = [1482,623.7,504.4,0,-360];
scatter(M_book,P_book,'filled','DisplayName',' McCormak Example 10.2')
plot(Mn/12000,P/1000,'DisplayName','Cook Algorithm')
xlabel('Mn (k-ft)')
ylabel('P (k)')
legend('Location','eastoutside')
grid on
box on