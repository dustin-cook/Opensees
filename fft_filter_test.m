% test FFT
clear all
close all
clc

% Inputs                   
dt = 1/100;             % Sampling period       
L = 2500;             % Length of signal
t = (0:L-1)*dt;        % Time vector
low_freq = 1;
high_freq = 15;

% Define the signal
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
accel = load([pwd '\ground_motions\ICSB_1979\chan_11_accel_short.tcl']);

% hold on
% for i = 1:99
%     plot(t,accel*(1-i/100),'color',[1,1,1]*i/100)
% end
% xlim([4,14])
% Noise
% X = S + 2*randn(size(t));

% Filter With FFT
[ X_filt ] = fn_fft_accel_filter( accel', dt, t, high_freq, low_freq );