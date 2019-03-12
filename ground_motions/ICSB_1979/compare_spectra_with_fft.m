clear all
close all
clc

%% Load in data
% Time signal
gm_NS_bld = dlmread('chan_11_accel.tcl','\n');
gm_EW_bld = dlmread('chan_13_accel.tcl','\n');
gm_NS_grd = dlmread('gm_ns_free_field.tcl','\n');
gm_EW_grd = dlmread('gm_ew_free_field.tcl','\n');

% Spectra
spectra_NS_bld = readtable('spectra_chan_11_accel.csv','ReadVariableNames',true);
spectra_EW_bld = readtable('spectra_chan_13_accel.csv','ReadVariableNames',true);
spectra_NS_grd = readtable('spectra_ns.csv','ReadVariableNames',true);
spectra_EW_grd = readtable('spectra_ew.csv','ReadVariableNames',true);

%% Define Inputs
n = size(gm_EW_grd,2);
dt = 0.01;
time = (0:length(gm_EW_grd)-1)*dt;
Fs = 1/dt;
freq = Fs*(0:(n/2))/n;
period = 1./freq;

%% use FFT to convert signal to frequency domain
fft_spectra_raw = fft(gm_EW_grd);
P2 = abs(fft_spectra_raw/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P1_scale = P1*max(abs(gm_EW_grd))/P1(1);

hold on
plot(period,P1)
plot(spectra_EW_grd.period,spectra_EW_grd.psa_1)
xlim([0,5])