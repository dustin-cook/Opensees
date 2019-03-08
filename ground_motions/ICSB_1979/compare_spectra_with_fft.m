clear all
close all
clc

%% Load in data
NS_building = dlmread('chan_11_accel.tcl','\n');
EW_building = dlmread('chan_13_accel.tcl','\n');
NS_ground = dlmread('gm_ns_free_field.tcl','\n');
EW_ground = dlmread('gm_ew_free_field.tcl','\n');

n = size(NS_ground,2)/2;
time = (1:length(EW_ground))*0.01;
freq = 1 ./ time;%0:79/(2*n*dt);


%% use FFT to convert signal to frequency domain
spectra_raw = fft(EW_ground);
spectra = abs(real(spectra_raw))/n;

plot(freq,spectra)