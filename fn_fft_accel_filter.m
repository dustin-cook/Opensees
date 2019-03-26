function [ X_filt ] = fn_fft_accel_filter( X, dt, t_vec, high_freq, low_freq )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Initial Variables
Fs = 1/dt;
L = length(X);
f_vec = Fs*(0:(L/2))/L;

%% Convert to Freq Domain
% Compute the Fourier transform of the signal.
Y = fft(X);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%% Filter in feqency Domain
% Low Pass Fitler of the signal
f_new = f_vec(2:end);
Y_1 = Y(1:length(Y)/2);
Y_2 = Y(length(Y)/2+1:end);
% Low pass
Y_1(f_new>high_freq) = 0;
Y_2(f_new>high_freq) = 0;
% High pass
Y_1(f_new<low_freq) = 0;
Y_2(f_new<low_freq) = 0;

Y_filt = [Y_1,Y_2];

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum 
% P1 based on P2 and the even-valued signal length L for the filtered data.
P2_filt = abs(Y_filt/L);
P1_filt = P2_filt(1:L/2+1);
P1_filt(2:end-1) = 2*P1_filt(2:end-1);

%% Convert back to time domaun
% Inverse FFT
X_filt = real(ifft(Y_filt));

%% Plotters
% % Time Domain
% figure
% subplot(2,1,1)
% hold on
% plot(t_vec,X)
% plot(t_vec,X_filt)
% title('Time Domain')
% xlabel('Time')
% ylabel('X(t)')
% 
% % Frequency Domain
% subplot(2,1,2)
% hold on
% plot(f_vec,P1)
% plot(f_vec,P1_filt)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% hold off

end

