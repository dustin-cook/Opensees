function [ acc_com ] = fn_balaji_accel_filter( acc_raw, dt, Cut_Freq, order, f_type)
% Description: Low-Pass Butterworth Filter

% Created By: Balaji Paramasivam
% Modified By: Dustin Cook
% Date Modified: 3/25/19

% Inputs 
% acc_raw - Acceleration matrix. First column should be time and subsequent
%           columns are acceleration or data you want to filter.
% Cut_Freq - cut-off frequency (hertz)
% order - Order of Filterring (default to 5th order if unsure (ie int 5))
% f_type - Filter Type: 'low' | 'bandpass' | 'high' | 'stop'

% Outputs
% acc_com - Filtered Acceleration Matrix

%% Begin Method
% Basic Time History Properties
[length_eq,cols] = size(acc_raw);
n = 2^nextpow2(length_eq);

% Pad Time History with Zeros
add_points = n-length_eq;
add_time = acc_raw(end,1)*ones(add_points,1);
add_time = add_time+((1:1:add_points)*dt)';
acc_raw_1 = zeros(n,cols);
acc_raw_1(:,1) = vertcat(acc_raw(:,1),add_time);
for i = 2:cols
    acc_raw_1((1:length_eq),i) = acc_raw((1:length_eq),i);
end

%% Construct and apply a n-order, acausal, butterworth, band-pass filter to acceleration data 
for kk = 2:cols
    % Compute the Nyquist Frequency (Maximum frequency that can be constructed)
    F_Nyq=1/2/dt; 

    % Calculate the normalized cuttoff frequency 
    Wn = Cut_Freq/F_Nyq; 
    
    % Designs an order 5 bandpass digital Butterworth filter with normalized cutoff frequency Wn 
    [b1, length_eq] = butter(order, Wn, f_type);
    filtered_time_hist(:,1) = acc_raw_1(:,1); 
    filtered_time_hist(:,kk) = filtfilt(b1,length_eq,acc_raw_1(:,kk));
end

[kk1,kk2] = size(acc_raw);
time = filtered_time_hist((1:kk1),1);
acc_ff_1 = filtered_time_hist((1:kk1),(2:(kk2)));
acc_com = horzcat(time,acc_ff_1);

end % Function