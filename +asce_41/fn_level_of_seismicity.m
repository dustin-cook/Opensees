function [ seismicity ] = fn_level_of_seismicity( spectra )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Sds = interp1(spectra.period,spectra.psa_5,0.2);
Sd1 = interp1(spectra.period,spectra.psa_5,1);


if Sds >= 0.5 || Sd1 >=0.2
    seismicity = 'high';
elseif Sds >= 0.33 || Sd1 >=0.133
    seismicity = 'moderate';
elseif Sds >= 0.167 || Sd1 >=0.067
    seismicity = 'low';
else
    seismicity = 'very low';
end
end

