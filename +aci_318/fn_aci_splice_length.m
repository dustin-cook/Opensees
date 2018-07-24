function [ l_s ] = fn_aci_splice_length( fy, fc, d_b, l_dt_raw )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Does not consider lap splice lengths of differnt size bars
% 1. Assumes As provided / As req < 2.0

%% Calculate Splice Length Per ACI 318-14
if d_b > 1.41
    error('Bar diameter too large (greater than # 11) for splice calculation. Section 25.5.1.1')
end

% Table 25.5.2.1, assuming As provided is less than 2*As required
l_st = max([1.3*l_dt_raw,12]);

% Compresson lap splice, Section 12.5.5.1
if fy >= 60000
    l_sc = max([0.0005*fy*d_b,12]);
else
    l_sc = max([(0.0009*fy-24)*d_b,12]);
end

if fc < 3000
    l_sc = l_sc*(4/3);
end

%% Total Splice Length
l_s = max([l_st,l_sc]);

end

