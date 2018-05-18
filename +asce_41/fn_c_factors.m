function [ c1, c2 ] = fn_c_factors( site_class, num_stories, T, haz_class, DCR_max )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

if strcmp(site_class,'A') || strcmp(site_class,'B')
    a = 130;
elseif strcmp(site_class,'C')
    a = 90;
else
    a = 60;
end

% Calculate Cm
[ c_m ] = fn_cm( num_stories, T, haz_class );

u_strength = max([DCR_max*c_m/1.5,1]);

c1 = 1 + (u_strength-1)/(a*T^2);
c2 = 1 + (1/800)*((u_strength-1)/T)^2;


end

