function [ row_bal ] = fn_balanced_moment( fc, fy )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Based on EQ B-1 of ACI 318-11
% 3. Assumes Concrete at 0.002 strain and Steel at ey

%% Calculate total balanced reinforcemnt ratio accoding to EQ B-1 of ACI 318-11
row_bal = (0.85*0.85*fc/fy)*(87000/(87000+fy));

end

